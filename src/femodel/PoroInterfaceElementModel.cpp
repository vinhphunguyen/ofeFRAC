/*
 * 
 *  Copyright (C) 2016 Monash University. All rights reserved.
 *  
 *  This class implements a model of two dimensional poro-cohesive interface elements. 
 *  To be used with the HygroMechanicalModel for the bulk (fully saturated medium)
 *  
 *  - bottom and upper faces: displacement and pore pressure dofs.
 *  - mid-plane: fracutring fluid pressure 
 *  The so-called triple-noded interface elements.
 *  Different to conventional interface elements in Nint, which is 
 *  now [Nint -Nint] i.e. a matrix of 2x8 for linear interface elements.
 *
 *  Author: V.P. Nguyen, phu.nguyen@monash.edu
 *  Date: 27 September 2016
 *
 * Developement status:
 *
 * 19 October 2016: 
 *  - done with the constraint of pressure fields.
 *
 */

#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <jem/base/System.h>
#include <jem/base/Array.h>
#include <jem/util/Properties.h>
#include <jem/util/StringUtils.h>
#include <jem/util/Event.h>
#include <jem/io/FileWriter.h>
#include <jem/io/PrintWriter.h>

#include <jive/util/XDofSpace.h>
#include <jive/util/Printer.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/fem/Globdat.h>
#include <jive/fem/NodeGroup.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/geom/BShapeFactory.h>
#include <jive/geom/IShapeFactory.h>
#include <jive/geom/InterfaceShape.h>
#include <jive/geom/InternalShape.h>
#include <jive/geom/BoundaryShape.h>
#include <jive/geom/ParametricEdge.h>
#include <jive/geom/StdLine.h>
#include <jive/mp/ItemMask.h>


#include "PoroInterfaceElementModel.h"
#include "util/utilities.h"


using std::cout;

using namespace jem;
using jem::newInstance;
using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jem::io::PrintWriter;
using jem::io::FileWriter;
using jive::model::StateVector;
using jive::util::Printer;
using jive::mp::ItemMask;
using jive::Cubix;
using jive::fem::NodeGroup;


typedef MatmulChain<double,1>      MChain1;
typedef MatmulChain<double,2>      MChain2;
typedef MatmulChain<double,3>      MChain3;
typedef MatmulChain<double,4>      MChain4;


// -----------------------------------------------------------------
//   static data
// -----------------------------------------------------------------

const char* PoroInterfaceElementModel::SHAPE_PROP     = "shape";
const char* PoroInterfaceElementModel::MATERIAL_PROP  = "material";
const char* PoroInterfaceElementModel::MESH_PROP      = "meshFile";
const char* PoroInterfaceElementModel::THICKNESS_PROP = "thickness";
const char* PoroInterfaceElementModel::DISP_NODES_PROP= "dNodeGroup";
const char* PoroInterfaceElementModel::DIRICHLET_NODES_PROP= "dirNodeGroup";
const char* PoroInterfaceElementModel::BIOT_COEF_PROP    = "alpha";
const char* PoroInterfaceElementModel::VISCOSITY_PROP    = "mu";
const char* PoroInterfaceElementModel::DTIME_PROP        = "dtime";

// -----------------------------------------------------------------
//   constructor
// -----------------------------------------------------------------

PoroInterfaceElementModel::PoroInterfaceElementModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat ) : Super(name)
{
  using jem:: util::Event;
  using jive::util::joinNames;
  using namespace jive::geom;
  using jive::StringVector;
  using jive::fem::NodeGroup;

  const String context = getContext();

  Properties    myProps = props.getProps   ( myName_ );
  Properties    myConf  = conf .makeProps  ( myName_ );
  Properties    shProps = myProps.getProps ( SHAPE_PROP );
 
  // build nodes_ and dofs_ and cons_ 
  
  egroup_     = ElementGroup::get ( myConf, myProps, globdat, context );
  elems_      = egroup_.getElements ( );
  nodes_      = elems_.getNodes     ( );
  ielemCount_ = egroup_.size        ( );
  ielems_.resize ( ielemCount_ );
  ielems_     = egroup_.getIndices  ( );
  
  //nodes_   = NodeSet::find    ( globdat                   );
  dofs_       = XDofSpace::get   ( nodes_.getData(), globdat );
  cons_       = Constraints::get ( dofs_, globdat            );
  
  // dofs
  // displacement and pore pressure dofs already defined in HygroMechanicalModel
  // so just retrieve them

  rank_      = nodes_.rank   ( );
  
  aaDofTypes_.resize ( rank_ );        // displacement dof types
  pwDofTypes_.resize ( 1 );            // pore pressure dof type
  pfDofTypes_.resize ( 1 );            // fracturing pressure dof type

  aaDofTypes_[0] = 0;
  aaDofTypes_[1] = 1;
  pwDofTypes_[0] = 2;
    
  // have to define new dof type for the fracturing fluid pressure 

  pfDofTypes_[0] = dofs_->addType ( "dpf" );

  // PARALLEL: filter our element indices array ielems_ so that each processor process
  // unique elements.

  idx_t i,j;

  Ref<ItemMask> mask = ItemMask::get ( elems_.getData(), globdat  );
  
  for ( i = j = 0; i < ielemCount_; i++ )
  {
    idx_t  ielem = ielems_[i];

    if ( mask->getValue ( ielem ) )
    {
        ielems_[j++] = ielem;
    }
  }

  ielems_.reshape ( j );
  ielemCount_ = j;

  // read interface mesh and build ielems_
  
  StringVector mats;

  myProps.get ( mats, MATERIAL_PROP );
  myConf. set ( MATERIAL_PROP, mats );

  const int matCount = mats.size ();

  materials_.  resize ( matCount   );

  std::vector<int> nul; int temp, temp1; IntMatrix temp2;
  std::map<int,std::vector<int> > matMap; 

  readInterfaceMesh ( temp, temp1, temp2, matMap, dupNodes_,
                      myProps, myConf, nul );

  
  myConf.erase ( "bulk1s");
  myConf.erase ( "bulk2s");


  // shape + integration rule
  
  Ref<BShape> bShape;

  try                                // for standard interface elements
  {
    bShape = BShapeFactory::newInstance
     ( joinNames (myName_, SHAPE_PROP ), conf, props );
  }
  catch ( const jem::Exception& ex ) // for NURBS interface elements
  {
    myConf.set ( SHAPE_PROP, shProps );

    Properties  shapeProps = myProps.getProps ( SHAPE_PROP );
    String      shapeName;

    shapeProps.get( shapeName, "type" );
    
    if      ( shapeName == "Line3B" )
    {
      Ref<StdShape> sshape = newInstance<StdLine3B>();

      bShape = newInstance<ParametricEdge> (
                    "edge",
                    StdLine::getNewtonCotesScheme (3),
                    //StdLine::getGaussScheme (3),
                    sshape );
    }
    //else if ( shapeName == "Line4B" )
    //{
    //  Ref<StdShape> sshape = newInstance<StdLine4B>();

    //  bShape = newInstance<ParametricEdge> (
    //                "edge",
    //               // StdLine::getNewtonCotesScheme (4),
    //                StdLine::getGaussScheme (6),
    //                sshape );
    //}
    //else if ( shapeName == "Line5B" )
    //{
    //  Ref<StdShape> sshape = newInstance<StdLine5B>();

    //  bShape = newInstance<ParametricEdge> (
    //                "edge",
    //                StdLine::getNewtonCotesScheme (5),
    //                sshape );
    //}
    else
    {
       throw IllegalInputException (
         context, String::format ( "unsupported interface shape object"));
    }

    if ( bShape == NIL ) throw;
  }

  // finally create the interface shape object
  // FShape is short for InterfaceShape (jive/geom/InterfaceShape.h)

  shape_   = newInstance<FShape> ( bShape );

  // build StdShape corresponds to bShape
  // boundary shape = "BLine2" => InternalShape corresponds to it is "Line2"

  Properties  shapeProps = myProps.getProps ( SHAPE_PROP );
  String      type;
  shapeProps.get   ( type, "type" );
  String      nType      ( type.begin()+1, type.end() );
  shapeProps.set ( "type", nType );

  // sshape_: used to compute derivatives of shape functions of the flow elements

  sshape_ = IShapeFactory::newInstance ( joinNames ( myName_, SHAPE_PROP ), conf, props );

  // building the mapping matrix ipMap_

  nodeCount_  = shape_->nodeCount ( ) ;
  nodeCountF_ = nodeCount_ + nodeCount_/2;
  ipCount_    = shape_->ipointCount ();
  dofCount_   = nodeCount_*2;            // number of displacement dofs
  nodeCount2_ = nodeCount_/2;
  
  // write some nodal outputs to files
  // for output for every nodes, using the jive built-in Output module
  
  Properties outputProps, wellboreProps;

  write_ = false;
  
  if ( myProps.find ( outputProps ,  "output" ) )
  {
    Properties outputConf = myConf. makeProps ( "output" );

    outputProps.get  ( interval_ , "interval"          );
    outputProps.get  ( fileName_ , "file"              );

    outputConf .set  ( "interval",  interval_          );
    outputConf .set  ( "file",      fileName_          );

    write_    = true;
  }
   
  myConf .set  ( "writeOut",       write_        );
  
  
  // this is to get the flow nodes and add the dofs associated to those nodes
  // also compute the s-coordinates for the interface elements, required for H_long matrix
  // also find the upper/lower nodes for output.

  //std::cout << "gggg\n";
  IdxVector inodesF  ( nodeCountF_   );
  IdxVector inodesA  ( nodeCountF_/3 );
  IdxVector inodesB  ( nodeCountF_/3 );
  IdxVector inodesM  ( nodeCountF_/3 );
  IdxVector inodesAF ( ielemCount_ * inodesM.size() );
  IdxVector inodesBF ( ielemCount_ * inodesM.size() );
  IdxVector inodesMF ( ielemCount_ * inodesM.size() );
  idx_t kk(0), ielem;

  coords1D_.resize ( nodeCount2_, ielemCount_ );
  Matrix    coords ( rank_, nodeCount2_       );
  double    length;

  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    ielem = ielems_[ie];
    elems_.getElemNodes  ( inodesF, ielem  );

    inodesA = inodesF[slice(BEGIN,nodeCount2_)      ];   // nodes of top face
    inodesB = inodesF[slice(nodeCount2_,nodeCount_) ];   // nodes of bottom face
    inodesM = inodesF[slice(nodeCount_,END)         ];   // nodes of middle face i.e. flow element

    for ( idx_t j = 0; j < inodesM.size(); j++, kk++ )
    {
      inodesAF[kk] = inodesA[j]; 
      inodesBF[kk] = inodesB[j]; 
      inodesMF[kk] = inodesM[j]; 
    }
    
    nodes_.getSomeCoords ( coords, inodesA );
    
    length = sqrt ( (coords(0,0)-coords(0,nodeCount2_-1))*(coords(0,0)-coords(0,nodeCount2_-1))  + 
                    (coords(1,0)-coords(1,nodeCount2_-1))*(coords(1,0)-coords(1,nodeCount2_-1))  );

    if ( nodeCount2_ == 2 )
    {
      coords1D_(0,ie) = 0.;
      coords1D_(1,ie) = length;
    }
    else
    {
      coords1D_(0,ie) = 0.;
      coords1D_(1,ie) = 0.5*length;
      coords1D_(2,ie) = length;
    }
  }

  // the following works for linear and quadratic elements 

  //idx_t nnode = ielemCount_ * ( nodeCount2_ - 1 ) + 1;

  //std::cout << "fgggg\n";
  //inodesA_ = jive::util::makeUnique ( inodesAF );
  //inodesB_ = jive::util::makeUnique ( inodesBF );
  //inodesM_ = jive::util::makeUnique ( inodesMF );

  dofs_->addDofs ( jive::util::makeUnique ( inodesMF ) , pfDofTypes_);
  
  if ( write_ )
  {
  
    inodesA_.resize ( inodesAF.size () );
    inodesB_.resize ( inodesBF.size () );
    inodesM_.resize ( inodesAF.size () );

    inodesA_ = inodesAF;
    inodesB_ = inodesBF;
    inodesM_ = inodesMF;

    System::out() << inodesA_ << "\n";
    System::out() << inodesB_ << "\n";
    System::out() << inodesM_ << "\n";
  }
  
  writeWellbore_ = false;

  if ( myProps.find ( wellboreProps, "wellboreOut" ) )
  {
    String file;

    Properties wellboreConf = myConf. makeProps ( "wellboreOut" );

    wellboreProps.get  ( wellboreNodes_, "nodes" );
    wellboreProps.get  ( file,           "file"  );

    wellboreConf.set ( "nodes", wellboreNodes_ );
    wellboreConf.set ( "file",  file           );
  
    wellboreUDofA_.resize ( rank_ );
    wellboreUDofB_.resize ( rank_ );
    wellborePDofA_.resize ( 1 );
    wellborePDofB_.resize ( 1 );
    wellborePDofM_.resize ( 1 );
    
    dofs_->getDofIndices ( wellboreUDofA_,  nodes_.findNode(wellboreNodes_[0]),  aaDofTypes_ );
    dofs_->getDofIndices ( wellboreUDofB_,  nodes_.findNode(wellboreNodes_[1]),  aaDofTypes_ );
    dofs_->getDofIndices ( wellborePDofA_,  nodes_.findNode(wellboreNodes_[0]),  pwDofTypes_ );
    dofs_->getDofIndices ( wellborePDofB_,  nodes_.findNode(wellboreNodes_[1]),  pwDofTypes_ );
    dofs_->getDofIndices ( wellborePDofM_,  nodes_.findNode(wellboreNodes_[2]),  pfDofTypes_ );

    writeWellbore_ = true;

    wellboreOut_   = newInstance<PrintWriter> ( newInstance<FileWriter> ( file ));
  
    damagedElems_.resize ( ielemCount_ );
    damagedElems_ = 0;
  }
  
  //System::out () << inodesFL << "\n";
  
  System::out () << "Number of interface elements       ..... " << ielemCount_ << "\n"
                 << "Number of all nodes per element    ..... " << nodeCountF_ << "\n" 
                 << "Number of solid node per element   ..... " << nodeCount_ 
		 <<"\n";
  
  // materials

  Ref<CohesiveMaterial> mat;
  String                matName;
  std::vector<int>      ielems;

  for ( int iMat = 0; iMat < matCount; iMat++ )
  {
    matName = mats[iMat];
    mat     = newCohesiveMaterial ( matName, myConf, myProps, globdat );

    Properties matProp = myProps.getProps  ( matName );
    Properties matConf = myConf. makeProps ( matName );
    
    mat->configure ( matProp, globdat );
    mat->getConfig ( matConf, globdat );

    // get the elements of this material

    ElementGroup  egroup    = ElementGroup::get ( matConf, matProp, globdat, context );
    IdxVector     eIndices  = egroup.getIndices ();

    idx_t elemCount = eIndices.size ();
    idx_t ielem;
    //System::out () << "Number of interface elements ..... " << elemCount << "\n";
  
    //System::out () << "Number of interface elements ..... " << elemCount << "\n";

    for ( i = j = 0; i < elemCount; i++ )
    {
        ielem   = eIndices[i];

        if ( mask->getValue ( ielem ) )
        {
          elemMatMap_[ielem] = iMat;
          j++;
        }
    }

    mat->allocPoints ( j * ipCount_ );
    materials_[iMat] = mat;

  }
    //System::out () <<  elemMatMap_ << "\n";
  
  // initalize the mapping between integration points  and the material points

  initializeIPMPMap_ ();

  activeElems_.resize(ielemCount_);
  activeElems_ = 1;
    
  // compute matrix of shape functions

  getFShapeFuncs_  = getFShapeFunc  ( rank_ );
  
  globdat.get ( dtime_, "dtime" );

  thickness_ = 1.0;
  loc2Glob_  = 1.0;
  alpha_     = 1;
  w0_        = 0.;
  wf_        = 0;

  first_        = true;
  isConstraint_ = false;
  isConstant_   = false;

  if ( myProps.find ( permea_, "kl" ) )
  {
    isConstant_ = true;
    myConf. set  ( "kl", permea_ );
    myConf. set  ( "cstPermeability", isConstant_ );
  }

}
     

PoroInterfaceElementModel::~PoroInterfaceElementModel()
{
}

// -----------------------------------------------------------------
//   configure
// -----------------------------------------------------------------

void PoroInterfaceElementModel::configure

    ( const Properties&       props,
      const Properties&       globdat )
{
  Properties  myProps  = props  .findProps ( myName_ );
  
  myProps.find  ( w0_ , "w0" );
  
  if ( myProps.find ( wf_ , "wf" ) )
  {
    wf_  = 1.0 / wf_;
  }
  
  myProps.find ( thickness_, THICKNESS_PROP );
  myProps.find ( loc2Glob_, "rotation" );
  myProps.find ( alpha_ , BIOT_COEF_PROP );
  myProps.get  ( mu_ , VISCOSITY_PROP );
  myProps.get  ( kt_ , "kt" );
  myProps.find ( isConstraint_ ,  "pressureConstraint" );
  
  mu_  = 1.0 / (12.0*mu_);
}

// -----------------------------------------------------------------
//   getConfig
// -----------------------------------------------------------------


void  PoroInterfaceElementModel::getConfig

    ( const Properties&       conf )             const
{
  Properties  myConf  = conf  .makeProps ( myName_ );

  myConf. set   ( "w0", w0_  );
  myConf. set  ( "kt", kt_  );
  myConf. set  ( VISCOSITY_PROP, mu_  );
  myConf. set  ( BIOT_COEF_PROP, alpha_  );
  myConf. set  ( THICKNESS_PROP, thickness_ );
  myConf .set  ( "pressureConstraint", isConstraint_   );
}    

// -----------------------------------------------------------------
//   constraintInternalNodes_
// -----------------------------------------------------------------

void PoroInterfaceElementModel::constraintInternalNodes_

    ( const Properties&       props,
      const Properties&       globdat )
{
  System::out() << "PoroInterfaceElementModel: constraint internal nodes...\n";

  Assignable<NodeSet>   nodes    = NodeSet::find ( globdat );
  Assignable<NodeGroup> intNodes = NodeGroup::find ( "internal", nodes, globdat );

  IdxVector  inodes;

  idx_t nodeCount = intNodes.size();

  //System::out() << " there are " << nodeCount << " internal nodes\n";
  
  inodes.resize ( nodeCount );
  inodes = intNodes.getIndices ();

  idx_t in, idof;

  for ( idx_t i = 0; i < nodeCount; i++ )
  {
    in   = inodes[i];
    idof = dofs_->findDofIndex ( in, 2 );
  
    //System::out() << idof << ", " << in << " \n";

    // find nodes with p_w dof \ie not fracture fluid nodes

    if ( idof > 0 )
    {
      if ( cons_->isSlaveDof ( idof ) ) 
      {
        continue;
      }

      cons_->addConstraint ( idof, 0. );
    }
  }

  //cons_->printTo ( Printer::get() );

  System::out() << "PoroInterfaceElementModel: constraint internal nodes...done\n";
}

// -----------------------------------------------------------------
//   removeIntConstraints_
// -----------------------------------------------------------------

void PoroInterfaceElementModel::removeIntConstraints_

    ( const Properties&       globdat )
{
  System::out() << "PoroInterfaceElementModel: remove internal constraints...\n";

  Assignable<NodeSet>   nodes    = NodeSet::find ( globdat );
  Assignable<NodeGroup> intNodes = NodeGroup::find ( "internal", nodes, globdat );

  IdxVector  inodes;

  idx_t nodeCount = intNodes.size();

  //System::out() << " there are " << nodeCount << " internal nodes\n";
  
  inodes.resize ( nodeCount );
  inodes = intNodes.getIndices ();

  idx_t in, pwDof, pfDof;

  for ( idx_t i = 0; i < nodeCount; i++ )
  {
    in    = inodes[i];
    pwDof = dofs_->findDofIndex ( in, 2 );
    pfDof = dofs_->findDofIndex ( in, 3 );
  
    // find nodes with p_w dof \ie not fracture fluid nodes

    if ( pwDof > 0 ) cons_->eraseConstraint ( pwDof );
    if ( pfDof > 0 ) cons_->eraseConstraint ( pfDof );
  }

  System::out() << "PoroInterfaceElementModel:  remove internal constraints...done\n";
}

    
// -----------------------------------------------------------------
//   takeAction
// -----------------------------------------------------------------    

bool PoroInterfaceElementModel::takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat )
{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  if ( action == Actions::INIT )
  {
    if ( isConstraint_ ) initPressureConstraints_ (  );

    constraintInternalNodes_ ( params, globdat );

    return true;
  }

  if ( action == Actions::GET_EXT_VECTOR )
  {
    return true;
  }

  if ( action == Actions::COMMIT )
  {
    for ( int im = 0; im < materials_.size(); im++ )
    {
      materials_[im]->commit ();
    }

    if (write_)         writeNodalOutput_    ( globdat );
    if (writeWellbore_) writeWellboreOutput_ ( globdat );

    if ( first_ )
    {
      removeIntConstraints_ ( globdat );
      if ( isConstraint_ ) initPressureConstraints_ (  );
    }

    first_ = false;

    return true;
  }

  if ( action == "CHECK_COMMIT" )
  {
    //checkCommit_   ( globdat, params  );

    return true;
  }
  
  if ( action == Actions::GET_MATRIX0 )
  {
    Ref<MatrixBuilder>  mbuilder;

    Vector  disp, disp0;
    Vector  force;

    // Get the current displacements.

    StateVector::get    ( disp,  dofs_, globdat );
    StateVector::getOld ( disp0, dofs_, globdat );

    // Get the matrix builder and the internal force vector.

    params.get ( mbuilder, ActionParams::MATRIX0    );
    params.get ( force,    ActionParams::INT_VECTOR );

    getMatrix_ ( *mbuilder, force, disp0, disp );

    // store the internal force in globdat for use later

    globdat.set ( ActionParams::INT_VECTOR, force );

    return true;
  }

  
  if ( action == Actions::GET_INT_VECTOR )
  {
    Vector  disp, disp0;
    Vector  force;

    StateVector::get    ( disp,  dofs_, globdat );
    StateVector::getOld ( disp0, dofs_, globdat );

    params.get ( force,    ActionParams::INT_VECTOR );

    getIntForce_ ( force, disp0, disp );

    // store the internal force in globdat for use later

    globdat.set ( ActionParams::INT_VECTOR, force );

    return true;
  }

  return false;
}

//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void PoroInterfaceElementModel::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   disp0,
    const Vector&   disp )

{
  Matrix      coords     ( rank_, nodeCount_    );
  Matrix      coordsM    ( 1,     nodeCount2_   );
  Matrix      coheStiff  ( rank_,     rank_     ); // T matrix relates traction and rate of disp jump
  Matrix      Q          ( rank_,     rank_     ); // transformation matrix
  Matrix      normals    ( rank_, nodeCount2_   );
  Matrix      intScheme  ( rank_, ipCount_      );

  Cubix       xgrads     (1, nodeCount2_, ipCount_ );

  Vector      fintU       ( dofCount_  );
  Vector      fintP       ( nodeCount2_ );

  Matrix      Kuu        ( dofCount_,  dofCount_  );   // stiffness tangent matrix
  Matrix      Qup        ( dofCount_,  nodeCount2_ );   // 
  Matrix      Cpp        ( nodeCount2_, nodeCount2_ );   // Compressibility matrix
  Matrix      Hpp        ( nodeCount2_, nodeCount2_ );   // longitudinal permeability matrix
  Matrix      Htp        ( nodeCount2_, nodeCount2_ );   // transversal permeability matrix
  Matrix      CH         ( nodeCount2_, nodeCount2_ );   // C + H 
  
  Vector      elemDisp   ( dofCount_  ); // nodal displacement of upper face
  Vector      elemPresM  ( nodeCount2_); // nodal pressure of mid face
  Vector      elemPresA  ( nodeCount2_); // nodal pressure of top face
  Vector      elemPresB  ( nodeCount2_); // nodal pressure of top face

  Vector      elemDisp0   ( dofCount_  ); // nodal displacement of upper face
  Vector      elemPresM0  ( nodeCount2_); // nodal pressyre  of mid face
  Vector      elemPresA0  ( nodeCount2_); // nodal pressure of top face
  Vector      elemPresB0  ( nodeCount2_); // nodal pressure of top face
  Vector      du          ( dofCount_  ); 
  Vector      dP          ( nodeCount2_); // nodal pressyre  of mid face

  Vector      jump        ( rank_ );      // displacement jump at 1 integration point  
  Vector      jump0       ( rank_ );      // displacement jump at 1 integration point  
  Vector      djump       ( rank_ );      // displacement jump at 1 integration point  
  Vector      elemDJump   ( dofCount_ ); // nodal displacement jump 
  Vector      elemDJump0  ( dofCount_ ); // nodal displacement jump 
  Vector      traction    ( rank_ );

  IdxVector   inodesF    ( nodeCountF_ ); // inodes of entire element (6 for linear elements)
  IdxVector   inodes     ( nodeCount_  ); // inodes of standard interface element (4 for linear elements)
  IdxVector   inodesM    ( nodeCount2_ ); // inodes of mid face (or flow element)
  IdxVector   inodesA    ( nodeCount2_ ); // inodes of upper nodes
  IdxVector   inodesB    ( nodeCount2_ ); // inodes of lower nodes

  IdxVector   iaadofs    ( dofCount_   ); // displacement dofs
  IdxVector   ipfdofsM   ( nodeCount2_ ); // fracturing fluid dof (flow element)
  IdxVector   ipwdofsA   ( nodeCount2_ ); // flow fluid dof of top nodes
  IdxVector   ipwdofsB   ( nodeCount2_ ); // flow fluid dof of bottom nodes
  
  Vector      ipWeights  ( ipCount_    );
  Vector      ipWeightsM ( ipCount_    );
  Vector      Nf         ( nodeCount2_ );
  Vector      m          ( 2           );

  m[0] = 1.; m[1] = 0.;

  Matrix      N          ( rank_, nodeCount_ * 2 );
  Matrix      Nt         = N.transpose ( ); 

  MChain2     mc2;
  MChain3     mc3;

  int         ielem, iMat, cpoint;

  Ref<CohesiveMaterial> material;

  std::vector<int>      mpoints;

  double                wi, a, b;

  intScheme = shape_->getIntegrationScheme ();

  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    ielem = ielems_[ie];
    elems_.getElemNodes  ( inodesF, ielem  );

    inodes  = inodesF[slice(BEGIN,nodeCount_)];   // nodes of upper/lower faces
    inodesA = inodes[slice(BEGIN,nodeCount2_)];   // nodes of upper face
    inodesB = inodes[slice(nodeCount2_,END)  ];   // nodes of lower face
    inodesM = inodesF[slice(nodeCount_,END)  ];   // nodes of middle face i.e. flow element
                
    nodes_.getSomeCoords ( coords, inodes );
    
    /*
    length = sqrt ( (coords(0,0)-coords(0,nodeCount2_-1))*(coords(0,0)-coords(0,nodeCount2_-1))  + 
                    (coords(1,0)-coords(1,nodeCount2_-1))*(coords(1,0)-coords(1,nodeCount2_-1))  );

    // only works for linear elements!!!
    coordsM(0,0) = 0.;
    coordsM(0,1) = length;

    System::out() << inodesA << "\n";
    System::out() << inodesB << "\n";    
    System::out() << inodesM << "\n";
    System::out() << length << "\n";
    */

    coordsM(0,ALL) = coords1D_( ALL, ie );

    // get the dofs associated with nodes of each element
    // three dofs: aa(disp), pw(pore pressure) and pf(fracturing fluid pressure)

    dofs_->getDofIndices ( iaadofs,  inodes,  aaDofTypes_ );
    dofs_->getDofIndices ( ipfdofsM, inodesM, pfDofTypes_ );
    dofs_->getDofIndices ( ipwdofsA, inodesA, pwDofTypes_ );
    dofs_->getDofIndices ( ipwdofsB, inodesB, pwDofTypes_ );
    
    // get nodal values

    elemDisp      = select ( disp,  iaadofs  );  // displacements
    elemPresA     = select ( disp,  ipwdofsA );  // pressure of upper  nodes
    elemPresB     = select ( disp,  ipwdofsB );  // pressure of lower nodes
    elemPresM     = select ( disp,  ipfdofsM );  // pressure of mid nodes

    elemDisp0     = select ( disp0, iaadofs );  // same, converged values
    elemPresA0    = select ( disp0, ipwdofsA ); 
    elemPresB0    = select ( disp0, ipwdofsB ); 
    elemPresM0    = select ( disp0, ipfdofsM ); 

    // get the integration points and weights on the cohesive surface

    Matrix functions = shape_->getShapeFunctions     ();
    shape_->getIntegrationWeights ( ipWeights, coords );
    ipWeights *= thickness_;

    //System::out() << ipWeights << "\n";

    getTransformationMatrix_ ( Q, coords(ALL,slice(BEGIN,nodeCount2_)) );

    //mpoints  = materialMap_[ie];
    
    iMat     = elemMatMap_[ielem]; 
    material = materials_[iMat];
    
    //System::out() << "iMat: " << iMat << "\n";
    //System::out() << "ielem: " << ielem << "\n";

    //shape_->getNormals ( normals, ipWeights,  coords );
    //System::out() << "normals: " << normals << "\n";
    
    sshape_->getShapeGradients ( xgrads, ipWeightsM, coordsM );
    
    /*
    System::out() << "coordsM    "  << coordsM    << "\n";
    System::out() << "ipWeights  "  << ipWeights  << "\n";
    System::out() << "ipWeightsM "  << ipWeightsM << "\n";
    //System::out() << xgrads << "\n";
    */
    
    // Initialization the internal forces
    // and the tangent stiffness matrix
    
    fintU = 0.0; fintP = 0.0;

    Kuu   = 0.0; Qup   = 0.0; Cpp   = 0.0; 
    Hpp   = 0.0; Htp   = 0.0;    

    // loop over these integration points

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // compute N-matrix from shape functions

      Nf = functions(ALL,ip);
      getFShapeFuncs_ ( N, Nf);
      //System::out() << "N: " << N << "\n";

      // compute the displacement jump at integration point ip
      // using interpolation with shape functions

      jump   = matmul ( N, elemDisp  );
      jump0  = matmul ( N, elemDisp0 );
      djump  = jump - jump0;
    
      //System::out() << "jumAB: " << jump << "\n";

      // get the correct material point of this int point

      cpoint = ipMpMap_ (ielem,ip); //System::out() << "cpoint: " << cpoint << "\n\n";

      // compute traction and stiffness in cohesive material

      jump0 = matmul ( Q, jump0 ); // from global coords. to local coords.
      djump = matmul ( Q, djump );

      material->update ( traction, coheStiff, jump0, djump, cpoint );

      // transform the cohesive tangent matrix to global system

      coheStiff   = mc3.matmul ( Q, coheStiff, Q );

      // now, compute the cohesive matrix and internal force

      a          = Q(0,0)*jump[0] + Q(0,1)*jump[1];
      //System::out() << "[[u]]: " << a << "\n";
      if ( a < 0 ) a = 0.;
      b          = a*a*a*mu_;

      if (isConstant_) b = permea_;

      wi         = ipWeights[ip];

      //System::out() << "rb,wa: " << b << "\n";
      //System::out() << "b: " << b << "\n";
      //b = 0. ? b < 0. : b;

      fintU += wi * mc2.matmul ( Nt, Q, traction  );
      Kuu   += wi * mc3.matmul ( Nt, coheStiff, N );
      Qup   -= wi * alpha_  * mc3.matmul ( Nt, Q, matmul ( m, Nf ));
      Cpp   -= wi * wf_     * matmul ( Nf, Nf );
      Htp   += wi * kt_     * matmul ( Nf, Nf );
      Hpp   -= wi * b       * matmul ( xgrads(0,ALL,ip), xgrads(0,ALL,ip) );
    }
    
    /*
    System::out() << Kuu << "\n";
    System::out() << Qup << "\n";
    System::out() << Cpp << "\n";
    System::out() << Hpp << "\n";
    */
    
    // Assembly ...
    // Refer to Equation (68) in my note on the topic
    
    Matrix QupT = Qup.transpose ( );

    Htp *= dtime_;
    //System::out() << Htp << "\n";

    CH   = Cpp + dtime_ * Hpp - 2.0*Htp;

    du    = elemDisp  - elemDisp0;
    dP    = elemPresM - elemPresM0;

    fintU += matmul ( Qup,  elemPresM  );
    fintP += matmul ( QupT, du )  + matmul ( Cpp, dP ) + dtime_ * matmul ( Hpp, elemPresM ) 
          - 2. * matmul ( Htp, elemPresM ) 
          +      matmul ( Htp, elemPresA ) 
          +      matmul ( Htp, elemPresB );


    // assembly of forces due to cohesive traction

    select ( force, iaadofs  ) += fintU;
    select ( force, ipfdofsM ) += fintP;
  
    //System::out() << "fint: " << force << "\n";

    // assembly of stiffness due to cohesive traction (elemDispJump)

    mbuilder.addBlock ( iaadofs,  iaadofs,  Kuu  );
    mbuilder.addBlock ( iaadofs,  ipfdofsM, Qup  );
    mbuilder.addBlock ( ipfdofsM, iaadofs,  QupT );
    mbuilder.addBlock ( ipfdofsM, ipfdofsM, CH   );
    mbuilder.addBlock ( ipfdofsM, ipwdofsA, Htp  );
    mbuilder.addBlock ( ipfdofsM, ipwdofsB, Htp  );
    mbuilder.addBlock ( ipwdofsA, ipfdofsM, Htp  );
    mbuilder.addBlock ( ipwdofsB, ipfdofsM, Htp  );
  }
}

//-----------------------------------------------------------------------
//   getIntForce_
//-----------------------------------------------------------------------


void PoroInterfaceElementModel::getIntForce_

  ( const Vector&   force,
    const Vector&   disp0,
    const Vector&   disp )

{
}

// ---------------------------------------------------------------------
//   getTransformationMatrix_
// ---------------------------------------------------------------------
// Attention with the direction of the local coordinate system
// It must be consistent with the way the disp jump is computed
// disp jump = disp of upper face - disp of lower face
// => the normal direction points from lower face to upper face


void PoroInterfaceElementModel::getTransformationMatrix_

  ( const Matrix& Q,
    const Matrix& coord )
{
  double x0 = coord(0,1) - coord(0,0);
  double y0 = coord(1,1) - coord(1,0);

  double alpha     = ::atan2 ( y0, x0 );
  double sinAlpha  = ::sin (alpha);
  double cosAlpha  = ::cos (alpha);

  Q(0,0) = - sinAlpha; Q(0,1) = cosAlpha;
  Q(1,0) =   cosAlpha; Q(1,1) = sinAlpha;

  Q *= loc2Glob_;
  //System::out() << Q << "\n\n";
}

// ---------------------------------------------------------------------
//   initPressureConstraints_
// ---------------------------------------------------------------------

void PoroInterfaceElementModel::initPressureConstraints_ ()
{
  const int            nodeCount = dupNodes_.size ();

  int                  pwDof;
  int                  pwDofM;     // master pore pressure dof
  int                  inCount;
  int                  inter; // 0: standard node
                              // 1: interface node

  int pfType  = pfDofTypes_[0];
  int pwType  = pwDofTypes_[0];

  std::vector<int>     nodes;

  for ( int in = 0; in < nodeCount; in++ )
  {
    nodes   = dupNodes_[in];
    inCount = nodes.size ();
    inter   = nodes[inCount-1];

    if ( inCount == 2 ) continue;

    pwDofM = dofs_->getDofIndex ( nodes[0], pwType ) ; 

    for ( int jn = 1; jn < inCount-2; jn++ )
    {
      //System::out() << nodes[0] << ", " << nodes[jn] << ", " << nodes[inCount-2]<<  "\n";
	  pwDof = dofs_->getDofIndex ( nodes[jn], pwType ) ; 
      cons_->addConstraint ( pwDof, pwDofM, 1. ); 
    }
	 
    pwDof = dofs_->getDofIndex ( nodes[inCount-2], pfType ) ; 
    cons_->addConstraint ( pwDof, pwDofM, 1. ); 
  }

  //cons_->printTo ( Printer::get() );

  System::out() << "finished the constraint of bulk interface dofs...\n";
}

/*
 * Old implementation, element based
 * which might not be working for complex cases.
 
void PoroInterfaceElementModel::initPressureConstraints_ ()
{

  cout << "initializing the constraint of bulk interface dofs...\n";

  IdxVector   inodesF    ( nodeCountF_ ); // inodes of entire element (6 for linear elements)
  IdxVector   inodes     ( nodeCount_  ); // inodes of standard interface element (4 for linear elements)
  IdxVector   inodesM    ( nodeCount2_ ); // inodes of mid face (or flow element)
  IdxVector   inodesA    ( nodeCount2_ ); // inodes of upper nodes
  IdxVector   inodesB    ( nodeCount2_ ); // inodes of lower nodes

  IdxVector   ipfdofsM   ( nodeCount2_ ); // fracturing fluid dof (flow element)
  IdxVector   ipwdofsA   ( nodeCount2_ ); // flow fluid dof of top nodes
  IdxVector   ipwdofsB   ( nodeCount2_ ); // flow fluid dof of bottom nodes
  
  idx_t       ielem, idofA, idofB, idofM;

  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    ielem = ielems_[ie];
    elems_.getElemNodes  ( inodesF, ielem  );

    inodes  = inodesF[slice(BEGIN,nodeCount_)];   // nodes of upper/lower faces
    inodesA = inodes[slice(BEGIN,nodeCount2_)];   // nodes of upper face
    inodesB = inodes[slice(nodeCount2_,END)  ];   // nodes of lower face
    inodesM = inodesF[slice(nodeCount_,END)  ];   // nodes of middle face i.e. flow element
                
    // get the dofs associated with nodes of each element
    // three dofs: pw(pore pressure) and pf(fracturing fluid pressure)

    dofs_->getDofIndices ( ipfdofsM, inodesM, pfDofTypes_ );
    dofs_->getDofIndices ( ipwdofsA, inodesA, pwDofTypes_ );
    dofs_->getDofIndices ( ipwdofsB, inodesB, pwDofTypes_ );

    for ( idx_t i = 0; i < nodeCount2_; i++ )
    {
      idofA  = ipwdofsA[i];    
      idofB  = ipwdofsB[i];    
      idofM  = ipfdofsM[i];    

      if ( cons_->isSlaveDof ( idofB ) == false ) cons_->addConstraint ( idofB, idofA, 1. ); 
      if ( cons_->isSlaveDof ( idofM ) == false ) cons_->addConstraint ( idofM, idofA, 1. ); 
    }
  }
  
  //cons_->printTo ( Printer::get() );  Printer::flush ();
}
*/

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------


// called after convergence to check if there is any inactive interface
// becoming active. Procedure is:
// For all inactive interface elements ie:
//   for all integration points ip
//     - compute traction at ip, t
//     - check activation criterion for t, if yes
//     - activate ie by removing constraints
//     - if no, go to next ipoint
// Question:
// How many interface elements should be allowed to be activated
// during one time step?


void PoroInterfaceElementModel::checkCommit_

  ( const Properties& globdat,
    const Properties& params )
{
  using jive::model::ActionParams;

  cout << "Checking activation of interface elements ...\n";
}

//-----------------------------------------------------------------------
//   initializeIPMPMap_
//-----------------------------------------------------------------------

void  PoroInterfaceElementModel::initializeIPMPMap_ ( )

{
        int   ipoint;

  IntVector   matMap ( materials_.size ( ) );

  matMap = 0;

  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    int  ielem = ielems_[ie];

    // get correct material ID

    int matID = elemMatMap_[ielem]; 
    ipoint    = matMap[matID];     // System::out() << "ielem: " << ielem<< "\n\n";

    // loop over integration points 

    for ( idx_t ip = 0; ip < ipCount_; ip++, ipoint++, matMap[matID]++  )
    {
      ipMpMap_ (ielem,ip) = ipoint; //System::out() << "ipMpMap_: " << ipoint << "\n\n";
    }

  }
}

//-----------------------------------------------------------------------
//   computeForceForInjection_
//-----------------------------------------------------------------------

void  PoroInterfaceElementModel::computeForceForInjection_ ( const Vector& fext )

{
  int idof = dofs_->getDofIndex (  q0Node_, pfDofTypes_[0] );

  double force = -dtime_ * Q0_;

    
  fext[idof] += force;

  //System::out() << "fext: " << fext << "\n";
}

//-----------------------------------------------------------------------
//   writeNodalOutput_
//-----------------------------------------------------------------------
// write to files with the following format
// x-coordTop y-coordTop pressureTop x-coordBot y-coordB pressureBot x-coordM y-coordM pressureMid
// x-coordTop y-coordTop pressureTop x-coordBot y-coordB pressureBot x-coordM y-coordM pressureMid
// x-coordTop y-coordTop pressureTop x-coordBot y-coordB pressureBot x-coordM y-coordM pressureMid
// etc.

void PoroInterfaceElementModel::writeNodalOutput_

   ( const Properties&  globdat )
{
  using jive::model::StateVector;
  
  const String   context = getContext ();
        int      step;

  globdat.get ( step, Globdat::TIME_STEP );


  if ( ( step % interval_ ) != 0 ) return;
  
  System::out() << "Writing output for crack nodes...\n";

  Ref<PrintWriter> out   = newInstance<PrintWriter> (newInstance<FileWriter> 
		                    ( fileName_ + String(step)  + ".dat" ));

  // get the state vector

  Vector       state;
  StateVector::get ( state, dofs_, globdat );

  // get the node coordinates 

  const idx_t  nodeCount = inodesA_.size();

  //System::out() << nodeCount;

  Matrix    coordsA     ( rank_, nodeCount );
  Matrix    coordsB     ( rank_, nodeCount );
  Matrix    coordsM     ( rank_, nodeCount );
  
  nodes_.getSomeCoords ( coordsA, inodesA_ );
  nodes_.getSomeCoords ( coordsB, inodesB_ );
  nodes_.getSomeCoords ( coordsM, inodesM_ );

  idx_t dof1, dof2, dof3, axdof1, aydof1, axdof2, aydof2, idA, idB, idM;
    
  idx_t  pwType = pwDofTypes_[0];
  idx_t  pfType = pfDofTypes_[0];

  idx_t  axType = aaDofTypes_[0];
  idx_t  ayType = aaDofTypes_[1];

  for ( idx_t i = 0; i < nodeCount; i++ ) 
  {
    idA   = inodesA_[i];
    idB   = inodesB_[i];
    idM   = inodesM_[i];
    
    axdof1 = dofs_->findDofIndex ( idA, axType );
    aydof1 = dofs_->findDofIndex ( idA, ayType );
    
    axdof2 = dofs_->findDofIndex ( idB, axType );
    aydof2 = dofs_->findDofIndex ( idB, ayType );

    dof1   = dofs_->findDofIndex ( idA, pwType );
    dof2   = dofs_->findDofIndex ( idB, pwType );
    dof3   = dofs_->findDofIndex ( idM, pfType );

    *out << String::format( "%12.3E %12.3E %6.3E %6.3E  " , coordsA(0,i), coordsA(1,i), state[axdof1], state[aydof1] );
    *out << String::format( "%12.8E ", state[dof1] );
    
    *out << String::format( "%12.3E %12.3E %6.3E %6.3E  " , coordsB(0,i), coordsB(1,i), state[axdof2], state[aydof2] );
    *out << String::format( "%12.8E ", state[dof2] );
   
    *out << String::format( "%12.3E %12.3E  "       , coordsM(0,i), coordsM(1,i) );
    *out << String::format( "%12.8E ", state[dof3] );

     out->printLine();
     out->flush();
  }
}


//-----------------------------------------------------------------------
//   writeWellboreOutput_
//-----------------------------------------------------------------------

void PoroInterfaceElementModel::writeWellboreOutput_

   ( const Properties&  globdat )
{
  // get the state vector

  Vector       state;
  StateVector::get ( state, dofs_, globdat );

  Vector dispA(2), dispB(2);

  dispA      = select ( state,  wellboreUDofA_ );
  dispB      = select ( state,  wellboreUDofB_ );

  double pA  = state[wellborePDofA_[0]];
  double pB  = state[wellborePDofB_[0]];
  double pM  = state[wellborePDofM_[0]];

  // find the tip => crack length 

  
  IdxVector   inodesF    ( nodeCountF_ ); // inodes of entire element (6 for linear elements)
  IdxVector   inodes     ( nodeCount2_ ); // inodes of mid face (or flow element)
  Matrix      coords     ( rank_, nodeCount2_ );
  Matrix      gcoords    ( rank_, ipCount_    );
  idx_t       ielem, iMat, cpoint, ig(0);
  double      dam(0.), xmax(0.), x1, x2;

  Ref<CohesiveMaterial> material;

  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    ielem = ielems_[ie];
    elems_.getElemNodes  ( inodesF, ielem  );
    inodes  = inodesF[slice(BEGIN,nodeCount2_)];   // nodes of upper/lower faces

    nodes_.getSomeCoords ( coords, inodes );

    //System::out() << coords << "\n";
    x1       = coords(0,0);
    x2       = coords(0,nodeCount2_-1);
    x2       = max (x1, x2);
    
    iMat     = elemMatMap_[ielem]; 
    material = materials_[iMat];
   
    //shape_->getGlobalIntegrationPoints ( gcoords, coords );

    // loop over these integration points

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      cpoint = ipMpMap_ (ielem,ip); //System::out() << "cpoint: " << cpoint << "\n\n";
      //dam    = material->giveHistory ( cpoint );

      //if (dam > 0. )
      if ( material->justTractionFree ( cpoint ) )
      {
        ig++;
      }
    }
    if ( ig == ipCount_ ){
      xmax = max (x2, xmax );
    }
      
    ig = 0;
  }
    
  *wellboreOut_ << String::format( "%12.8E  ", dispA[1] - dispB[1]     );
  *wellboreOut_ << String::format( "%12.8E %12.8E %12.8E ", pA, pB, pM );
  *wellboreOut_ << String::format( "%12.8E  ", xmax );

  wellboreOut_->printLine();
  wellboreOut_->flush();
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newInterfaceModel
//-----------------------------------------------------------------------


static Ref<Model>     newPoroInterfaceElementModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<PoroInterfaceElementModel> ( name, conf, props, globdat );
}

//-----------------------------------------------------------------------
//   declareInterfaceElementModel
//-----------------------------------------------------------------------

void declarePoroInterfaceElementModel ()
{
  jive::model::ModelFactory::declare ( "PoroInterface", & newPoroInterfaceElementModel );
}



