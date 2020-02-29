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

#include <jem/base/Float.h>
#include <jem/base/System.h>
#include <jem/base/Array.h>
#include <jem/util/Properties.h>
#include <jem/util/StringUtils.h>
#include <jem/util/Event.h>
#include <jem/io/FileWriter.h>
#include <jem/io/PrintWriter.h>

#include <jive/util/Database.h>
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


#include "RigidPoroInterfaceElementModel.h"
#include "util/utilities.h"
#include "material/XCohesiveMat.h"
#include "material/TuronXCohesiveMat.h"


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
using jive::util::Database;
using jive::util::IntDBColumn;
using jive::mp::ItemMask;



typedef MatmulChain<double,1>      MChain1;
typedef MatmulChain<double,2>      MChain2;
typedef MatmulChain<double,3>      MChain3;
typedef MatmulChain<double,4>      MChain4;


// -----------------------------------------------------------------
//   static data
// -----------------------------------------------------------------

const char* RGPoroInterfaceElementModel::SHAPE_PROP     = "shape";
const char* RGPoroInterfaceElementModel::MATERIAL_PROP  = "material";
const char* RGPoroInterfaceElementModel::MESH_PROP      = "meshFile";
const char* RGPoroInterfaceElementModel::THICKNESS_PROP = "thickness";
const char* RGPoroInterfaceElementModel::DISP_NODES_PROP= "dNodeGroup";
const char* RGPoroInterfaceElementModel::DIRICHLET_NODES_PROP= "dirNodeGroup";
const char* RGPoroInterfaceElementModel::BIOT_COEF_PROP    = "alpha";
const char* RGPoroInterfaceElementModel::VISCOSITY_PROP    = "mu";
const char* RGPoroInterfaceElementModel::DTIME_PROP        = "dtime";

// -----------------------------------------------------------------
//   constructor
// -----------------------------------------------------------------

RGPoroInterfaceElementModel::RGPoroInterfaceElementModel

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
  
  allElems_   = ElemSet::find    ( globdat                   );
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
  
  // this is to get the flow nodes and add the dofs associated to those nodes
  // also find the upper/lower nodes for output.

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

  idx_t nnode = ielemCount_ * ( nodeCount2_ - 1 ) + 1;

  inodesA_.resize ( nnode );
  inodesB_.resize ( nnode );
  inodesM_.resize ( nnode );
  
  inodesA_ = jive::util::makeUnique ( inodesAF );
  inodesB_ = jive::util::makeUnique ( inodesBF );
  inodesM_ = jive::util::makeUnique ( inodesMF );

  std::cout << "fffffff\n";

  dofs_->addDofs ( jive::util::makeUnique ( inodesMF ) , pfDofTypes_);
  
  //System::out () << inodesA_ << "\n";
  
  System::out () << "Number of interface elements       ..... " << ielemCount_ << "\n"
                 << "Number of all nodes per element    ..... " << nodeCountF_ << "\n" 
                 << "Number of solid node per element   ..... " << nodeCount_ 
		 <<"\n";
  
  // initialize ipStatus_ (store beta: 0 for elastic, 1 for damaged)

  ipStatus_.resize ( ielemCount_, ipCount_);
  ipStatus_ = 0;
  
  // materials

  Ref<XCohesiveMat>     mat;
  String                matName;
  std::vector<int>      ielems;

  for ( int iMat = 0; iMat < matCount; iMat++ )
  {
    matName = mats[iMat];
    mat     = newXCohesiveMat ( matName, myConf, myProps, globdat );

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

  // compute matrix of shape functions

  getFShapeFuncs_  = getFShapeFunc  ( rank_ );
  getShapeGrads_   = getShapeGradsFunc ( rank_ );

  // get thickness

  thickness_ = 1.0;

  myProps.find ( thickness_, THICKNESS_PROP );
  myConf. set  ( THICKNESS_PROP, thickness_ );

  tractionFree_  = false;
  initConsDone1_ = false;
  initConsDone2_ = false;

  String dispName;

  if ( myProps.find ( dispName , DISP_NODES_PROP ) )
  {
    myConf. set  ( DISP_NODES_PROP, dispName  );

    Assignable<NodeGroup> group = 
      
	      NodeGroup::find ( dispName, nodes_, globdat );

    if ( group == NIL )
    {
      System::err() << "InterfaceElementModel : node group == NIL!\n";
    }

    dispNodes_ . resize ( group.size () );
    dispNodes_ = group.getIndices ();
  }

  if ( myProps.find ( dispName , DIRICHLET_NODES_PROP ) )
  {
    myConf. set  ( DIRICHLET_NODES_PROP, dispName  );

    Assignable<NodeGroup> group = 
      
	      NodeGroup::find ( dispName, nodes_, globdat );

    if ( group == NIL )
    {
      System::err() << "InterfaceElementModel : node group == NIL!\n";
    }

    diriNodes_ . resize ( group.size () );
    diriNodes_ = group.getIndices ();
  }

  loc2Glob_ = 1.0;
  myProps.find ( loc2Glob_, "rotation" );

  alpha_ = 1;
  myProps.find ( alpha_ , BIOT_COEF_PROP );
  myConf. set  ( BIOT_COEF_PROP, alpha_  );
  
  globdat.get ( dtime_, "dtime" );
  
  myProps.get  ( mu_ , VISCOSITY_PROP );
  myConf. set  ( VISCOSITY_PROP, mu_  );
  
  w0_ = 0.;
  myProps.find  ( w0_ , "w0" );
  myConf. set   ( "w0", w0_  );
  
  Properties q0Props;

  if ( myProps.find ( q0Props, "injection" ) )
  {
    Properties q0Conf = myConf. makeProps ( "injection" );
    q0Props.get  ( Q0_ , "rate" );
    q0Conf. set  ( "rate", Q0_ );
    q0Props.get  ( q0Node_, "node" );
    q0Conf. set  ( "node", q0Node_ );

    q0Node_ = nodes_.findNode( q0Node_ );
  }
  
  // usually 1/wf=0
  wf_ = 0;

  if ( myProps.find ( wf_ , "wf" ) )
  {
    myConf. set  ( "wf", wf_  );
    wf_  = 1.0 / wf_;
  }
  
  myProps.get  ( kt_ , "kt" );
  myConf. set  ( "kt", kt_  );

  mu_  = 1.0 / (12.0*mu_);

  isConstant_ = false;
  myConf. set  ( "cstPermeability", isConstant_ );

  if ( myProps.find ( permea_, "kl" ) )
  {
    isConstant_ = true;
    myConf. set  ( "kl", permea_ );
    myConf. set  ( "cstPermeability", isConstant_ );
  }

  // write some nodal outputs to files
  // for output for every nodes, using the jive built-in Output module

  write_ = false;
  if ( myProps.find ( fileName_ ,  "outputFile" ) )
  {
    myProps.get  ( interval_ , "outputInterval"    );
    myConf .set  ( "outputInterval", interval_     );
    myConf .set  ( "outputFile",     fileName_     );

    write_    = true;
  }
   myConf .set  ( "writeOut",       write_        );
 
  isPrsConstraint_ = false;
  isDspConstraint_ = true;

  myProps.find ( isPrsConstraint_ ,  "pressureConstraint" );
  myConf .set  ( "pressureConstraint", isPrsConstraint_   );
  
  myProps.find ( isDspConstraint_ ,  "dispConstraint" );
  myConf .set  ( "dispConstraint", isDspConstraint_   );

  Properties wellboreProps;

  if ( myProps.find ( wellboreProps, "wellboreOut" ) )
  {
    String file;

    Properties wellboreConf = myConf. makeProps ( "wellboreOut" );
    wellboreProps.get  ( wellboreNodes_, "nodes" );
    wellboreProps.get  ( file,           "file"  );
  
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
  
   myProps.find ( tractionFree_, "tractionFree" );
  
  // read element database for the bulk neibouring elements

  IntDBColumn*   dbcol = NULL;
  Ref<Database>  dbase = Database::get       ( "neighbours", allElems_.getData(), globdat, context );

  dbcol  = dbase->getIntColumn ( "connect" );
  ibulkElems_.resize(ielemCount_,2);

  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    idx_t  ielem = ielems_[ie];

    Tuple<idx_t,2>  ibulk;
    idx_t           jelem, kelem;

    dbcol->getValues ( ibulk, ielem );

    jelem             = allElems_.findElement ( ibulk[0] );
    kelem             = allElems_.findElement ( ibulk[1] );
    ibulkElems_(ie,0) = jelem;
    ibulkElems_(ie,1) = kelem;

    if ( ( jelem < 0 ) || ( kelem < 0 ) || ( jelem == kelem ) ){
      System::err() << ielem << "\n";
      System::err() << ibulk << "\n";
    }
  }
  
  // shape of volumetric elements

  ishape_ = IShapeFactory::newInstance (joinNames(myName_, "ishape"), conf, props);
  inodeCount_ = ishape_->nodeCount ();
  idofCount_  = inodeCount_ * 2;

  
  // get materials of volumetric elements stored in globdat

  //globdat.get( belemMatMap_, "elemMatMap" );

  damagedElems_.resize(ielemCount_);
  damagedElems_ = 1;
  
  checkInterval_ = 1;

  myProps.find ( checkInterval_, "checkInterval" );
  myConf. set  ( "checkInterval", checkInterval_ );

  resolve_ = true;

  myProps.find ( resolve_, "resolve" );
  myConf. set  ( "resolve", resolve_ );
}
     

RGPoroInterfaceElementModel::~RGPoroInterfaceElementModel()
{
}

// -----------------------------------------------------------------
//   configure
// -----------------------------------------------------------------

void RGPoroInterfaceElementModel::configure

    ( const Properties&       props,
      const Properties&       globdat )
{
}

// -----------------------------------------------------------------
//   getConfig
// -----------------------------------------------------------------


void  RGPoroInterfaceElementModel::getConfig

    ( const Properties&       conf )             const
{
}    
    
// -----------------------------------------------------------------
//   takeAction
// -----------------------------------------------------------------    

bool RGPoroInterfaceElementModel::takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat )
{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  if ( action == Actions::INIT )
  {
     if ( isPrsConstraint_ ) initPressureConstraints_ (  );
     
     if ( isDspConstraint_) initDisplacementConstraints_ (  );

     zeroPressure_ ( );

    return true;
  }

  if ( action == Actions::GET_EXT_VECTOR )
  {
    //updateBcConstraints_      ();
    //constraintDirichletNodes_ ();  // do once only
    
    // Get the external force vector.
    Vector  force;
    params.get ( force,    ActionParams::EXT_VECTOR );
    computeForceForInjection_(force);

    return true;
  }

  if ( action == Actions::COMMIT )
  {

    for ( int im = 0; im < materials_.size(); im++ )
    {
      materials_[im]->commit ();
    }

    System::out() << myName_ << " is ,w commit\n";
    //System::out() << "Materials commited\n";

    if (write_)         writeNodalOutput_    ( globdat );
    if (writeWellbore_) writeWellboreOutput_ ( globdat );
    
    return true;
  }

  //if ( action == Actions::CHECK_COMMIT )
  if ( action == "CHECKCOMMIT" )
  {
    Vector  disp;
    StateVector::get ( disp, dofs_, globdat );

    System::out() << myName_ << " is checking commit\n";
    
    checkCommit_   ( globdat, params, disp  );

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



void RGPoroInterfaceElementModel::getMatrix_

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

  Ref<XCohesiveMat>     material;

  std::vector<int>      mpoints;

  double                wi, a, b, length;

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

    //System::out() << "jumAB: " << elemDispJump << "\n";

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

      a = 0.; traction=0.; coheStiff=0.;
      double alpha=0.;

      if ( ipStatus_(ie,ip) )
      {
        material->update ( traction, coheStiff, jump0, djump, cpoint );
        coheStiff   = mc3.matmul ( Q, coheStiff, Q ); // transform the cohesive tangent matrix to global system

        a          = Q(0,0)*jump[0] + Q(0,1)*jump[1];
        alpha =alpha_;

        /*
        if ( jem::Float::isNaN( sum(  traction  ) ) )
        {
        System::out() << "traction: " << traction << "\n";
        System::out() << "tangent: " << coheStiff << "\n";
        }
        */
      if ( a < 0 ) a = 0.;
      b          = a*a*a*mu_;

      if (isConstant_) b = permea_;

      wi         = ipWeights[ip];

      //System::out() << "rb,wa: " << b << "\n";
      //System::out() << "b: " << b << "\n";
      //b = 0. ? b < 0. : b;

      fintU += wi * mc2.matmul ( Nt, Q, traction  );
      Kuu   += wi * mc3.matmul ( Nt, coheStiff, N );
      Qup   -= wi * alpha  * mc3.matmul ( Nt, Q, matmul ( m, Nf ));
      Cpp   -= wi * wf_     * matmul ( Nf, Nf );
      Htp   += wi * kt_     * matmul ( Nf, Nf );
      Hpp   -= wi * b       * matmul ( xgrads(0,ALL,ip), xgrads(0,ALL,ip) );
      }

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


void RGPoroInterfaceElementModel::getIntForce_

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


void RGPoroInterfaceElementModel::getTransformationMatrix_

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
//   initDisplacementConstraints_
// ---------------------------------------------------------------------


void RGPoroInterfaceElementModel::initDisplacementConstraints_ ()
{

  cout << "initializing the displacement constraint...\n";

  IdxVector   inodesF    ( nodeCountF_ ); // inodes of entire element (6 for linear elements)
  IdxVector   inodes     ( nodeCount_  ); // inodes of standard interface element (4 for linear elements)
  IdxVector   inodesA    ( nodeCount2_ ); // inodes of upper nodes
  IdxVector   inodesB    ( nodeCount2_ ); // inodes of lower nodes

  IdxVector   iaadofsA   ( nodeCount_ ); // flow fluid dof of top nodes
  IdxVector   iaadofsB   ( nodeCount_ ); // flow fluid dof of bottom nodes
  
  idx_t       ielem, idofA, idofB;

  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    ielem = ielems_[ie];
    elems_.getElemNodes  ( inodesF, ielem  );

    inodes  = inodesF[slice(BEGIN,nodeCount_)];   // nodes of upper/lower faces
    inodesA = inodes[slice(BEGIN,nodeCount2_)];   // nodes of upper face
    inodesB = inodes[slice(nodeCount2_,END)  ];   // nodes of lower face
                
    // get the dofs associated with nodes of each element

    dofs_->getDofIndices ( iaadofsA, inodesA, aaDofTypes_ );
    dofs_->getDofIndices ( iaadofsB, inodesB, aaDofTypes_ );

    for ( idx_t i = 0; i < nodeCount_; i++ )
    {
      idofA  = iaadofsA[i];    
      idofB  = iaadofsB[i];    

      if ( cons_->isSlaveDof ( idofB ) == false ) cons_->addConstraint ( idofB, idofA, 1. ); 
    }
  }
  
  cout << "initializing the displacement dofs...done\n";
  //cons_->printTo ( Printer::get() );  Printer::flush ();
}


// ---------------------------------------------------------------------
//   removeDisplacementConstraints
// ---------------------------------------------------------------------


void RGPoroInterfaceElementModel::removeDisplacementConstraints_ ( int ie )
{

  IdxVector   inodesF    ( nodeCountF_ ); // inodes of entire element (6 for linear elements)
  IdxVector   inodes     ( nodeCount_  ); // inodes of standard interface element (4 for linear elements)
  IdxVector   inodesA    ( nodeCount2_ ); // inodes of upper nodes
  IdxVector   inodesB    ( nodeCount2_ ); // inodes of lower nodes

  IdxVector   iaadofsA   ( nodeCount_  ); // displacement dof of top nodes
  IdxVector   iaadofsB   ( nodeCount_  ); // displacement dof of bottom nodes
  IdxVector   ipwdofsA   ( nodeCount2_ ); // flow fluid dof of top nodes
  IdxVector   ipwdofsB   ( nodeCount2_ ); // flow fluid dof of bottom nodes
  
  idx_t       ielem, idofA, idofB;

  // Iterate over all elements assigned to this model.

    ielem = ielems_[ie];
    elems_.getElemNodes  ( inodesF, ielem  );

    inodes  = inodesF[slice(BEGIN,nodeCount_)];   // nodes of upper/lower faces
    inodesA = inodes[slice(BEGIN,nodeCount2_)];   // nodes of upper face
    inodesB = inodes[slice(nodeCount2_,END)  ];   // nodes of lower face
                
    // get the dofs associated with nodes of each element

    dofs_->getDofIndices ( iaadofsA, inodesA, aaDofTypes_ );
    dofs_->getDofIndices ( iaadofsB, inodesB, aaDofTypes_ );
    dofs_->getDofIndices ( ipwdofsA, inodesA, pwDofTypes_ );
    dofs_->getDofIndices ( ipwdofsB, inodesB, pwDofTypes_ );

    for ( idx_t i = 0; i < nodeCount_; i++ )
    {
      idofA  = iaadofsA[i];    
      idofB  = iaadofsB[i];    

      if ( cons_->isSlaveDof ( idofB ) == true ) cons_->eraseConstraint ( idofB ); 
    }
      
    if ( cons_->isMasterDof ( ipwdofsA[0] ) ) cons_->unsetRvalue ( ipwdofsA[0] ); 
    if ( cons_->isMasterDof ( ipwdofsB[0] ) ) cons_->unsetRvalue ( ipwdofsB[0] ); 
  
  //cons_->printTo ( Printer::get() );  Printer::flush ();
}

void RGPoroInterfaceElementModel::initPressureConstraints_ ()
{
  const int            nodeCount = dupNodes_.size ();

  int                  pwDof, ydof;
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
	   pwDof = dofs_->getDofIndex ( nodes[jn], pwType ) ; 
       cons_->addConstraint ( pwDof, pwDofM, 1. ); 
    }
	 
    pwDof = dofs_->getDofIndex ( nodes[inCount-2], pfType ) ; 
    cons_->addConstraint ( pwDof, pwDofM, 1. ); 
  }
  
  cout << "initializing the pressure dofs...done\n";
}


// ---------------------------------------------------------------------
//   initPressureConstraints_
// ---------------------------------------------------------------------

/*
void RGPoroInterfaceElementModel::initPressureConstraints_ ()
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

// ---------------------------------------------------------------------
//   constraintDirichletNodes_
// ---------------------------------------------------------------------

void RGPoroInterfaceElementModel::constraintDirichletNodes_ ()
{
}


//-----------------------------------------------------------------------
//   updateBcConstraint_
//-----------------------------------------------------------------------

// NOTE: this has to be called every time steps!!!
// There are constraints between duplicated nodes and Dirichlet nodes
// When the prescribed displacement of Dirichlet nodes change, it is 
// necessary to update their slave dofs so that they're still tied
// together. This is the case for displacement control only.

void RGPoroInterfaceElementModel::updateBcConstraints_ ()
{
}


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


void RGPoroInterfaceElementModel::checkCommit_

  ( const Properties& globdat,
    const Properties& params,
    const Vector&     disp )
{
  using jive::model::ActionParams;

  int      step;

  globdat.get ( step, Globdat::TIME_STEP );

  Matrix     eTractions ( rank_, ipCount_              );     
  Cubix      tractions  ( ielemCount_, ipCount_, rank_ );     
  Vector     failure    ( ielemCount_ );

  tractions = 0.;
  failure   = 0.;

  int ig(0), critElem, ielem, iMat, cpoint;
  double maxT(0.), maxTe;

  // do for every checkInterval_ step

  if ( ( step % checkInterval_ ) == 0 )
  {
    for ( int ie = 0; ie < ielemCount_; ie++ )
    {
      // skip elements already cracked at all Gauss points

      if  ( testany ( ipStatus_(ie,ALL) ) ) continue;

      maxTe = evalElementFailure_ ( ie, ig, disp, globdat, eTractions );

      tractions ( ie, ALL, ALL ) = eTractions;

      if ( maxTe > maxT )
      {
        maxT     = maxTe;
        critElem = ie;
      }
    }
  }

  if ( maxT > 1. )       
  {
     ielem             = ielems_[critElem];
     iMat              = elemMatMap_[ielem];
     Ref<XCohesiveMat> material = materials_[iMat];

     for ( int ip = 0; ip < ipCount_; ip++ )
     {
       cpoint = ipMpMap_ (ielem,ip); //System::out() << "cpoint: " << cpoint << "\n\n";
       dynamicCast<TuronXCohesiveMat>(material)->initShift ( cpoint, tractions(critElem, ip, ALL) );
       ipStatus_(critElem,ip) = 1;
       ig++;
     }
     //System::out() << tractions <<"\n";
     removeDisplacementConstraints_ ( critElem );
  
     if ( resolve_ )
     {
       params.set ( "accept", false );
     }
  }

}

//-----------------------------------------------------------------------
//   initializeIPMPMap_
//-----------------------------------------------------------------------

void  RGPoroInterfaceElementModel::initializeIPMPMap_ ( )

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

void  RGPoroInterfaceElementModel::computeForceForInjection_ ( const Vector& fext )

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

void RGPoroInterfaceElementModel::writeNodalOutput_

   ( const Properties&  globdat )
{
  using jive::model::StateVector;
  const String   context = getContext ();
        int      step;

  globdat.get ( step, Globdat::TIME_STEP );

  if ( ( step % interval_ ) != 0 ) return;

  Ref<PrintWriter> out   = newInstance<PrintWriter> (newInstance<FileWriter> 
		                    ( fileName_ + String(step)  + ".dat" ));

  // get the state vector

  Vector       state;
  StateVector::get ( state, dofs_, globdat );

  // get the node coordinates 

  int  nodeCount = inodesA_.size();
  Matrix    coordsA     ( rank_, nodeCount );
  Matrix    coordsB     ( rank_, nodeCount );
  Matrix    coordsM     ( rank_, nodeCount );
  
  nodes_.getSomeCoords ( coordsA, inodesA_ );
  nodes_.getSomeCoords ( coordsB, inodesB_ );
  nodes_.getSomeCoords ( coordsM, inodesM_ );

  int dof1, dof2, dof3, idA, idB, idM;

  for ( int i = 0; i < nodeCount; i++ ) 
  {
    idA   = inodesA_[i];
    idB   = inodesB_[i];
    idM   = inodesM_[i];

    dof1 = dofs_->findDofIndex ( idA, pwDofTypes_[0] );
    dof2 = dofs_->findDofIndex ( idB, pwDofTypes_[0] );
    dof3 = dofs_->findDofIndex ( idM, pfDofTypes_[0] );

    *out << String::format( "%12.3E %12.3E  "       , coordsA(0,i), coordsA(1,i) );
    *out << String::format( "%12.8E ", state[dof1] );
    
    *out << String::format( "%12.3E %12.3E  "       , coordsB(0,i), coordsB(1,i) );
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

void RGPoroInterfaceElementModel::writeWellboreOutput_

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
  idx_t       ielem, iMat, cpoint;
  double      dam(0.), xmax(0.), x1, x2;

  Ref<XCohesiveMat> material;

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

    // loop over these integration points

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      cpoint = ipMpMap_ (ielem,ip); //System::out() << "cpoint: " << cpoint << "\n\n";
      dam    = material->giveHistory ( cpoint );

      if (dam > 0. )
      {
        xmax              = max (x2, xmax );
        //System::out() << x1 << " " << x2 << " " << xmax << "\n";
        continue; // only one gauss point damaged => element damaged
      }
    }
  }
    
  *wellboreOut_ << String::format( "%12.8E  ", dispA[1] - dispB[1]     );
  *wellboreOut_ << String::format( "%12.8E %12.8E %12.8E ", pA, pB, pM );
  *wellboreOut_ << String::format( "%12.8E  ", xmax );

  wellboreOut_->printLine();
  wellboreOut_->flush();
}

//--------------------------------------------------------------------
//     evalElementFailure
//--------------------------------------------------------------------

double  RGPoroInterfaceElementModel::evalElementFailure_

( int               ie,
  int&              ig,
  const Vector&     disp,
  const Properties& globdat,
  const Matrix&     tractions )
{
   IdxVector          inodesF    ( nodeCountF_ ); // inodes of entire element (6 for linear elements)
   IdxVector          inodes     ( nodeCount_  );
   IdxVector          inodesB1   ( inodeCount_ );
   IdxVector          inodesB2   ( inodeCount_ );
   IdxVector          idofsB1    ( idofCount_ );
   IdxVector          idofsB2    ( idofCount_ );
   Vector             elemDisp1  ( idofCount_ );
   Vector             elemDisp2  ( idofCount_ );
   Vector             stress     ( 4 );
   Vector             traction   ( 2 );
   Vector             xi1        ( 2 );
   Vector             xi2        ( 2 );
   Vector             normal     ( 2 );
   Vector             ipWeights (ipCount_);

   Matrix             coordsB1   ( rank_, inodeCount_   );
   Matrix             coordsB2   ( rank_, inodeCount_   );
   Matrix             coords     ( rank_,  nodeCount_   );
   Matrix             gcoords    ( rank_, ipCount_      );
   Matrix             normals    ( rank_, ipCount_      );
   Matrix             Q          ( rank_, rank_         );

   Matrix             grads1     ( rank_, inodeCount_ );
   Matrix             grads2     ( rank_, inodeCount_ );
   Matrix             B1         ( 4, idofCount_ );
   Matrix             B2         ( 4, idofCount_ );
   Matrix             C1(4,4), C2(4,4);

   double             ipFailure;
   idx_t              cpoint;

   MChain2            mc2;

   IdxVector   inodesA    ( nodeCount2_  ); // inodes of upper face
   IdxVector   inodesB    ( nodeCount2_  ); // inodes of lower face
   IdxVector   idofsA     ( nodeCount_   ); // dofs of upper face
   IdxVector   idofsB     ( nodeCount_   ); // dofs of lower face

   int                count(0);

   //System::out() << "bulk1 " << "\n";
   int  ielem = ielems_[ie];
   allElems_.getElemNodes  ( inodesF, ielem  );
   inodes  = inodesF[slice(BEGIN,nodeCount_)];   // nodes of upper/lower faces
   nodes_.getSomeCoords ( coords, inodes );

   //System::out() << "bulk1 " << "\n";
   // to get the averaged stress


   idx_t bulk1         = ibulkElems_(ie,0);
   idx_t bulk2         = ibulkElems_(ie,1);

   //System::out() << "bulk1 " << bulk1 << "\n";
   //System::out() << "bulk2 " << bulk2 << "\n";

   allElems_.getElemNodes ( inodesB1, bulk1 );
   allElems_.getElemNodes ( inodesB2, bulk2 );

   dofs_->getDofIndices ( idofsB1, inodesB1, aaDofTypes_ );
   dofs_->getDofIndices ( idofsB2, inodesB2, aaDofTypes_ );

   elemDisp1     = select ( disp, idofsB1 );
   elemDisp2     = select ( disp, idofsB2 );

   nodes_.getSomeCoords ( coordsB1, inodesB1 );
   nodes_.getSomeCoords ( coordsB2, inodesB2 );

   // get the global coords of the interface GPs
   // used to inverse map back to the bulk parent domains

   shape_->getGlobalIntegrationPoints ( gcoords, coords );

   idx_t  iMat = 0;
   Ref<Material> bmat;
   globdat.get( bmat, String::format("mat%d",iMat) );
   bmat->update (stress,C1,stress,0);
   C2 = C1; // @todo: assume 1 material only

   //shape_->getIntegrationWeights ( ipWeights, coords );
   shape_->getNormals ( normals, ipWeights,  coords );

   iMat     = elemMatMap_[ielem];
   Ref<XCohesiveMat> material = materials_[iMat];

   double gamma = 0.5;

   Matrix functions = shape_->getShapeFunctions     ();

   double maxT(0.);

   // loop over these integration points

   for ( int ip = 0; ip < ipCount_; ip++ )
   {
     // already checked where this is called
     // if ( ipStatus_(ie,ip) == 1 ) continue; // already damaged

      ishape_->getLocalPoint ( xi1, gcoords(ALL,ip), 1e-6, coordsB1 );
      ishape_->getLocalPoint ( xi2, gcoords(ALL,ip), 1e-6, coordsB2 );

      ishape_->evalShapeGradients ( grads1, xi1, coordsB1 );
      ishape_->evalShapeGradients ( grads2, xi2, coordsB2 );

      // from bulk shape functions and derivatives, compute N and B matrices

      getShapeGrads_ ( B1, grads1  );
      getShapeGrads_ ( B2, grads2  );

      //System::out() << "sigma1: " <<  gamma*mc2.matmul(C1,B1,elemDisp1) << "\n";
      //System::out() << "sigma2: " <<  gamma*mc2.matmul(C2,B2,elemDisp2) << "\n";

      stress  =     gamma*mc2.matmul(C1,B1,elemDisp1);
      stress += (1-gamma)*mc2.matmul(C2,B2,elemDisp2);
      normal       = -normals(ALL,ip);

      double c2 = normal[0]*normal[0];
      double s2 = normal[1]*normal[1];
      double cs = normal[0]*normal[1];

      traction[0] =   c2 * stress[0] + s2 * stress[1] +    2 * cs   * stress[3];
      traction[1] = - cs * stress[0] + cs * stress[1] + ( c2 - s2 ) * stress[3];

      tractions(ip,ALL) = traction;

      cpoint = ipMpMap_ (ielem,ip); //System::out() << "cpoint: " << cpoint << "\n\n";
      ipFailure = material->evalFailure (traction,cpoint);
      
      if (ipFailure > 1.) 
      {
        count++;
        maxT = max ( maxT, ipFailure );
      }
   }

   return maxT;
}

void RGPoroInterfaceElementModel::zeroPressure_ ()
{
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

    for ( idx_t i = 0; i < nodeCount2_-1; i++ )
    {
      idofA  = ipwdofsA[i];    
      idofB  = ipwdofsB[i];    
      idofM  = ipfdofsM[i];    

      if ( cons_->isMasterDof ( idofB ) ) cons_->addConstraint ( idofB, 0. ); 
      if ( cons_->isMasterDof ( idofA ) ) cons_->addConstraint ( idofA, 0. ); 
    }
  }
  
  //cons_->printTo ( Printer::get() );  Printer::flush ();
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newInterfaceModel
//-----------------------------------------------------------------------


static Ref<Model>     newRGPoroInterfaceElementModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<RGPoroInterfaceElementModel> ( name, conf, props, globdat );
}

//-----------------------------------------------------------------------
//   declareInterfaceElementModel
//-----------------------------------------------------------------------

void declareRGPoroInterfaceElementModel ()
{
  jive::model::ModelFactory::declare ( "RGPoroInterface", & newRGPoroInterfaceElementModel );
}



