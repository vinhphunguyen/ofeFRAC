/*
 * 
 *  Copyright (C) 2016 Monash University. All rights reserved.
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


#include "PerforationInterfaceModel.h"
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


typedef MatmulChain<double,1>      MChain1;
typedef MatmulChain<double,2>      MChain2;
typedef MatmulChain<double,3>      MChain3;
typedef MatmulChain<double,4>      MChain4;


// -----------------------------------------------------------------
//   static data
// -----------------------------------------------------------------

const char* PerforationInterfaceModel::SHAPE_PROP     = "shape";
const char* PerforationInterfaceModel::MATERIAL_PROP  = "material";
const char* PerforationInterfaceModel::MESH_PROP      = "meshFile";
const char* PerforationInterfaceModel::THICKNESS_PROP = "thickness";
const char* PerforationInterfaceModel::BIOT_COEF_PROP    = "alpha";
const char* PerforationInterfaceModel::VISCOSITY_PROP    = "mu";

// -----------------------------------------------------------------
//   constructor
// -----------------------------------------------------------------

PerforationInterfaceModel::PerforationInterfaceModel

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
  
  /*
  inodesA_ = jive::util::makeUnique ( inodesAF );
  inodesB_ = jive::util::makeUnique ( inodesBF );
  inodesM_ = jive::util::makeUnique ( inodesMF );
  */

  dofs_->addDofs ( jive::util::makeUnique ( inodesMF ) , pfDofTypes_);
  
  //System::out () << inodesA_ << "\n";
  
  System::out () << "Number of interface elements       ..... " << ielemCount_ << "\n"
                 << "Number of all nodes per element    ..... " << nodeCountF_ << "\n" 
                 << "Number of solid node per element   ..... " << nodeCount_ 
		 <<"\n";
  
  // compute matrix of shape functions

  getFShapeFuncs_  = getFShapeFunc  ( rank_ );

  // get thickness

  loc2Glob_  = 1.0;
  thickness_ = 1.0;

  myProps.find ( thickness_, THICKNESS_PROP );
  myConf. set  ( THICKNESS_PROP, thickness_ );

  alpha_ = 1;
  myProps.find ( alpha_ , BIOT_COEF_PROP );
  myConf. set  ( BIOT_COEF_PROP, alpha_  );
  
  globdat.get ( dtime_, "dtime" );
  
  myProps.get  ( mu_ , VISCOSITY_PROP );
  myConf. set  ( VISCOSITY_PROP, mu_  );
  
  Properties q0Props;

  if ( myProps.find ( q0Props, "injection" ) )
  {
    Properties q0Conf = myConf. makeProps ( "injection" );

    IdxVector q0Nodes;

    q0Props.get  ( Q0_ , "rate" );
    q0Conf. set  ( "rate", Q0_ );
    q0Props.get  ( q0Nodes, "node" );
    q0Conf. set  ( "node", q0Nodes );
    q0Dofs_.resize ( q0Nodes.size());

    for (idx_t i = 0; i < q0Nodes.size(); i++ )
    {
      idx_t qnode = nodes_.findNode( q0Nodes[i] );
      q0Dofs_[i]  = dofs_->getDofIndex (  qnode, pfDofTypes_[0] );
     }
  }
  
  // usually 1/wf=0
  wf_ = 0;

  if ( myProps.find ( wf_ , "wf" ) )
  {
    myConf. set  ( "wf", wf_  );
    wf_  = 1.0 / wf_;
  }
  
  mu_  = 1.0 / (12.0*mu_);

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
  }
  
  std::vector<int> nul; int temp, temp1; IntMatrix temp2;
  std::map<int,std::vector<int> > matMap; 

  readInterfaceMesh ( temp, temp1, temp2, matMap, dupNodes_,
                      myProps, myConf, nul );

  first_ = true;
  
  myConf.erase ( "bulk1s");
  myConf.erase ( "bulk2s");
  
  isConstraint_ = false;

  myProps.find ( isConstraint_ ,  "pressureConstraint" );
  myConf .set  ( "pressureConstraint", isConstraint_   );
}
     

PerforationInterfaceModel::~PerforationInterfaceModel()
{
}

// -----------------------------------------------------------------
//   configure
// -----------------------------------------------------------------

void PerforationInterfaceModel::configure

    ( const Properties&       props,
      const Properties&       globdat )
{
}

// -----------------------------------------------------------------
//   getConfig
// -----------------------------------------------------------------


void  PerforationInterfaceModel::getConfig

    ( const Properties&       conf )             const
{
}    
    
// -----------------------------------------------------------------
//   takeAction
// -----------------------------------------------------------------    

bool PerforationInterfaceModel::takeAction

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
    
    initDispConstraints_ (  );
    
    return true;
  }

  if ( action == Actions::GET_EXT_VECTOR )
  {
    // Get the external force vector.
    
    Vector  force;
    params.get ( force,    ActionParams::EXT_VECTOR );
    
    computeForceForInjection_ ( force );

    return true;
  }

  if ( action == Actions::COMMIT )
  {
    if (writeWellbore_) writeWellboreOutput_ ( globdat );
    
    if ( first_ )
    {
      removeDispConstraints_ ( );
    }

    first_ = false;

    return true;
  }

  if ( action == Actions::CHECK_COMMIT )
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



void PerforationInterfaceModel::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   disp0,
    const Vector&   disp )

{
  Matrix      coords     ( rank_, nodeCount_    );
  Matrix      coordsM    ( 1,     nodeCount2_   );
  Matrix      Q          ( rank_,     rank_     ); // transformation matrix
  Matrix      normals    ( rank_, nodeCount2_   );
  Matrix      intScheme  ( rank_, ipCount_      );

  Cubix       xgrads     (1, nodeCount2_, ipCount_ );

  Vector      fintU       ( dofCount_  );
  Vector      fintP       ( nodeCount2_ );

  Matrix      Qup        ( dofCount_,  nodeCount2_ );   // 
  Matrix      Cpp        ( nodeCount2_, nodeCount2_ );   // Compressibility matrix
  Matrix      Hpp        ( nodeCount2_, nodeCount2_ );   // longitudinal permeability matrix
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
  Vector      elemDJump   ( dofCount_ ); // nodal displacement jump 
  Vector      elemDJump0  ( dofCount_ ); // nodal displacement jump 

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

  int         ielem, iMat;

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

   Qup   = 0.0; Cpp   = 0.0; Hpp   = 0.0; 

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

      a          = Q(0,0)*jump[0] + Q(0,1)*jump[1];
      /*
      System::out() << "[[u]]: " << a << "\n";
      System::out() << "[[u]]: " << jump << "\n";
      System::out() << "[[u]]: " << Q << "\n";
      */
      if ( a < 0 ) a = 0.;
      b          = a*a*a*mu_;

      wi         = ipWeights[ip];

      //System::out() << "rb,wa: " << b << "\n";
      //System::out() << "b: " << b << "\n";
      //b = 0. ? b < 0. : b;

      Qup   -= wi * alpha_  * mc3.matmul ( Nt, Q, matmul ( m, Nf ));
      Cpp   -= wi * wf_     * matmul ( Nf, Nf );
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

    //System::out() << Htp << "\n";

    CH   = Cpp + dtime_ * Hpp;

    du    = elemDisp  - elemDisp0;
    dP    = elemPresM - elemPresM0;

    fintU += matmul ( Qup,  elemPresM  );
    fintP += matmul ( QupT, du )  + matmul ( Cpp, dP ) + dtime_ * matmul ( Hpp, elemPresM ) ;


    // assembly of forces 

    select ( force, iaadofs  ) += fintU;
    select ( force, ipfdofsM ) += fintP;
  
    //System::out() << "fint: " << force << "\n";

    // assembly of stiffness due to cohesive traction (elemDispJump)

    mbuilder.addBlock ( iaadofs,  ipfdofsM, Qup  );
    mbuilder.addBlock ( ipfdofsM, iaadofs,  QupT );
    mbuilder.addBlock ( ipfdofsM, ipfdofsM, CH   );
  }
}

//-----------------------------------------------------------------------
//   getIntForce_
//-----------------------------------------------------------------------


void PerforationInterfaceModel::getIntForce_

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


void PerforationInterfaceModel::getTransformationMatrix_

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

//-----------------------------------------------------------------------
//   computeForceForInjection_
//-----------------------------------------------------------------------

void  PerforationInterfaceModel::computeForceForInjection_ ( const Vector& fext )

{
  double force(0.);
  
  if ( first_ == false ) force = -dtime_ * Q0_;

  select(fext, q0Dofs_) += force;

  //first_ = false;

  //System::out() << "fext: " << fext << "\n";
}

// ---------------------------------------------------------------------
//   initPressureConstraints_
// ---------------------------------------------------------------------

void PerforationInterfaceModel::initPressureConstraints_ ()
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
	   pwDof = dofs_->getDofIndex ( nodes[jn], pwType ) ; 
       cons_->addConstraint ( pwDof, pwDofM, 1. ); 
    }
	 
    pwDof = dofs_->getDofIndex ( nodes[inCount-2], pfType ) ; 
    cons_->addConstraint ( pwDof, pwDofM, 1. ); 
  }
}

//-----------------------------------------------------------------------
//   initDispConstraints_
//-----------------------------------------------------------------------

void PerforationInterfaceModel::initDispConstraints_ ()
{

  IdxVector   inodesF    ( nodeCountF_ ); // inodes of entire element (6 for linear elements)
  IdxVector   inodes     ( nodeCount_  ); // inodes of standard interface element (4 for linear elements)
  IdxVector   inodesA    ( nodeCount2_ ); // inodes of upper nodes
  IdxVector   inodesB    ( nodeCount2_ ); // inodes of lower nodes

  idx_t       ielem, inodeA, inodeB, idofAx, idofAy, idofBx, idofBy;
  
  idx_t       xtype  = aaDofTypes_[0];
  idx_t       ytype  = aaDofTypes_[1];

  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    ielem = ielems_[ie];
    elems_.getElemNodes  ( inodesF, ielem  );

    inodes  = inodesF[slice(BEGIN,nodeCount_)];   // nodes of upper/lower faces
    inodesA = inodes[slice(BEGIN,nodeCount2_)];   // nodes of upper face
    inodesB = inodes[slice(nodeCount2_,END)  ];   // nodes of lower face

    for ( idx_t i = 0; i < nodeCount2_; i++ )
    {
      inodeA = inodesA[i];
      inodeB = inodesB[i];

      //System::out() << inodeA << " " << inodeB << "\n";
    
      idofAx = dofs_->getDofIndex ( inodeA, xtype );
      idofAy = dofs_->getDofIndex ( inodeA, ytype );

      idofBx = dofs_->getDofIndex ( inodeB, xtype );
      idofBy = dofs_->getDofIndex ( inodeB, ytype );

      cons_->addConstraint ( idofAx, idofBx, 1. ); 
      cons_->addConstraint ( idofAy, idofBy, 1. ); 
    }
  }
  
  //cons_->printTo ( Printer::get() );  Printer::flush ();
}

//-----------------------------------------------------------------------
//   removeDispConstraints_
//-----------------------------------------------------------------------

void PerforationInterfaceModel::removeDispConstraints_ ()
{

  IdxVector   inodesF    ( nodeCountF_ ); // inodes of entire element (6 for linear elements)
  IdxVector   inodes     ( nodeCount_  ); // inodes of standard interface element (4 for linear elements)
  IdxVector   inodesA    ( nodeCount2_ ); // inodes of upper nodes
  IdxVector   inodesB    ( nodeCount2_ ); // inodes of lower nodes

  idx_t       ielem, inodeA, inodeB, idofAx, idofAy, idofBx, idofBy;
  
  idx_t       xtype  = aaDofTypes_[0];
  idx_t       ytype  = aaDofTypes_[1];


  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    ielem = ielems_[ie];
    elems_.getElemNodes  ( inodesF, ielem  );

    inodes  = inodesF[slice(BEGIN,nodeCount_)];   // nodes of upper/lower faces
    inodesA = inodes[slice(BEGIN,nodeCount2_)];   // nodes of upper face
    inodesB = inodes[slice(nodeCount2_,END)  ];   // nodes of lower face

    //System::out() << inodes << "\n";
    //System::out() << nodeCount2_ << "\n";
    //System::out() << inodesA << "\n";
    //System::out() << inodesB << "\n";

    for ( idx_t i = 0; i < nodeCount2_; i++ )
    {
      inodeA = inodesA[i];
      inodeB = inodesB[i];
    
      idofAx = dofs_->getDofIndex ( inodeA, xtype );
      idofAy = dofs_->getDofIndex ( inodeA, ytype );

      idofBx = dofs_->getDofIndex ( inodeB, xtype );
      idofBy = dofs_->getDofIndex ( inodeB, ytype );

      cons_->eraseConstraint ( idofAx ); 
      cons_->eraseConstraint ( idofAy ); 
      
      cons_->eraseConstraint ( idofBx ); 
      cons_->eraseConstraint ( idofBy ); 
    }
  }
  
  //cons_->printTo ( Printer::get() );  Printer::flush ();
}

//-----------------------------------------------------------------------
//   writeWellboreOutput_
//-----------------------------------------------------------------------

void PerforationInterfaceModel::writeWellboreOutput_

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

    // loop over these integration points

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
    }
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


static Ref<Model>     newPerforationInterfaceModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<PerforationInterfaceModel> ( name, conf, props, globdat );
}

//-----------------------------------------------------------------------
//   declareInterfaceElementModel
//-----------------------------------------------------------------------

void declarePerforationInterfaceModel ()
{
  jive::model::ModelFactory::declare ( "PerforationInterface", & newPerforationInterfaceModel );
}



