/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implements a model of two dimensional interface elements. 
 *  Both linear 4 node and quadratic 6 node interface element are supported.
 *  Both Gauss and Newton-Cotes integration schemes are implemented.
 *  
 *
 *  The overall procedure of using interface elements is 
 *
 *   (1) Build a standard FE mesh
 *   (2) Using a small program to separate this mesh and insert interface
 *       elements in between interelement boundaries. The modified mesh
 *       is saved to a file named *-solid.mesh and the interface elements
 *       are written to a file named *-interface.mesh.
 *   (3) The problem is then modelled using two models:
 *
 *      (a) Continuum model using the *-solid.mesh
 *      (b) Interface model (this class) using the *-interface.mesh
 *  
 *
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 5 July 2012
 *
 */


#include <iostream>
#include <stdio.h>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <jem/base/System.h>
#include <jem/base/Array.h>
#include <jem/util/Properties.h>
#include <jem/util/StringUtils.h>
#include <jem/util/Event.h>

#include <jive/util/XDofSpace.h>
#include <jive/util/Printer.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/fem/Globdat.h>
#include <jive/fem/NodeGroup.h>
#include <jive/geom/InternalShape.h>
#include <jive/geom/BoundaryShape.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/geom/BShapeFactory.h>
#include <jive/geom/InterfaceShape.h>


#include "DiscreteInterfaceModel.h"
#include "util/utilities.h"


using std::cout;

using namespace jem;
using jem::newInstance;
using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jive::model::StateVector;
using jive::util::Printer;
using jive::IdxVector;


typedef MatmulChain<double,1>      MChain1;
typedef MatmulChain<double,2>      MChain2;
typedef MatmulChain<double,3>      MChain3;
typedef MatmulChain<double,4>      MChain4;


// -----------------------------------------------------------------
//   static data
// -----------------------------------------------------------------

const char* DiscreteInterfaceElementModel::SHAPE_PROP     = "shape";
const char* DiscreteInterfaceElementModel::MATERIAL_PROP  = "material";
const char* DiscreteInterfaceElementModel::MESH_PROP      = "meshFile";
const char* DiscreteInterfaceElementModel::THICKNESS_PROP = "thickness";
const char* DiscreteInterfaceElementModel::DISP_NODES_PROP= "dNodeGroup";
const char* DiscreteInterfaceElementModel::DIRICHLET_NODES_PROP= "dirNodeGroup";

// -----------------------------------------------------------------
//   constructor
// -----------------------------------------------------------------

DiscreteInterfaceElementModel::DiscreteInterfaceElementModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat ) : Super(name)
{
  using jem:: util::Event;
  using jive::util::joinNames;
  using jive::geom::BShapeFactory;
  using jive::geom::BShape;
  using jive::StringVector;
  using jive::fem::NodeGroup;

  Properties    myProps = props.getProps  ( myName_ );
  Properties    myConf  = conf .makeProps ( myName_ );
 
  // build nodes_ and dofs_ and cons_ 
  
  nodes_   = NodeSet::find    ( globdat                   );
  dofs_    = XDofSpace::get   ( nodes_.getData(), globdat );
  cons_    = Constraints::get ( dofs_, globdat            );

  // read interface mesh and build ielems_
  
  std::map<int,std::vector<int> > matMap; 
  
  rank_      = nodes_.rank   ( );
  
  dofTypes_.resize ( rank_ );
  dofTypes_[0] = 0;
  dofTypes_[1] = 1;

  StringVector mats;

  myProps.get ( mats, MATERIAL_PROP );
  myConf. set ( MATERIAL_PROP, mats );

  const int matCount = mats.size ();

  materials_.  resize ( matCount   );

  std::vector<int> nul;

  readDiscreteInterfaceMesh ( elemCount_, nodeCount_, ielems_, matMap, dupNodes_,
                              myProps, myConf, nul );

  ielems_ -= 1;
  System::out() << ielems_ << "\n\n";

  nodeCount_ = 2;

  System::out () << "Number of interface elements ..... " << elemCount_ << "\n"
                 << "Number of node per element   ..... " << nodeCount_ 
		 <<"\n";
  // materials

  Ref<CohesiveMaterial> mat;
  String                matName;
  String                groupName;
  std::vector<int>      ielems;
  int                   ieCount;

  materialMap_.resize ( elemCount_ );

  for ( int im = 0; im < matCount; im++ )
  {
    matName = mats[im];
    mat     = newCohesiveMaterial ( matName, myConf, myProps, globdat );

    Properties matProp = myProps.getProps  ( matName );
    Properties matConf = myConf. makeProps ( matName );
    
    mat->configure ( matProp, globdat );
    mat->getConfig ( matConf, globdat );

    ielems  = matMap[im];
    ieCount = ielems.size ( );

    //int  id = 0;

    for ( int ie = 0; ie < ieCount; ie++ )
    {
      materialMap_[ielems[ie]].push_back ( im );
    }

    mat->allocPoints ( ieCount );

    materials_[im] = mat;
  }

  activeElems_.resize(elemCount_);
  activeElems_ = 1;
    
  // compute matrix of shape functions
  //getShapeFuncs_  = getShapeFunc  ( rank_ );

  // get thickness

  thickness_ = 1.0;

  myProps.find ( thickness_, THICKNESS_PROP );
  myConf. set  ( THICKNESS_PROP, thickness_ );

  initConsDone_ = false;

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

}
     

DiscreteInterfaceElementModel::~DiscreteInterfaceElementModel()
{
}

// -----------------------------------------------------------------
//   configure
// -----------------------------------------------------------------

void DiscreteInterfaceElementModel::configure

    ( const Properties&       props,
      const Properties&       globdat )
{
}

// -----------------------------------------------------------------
//   getConfig
// -----------------------------------------------------------------


void  DiscreteInterfaceElementModel::getConfig

    ( const Properties&       conf )             const
{
}    
    
// -----------------------------------------------------------------
//   takeAction
// -----------------------------------------------------------------    

bool DiscreteInterfaceElementModel::takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat )
{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  if ( action == Actions::INIT )
  {
    initConstraints_ (  );

    return true;
  }

  if ( action == Actions::GET_EXT_VECTOR )
  {
    updateBcConstraints_      ();
    constraintDirichletNodes_ ();  // do once only

    return true;
  }

  if ( action == Actions::COMMIT )
  {
    for ( int im = 0; im < materials_.size(); im++ )
    {
      materials_[im]->commit ();
    }

    return true;
  }

  if ( action == "CHECK_COMMIT" )
  {
    checkCommit_   ( globdat, params  );

    return true;
  }
  
  if ( action == Actions::GET_MATRIX0 )
  {
    Ref<MatrixBuilder>  mbuilder;

    Vector  disp;
    Vector  force;

    // Get the current displacements.

    StateVector::get ( disp, dofs_, globdat );

    // Get the matrix builder and the internal force vector.

    params.get ( mbuilder, ActionParams::MATRIX0    );
    params.get ( force,    ActionParams::INT_VECTOR );

    getMatrix_ ( *mbuilder, force, disp );

    // store the internal force in globdat for use later

    globdat.set ( ActionParams::INT_VECTOR, force );

    return true;
  }

  return false;
}

//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void DiscreteInterfaceElementModel::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   disp )

{
  Matrix      coords     ( rank_, nodeCount_    );
  
  Matrix      coheStiff  ( 4, 4 ); // T matrix relates traction and rate of disp jump
  Matrix      elemMat    ( 4, 4 );
  Vector      elemForce  ( 4 );
  Matrix      Q          ( rank_,   rank_     ); // transformation matrix
  Matrix      R          ( 4, 4 );
  Matrix      Rt         = R.transpose ( ); 
  
  Vector      disp1  ( 2 ); // nodal displacement of node 1
  Vector      disp2  ( 2 ); // nodal displacement of node 2
  Vector      disp1Prime(2), disp2Prime(2); // local displacement of node 1 and 2.
  Vector      jump       ( 4 );          // u1' and u2' in one vector
  Vector      xx         ( 4 );          
  Vector      stiffness  ( rank_ );

  IdxVector   inodes     ( 2 ); // inodes of both faces
  IdxVector   idofs1     ( 2 ); // dofs of first node of the spring
  IdxVector   idofs2     ( 2 ); // dofs of second node of the spring
  IdxVector   idofs      ( 4 ); // dofs of two nodes of the spring
  
  int         node1, node2;                  // node1 and node2 of the spring

  MChain2     mc2;
  MChain3     mc3;

  double      kn, kt; // normal and tangential stiffness 
  double      f1x, f1y;

  Ref<CohesiveMaterial> material;

  std::vector<int>      mpoints;

  R = 0.;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < elemCount_; ie++ )
  {
    // if this element is inactive, then discard it 

    if ( activeElems_[ie] == 0 ) continue;

    // Get the element connectivity 
    // nodal coordinates and nodal DOFs of element ie.

    inodes  = ielems_(ie,ALL);
    node1   = inodes[0];
    node2   = inodes[1];
                
    nodes_.getSomeCoords ( coords, inodes );
    
    System::out() << inodes << "\n";
    System::out() << coords << "\n";
   
    // get the dofs associated with nodes of each element

    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );
    dofs_->getDofIndices ( idofs1, node1,  dofTypes_ );
    dofs_->getDofIndices ( idofs2, node2,  dofTypes_ );
    
    // get nodal displacements of two elements 
    // need both at the same time to compute the displacement jump

    disp1     = select ( disp, idofs1 );
    disp2     = select ( disp, idofs2 ); 

    // transform to local coordinate system

    getTransformationMatrix_ ( Q, coords );

    R(slice(0,2),slice(0,2)) = Q;
    R(slice(2,4),slice(2,4)) = Q;
    //R(0,0) = Q(0,0); R(0,1) = Q(0,1);
    //R(1,0) = Q(1,0); R(1,1) = Q(1,1);
    //R(2,2) = Q(0,0); R(2,3) = Q(0,1);
    //R(3,2) = Q(1,0); R(3,3) = Q(1,1);
    
    disp1Prime = matmul ( Q, disp1 );
    disp2Prime = matmul ( Q, disp2 );

    jump[slice(0,2)] = disp1Prime;
    jump[slice(2,4)] = disp2Prime;

    System::out() << "dispA: " << jump << "\n";

    // get the integration points and weights on the cohesive surface

    mpoints  = materialMap_[ie];
    material = materials_[mpoints[0]];

    // get the correct material point of this int point

    material->update ( stiffness, coheStiff, jump, ie );

    kn = stiffness[0];
    kt = stiffness[1];

    f1x =  kn * ( disp1Prime[0] - disp2Prime[0] );
    f1y =  kt * ( disp1Prime[1] - disp2Prime[1] );

    xx[0] =  f1x;
    xx[1] =  f1y;
    xx[2] = -f1x;
    xx[3] = -f1y;

    //System::out() << "xx: " << kn << " " << f1x << " " << "\n\n";
    printf (" xx %5.12f f1x %5.12f \n", kn, f1x );
    //System::out() << "stiff: " << coheStiff << "\n\n";

    // compute internal force and stiffness matrix
    
    elemMat   = mc3.matmul ( Rt, coheStiff, R );
    elemForce = matmul     ( Rt, xx );

    // assembly of forces due to cohesive traction
    // assembly of stiffness due to cohesive traction (elemDispJump)

    select ( force, idofs ) += elemForce;
    mbuilder.addBlock ( idofs, idofs, elemMat );
  }
}

// ---------------------------------------------------------------------
//   getTransformationMatrix_
// ---------------------------------------------------------------------

void DiscreteInterfaceElementModel::getTransformationMatrix_

  ( const Matrix& Q,
    const Matrix& coord )
{
  double x0 = coord(0,1) - coord(0,0);
  double y0 = coord(1,1) - coord(1,0);

  double alpha     = ::atan2 ( y0, x0 );

  alpha += 1.5707963267948966;

  double sinAlpha  = ::sin (alpha);
  double cosAlpha  = ::cos (alpha);

  Q(1,0) = - sinAlpha; Q(1,1) = cosAlpha;
  Q(0,0) =   cosAlpha; Q(0,1) = sinAlpha;

  //Q *= -1.;
  //System::out() << Q << "\n\n";
}

// ---------------------------------------------------------------------
//   initConstraints_
// ---------------------------------------------------------------------

// The aim of this function is to constraint the bulk interface elements.
// Doing this contraint rather than using high dummy stiffness has two advantages
// namely (1) no ill conditioned matrix, (2) reduced number of dofs.
// This constraint will be removed when the traction exceeds the tensile strength.
// Used only when there are more than 1 materials.
// Attention:
//  pay attention to node with imposed displacements. 
//  this function is always executed after Dirichlet BCs is applied.

void DiscreteInterfaceElementModel::initConstraints_ ()
{
  cons_->printTo ( Printer::get() );  Printer::flush ();

  cout << "initializing the constraint of bulk interface dofs...\n";

  if ( materials_.size () == 1 ) return;

  // loop over elements, do only for bulk elements
  // which are elements with material = 0

  const int            nodeCount = dupNodes_.size ();

  int                  xdof, ydof;
  int                  xdofM, ydofM;
  int                  inCount;
  int                  inter; // 0: standard node
                              // 1: interface node

  int xtype  = dofTypes_[0];
  int ytype  = dofTypes_[1];

  std::vector<int>     nodes;

  for ( int in = 0; in < nodeCount; in++ )
  {
    nodes   = dupNodes_[in];
    inCount = nodes.size ();
    inter   = nodes[inCount-1];

    if ( inCount == 2 ) continue;

    if ( inter == 0 )  // standard bulk node
    {
      xdofM = dofs_->getDofIndex ( nodes[0], xtype ) ; 
      ydofM = dofs_->getDofIndex ( nodes[0], ytype ) ; 

      for ( int jn = 1; jn < inCount-1; jn++ )
      {
	xdof = dofs_->getDofIndex ( nodes[jn], xtype ) ; 
	ydof = dofs_->getDofIndex ( nodes[jn], ytype ) ; 
	
	addConstraint_ ( xdof, ydof, xdofM, ydofM );
      }
    }
    else  // interfacial node
    {
      xdofM = dofs_->getDofIndex ( nodes[1], xtype ) ; 
      ydofM = dofs_->getDofIndex ( nodes[1], ytype ) ; 

      for ( int jn = 2; jn < inCount-1; jn++ )
      {
	xdof = dofs_->getDofIndex ( nodes[jn], xtype ) ; 
	ydof = dofs_->getDofIndex ( nodes[jn], ytype ) ; 
	
	addConstraint_ ( xdof, ydof, xdofM, ydofM );
      }
    }
  }

  cons_->compress ();
 
  // This is for debugging only.

  System::out() << "after \n";
  cons_->printTo ( Printer::get() );  Printer::flush ();
}

// ---------------------------------------------------------------------
//   constraintDirichletNodes_
// ---------------------------------------------------------------------

void DiscreteInterfaceElementModel::constraintDirichletNodes_ ()
{
  if ( diriNodes_.size () == 0 ) return;

  if ( initConsDone_ ) return;

  cout << "Constraint Dirichlet nodes...\n";

  const int    diriNodesCount = diriNodes_.size ();

        int    index;
	int    xtype = dofTypes_[0];
	int    ytype = dofTypes_[1];
	int    xdof, ydof, xdofM, ydofM;
        int    inCount;

  std::vector<int> nodes;

  for ( int in = 0; in < diriNodesCount; in++ )
  {
    index   = diriNodes_[in]; 
    nodes   = dupNodes_[index];
    inCount = nodes.size ();

    if ( inCount == 2 ) continue;

    xdofM = dofs_->getDofIndex ( nodes[0], xtype ) ; 
    ydofM = dofs_->getDofIndex ( nodes[0], ytype ) ; 

    xdof  = dofs_->getDofIndex ( nodes[1], xtype ) ; 
    ydof  = dofs_->getDofIndex ( nodes[1], ytype ) ; 
      
    addConstraint_ ( xdof, ydof, xdofM, ydofM );
  }

  initConsDone_ = true;
}


//-----------------------------------------------------------------------
//   updateBcConstraint_
//-----------------------------------------------------------------------

// There are constraints between duplicated nodes and Dirichlet nodes
// When the prescribed displacement of Dirichlet nodes change, it is 
// necessary to update their slave dofs so that they're still tied
// together. This is the case for displacement control only.

void DiscreteInterfaceElementModel::updateBcConstraints_ ()
{
  if ( dispNodes_.size () == 0 ) return;

  cout << "updating constraints ...\n";

  const int    dispNodesCount = dispNodes_.size ();

        int    index;
	int    xtype = dofTypes_[0];
	int    ytype = dofTypes_[1];
	int    xdof, ydof, xdofM, ydofM;
        int    inCount, inter;

  std::vector<int> nodes;

  for ( int in = 0; in < dispNodesCount; in++ )
  {
    index   = dispNodes_[in];

    nodes   = dupNodes_[index];
    inCount = nodes.size ();
    inter   = nodes[inCount-1];

    if ( inCount == 2 ) continue;

  //  if ( inter == 0 )  // standard bulk node
    //{
      xdofM = dofs_->getDofIndex ( nodes[0], xtype ) ; 
      ydofM = dofs_->getDofIndex ( nodes[0], ytype ) ; 

      for ( int jn = 1; jn < inCount-1; jn++ )
      {
	      xdof = dofs_->getDofIndex ( nodes[jn], xtype ) ; 
	      ydof = dofs_->getDofIndex ( nodes[jn], ytype ) ; 
	
	      addConstraint_ ( xdof, ydof, xdofM, ydofM );
      }
   // }
  /*  else  // interfacial node
    {
      xdofM = dofs_->getDofIndex ( nodes[1], xtype ) ; 
      ydofM = dofs_->getDofIndex ( nodes[1], ytype ) ; 

      for ( int jn = 2; jn < inCount-1; jn++ )
      {
	xdof = dofs_->getDofIndex ( nodes[jn], xtype ) ; 
	ydof = dofs_->getDofIndex ( nodes[jn], ytype ) ; 
	
	addConstraint_ ( xdof, ydof, xdofM, ydofM );
      }
    }*/

  }
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


void DiscreteInterfaceElementModel::checkCommit_

  ( const Properties& globdat,
    const Properties& params )
{
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newInterfaceModel
//-----------------------------------------------------------------------


static Ref<Model>     newDiscreteInterfaceElementModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<DiscreteInterfaceElementModel> ( name, conf, props, globdat );
}

//-----------------------------------------------------------------------
//   declareInterfaceElementModel
//-----------------------------------------------------------------------

void declareDiscreteInterfaceElementModel ()
{
  jive::model::ModelFactory::declare ( "DiscreteInterface", & newDiscreteInterfaceElementModel );
}



