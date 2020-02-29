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
 *  Date: 5 July 2009
 *
 */


#include <iostream>
#include <map>
#include <iterator>


#include <jem/base/System.h>
#include <jem/base/Array.h>
#include <jem/util/StringUtils.h>
#include <jem/util/Event.h>
#include <jem/numeric/algebra/utilities.h>


#include <jive/util/XDofSpace.h>
#include <jive/util/Printer.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/fem/Globdat.h>
#include <jive/fem/NodeGroup.h>
#include <jive/geom/StdShape.h>
#include <jive/geom/BoundaryShape.h>
#include <jive/geom/BShapeFactory.h>
#include <jive/geom/StdShapeFactory.h>
#include <jive/geom/InterfaceShape.h>
#include <jive/geom/Geometries.h>
#include <jive/geom/error.h>
#include <jive/geom/BoundaryTriangle.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>


#include "InterfaceElement3DModel.h"
#include "util/utilities.h"

extern "C"
{
  #include <math.h>
}

using std::cout;

using namespace jem;
using jem::newInstance;
using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jive::model::StateVector;
using jive::util::Printer;
using jive::Cubix;

typedef MatmulChain<double,1>      MChain1;
typedef MatmulChain<double,2>      MChain2;
typedef MatmulChain<double,3>      MChain3;
typedef MatmulChain<double,4>      MChain4;


// -----------------------------------------------------------------
//   static data
// -----------------------------------------------------------------

const char* InterfaceElement3DModel::SHAPE_PROP     = "shape";
const char* InterfaceElement3DModel::MATERIAL_PROP  = "material";
const char* InterfaceElement3DModel::DISP_NODES_PROP= "dNodeGroup";
const char* InterfaceElement3DModel::DIRICHLET_NODES_PROP= "dirNodeGroup";
const char* InterfaceElement3DModel::LARGE_DISP_PROP= "largeDisp";

// -----------------------------------------------------------------
//   constructor
// -----------------------------------------------------------------

InterfaceElement3DModel::InterfaceElement3DModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat ) : Super(name)
{
  using jem::util::Event;
  using jive::util::joinNames;
  using jive::geom::StdShapeFactory;
  using jive::geom::BShapeFactory;
  using jive::geom::BShape;
  using jive::geom::Geometries;
  using jive::StringVector;
  using jive::fem::NodeGroup;

  Properties    myProps = props.getProps   ( myName_    );
  Properties    myConf  = conf .makeProps  ( myName_    );
  Properties    shProps = myProps.getProps ( SHAPE_PROP );

  // build nodes_ and dofs_ and cons_ 
  
  nodes_   = NodeSet::find    ( globdat                   );
  dofs_    = XDofSpace::get   ( nodes_.getData(), globdat );
  cons_    = Constraints::get ( dofs_, globdat            );

  // read interface mesh and build ielems_
  
  std::map<int,std::vector<int> > matMap; 
  
  dofTypes_.resize ( 3 );

  dofTypes_[0] = 0;
  dofTypes_[1] = 1;
  dofTypes_[2] = 2;

  StringVector mats;

  myProps.get ( mats, MATERIAL_PROP );
  myConf. set ( MATERIAL_PROP, mats );

  const idx_t matCount = mats.size ();

  materials_.  resize ( matCount   );

  readInterfaceMesh ( elemCount_, nodeCount_, ielems_, matMap, dupNodes_,
                      myProps, myConf, oppositeVertices_ );

  std::for_each(oppositeVertices_.begin(), oppositeVertices_.end(), [](int& d) { d-=1;});

  dofCount_   = 3 * nodeCount_;
  dofCount2_  = dofCount_  / 2;
  nodeCount2_ = nodeCount_ / 2;

  // shape + integration rule
  
  Ref<BShape> bShape;

  // Since Newton-Cotes for triangle is not implemented in
  // jem-jive, it is supported here by the user

  try
  {
    bShape = BShapeFactory::newInstance
     ( joinNames (myName_, SHAPE_PROP ), conf, props );
  }
  catch ( const jem::Exception& ex )
  {
    myConf.set ( SHAPE_PROP, shProps );

    bShape = makeBTriangleNC ( shProps );

    if ( bShape == NIL ) throw;
  }

  shape_   = newInstance<FShape> ( bShape );

  // build the StdShape corresponds to the bShape

  const String shapeGeom = bShape->getGeometry ();
        String type;

  if      ( nodeCount2_ == 3 || nodeCount2_ == 4 ) // triangle 3 or quad 4 
  {
    type = "Linear";
  }
  else if ( nodeCount2_ == 6 || nodeCount2_ == 8) // triangle 6 or quad 8 
  {
    type = "Quadratic";
  }
  else
  {
    type = "BiQuadratic";
  }

  sshape_ = StdShapeFactory::newInstance ( type, shapeGeom, conf, props );

  // building the mapping matrix ipMap_

  ipCount_ = shape_->ipointCount ();
  
  // materials

  Ref<CohesiveMaterial> mat;
  String                matName;
  String                groupName;
  std::vector<int>      ielems;
  int                   ieCount;
  
  if (matCount != matMap.size() )
  {
    throw IllegalInputException ( "Interface3DModel: wrong material number", "wrong material");
  }

  materialMap_.resize ( elemCount_ );


  for ( int im = 0; im < matCount; im++ )
  {
    matName = mats[im];
    mat     = newCohesiveMaterial ( matName, myConf, myProps, globdat );

    Properties matProp = myProps.getProps  ( matName );
    Properties matConf = myConf. makeProps ( matName );
    
    mat->configure ( matProp, globdat );
    mat->getConfig ( matConf, globdat );
    
    int index;  
    matProp.get ( index, "index" );
    matConf.set ( "index", index );

    ielems  = matMap[index];
    ieCount = ielems.size ( );

    int  id = 0;

    for ( int ie = 0; ie < ieCount; ie++ )
    {
      materialMap_[ielems[ie]].push_back ( im );

      for ( int ip = 0; ip < ipCount_; ip++, id++ )
      {
        materialMap_[ielems[ie]].push_back ( id );
      }
    }

    mat->allocPoints ( ieCount * ipCount_ );

    materials_[im] = mat;
  }
    
  activeElems_.resize(elemCount_);
  activeElems_ = 1;

  // compute matrix of shape functions

  getShapeFuncs_  = getShapeFunc  ( 3 );
  getUpdateCoordFuncs_  = getUpdateCoordFunc  ( 3 );

  initConsDone_ = false;

  String dispName;

  if ( myProps.find ( dispName , DISP_NODES_PROP ) )
  {
    myConf. set  ( DISP_NODES_PROP, dispName  );

    Assignable<NodeGroup> group = 
      
	      NodeGroup::find ( dispName, nodes_, globdat );

    if ( group == NIL)
    {
      System::err() << "InterfaceElement3DModel : node group == NIL!\n";
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
      System::err() << "InterfaceElement3DModel : node group == NIL!\n";
    }

    diriNodes_ . resize ( group.size () );
    diriNodes_ = group.getIndices ();
  }

  largeDisp_  = false;

  myProps.find ( largeDisp_, LARGE_DISP_PROP );
  myConf. set  ( LARGE_DISP_PROP, largeDisp_ );

  System::out () << "Number of interface elements ..... " << elemCount_ << "\n"
                 << "Number of node per element   ..... " << nodeCount_ << "\n"
                 << "Number of intergration point ..... " << ipCount_   <<"\n";
}
     

InterfaceElement3DModel::~InterfaceElement3DModel()
{
}

// -----------------------------------------------------------------
//   configure
// -----------------------------------------------------------------

void InterfaceElement3DModel::configure

    ( const Properties&       props,
      const Properties&       globdat )
{
}

// -----------------------------------------------------------------
//   getConfig
// -----------------------------------------------------------------


void  InterfaceElement3DModel::getConfig

    ( const Properties&       conf )             const
{
}    
    
// -----------------------------------------------------------------
//   takeAction
// -----------------------------------------------------------------    

bool InterfaceElement3DModel::takeAction

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
  
  if ( action == Actions::GET_INT_VECTOR )
  {

    Vector  disp;
    Vector  force;

    // Get the current displacements.

    StateVector::get ( disp, dofs_, globdat );

    // Get the matrix builder and the internal force vector.

    params.get ( force,    ActionParams::INT_VECTOR );

    getIntForce_ ( force, disp );

    // store the internal force in globdat for use later

    globdat.set ( ActionParams::INT_VECTOR, force );

    return true;
  }

  return false;
}

//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void InterfaceElement3DModel::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   disp )

{
  Matrix	coords     ( 3, nodeCount_ ); // coords in the undeformed configuration
  Matrix	uCoords    ( 3, nodeCount_ ); // coords in the deformed configuration
  Matrix	hCoords    ( 3, nodeCount2_ ); // coords in the undeformed config of 1 surface
  Matrix	uCoordsA   ( 3, nodeCount2_ ); // coords in the deformed config of surface +
  Matrix	uCoordsB   ( 3, nodeCount2_ ); // coords in the deformed config of surface -
  
  Matrix	coheStiff  ( 3, 3 ); // T matrix relates traction and rate of disp jump
  Matrix	Q          ( 3, 3 ); // transformation matrix
  Matrix	Qt         = Q.transpose ( ); 

  Matrix	elemMat    ( nodeCount2_ * 3, nodeCount2_ * 3 );
  Vector	elemForce  ( nodeCount2_ * 3 );
  
  Vector	elemDispA   ( dofCount2_ ); // nodal displacement of upper face
  Vector	elemDispB   ( dofCount2_ ); // nodal displacement of lower face
  Vector	elemDispJump( dofCount2_ ); // nodal displacement jump 

  Vector	traction   ( 3 );          // traction at 1 integration point 
  Vector        jump       ( 3 );          // displacement jump at 1 integration point  

  IdxVector	inodes     ( nodeCount_  ); // inodes of both faces
  IdxVector	inodesA    ( nodeCount2_ ); // inodes of upper face
  IdxVector	inodesB    ( nodeCount2_ ); // inodes of lower face
  IdxVector	idofsA     ( dofCount2_  ); // dofs of upper face
  IdxVector	idofsB     ( dofCount2_  ); // dofs of lower face
  
  Vector	ipWeights  ( ipCount_    );
  Vector	ff         ( nodeCount2_ );

  Matrix	N          ( 3, nodeCount2_ * 3 );
  Matrix	Nt         = N.transpose ( ); 

  Matrix        grads      ( 2, nodeCount2_ ); // gradients of shape functions
                                               // used to compute the transformation matrix Q
  MChain2	mc2; MChain1 mc1;
  MChain3	mc3;

  int		cpoint;

  Ref<CohesiveMaterial> material;

  std::vector<int>      mpoints;

  double                wi;

  Vector        center(2); center[0] = center[1] = 0.;
  Vector        null(2);

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < elemCount_; ie++ )
  {
    // if this element is inactive, then discard it 

    if ( activeElems_[ie] == 0 ) continue;

    // Get the element connectivity 
    // nodal coordinates and nodal DOFs of element ie.

    inodes  = ielems_(ie,ALL); 
    inodesA = inodes[slice(BEGIN,nodeCount2_)]; 
    inodesB = inodes[slice(nodeCount2_,END)  ];
    
    // get coords in undeformed config

    nodes_.getSomeCoords ( coords, inodes );

    hCoords = coords(ALL,slice(nodeCount2_,END));
    
    // get the dofs associated with nodes of each element

    dofs_->getDofIndices ( idofsA, inodesA, dofTypes_ );
    dofs_->getDofIndices ( idofsB, inodesB, dofTypes_ );
    
    // get nodal displacements of two elements 
    // need both at the same time to compute the displacement jump

    elemDispA     = select ( disp, idofsA );
    elemDispB     = select ( disp, idofsB ); 

    elemDispJump  = elemDispA - elemDispB; 

    // if large displacement is needed,

    if ( largeDisp_ )
    {
      getUpdateCoordFuncs_   ( uCoordsA, hCoords, elemDispA );
      getUpdateCoordFuncs_   ( uCoordsB, hCoords, elemDispB );

      coords(ALL,slice(0,nodeCount2_))   = uCoordsA;
      coords(ALL,slice(nodeCount2_,END)) = uCoordsB;
    }

    //System::out() << "disp A " << elemDispA << "\n";
    //System::out() << "disp B " << elemDispB << "\n";
    //System::out() << "B-A    " << elemDispJump << "\n";

    // get the integration points and weights on the cohesive surface

    Matrix functions = shape_->getShapeFunctions     ();
                       shape_->getIntegrationWeights ( ipWeights, coords );

    Matrix ischeme   = shape_->getIntegrationScheme ();		       

    mpoints  = materialMap_[ie];
    material = materials_[mpoints[0]];

    sshape_->evalShapeGradients ( ff, grads, center );

    getTransformationMatrix_ ( Q, oppositeVertices_[ie], grads, coords );

    elemMat    = 0.0;
    elemForce  = 0.0;

    // loop over these integration points

    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      // compute N-matrix from shape functions

      getShapeFuncs_ ( N, functions(ALL,ip) );

      // compute the displacement jump at integration point ip
      // using interpolation with shape functions

      jump[0] = dot ( functions(ALL,ip),  elemDispJump[slice(0,END,3)] );
      jump[1] = dot ( functions(ALL,ip),  elemDispJump[slice(1,END,3)] );
      jump[2] = dot ( functions(ALL,ip),  elemDispJump[slice(2,END,3)] );

      // get the correct material point of this int point

      cpoint = mpoints[1+ip]; 

      //sshape_->evalShapeGradients ( grads, ischeme(slice(1,2),ip), jump );
      //sshape_->evalShapeGradients ( grads, center, center );

      //getTransformationMatrix_ ( Q, oppositeVertices_[ie], grads,
//	                         coords(ALL,slice(nodeCount2_,END)) );

      // compute traction and stiffness in cohesive material

      jump = matmul ( Qt, jump ); 
      material->update ( traction, coheStiff, jump, cpoint );

      // transform the cohesive tangent matrix to global system

      coheStiff   = mc3.matmul ( Q, coheStiff, Qt ); 

      //System::out() << "global tangent:  " << coheStiff << "\n";
      //System::out() << "global traction: " << mc1.matmul ( Q, traction) << "\n\n";

      // now, compute the cohesive matrix and internal force
      
      wi         = ipWeights[ip];
      elemForce += wi * mc2.matmul ( Nt, Q, traction  );
      elemMat   += wi * mc3.matmul ( Nt, coheStiff, N );
    }

    // assembly of forces due to cohesive traction

    select ( force, idofsA ) += elemForce;
    select ( force, idofsB ) -= elemForce;

    // assembly of stiffness due to cohesive traction (elemDispJump)

    mbuilder.addBlock ( idofsA, idofsA, elemMat );
    mbuilder.addBlock ( idofsB, idofsB, elemMat );

    elemMat = -elemMat;

    mbuilder.addBlock ( idofsA, idofsB, elemMat );
    mbuilder.addBlock ( idofsB, idofsA, elemMat );
  }
}

//-----------------------------------------------------------------------
//   getIntForce_
//-----------------------------------------------------------------------


void InterfaceElement3DModel::getIntForce_

  ( const Vector&   force,
    const Vector&   disp )

{
  Matrix	coords     ( 3, nodeCount_ ); // coords in the undeformed configuration
  Matrix	uCoords    ( 3, nodeCount_ ); // coords in the deformed configuration
  Matrix	hCoords    ( 3, nodeCount2_ ); // coords in the undeformed config of 1 surface
  Matrix	uCoordsA   ( 3, nodeCount2_ ); // coords in the deformed config of surface +
  Matrix	uCoordsB   ( 3, nodeCount2_ ); // coords in the deformed config of surface -
  
  Matrix	coheStiff  ( 3, 3 ); // T matrix relates traction and rate of disp jump
  Matrix	Q          ( 3, 3 ); // transformation matrix
  Matrix	Qt         = Q.transpose ( ); 

  Vector	elemForce  ( nodeCount2_ * 3 );
  
  Vector	elemDispA   ( dofCount2_ ); // nodal displacement of upper face
  Vector	elemDispB   ( dofCount2_ ); // nodal displacement of lower face
  Vector	elemDispJump( dofCount2_ ); // nodal displacement jump 

  Vector	traction   ( 3 );          // traction at 1 integration point 
  Vector        jump       ( 3 );          // displacement jump at 1 integration point  

  IdxVector	inodes     ( nodeCount_  ); // inodes of both faces
  IdxVector	inodesA    ( nodeCount2_ ); // inodes of upper face
  IdxVector	inodesB    ( nodeCount2_ ); // inodes of lower face
  IdxVector	idofsA     ( dofCount2_  ); // dofs of upper face
  IdxVector	idofsB     ( dofCount2_  ); // dofs of lower face
  
  Vector	ipWeights  ( ipCount_    );
  Vector	ff         ( nodeCount2_ );

  Matrix	N          ( 3, nodeCount2_ * 3 );
  Matrix	Nt         = N.transpose ( ); 

  Matrix        grads      ( 2, nodeCount2_ ); // gradients of shape functions
                                               // used to compute the transformation matrix Q
  MChain2	mc2; MChain1 mc1;
  MChain3	mc3;

  int		cpoint;

  Ref<CohesiveMaterial> material;

  std::vector<int>      mpoints;

  double                wi;

  Vector        center(2); center[0] = center[1] = 0.;
  Vector        null(2);

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < elemCount_; ie++ )
  {
    // if this element is inactive, then discard it 

    if ( activeElems_[ie] == 0 ) continue;

    // Get the element connectivity 
    // nodal coordinates and nodal DOFs of element ie.

    inodes  = ielems_(ie,ALL); 
    inodesA = inodes[slice(BEGIN,nodeCount2_)]; 
    inodesB = inodes[slice(nodeCount2_,END)  ];
    
    // get coords in undeformed config

    nodes_.getSomeCoords ( coords, inodes );

    hCoords = coords(ALL,slice(nodeCount2_,END));
    
    // get the dofs associated with nodes of each element

    dofs_->getDofIndices ( idofsA, inodesA, dofTypes_ );
    dofs_->getDofIndices ( idofsB, inodesB, dofTypes_ );
    
    // get nodal displacements of two elements 
    // need both at the same time to compute the displacement jump

    elemDispA     = select ( disp, idofsA );
    elemDispB     = select ( disp, idofsB ); 

    elemDispJump  = elemDispA - elemDispB; 

    // if large displacement is needed,

    if ( largeDisp_ )
    {
      getUpdateCoordFuncs_   ( uCoordsA, hCoords, elemDispA );
      getUpdateCoordFuncs_   ( uCoordsB, hCoords, elemDispB );

      coords(ALL,slice(0,nodeCount2_))   = uCoordsA;
      coords(ALL,slice(nodeCount2_,END)) = uCoordsB;
    }

    //System::out() << "disp A " << elemDispA << "\n";
    //System::out() << "disp B " << elemDispB << "\n";
    //System::out() << "B-A    " << elemDispJump << "\n";

    // get the integration points and weights on the cohesive surface

    Matrix functions = shape_->getShapeFunctions     ();
                       shape_->getIntegrationWeights ( ipWeights, coords );

    Matrix ischeme   = shape_->getIntegrationScheme ();		       

    mpoints  = materialMap_[ie];
    material = materials_[mpoints[0]];

    sshape_->evalShapeGradients ( ff, grads, center );

    getTransformationMatrix_ ( Q, oppositeVertices_[ie], grads, coords );

    elemForce  = 0.0;

    // loop over these integration points

    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      // compute N-matrix from shape functions

      getShapeFuncs_ ( N, functions(ALL,ip) );

      // compute the displacement jump at integration point ip
      // using interpolation with shape functions

      jump[0] = dot ( functions(ALL,ip),  elemDispJump[slice(0,END,3)] );
      jump[1] = dot ( functions(ALL,ip),  elemDispJump[slice(1,END,3)] );
      jump[2] = dot ( functions(ALL,ip),  elemDispJump[slice(2,END,3)] );

      // get the correct material point of this int point

      cpoint = mpoints[1+ip]; 

      //sshape_->evalShapeGradients ( grads, ischeme(slice(1,2),ip), jump );
      //sshape_->evalShapeGradients ( grads, center, center );

      //getTransformationMatrix_ ( Q, oppositeVertices_[ie], grads,
//	                         coords(ALL,slice(nodeCount2_,END)) );

      // compute traction and stiffness in cohesive material

      jump = matmul ( Qt, jump ); 
      material->update ( traction, coheStiff, jump, cpoint );

      // transform the cohesive tangent matrix to global system

      coheStiff   = mc3.matmul ( Q, coheStiff, Qt ); 

      //System::out() << "global tangent:  " << coheStiff << "\n";
      //System::out() << "global traction: " << mc1.matmul ( Q, traction) << "\n\n";

      // now, compute the cohesive matrix and internal force
      
      wi         = ipWeights[ip];
      elemForce += wi * mc2.matmul ( Nt, Q, traction  );
    }

    // assembly of forces due to cohesive traction

    select ( force, idofsA ) += elemForce;
    select ( force, idofsB ) -= elemForce;
  }
}

// ---------------------------------------------------------------------
//   getTransformationMatrix_
// ---------------------------------------------------------------------
// Attention with the direction of the local coordinate system
// It must be consistent with the way the disp jump is computed
// disp jump = disp of upper face - disp of lower face
// => the normal direction points from lower face to upper face


void InterfaceElement3DModel::getTransformationMatrix_

  ( const Matrix& Q,
          int     oppVertex,
    const Matrix& xgrads,
    const Matrix& coords )
{
  using jem::numeric::crossProduct;
  using jem::numeric::dotProduct;
  using namespace jive::geom;

  Tuple<double,3>     n, s1, s2;
  Tuple<double,3>     p0, p1, v0;

  Vector              p(3);
  double              xg0, xg1;

  nodes_.getNodeCoords ( p, oppVertex );

  p0[0]    = p[0];        p0[1] = p[1];        p0[2] = p[2]; 
  p1[0]    = coords(0,0); p1[1] = coords(1,0); p1[2] = coords(2,0);

  v0       = p0 - p1;

  s1 = 0.0;
  s2 = 0.0;

  Matrix  midCoord(3,nodeCount2_);

  midCoord = 0.5 * ( coords(ALL,slice(0,nodeCount2_)) + 
	             coords(ALL,slice(nodeCount2_,END)) );

  for ( int i = 0; i < nodeCount2_; i++ )
  {
    xg0   = xgrads(0,i);
    xg1   = xgrads(1,i);

    s1[0] += midCoord(0,i) * xg0;
    s1[1] += midCoord(1,i) * xg0;
    s1[2] += midCoord(2,i) * xg0;

    s2[0] += midCoord(0,i) * xg1;
    s2[1] += midCoord(1,i) * xg1;
    s2[2] += midCoord(2,i) * xg1;
  }

  n        = crossProduct ( s1, s2 );
  double a = ::sqrt ( dotProduct( n, n   ) );
  double b = ::sqrt ( dotProduct( s1, s1 ) );

  if ( jem::isTiny( a ) || jem::isTiny( b ) )
  {
    zeroVectorError ( "normal vector", "normal" );
  }

  // normalize n and s1

  n   = (1.0 / a) * n;
  s1  = (1.0 / b) * s1;

  // flip direction if n points outward

  if ( dotProduct ( n, v0 ) < 0 )
  {
    n *= -1.;
  }

  // make s2 orthogonal to n and s1

  s2  = crossProduct ( n, s1 );

  // build the transformation matrix Q
  
  Q(0,0) = n[0]; Q(0,1) = s1[0]; Q(0,2) = s2[0];
  Q(1,0) = n[1]; Q(1,1) = s1[1]; Q(1,2) = s2[1];
  Q(2,0) = n[2]; Q(2,1) = s1[2]; Q(2,2) = s2[2]; 
  
  //System::out() << Q << "\n";
}

// ---------------------------------------------------------------------
//   getTransformationMatrix_
// ---------------------------------------------------------------------

void InterfaceElement3DModel::getTransformationMatrix_

  ( const Matrix& Q,
          int     oppVertex,
    const Matrix& coord )
{
  using jem::numeric::crossProduct;
  using jem::numeric::dotProduct;
  using namespace jive::geom;

  Vector              p(3);
  Tuple<double,3>     n;
  Tuple<double,3>     p0, p1, p2, p3;
  Tuple<double,3>     v0, v1, v2;

  nodes_.getNodeCoords ( p, oppVertex );

  p0[0]    = p[0];       p0[1] = p[1];       p0[2] = p[2]; 

  p1[0]    = coord(0,0); p1[1] = coord(1,0); p1[2] = coord(2,0);
  p2[0]    = coord(0,1); p2[1] = coord(1,1); p2[2] = coord(2,1);
  p3[0]    = coord(0,2); p3[1] = coord(1,2); p3[2] = coord(2,2);

  v0       = p0 - p1;
  v1       = p2 - p1;
  v2       = p3 - p1;

  n        = crossProduct ( v1, v2 );
  double a = ::sqrt ( dotProduct( n, n   ) );
  double b = ::sqrt ( dotProduct( v1, v1 ) );

  if ( jem::isTiny( a ) || jem::isTiny( b ) )
  {
    zeroVectorError ( "normal vector", "normal" );
  }

  n   = (1.0 / a) * n;
  v1  = (1.0 / b) * v1;

  // flip the normal so that it points inward

  if ( dotProduct( n, v0) < 0. )
  {
    n *= -1.;
  }

  v2  = crossProduct ( n, v1 );

  Q(0,0) = n[0]; Q(0,1) = v1[0]; Q(0,2) = v2[0];
  Q(1,0) = n[1]; Q(1,1) = v1[1]; Q(1,2) = v2[1];
  Q(2,0) = n[2]; Q(2,1) = v1[2]; Q(2,2) = v2[2]; 
  
  //System::out() << Q << "\n";
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

void InterfaceElement3DModel::initConstraints_ ()
{
  cout << "initializing the constraint of bulk interface dofs...\n";

  if ( materials_.size () == 1 ) return;

  // loop over elements, do only for bulk elements
  // which are elements with material = 0

  const int            nodeCount = dupNodes_.size ();

  int                  xdof,  ydof, zdof;
  int                  xdofM, ydofM, zdofM;
  int                  inCount;
  int                  inter; // 0: standard node
                              // 1: interface node

  int xtype  = dofTypes_[0];
  int ytype  = dofTypes_[1];
  int ztype  = dofTypes_[2];

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
      zdofM = dofs_->getDofIndex ( nodes[0], ztype ) ; 

      for ( int jn = 1; jn < inCount-1; jn++ )
      {
	xdof = dofs_->getDofIndex ( nodes[jn], xtype ) ; 
	ydof = dofs_->getDofIndex ( nodes[jn], ytype ) ; 
	zdof = dofs_->getDofIndex ( nodes[jn], ztype ) ; 
	
	addConstraint_ ( xdof,  ydof,  zdof, xdofM, ydofM, zdofM );
      }
    }
    else  // interfacial node
    {
      xdofM = dofs_->getDofIndex ( nodes[1], xtype ) ; 
      ydofM = dofs_->getDofIndex ( nodes[1], ytype ) ; 
      zdofM = dofs_->getDofIndex ( nodes[1], ztype ) ; 

      for ( int jn = 2; jn < inCount-1; jn++ )
      {
	xdof = dofs_->getDofIndex ( nodes[jn], xtype ) ; 
	ydof = dofs_->getDofIndex ( nodes[jn], ytype ) ; 
	zdof = dofs_->getDofIndex ( nodes[jn], ztype ) ; 
	
	addConstraint_ ( xdof,  ydof,  zdof, xdofM, ydofM, zdofM );
      }
    }
  }

  cons_->compress ();
}

// ---------------------------------------------------------------------
//   constraintDirichletNodes_
// ---------------------------------------------------------------------

void InterfaceElement3DModel::constraintDirichletNodes_ ()
{
  if ( diriNodes_.size () == 0 ) return;

  if ( initConsDone_ ) return;

  cout << "Constraint Dirichlet nodes...\n";

  const int    diriNodesCount = diriNodes_.size ();

        int    index;

  int    xtype = dofTypes_[0];
  int    ytype = dofTypes_[1];
  int    ztype = dofTypes_[2];

  int    xdof, ydof, zdof, xdofM, ydofM, zdofM;
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
    zdofM = dofs_->getDofIndex ( nodes[0], ztype ) ; 
    
    for ( int jn = 1; jn < inCount-1; jn++ )
    {
      xdof = dofs_->getDofIndex ( nodes[jn], xtype ) ; 
      ydof = dofs_->getDofIndex ( nodes[jn], ytype ) ; 
      zdof = dofs_->getDofIndex ( nodes[jn], ztype ) ; 
      
      addConstraint_ ( xdof, ydof, zdof, xdofM, ydofM, zdofM );
    }
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

void InterfaceElement3DModel::updateBcConstraints_ ()
{
  if ( dispNodes_.size () == 0 ) return;

  cout << "updating constraints ...\n";

  const int    dispNodesCount = dispNodes_.size ();

        int    index;

	int    xtype = dofTypes_[0];
	int    ytype = dofTypes_[1];
	int    ztype = dofTypes_[2];

	int    xdof, ydof, zdof, xdofM, ydofM, zdofM;
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
      zdofM = dofs_->getDofIndex ( nodes[0], ztype ) ; 

      for ( int jn = 1; jn < inCount-1; jn++ )
      {
	xdof = dofs_->getDofIndex ( nodes[jn], xtype ) ; 
	ydof = dofs_->getDofIndex ( nodes[jn], ytype ) ; 
	zdof = dofs_->getDofIndex ( nodes[jn], ztype ) ; 
	
	addConstraint_ ( xdof, ydof, zdof, xdofM, ydofM, zdofM );
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


void InterfaceElement3DModel::checkCommit_

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


static Ref<Model>     newInterfaceElement3DModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<InterfaceElement3DModel> ( name, conf, props, globdat );
}

//-----------------------------------------------------------------------
//   declareInterfaceElement3DModel
//-----------------------------------------------------------------------

void declareInterfaceElement3DModel ()
{
  jive::model::ModelFactory::declare ( "Interface3D", & newInterfaceElement3DModel );
}

//-----------------------------------------------------------------------
//   makeBTriangleNC
//-----------------------------------------------------------------------
// From Frans van der Meer

Ref<BShape> makeBTriangleNC ( const Properties& props )
{
  using jive::Matrix;
  using jive::geom::BoundaryTriangle3;
  using jive::geom::BoundaryTriangle6;

  String     type;
  String     intScheme;
  Matrix     ischeme;

  props.find (      type,      "type" );
  props.find ( intScheme, "intScheme" );

  // Set up the Newton-Cotes scheme for the interior of the
  // triangle. The first matrix row contains the weights and the other
  // two rows contain the coordinates. Both the weights and
  // coordinates must be specified in the local coordinate system of a
  // standard triangle. The weights must add up to 0.5 (the area of a
  // standard triangle).

  if ( intScheme == "NewtonCotes3" )
  {
    ischeme.resize ( 3, 3 );

    ischeme(0,ALL) = 1.0 / 6.0;

    ischeme(1,0) = 0.0; ischeme(2,0) = 0.0;
    ischeme(1,1) = 1.0; ischeme(2,1) = 0.0;
    ischeme(1,2) = 0.0; ischeme(2,2) = 1.0;
  }
  else if ( intScheme == "NewtonCotes6" )   
  {
    ischeme.resize ( 3, 6 );

    // this is not correct !!!

    //ischeme(0,slice(0,END,2)) = 1.0 / 24.0;
    //ischeme(0,slice(1,END,2)) = 3.0 / 24.0;

    // this weight is taken from Sylvester's 1970 paper

    ischeme(0,slice(0,END,2)) = 0.;
    ischeme(0,slice(1,END,2)) = 1.0 / 6.0;

    ischeme(1,0) = 0.0; ischeme(2,0) = 0.0;
    ischeme(1,1) = 0.5; ischeme(2,1) = 0.0;
    ischeme(1,2) = 1.0; ischeme(2,2) = 0.0;
    ischeme(1,3) = 0.5; ischeme(2,3) = 0.5;
    ischeme(1,4) = 0.0; ischeme(2,4) = 1.0;
    ischeme(1,5) = 0.0; ischeme(2,5) = 0.5;
  }
  else
  {
    return NIL;
  }

  if      ( type == "BTriangle3" ) 
  {
    return BoundaryTriangle3::getShape ( "BTriangleNC3", ischeme );
  }
  else if ( type == "BTriangle6" )
  {
    return BoundaryTriangle6::getShape ( "BTriangleNC6", ischeme );
  }

  return NIL;
}

