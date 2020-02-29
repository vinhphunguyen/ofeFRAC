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
 *  NURBS elements with C0 continuity across boundaries are supported.
 *  Only quadratic and cubic NURBS are implemented however.
 *
 *  Shape object used for interface elements: InterfaceShape(Ref<BoundaryShape>)
 *  BoundaryShape functionalities will be enriched by InterfaceShape for example
 *  integration along the mid-line.
 *
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 5 July 2009
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
#include <jive/geom/InterfaceShape.h>
#include <jive/geom/InternalShape.h>
#include <jive/geom/BoundaryShape.h>
#include <jive/geom/ParametricEdge.h>
#include <jive/geom/StdLine.h>
#include <jive/mp/ItemMask.h>


#include "InterfaceElementModel.h"
#include "util/utilities.h"


using std::cout;

using namespace jem;
using jem::newInstance;
using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jive::model::StateVector;
using jive::util::Printer;
using jive::mp::ItemMask;


typedef MatmulChain<double,1>      MChain1;
typedef MatmulChain<double,2>      MChain2;
typedef MatmulChain<double,3>      MChain3;
typedef MatmulChain<double,4>      MChain4;


// -----------------------------------------------------------------
//   static data
// -----------------------------------------------------------------

const char* InterfaceElementModel::SHAPE_PROP     = "shape";
const char* InterfaceElementModel::MATERIAL_PROP  = "material";
const char* InterfaceElementModel::MESH_PROP      = "meshFile";
const char* InterfaceElementModel::THICKNESS_PROP = "thickness";
const char* InterfaceElementModel::DISP_NODES_PROP= "dNodeGroup";
const char* InterfaceElementModel::DIRICHLET_NODES_PROP= "dirNodeGroup";

// -----------------------------------------------------------------
//   constructor
// -----------------------------------------------------------------

InterfaceElementModel::InterfaceElementModel

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
  
  egroup_ = ElementGroup::get ( myConf, myProps, globdat, context );

  elems_  = egroup_.getElements ( );
  nodes_  = elems_.getNodes     ( );

  ielemCount_ = egroup_.size    ( );

  ielems_.resize ( ielemCount_ );
  ielems_ = egroup_.getIndices ();
  
  //nodes_   = NodeSet::find    ( globdat                   );
  dofs_    = XDofSpace::get   ( nodes_.getData(), globdat );
  cons_    = Constraints::get ( dofs_, globdat            );
  
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
  
  rank_      = nodes_.rank   ( );
  
  dofTypes_.resize ( rank_ );
  dofTypes_[0] = 0;
  dofTypes_[1] = 1;

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

  // building the mapping matrix ipMap_

  nodeCount_  = shape_->nodeCount ( ) ;
  ipCount_    = shape_->ipointCount ();
  dofCount_   = nodeCount_ ;
  nodeCount2_ = nodeCount_/2;
  
  System::out () << "Number of interface elements ..... " << ielemCount_ << "\n"
                 << "Number of node per element   ..... " << nodeCount_ 
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

  getShapeFuncs_  = getShapeFunc  ( rank_ );

  // get thickness

  thickness_ = 1.0;

  myProps.find ( thickness_, THICKNESS_PROP );
  myConf. set  ( THICKNESS_PROP, thickness_ );

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

}
     

InterfaceElementModel::~InterfaceElementModel()
{
}

// -----------------------------------------------------------------
//   configure
// -----------------------------------------------------------------

void InterfaceElementModel::configure

    ( const Properties&       props,
      const Properties&       globdat )
{
}

// -----------------------------------------------------------------
//   getConfig
// -----------------------------------------------------------------


void  InterfaceElementModel::getConfig

    ( const Properties&       conf )             const
{
}    
    
// -----------------------------------------------------------------
//   takeAction
// -----------------------------------------------------------------    

bool InterfaceElementModel::takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat )
{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  if ( action == Actions::INIT )
  {
    // does not work well, use Discontinuous Galerkin/CZM model instead
    //initConstraints_ (  );

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


void InterfaceElementModel::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   disp0,
    const Vector&   disp )

{
  Matrix      coords     ( rank_, nodeCount_    );
  
  Matrix      coheStiff  ( rank_,     rank_     ); // T matrix relates traction and rate of disp jump
  Matrix      Q          ( rank_,     rank_     ); // transformation matrix
  Matrix      elemMat    ( nodeCount_,nodeCount_ );
  Vector      elemForce  ( nodeCount_ );
  
  Vector      elemDispA   ( dofCount_  ); // nodal displacement of upper face
  Vector      elemDispA0  ( dofCount_  ); // nodal displacement of upper face
  Vector      elemDispB   ( dofCount_  ); // nodal displacement of lower face
  Vector      elemDispB0  ( dofCount_  ); // nodal displacement of lower face
  Vector      jump        ( rank_ );      // displacement jump at 1 integration point  
  Vector      jump0       ( rank_ );      // displacement jump at 1 integration point  
  Vector      djump       ( rank_ );      // displacement jump at 1 integration point  
  Vector      elemDJump   ( dofCount_ ); // nodal displacement jump 
  Vector      elemDJump0  ( dofCount_ ); // nodal displacement jump 
  Vector      traction    ( rank_ );

  IdxVector   inodes     ( nodeCount_  ); // inodes of both faces
  IdxVector   inodesA    ( nodeCount2_ ); // inodes of upper face
  IdxVector   inodesB    ( nodeCount2_ ); // inodes of lower face
  IdxVector   idofsA     ( dofCount_   ); // dofs of upper face
  IdxVector   idofsB     ( dofCount_   ); // dofs of lower face
  
  Vector      ipWeights  ( ipCount_   );

  Matrix      N          ( rank_, nodeCount_ );
  Matrix      Nt         = N.transpose ( ); 

  MChain2     mc2;
  MChain3     mc3;

  idx_t         iMat, cpoint;

  Ref<CohesiveMaterial> material;

  std::vector<int>      mpoints;

  double                wi;

  Matrix      normals    ( rank_, nodeCount2_   );
  // Iterate over all elements assigned to this model.
  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    // Get the global element index.

    idx_t  ielem = ielems_[ie];
    elems_.getElemNodes  ( inodes, ielem  );

    // Get the element connectivity 
    // nodal coordinates and nodal DOFs of element ie.

    inodesA = inodes[slice(BEGIN,nodeCount2_)]; 
    inodesB = inodes[slice(nodeCount2_,END)  ];
                
    nodes_.getSomeCoords ( coords, inodes );
    
    //System::out() << inodes << "\n";
    //System::out() << inodesA << "\n";
    //System::out() << inodesB << "\n";    
   
    // get the dofs associated with nodes of each element

    dofs_->getDofIndices ( idofsA, inodesA, dofTypes_ );
    dofs_->getDofIndices ( idofsB, inodesB, dofTypes_ );
    
    // get nodal displacements of two elements 
    // need both at the same time to compute the displacement jump

    elemDispA     = select ( disp,  idofsA );
    elemDispA0    = select ( disp0, idofsA );
    elemDispB     = select ( disp,  idofsB ); 
    elemDispB0    = select ( disp0, idofsB ); 

    // compute the relative displacement 
    // ATTENTION: 
    //   if jump = uA - uB => select ( force, idofsA ) += elemForce;
    //                        select ( force, idofsB ) -= elemForce;
    //   and vice versa                        

    elemDJump   = elemDispA  - elemDispB;
    elemDJump0  = elemDispA0 - elemDispB0;

    //System::out() << "dispA: " << elemDispA << "\n";
    //System::out() << "dispB: " << elemDispB << "\n";
    //System::out() << "jumAB: " << elemDispJump << "\n";

    // get the shape funcs, integration points and weights on the cohesive surface

    Matrix functions = shape_->getShapeFunctions     ();
    shape_->getIntegrationWeights ( ipWeights, coords );
    ipWeights *= thickness_;

    getTransformationMatrix_ ( Q, coords(ALL,slice(BEGIN,nodeCount2_)) );

    //mpoints  = materialMap_[ie];
    
    iMat     = elemMatMap_[ielem]; 
    material = materials_[iMat];
    
    //System::out() << "iMat: " << iMat << "\n";
    //System::out() << "ielem: " << ielem << "\n";

    elemMat    = 0.0;
    elemForce  = 0.0;
    
    //shape_->getNormals ( normals, ipWeights,  coords );
    //System::out() << "normals: " << normals << "\n";

    // loop over these integration points

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // compute N-matrix from shape functions

      getShapeFuncs_ ( N, functions(ALL,ip) );

      // compute the displacement jump at integration point ip
      // using interpolation with shape functions

      jump[0]  = dot ( functions(ALL,ip),  elemDJump[slice(0,END,rank_)] );
      jump[1]  = dot ( functions(ALL,ip),  elemDJump[slice(1,END,rank_)] );
      
      jump0[0] = dot ( functions(ALL,ip),  elemDJump0[slice(0,END,rank_)] );
      jump0[1] = dot ( functions(ALL,ip),  elemDJump0[slice(1,END,rank_)] );

      djump    = jump - jump0;
    
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


void InterfaceElementModel::getIntForce_

  ( const Vector&   force,
    const Vector&   disp0,
    const Vector&   disp )

{
  Matrix      coords     ( rank_, nodeCount_    );
  Matrix      coheStiff  ( rank_,     rank_     ); // T matrix relates traction and rate of disp jump
  Matrix      Q          ( rank_,     rank_     ); // transformation matrix

  Vector      elemForce  ( nodeCount_ );
  Vector      elemDispA  ( dofCount_  ); // nodal displacement of upper face
  Vector      elemDispB  ( dofCount_  ); // nodal displacement of lower face
  Vector      jump       ( rank_ );          // displacement jump at 1 integration point  
  Vector      elemDispJump( dofCount_ );          // nodal displacement jump 
  Vector      traction   ( rank_ );

  IdxVector   inodes     ( nodeCount_  ); // inodes of both faces
  IdxVector   inodesA    ( nodeCount2_ ); // inodes of upper face
  IdxVector   inodesB    ( nodeCount2_ ); // inodes of lower face
  IdxVector   idofsA     ( dofCount_   ); // dofs of upper face
  IdxVector   idofsB     ( dofCount_   ); // dofs of lower face
  
  Vector      ipWeights  ( ipCount_   );

  Matrix      N          ( rank_, nodeCount_ );
  Matrix      Nt         = N.transpose ( ); 

  MChain2     mc2;

  int         iMat, cpoint;

  Ref<CohesiveMaterial> material;

  std::vector<int>      mpoints;

  double                wi;

  Matrix      normals    ( rank_, nodeCount2_ );
  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount_; ie++ )
  {
    // Get the element connectivity 
    // nodal coordinates and nodal DOFs of element ie.

    int  ielem = ielems_[ie];
    elems_.getElemNodes  ( inodes, ielem  );
    inodesA = inodes[slice(BEGIN,nodeCount2_)]; 
    inodesB = inodes[slice(nodeCount2_,END)  ];
                
    nodes_.getSomeCoords ( coords, inodes );
    
    //System::out() << inodes << "\n";
    //System::out() << inodesA << "\n";
    //System::out() << inodesB << "\n";    
   
    // get the dofs associated with nodes of each element

    dofs_->getDofIndices ( idofsA, inodesA, dofTypes_ );
    dofs_->getDofIndices ( idofsB, inodesB, dofTypes_ );
    
    // get nodal displacements of two elements 
    // need both at the same time to compute the displacement jump

    elemDispA     = select ( disp, idofsA );
    elemDispB     = select ( disp, idofsB ); 

    // compute the relative displacement 
    // ATTENTION: 
    //   if jump = uA - uB => select ( force, idofsA ) += elemForce;
    //                        select ( force, idofsB ) -= elemForce;
    //   and vice versa                        

    elemDispJump  = elemDispA - elemDispB;

    //System::out() << "dispA: " << elemDispA << "\n";
    //System::out() << "dispB: " << elemDispB << "\n";
    //System::out() << "jumAB: " << elemDispJump << "\n";

    // get the integration points and weights on the cohesive surface

    Matrix functions = shape_->getShapeFunctions     ();
    shape_->getIntegrationWeights ( ipWeights, coords );
    ipWeights *= thickness_;

    getTransformationMatrix_ ( Q, coords(ALL,slice(BEGIN,nodeCount2_)) );

    iMat     = elemMatMap_[ielem]; 
    material = materials_[iMat];

    elemForce  = 0.0;
    
    //shape_->getNormals ( normals, ipWeights,  coords );
    //System::out() << "normals: " << normals << "\n";

    // loop over these integration points

    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      // compute N-matrix from shape functions

      getShapeFuncs_ ( N, functions(ALL,ip) );

      // compute the displacement jump at integration point ip
      // using interpolation with shape functions

      for ( int ir = 0; ir < rank_; ++ir )
      {
        jump[ir] = dot ( functions(ALL,ip),  elemDispJump[slice(ir,END,rank_)] );
      }
    
      //System::out() << "jumAB: " << jump << "\n";


      // get the correct material point of this int point

      cpoint = ipMpMap_ (ielem,ip); //System::out() << "cpoint: " << cpoint << "\n\n";

      // compute traction and stiffness in cohesive material

      jump = matmul ( Q, jump );

      material->update ( traction, coheStiff, jump, cpoint );

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


void InterfaceElementModel::getTransformationMatrix_

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

void InterfaceElementModel::initConstraints_ ()
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

void InterfaceElementModel::constraintDirichletNodes_ ()
{
  if ( diriNodes_.size () == 0 ) return;
  if ( initConsDone1_ ) return;

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
    //Rvalx = cons_->hasRvalue ( xdofM );
    //Rvaly = cons_->hasRvalue ( ydofM );

    for ( int jn = 1; jn < inCount-1; jn++ )
    {
      xdof = dofs_->getDofIndex ( nodes[jn], xtype ) ; 
      ydof = dofs_->getDofIndex ( nodes[jn], ytype ) ; 
      
      cons_->addConstraint ( xdof, xdofM, 1.0 );
      cons_->addConstraint ( ydof, ydofM, 1.0 );
    }
  }

  initConsDone1_ = true;
}


//-----------------------------------------------------------------------
//   updateBcConstraint_
//-----------------------------------------------------------------------

// NOTE: this has to be called every time steps!!!
// There are constraints between duplicated nodes and Dirichlet nodes
// When the prescribed displacement of Dirichlet nodes change, it is 
// necessary to update their slave dofs so that they're still tied
// together. This is the case for displacement control only.

void InterfaceElementModel::updateBcConstraints_ ()
{
  if ( dispNodes_.size () == 0 ) return;

  cout << "updating constraints ...\n";

  const int    dispNodesCount = dispNodes_.size ();

        int    index;
	int    xtype = dofTypes_[0];
	int    ytype = dofTypes_[1];
	int    xdof, ydof, xdofM, ydofM;
        int    inCount;

  std::vector<int> nodes;

  for ( int in = 0; in < dispNodesCount; in++ )
  {
    index   = dispNodes_[in];

    nodes   = dupNodes_[index];
    inCount = nodes.size ();
    //inter   = nodes[inCount-1];

    if ( inCount == 2 ) continue;

    xdofM = dofs_->getDofIndex ( nodes[0], xtype ) ; 
    ydofM = dofs_->getDofIndex ( nodes[0], ytype ) ; 
    //Rvalx = cons_->hasRvalue ( xdofM );
    //Rvaly = cons_->hasRvalue ( ydofM );

    for ( int jn = 1; jn < inCount-1; jn++ )
    {
      xdof = dofs_->getDofIndex ( nodes[jn], xtype ) ; 
      ydof = dofs_->getDofIndex ( nodes[jn], ytype ) ; 

      cons_->addConstraint ( xdof, xdofM, 1.0 );
      cons_->addConstraint ( ydof, ydofM, 1.0 );
    }
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


void InterfaceElementModel::checkCommit_

  ( const Properties& globdat,
    const Properties& params )
{
  using jive::model::ActionParams;

  cout << "Checking activation of interface elements ...\n";

  Ref<CohesiveMaterial>  mat = materials_[0];

  bool                   isActivated = false; 

  int                    n1, n2;

  int                    xtype  = dofTypes_[0];
  int                    ytype  = dofTypes_[1];

  int                    xdof, ydof, xdof1, ydof1, xdof2, ydof2;

  IdxVector              inodes ( nodeCount_ );

  Vector                 traction(2);
  Vector                 fext;
  Vector                 ipWeights (ipCount_);

  Matrix                 coords (rank_, nodeCount_);

  //params.get ( fext,    ActionParams::EXT_VECTOR ); //error fext not defined
  globdat.get ( fext,    ActionParams::INT_VECTOR ); //error fext not defined

  // loop over all interface elements

  for ( int ie = 0; ie < ielemCount_; ie++ )
  {
    // ignore already active one

    if ( activeElems_[ie] == 1 ) continue;

    int  ielem = ielems_[ie];
    elems_.getElemNodes  ( inodes, ielem  );

    nodes_. getSomeCoords         ( coords,    inodes );
    shape_->getIntegrationWeights ( ipWeights, coords );
    ipWeights *= thickness_;

    // loop over integration points

    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      n1 = inodes[ip];
      n2 = inodes[ip+2];

      xdof1 = dofs_->getDofIndex ( n1, xtype );
      ydof1 = dofs_->getDofIndex ( n1, ytype );

      xdof2 = dofs_->getDofIndex ( n2, xtype );
      ydof2 = dofs_->getDofIndex ( n2, ytype );

      // compute traction 

      traction[0] = 0.5 / ipWeights[ip] * ( fext[xdof1] + fext[xdof2] ) ;
      traction[1] = 0.5 / ipWeights[ip] * ( fext[ydof1] + fext[ydof2] ) ;

     // check activation 

      //isActivated = mat->evalFailure ( traction );

      if ( isActivated ) break;
    }

    // if activation needed
    // remove constraint so that nodes can move apart
    // and mark ie as active (=1)

    if ( isActivated )
    {
      for ( int in = 0; in < 2; in++ )
      {
	n1 = inodes[in];
	n2 = inodes[in+2];

        xdof = dofs_->getDofIndex ( n2, xtype );
        ydof = dofs_->getDofIndex ( n2, ytype );

	cons_->eraseConstraint ( xdof );
	cons_->eraseConstraint ( ydof );
      }

      activeElems_[ie] = 1;

      System::out() << "There is one interface activated\n";

      break; // only one element activated per time step
    }
  }
}

//-----------------------------------------------------------------------
//   initializeIPMPMap_
//-----------------------------------------------------------------------

void  InterfaceElementModel::initializeIPMPMap_ ( )

{
        int   ipoint;

  IdxVector   matMap ( materials_.size ( ) );

  matMap = 0;

  for ( int ie = 0; ie < ielemCount_; ie++ )
  {
    int  ielem = ielems_[ie];

    // get correct material ID

    int matID = elemMatMap_[ielem]; 
    ipoint    = matMap[matID];     // System::out() << "ielem: " << ielem<< "\n\n";

    // loop over integration points 

    for ( int ip = 0; ip < ipCount_; ip++, ipoint++, matMap[matID]++  )
    {
      ipMpMap_ (ielem,ip) = ipoint; //System::out() << "ipMpMap_: " << ipoint << "\n\n";
    }

  }

}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newInterfaceModel
//-----------------------------------------------------------------------


static Ref<Model>     newInterfaceElementModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<InterfaceElementModel> ( name, conf, props, globdat );
}

//-----------------------------------------------------------------------
//   declareInterfaceElementModel
//-----------------------------------------------------------------------

void declareInterfaceElementModel ()
{
  jive::model::ModelFactory::declare ( "Interface", & newInterfaceElementModel );
}



