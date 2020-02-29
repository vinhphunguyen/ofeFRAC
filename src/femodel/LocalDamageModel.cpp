/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the local damage model
 *  using finite elements. This local models is regularised
 *  by using the mesh adjusted softening law.
 *  Currently, only implemented for T3 and Q4 elements
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *
 */

#include <jem/base/Array.h>
#include <jem/base/Error.h>
#include <jem/base/System.h>
#include <jem/base/IllegalInputException.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/sparse/utilities.h>
#include <jem/numeric/func/Function.h>
#include <jem/util/Flex.h>
#include <jem/util/Properties.h>
#include <jive/util/utilities.h>
#include <jive/util/XTable.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/Assignable.h>
#include <jive/util/FuncUtils.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/geom/Geometries.h>
#include <jive/geom/InternalShape.h>
#include <jive/geom/IShapeFactory.h>
#include <jive/fem/ElementGroup.h>

#include "FemModels.h"
#include "util/utilities.h"
#include "material/DamageMaterial.h"


using namespace jem;

using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jem::numeric::Function;
using jem::util::SparseArray;
using jem::util::Flex;
using jem::util::Properties;
using jive::Vector;
using jive::IdxVector;
using jive::Matrix;
using jive::Cubix;
using jive::util::XTable;
using jive::util::XDofSpace;
using jive::util::Assignable;
using jive::algebra::MatrixBuilder;
using jive::model::Model;
using jive::geom::IShape;
using jive::fem::NodeSet;
using jive::fem::ElementSet;
using jive::fem::ElementGroup;


//=======================================================================
//   typedefs
//=======================================================================

// Some handy aliases to avoid some typing

typedef MatmulChain<double,1>  MChain1;
typedef MatmulChain<double,2>  MChain2;
typedef MatmulChain<double,3>  MChain3;

typedef ElementSet             ElemSet;
typedef ElementGroup           ElemGroup;



//=======================================================================
//   class LocalDamageModel
//=======================================================================


class LocalDamageModel : public Model
{
 public:

  typedef LocalDamageModel  Self;
  typedef Model             Super;

  static const char*        DOF_NAMES[3];

  static const char*        SHAPE_PROP;
  static const char*        MATERIAL_PROP;

  //! Constructor
                            LocalDamageModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  //! Configure itself by reading data from properties file

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  //! 

  virtual void              getConfig

    ( const Properties&       conf )             const;

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );


 protected:

  virtual                  ~LocalDamageModel       ();


 private:

  //! getMatrix_: compute the tangent stiffness and internal force vector
  //  Compute K and fint for the updated displacement disp

  void                      getMatrix_

    ( MatrixBuilder&          mbuilder,
      const Vector&           force,
      const Vector&           disp );

  //! compute the consistent mass matrix

  void                      getMatrix20_

    ( MatrixBuilder&          mbuilder );
 
  /*
   *  compute the lumped mass matrix
   */

  void                      getMatrix21_

    ( MatrixBuilder&          mbuilder );

  void                      checkCommit_

    ( const Properties&       params );

  // compute equivalent element size for all elements

  void                      computeEqvElemSize_();

  //! Make stress, strain and damage tables for visualisation

  bool                      getTable_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      getStress_

    ( XTable&                 table,
      const Vector&           weights,
      const Vector&           disp );

  void                      getHistory_

    ( XTable&                 table,
      const Vector&           weights );

 private:

  Assignable<ElemGroup>     egroup_;
  Assignable<ElemSet>       elems_;
  Assignable<NodeSet>       nodes_;

  int                       rank_;

  Ref<IShape>               shape_;

  Ref<XDofSpace>            dofs_;

  IdxVector                 dofTypes_;
  IdxVector                 ielems_;

  Vector                    eqvElemSize_; // equivalent element sizes

  // material used in damage model must be a Damage material

  Ref<DamageMaterial>       material_;

  ShapeGradsFunc            getShapeGrads_;
  ShapeFunc                 getShapeFuncs_;

  //bool                      updated_;

  int                       elemCount_;   // total number of elements 
  int                       nodeCount_;   // num of nodes per element
  int                       ipCount_;     // num of GP per element
  int                       strCount_;    // size of strain vector
  int                       dofCount_;    // num of DOFs per element
};


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  LocalDamageModel::DOF_NAMES[3]     = { "dx", "dy", "dz" };

const char*  LocalDamageModel::SHAPE_PROP       = "shape";
const char*  LocalDamageModel::MATERIAL_PROP    = "material";



//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


LocalDamageModel::LocalDamageModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super ( name )

{
  using jive::util::joinNames;
  using jive::util::FuncUtils;
  using jive::geom::IShapeFactory;

  Properties    myProps = props.getProps  ( myName_ );
  Properties    myConf  = conf .makeProps ( myName_ );

  const String  context = getContext ();

  // Get the element group assigned to this model.

  egroup_ = ElemGroup::get ( myConf, myProps, globdat, context );

  ielems_.resize ( egroup_.getIndices ().size() );

  ielems_ = egroup_.getIndices  ();
  elems_  = egroup_.getElements ();
  nodes_  = elems_ .getNodes    ();
  rank_   = nodes_ .rank        ();


  // Create an internal shape object for computing the element shape
  // functions.

  shape_ = IShapeFactory::newInstance (
    joinNames ( myName_, SHAPE_PROP ),
    conf,
    props
  );

  // Make sure that the rank of the shape matches the rank of the
  // mesh.

  if ( shape_->globalRank() != rank_ )
  {
    throw IllegalInputException (
      context,
      String::format (
        "shape has invalid rank: %d (should be %d)",
        shape_->globalRank (),
        rank_
      )
    );
  }

  // compute some constants ...

  elemCount_  = ielems_.size        ();
  nodeCount_  = shape_->nodeCount   ();
  ipCount_    = shape_->ipointCount ();
  strCount_   = STRAIN_COUNTS[rank_];
  dofCount_   = rank_ * nodeCount_;

  // Make sure that the number of spatial dimensions (the rank of the
  // mesh) is valid.

  if ( rank_ < 1 || rank_ > 3 )
  {
    throw IllegalInputException (
      context,
      String::format (
        "invalid node rank: %d (should be 1, 2 or 3)", rank_
      )
    );
  }

  // Make sure that each element has the same number of nodes as the
  // shape object.

  elems_.checkSomeElements ( context, ielems_, nodeCount_ );

  // Get the degree of freedom (DOF) space, and add the DOFs of this
  // model.

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  dofTypes_.resize ( rank_ );

  for ( int i = 0; i < rank_; i++ )
  {
    dofTypes_[i] = dofs_->addType ( DOF_NAMES[i] );
  }

  dofs_->addDofs (
    elems_.getUniqueNodesOf ( ielems_ ),
    dofTypes_
  );

  // Create a material model object.

  Ref<Material> tmpMat = newMaterial ( MATERIAL_PROP, myConf, myProps, globdat );
  material_            = dynamicCast<DamageMaterial> ( tmpMat );

  if ( material_ == NIL )
  {
    throw IllegalInputException (
      context,
      "material type is not a damage material"
    );
  }

  Properties  matProps = myProps.findProps ( MATERIAL_PROP );
  Properties  matConf  = myConf.makeProps  ( MATERIAL_PROP );

  material_->configure    ( matProps, globdat );
  material_->getConfig    ( matConf , globdat );
  material_->allocPoints  ( elemCount_ * ipCount_  );

  // Select the correct function for computing the B-matrix.

  getShapeGrads_ = getShapeGradsFunc ( rank_ );
  getShapeFuncs_ = getShapeFunc      ( rank_ );

  // compute the equivalent element sizes

  eqvElemSize_.resize ( elemCount_ );
  eqvElemSize_ = 0.0;

  computeEqvElemSize_ ();
}


LocalDamageModel::~LocalDamageModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void LocalDamageModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void LocalDamageModel::getConfig ( const Properties& conf ) const
{
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool LocalDamageModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  // This action will be called when the solver needs
  // the tangent stiffness matrix and the internal force vector

  if ( action == Actions::GET_MATRIX0 )
  {
    //System::out() << myName_ << ": " << "computing stiffness matrix ... \n";

    Ref<MatrixBuilder>  mbuilder;

    Vector  disp;
    Vector  force;

    // Get the current displacements.

    StateVector::get ( disp, dofs_, globdat );

    // Get the matrix builder and the internal force vector.

    params.get ( mbuilder, ActionParams::MATRIX0    );
    params.get ( force,    ActionParams::INT_VECTOR );

    getMatrix_ ( *mbuilder, force, disp );

    return true;
  }

  // This action will be called when time integrator, say Newmark
  // needs the mass matrix

  if ( action == Actions::GET_MATRIX2 )
  {
    Ref<MatrixBuilder> mbuilder;

    params.get ( mbuilder, ActionParams::MATRIX2 );

    getMatrix20_( *mbuilder );

    return true;
  }

  // This action is called when solver gets convergence
  // Most often, it is time when material updates history variables

  if ( action == Actions::COMMIT )
  {
    material_->commit ();

    return true;
  }

  if ( action == "CHECK_COMMIT" )
  {
    checkCommit_ ( params );
    return true;
  }

  // This action is called by the visualisation, say FemViewModule
  // to plot the stress, strains etc

  if ( action == Actions::GET_TABLE )
  {
    return getTable_ ( params, globdat );
  }

  return false;
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void LocalDamageModel::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   disp )

{
  Cubix       grads      ( rank_, nodeCount_, ipCount_ );
  Matrix      coords     ( rank_, nodeCount_ );

  Matrix      elemMat    ( dofCount_, dofCount_ );

  Vector      elemForce  ( dofCount_ );
  Vector      elemDisp   ( dofCount_ );

  Matrix      stiff      ( strCount_, strCount_ );

  Vector      stress     ( strCount_ );
  Vector      strain     ( strCount_ );

  Matrix      b          ( strCount_, dofCount_ );
  Matrix      bt         = b.transpose ();

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );

  Vector      ipWeights  ( ipCount_   );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;

  int         ipoint = 0;

  double      wip;

  Matrix   stiff0 = material_->giveElasticMatrix ( ); 

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < elemCount_; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems_[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem  );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    // Compute the spatial derivatives of the element shapoe
    // functions.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Get the displacements at the element nodes.

    elemDisp = select ( disp, idofs );

    // Assemble the element stiffness matrix.

    elemMat   = 0.0;
    elemForce = 0.0;

    double  h = eqvElemSize_[ie];

    for ( int ip = 0; ip < ipCount_; ip++, ipoint++ )
    {
      // Compute the B-matrix for this integration point.

      getShapeGrads_ ( b, grads(ALL,ALL,ip) );

      // Compute the strain vector of this integration point

      matmul ( strain, b, elemDisp );

      // Get the tangent stiffness matrix and the stress vector
      // from the Material given the current strain

      material_->update ( stress, stiff, strain, ipoint, h );

      // Compute the element force vector.

      wip = ipWeights[ip];

      elemForce += wip * mc1.matmul ( bt, stress );

      // Compute the stiffness matrix

      elemMat   += wip * mc3.matmul ( bt, stiff, b );

      // add the correction part (tangent stiffness matrix)

      if ( material_->isLoading ( ipoint ) )
      {
        double  kappa       = material_->giveHistory     ( ipoint   );
        double  dOmega      = material_->getdOmegadKappa ( kappa, h );
        Vector  dEpsBardEps = material_->getdEpsBardEps  ( strain   );	

        elemMat -= dOmega * wip * mc3.matmul 
                        ( bt, mc2.matmul ( stiff0, matmul ( strain, dEpsBardEps ) ), b );	
      }

    }

    // Add the element matrix to the global stiffness matrix.

    mbuilder.addBlock ( idofs, idofs, elemMat );

    // Add the element force vector to the global force vector.

    select ( force, idofs ) += elemForce;
  }
}

//-----------------------------------------------------------------------
//   getMatrix20_
//-----------------------------------------------------------------------

// compute the consistent mass matrix 

void LocalDamageModel::getMatrix20_

    ( MatrixBuilder&          mbuilder )
{
  Matrix      coords     ( rank_, nodeCount_ );

  Matrix      elemMat    ( dofCount_, dofCount_ );

  Matrix      R          ( rank_, rank_ );

  Matrix      sfuncs     = shape_->getShapeFunctions ();
  Matrix      N          ( rank_, rank_ * nodeCount_ );
  Matrix      Nt         = N.transpose ( ); 

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );

  Vector      ipWeights  ( ipCount_   );

  MChain3     mc3;

  double      rho = material_->giveRho ( );

  R = 0.0;

  for ( int i = 0; i < rank_ ; i++ )
  {
    R(i,i) = rho;
  }

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < elemCount_; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems_[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem    );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getIntegrationWeights ( ipWeights, coords );

    // Assemble the element matrix and the internal force vector.

    elemMat   = 0.0;

    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      // compute matrix of shape function N       

      getShapeFuncs_ ( N, sfuncs(ALL,ip) );

      // Add the contribution of this integration point.

      elemMat   += ipWeights[ip] * mc3.matmul ( Nt, R, N );
    }

    // Add to the global mass matrix.

    mbuilder.addBlock ( idofs, idofs, elemMat );

  }
}

//-----------------------------------------------------------------------
//   getMatrix21_
//-----------------------------------------------------------------------

// compute the lumped mass matrix 

void LocalDamageModel::getMatrix21_

    ( MatrixBuilder&          mbuilder )
{
  Matrix      coords     ( rank_, nodeCount_ );

  Matrix      elemMat    ( dofCount_, dofCount_ );

  Matrix      R          ( rank_, rank_ );

  Matrix      sfuncs     = shape_->getShapeFunctions ();
  Matrix      N          ( rank_, rank_ * nodeCount_ );
  Matrix      Nt         = N.transpose ( ); 

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );

  Vector      ipWeights  ( ipCount_   );

  MChain3     mc3;

  double      rho = material_->giveRho ( );

  R = 0.0;
 
  for ( int i = 0; i < rank_ ; i++ )
  {
    R(i,i) = rho;
  }

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < elemCount_; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems_[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem    );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getIntegrationWeights ( ipWeights, coords );

    // Assemble the element matrix and the internal force vector.

    elemMat   = 0.0;

    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      // compute matrix of shape function N       

      getShapeFuncs_ ( N, sfuncs(ALL,ip) );

      // Add the contribution of this integration point.

      elemMat   += ipWeights[ip] * mc3.matmul ( Nt, R, N );
    }

    // Add to the global mass matrix.

    mbuilder.addBlock ( idofs, idofs, elemMat );

  }
}



//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void LocalDamageModel::checkCommit_

  ( const Properties&  params )

{
  params.set ( "accept", true );
}

//-----------------------------------------------------------------------
//   computeEqvElemSize_
//-----------------------------------------------------------------------

void LocalDamageModel::computeEqvElemSize_ ()
{
  System::out() << "Calculating equivalent element size ...\n";

  using jive::geom::Geometries;

  const String  shapeGeom  = shape_->getGeometry ();

  // check what type of element

 if ( shapeGeom == Geometries::TRIANGLE &&
	    nodeCount_ == 3 )
 {
    double     a(0), b(0), c(0);
    IdxVector  elem_nodes (3); 
    Matrix     co(rank_, 3);

    //helper

    double s(0);
    double F(-1);

    for ( int ie = 0; ie < elemCount_; ie++ )
    {
      elems_.getElemNodes  ( elem_nodes, ie );
      nodes_.getSomeCoords ( co, elem_nodes );

      a = sqrt(pow((co(0,1) - co(0,0)),2) + pow((co(1,1) - co(1,0)),2));
      b = sqrt(pow((co(0,2) - co(0,1)),2) + pow((co(1,2) - co(1,1)),2));
      c = sqrt(pow((co(0,0) - co(0,2)),2) + pow((co(1,0) - co(1,2)),2));


      s = 0.5 * (a + b + c);
      F = sqrt ( s * (s-a) * (s-b) * (s-c) );

      double he = ( a * b * c ) / ( 4.0 * F );

      JEM_ASSERT ( he > 1e-08 );

      eqvElemSize_[ie] = he;
    }
  }
  else if ( shapeGeom == Geometries::SQUARE &&
	    nodeCount_ == 4 )
  {

    double    a(0), b(0), c(0), d(0);
    IdxVector elem_nodes ( 4 ); 
    Matrix    co         ( rank_, 4 );

    // helper

    double s(0);
    double F(-1);


    for ( int ie = 0; ie < elemCount_; ie++ )
    {
      elems_.getElemNodes  ( elem_nodes, ie );
      nodes_.getSomeCoords ( co, elem_nodes );

      a = sqrt(pow((co(0,1) - co(0,0)),2) + pow((co(1,1) - co(1,0)),2));
      b = sqrt(pow((co(0,2) - co(0,1)),2) + pow((co(1,2) - co(1,1)),2));
      c = sqrt(pow((co(0,3) - co(0,2)),2) + pow((co(1,3) - co(1,2)),2));
      d = sqrt(pow((co(0,0) - co(0,3)),2) + pow((co(1,0) - co(1,3)),2));

      s = 0.5 * (a + b + c + d);
      F = sqrt ( (s-a) * (s-b) * (s-c) * (s-d) );

      eqvElemSize_[ie] = ( sqrt ( (a*b + c*d) * (a*c + b*d) * (a*d + b*c) ) ) /
                         ( 4.0 * F );
    }
  }
  else
  {
     System::warn() << "LocalDamageModel::computeEqvElemSize()," 
                    << "Equivalent element size is not implemented and set to 1.\n";

     eqvElemSize_ = 1.0;
  }

  System::out() << "Calculating equivalent element size ...done!\n";
}


//-----------------------------------------------------------------------
//   getTable_
//-----------------------------------------------------------------------


bool LocalDamageModel::getTable_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  Ref<XTable>  table;
  Vector       weights;
  String       name;


  // Get the table, the name of the table, and the table row weights
  // from the action parameters.

  params.get ( table,   ActionParams::TABLE );
  params.get ( name,    ActionParams::TABLE_NAME );
  params.get ( weights, ActionParams::TABLE_WEIGHTS );


  // Stress value are computed in the nodes.

  Vector  disp;

  StateVector::get ( disp, dofs_, globdat );

  if ( name == "stress" &&
       table->getRowItems() == nodes_.getData() )
  {
    getStress_ ( *table, weights, disp );

    return true;
  }

  // History values are computed in the nodes.

  if ( name == "damage" &&
       table->getRowItems() == nodes_.getData() )
  {
    getHistory_ ( *table, weights );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   getStress_
//-----------------------------------------------------------------------

void LocalDamageModel::getStress_

  ( XTable&        table,
    const Vector&  weights,
    const Vector&  disp )

{
  Matrix     sfuncs     = shape_->getShapeFunctions ();

  Cubix      grads      ( rank_, nodeCount_, ipCount_ );
  Matrix     b          ( strCount_, dofCount_        );
  Matrix     coords     ( rank_, nodeCount_           );

  Matrix     ndStress   ( nodeCount_, strCount_ );
  Vector     ndWeights  ( nodeCount_ );

  Vector     stressIp   ( strCount_ );
  Vector     strainIp   ( strCount_ );

  Vector     elemDisp   ( dofCount_  );
  Vector     ipWeights  ( ipCount_   );

  IdxVector  idofs      ( dofCount_  );
  IdxVector  inodes     ( nodeCount_ );
  IdxVector  jcols      ( strCount_  );

  MChain1    mc1;

  int        ipoint;

  // Add the columns for the stress components to the table.

  switch ( strCount_ )
  {
  case 1:

    jcols[0] = table.addColumn ( "s_xx" );

    break;

  case 3:

    jcols[0] = table.addColumn ( "s_xx" );
    jcols[1] = table.addColumn ( "s_yy" );
    jcols[2] = table.addColumn ( "s_xy" );

    break;

  case 6:

    jcols[0] = table.addColumn ( "s_xx" );
    jcols[1] = table.addColumn ( "s_yy" );
    jcols[2] = table.addColumn ( "s_zz" );
    jcols[3] = table.addColumn ( "s_xy" );
    jcols[4] = table.addColumn ( "s_yz" );
    jcols[5] = table.addColumn ( "s_xz" );

    break;

  default:

    throw Error (
      JEM_FUNC,
      "unexpected number of stress components: " +
      String ( strCount_ )
    );
  }

  Matrix stiff0  = material_->giveElasticMatrix ();
  Vector omega   = material_->giveOmega         ();

  // Iterate over all elements assigned to this model.

  ipoint = 0;

  for ( int ie = 0; ie < elemCount_; ie++ )
  {
    int  ielem = ielems_[ie];

    elems_.getElemNodes ( inodes, ielem );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    shape_->getShapeGradients ( grads, ipWeights, coords );

    elemDisp = select ( disp, idofs );

    // Iterate over the integration points.

    ndStress  = 0.0;
    ndWeights = 0.0;

    for ( int ip = 0; ip < ipCount_; ip++, ipoint++ )
    {
      getShapeGrads_ ( b, grads(ALL,ALL,ip) );

      matmul         ( strainIp, b, elemDisp );

      // Extrapolate the integration point stresses to the nodes using
      // the transposed shape functions.

      stressIp   = mc1.matmul ( stiff0, strainIp );
      stressIp  *= 1.0 - omega[ipoint];

      ndStress  += matmul ( sfuncs(ALL,ip), stressIp );
      ndWeights += sfuncs(ALL,ip);
    }

    select ( weights, inodes ) += ndWeights;

    // Add the stresses to the table.

    table.addBlock ( inodes, jcols, ndStress );
  }
}


//-----------------------------------------------------------------------
//   getHistory_
//-----------------------------------------------------------------------


void LocalDamageModel::getHistory_

  ( XTable&          table,
    const Vector&    weights )
{
  Matrix     coords     ( rank_, nodeCount_ );
  Vector     ipWeights  ( ipCount_ );

  IdxVector  inodes     ( nodeCount_ );
  IdxVector  jcols      ( 1 );

  Matrix     ndHistory  ( nodeCount_, 1 );
  Vector     ndWeights  ( nodeCount_    );

  Matrix     sfuncs     = shape_->getShapeFunctions(); 

  // Add a column for the history parameter to the table.

  jcols[0] = table.addColumn ( "hist" );

  // Get the history from the material

  Vector  hist = material_->giveOmega(  );

  // Iterate over all elements assigned to this model.

  int ipoint = 0;

  weights = 1.0;

  for ( int ie = 0; ie < elemCount_; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems_[ie];

    // Get the element coordinates.

    elems_.getElemNodes  ( inodes, ielem  );

    ndHistory = 0.0;
    ndWeights = 0.0;

    for ( int ip = 0; ip < ipCount_; ip++, ipoint++ )
    {
      ndHistory(ALL,0) += abs ( sfuncs(ALL,ip) ) * hist[ipoint] ;
      ndWeights        += abs ( sfuncs(ALL,ip) );
    }

    //select ( weights, inodes ) += ndWeights;

    table.addBlock ( inodes, jcols, ndHistory );
  }
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newDamageModel
//-----------------------------------------------------------------------


static Ref<Model>     newDamageModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<LocalDamageModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareDamageModel
//-----------------------------------------------------------------------


void declareLocalDamageModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "LocalDamage", & newDamageModel );
}
