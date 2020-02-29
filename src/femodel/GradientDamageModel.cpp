/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the gradient enhanced damage model
 *  using finite elements. 
 *  The number of elements and the integration scheme/element
 *  can be changed during the simulation.
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 */


#include <jem/base/Array.h>
#include <jem/base/Error.h>
#include <jem/base/System.h>
#include <jem/base/IllegalInputException.h>
#include <jem/io/PrintWriter.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/sparse/utilities.h>
#include <jem/util/Properties.h>
#include <jem/util/Flex.h>
#include <jive/util/utilities.h>
#include <jive/util/XTable.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/Printer.h>
#include <jive/util/Assignable.h>
#include <jive/util/Constraints.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/geom/Geometries.h>
#include <jive/geom/InternalShape.h>
#include <jive/geom/IShapeFactory.h>
#include <jive/geom/StdSquare.h>
#include <jive/geom/StdCube.h>
#include <jive/geom/IShapeFactory.h>
#include <jive/geom/StdShapeFactory.h>
#include <jive/geom/StdBezierShape.h>
#include <jive/geom/ParametricArea.h>
#include <jive/geom/ParametricVolume.h>
#include <jive/geom/Geometries.h>
#include <jive/fem/ElementGroup.h>

#include "FemModels.h"
#include "util/utilities.h"
#include "material/DamageMaterial.h"
#include "material/Material.h"
#include "util/TbFiller.h"

using namespace jem;

using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jem::util::SparseArray;
using jem::util::Properties;
using jem::util::Flex;
using jive::Vector;
using jive::IdxVector;
using jive::IntMatrix;
using jive::Matrix;
using jive::Cubix;
using jive::util::XTable;
using jive::util::XDofSpace;
using jive::util::Assignable;
using jive::util::Constraints;
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
typedef MatmulChain<double,4>  MChain4;

typedef ElementSet             ElemSet;
typedef ElementGroup           ElemGroup;


//=======================================================================
//   class GradientModel
//=======================================================================


class GradientDamageModel : public Model
{
 public:

  typedef GradientDamageModel     Self;
  typedef Model             Super;

  static const char*        DISP_NAMES[3];

  static const char*        SHAPE_PROP;
  static const char*        MATERIAL_PROP;


                            GradientDamageModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       conf,
      const Properties&       globdat )  const;

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );


 protected:

  virtual                  ~GradientDamageModel  ();


 private:

  void                      getMatrix_

    ( MatrixBuilder&          mbuilder,
      const Vector&           force,
      const Vector&           state );

  void                      getMatrix2_

    ( MatrixBuilder&          mbuilder );

  void                      setConstraints_

    ( Constraints&            cons );

  bool                      getTable_

    ( const Properties&       params,
      const Properties&       globdat );

  /**
   * @brief      This function is called by getTable_.
   *
   * @param[out] table     The table
   * @param[in]  weights   The weights
   * @param[in]  contents  Filter for which data to be written to table
   * @param[in]  state     The nodal value vector
   */

  void                      getOutputData_

    ( Ref<XTable>             table,
      const Vector&           weights,
      const String&           contents,
      const Vector&           state );  


  /* params.set ( "accept", false ) => keep same load and resolve
   * params.set ( "accept", true )  => advance to new load step (normal case)
   */

  void                      checkCommit_

    ( const Properties&       params );

  /* Initialize the mapping between integration points
   * and material points 
   */

  void                      initializeIPMPMap_ ();


 private:

  Assignable<ElemGroup>     egroup_;
  Assignable<ElemSet>       elems_;
  Assignable<NodeSet>       nodes_;

  int                       rank_;

  Ref<IShape>               shape_;

  Ref<XDofSpace>            dofs_;
  IdxVector                 dofTypes_;
  IdxVector                 eqvTypes_;
  IdxVector                 dispTypes_;

  Matrix                    strain_;

  Ref<DamageMaterial>       material_;

  ShapeGradsFunc            getShapeGrads_;
  ShapeFunc                 getShapeFuncs_;

  /* 
   *  Mapping between integration points and material points  
   */

  SparseArray <int, 2>      ipMpMap_;

  /* 
   * Integer vector holds 1 or zero to indicate if element ID is 
   * active or not. Used to remove fully damaged element.
   */

  IdxVector                 isActive_;
};


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  GradientDamageModel::DISP_NAMES[3] = { "dx", "dy", "dz" };

const char*  GradientDamageModel::SHAPE_PROP    = "shape";
const char*  GradientDamageModel::MATERIAL_PROP = "material";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


GradientDamageModel::GradientDamageModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super ( name )

{
  using jive::util::joinNames;
  using jive::geom::IShapeFactory;

  Properties    myProps = props.getProps  ( myName_ );
  Properties    myConf  = conf .makeProps ( myName_ );

  const String  context = getContext ();


  // Get the element group assigned to this model.

  egroup_ = ElemGroup::get ( myConf, myProps, globdat, context );

  elems_  = egroup_.getElements ();
  nodes_  = elems_ .getNodes    ();
  rank_   = nodes_ .rank        ();

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

  // Make sure that each element has the same number of nodes as the
  // shape object.

  elems_.checkSomeElements (
    context,
    egroup_.getIndices (),
    shape_->nodeCount  ()
  );

  // Get the degree of freedom (DOF) space, and add the DOFs of this
  // model. Note that the first DOF types are the displacement DOFs,
  // and the last DOF type is the non-local equivalent strain DOF.

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  dofTypes_.resize ( rank_ + 1 );

  // Make a shallow copy of the first part of the dofTypes_ array.

  dispTypes_.ref ( dofTypes_[slice(BEGIN,rank_)] );

  for ( int i = 0; i < rank_; i++ )
  {
    dispTypes_[i] = dofs_->addType ( DISP_NAMES[i] );
  }

  // Make a shallow copy of the last part of the dofTypes_ array.

  eqvTypes_.ref ( dofTypes_[slice(rank_,END)] );

  eqvTypes_[0] = dofs_->addType ( "eqv" );

    
  dofs_->addDofs (
    elems_.getUniqueNodesOf ( egroup_.getIndices() ),
    dofTypes_
  );

  // Compute the total number of integration points.

  int  ipCount = shape_->ipointCount() * egroup_.size();

  // Allocate memory for the strains, including the non-local
  // equivalent strain.

  strain_.resize ( STRAIN_COUNTS[rank_] + 1, ipCount );

  strain_ = 0.0;

  // Create a material model object.
  // Material of gradien enhanced damage model MUST be a Damage material

  Ref<Material> tmpMat = newMaterial ( MATERIAL_PROP, myConf, myProps, globdat );
  material_            = dynamicCast<DamageMaterial> ( tmpMat );

  if ( material_ == NIL )
  {
    throw IllegalInputException (
      context,
      "material type is not a damage material"
    );
  }

  material_->allocPoints  ( ipCount );

  initializeIPMPMap_ ();

  // Select the correct function for computing the B-matrix.

  getShapeGrads_ = getShapeGradsFunc ( rank_ );

  IdxVector   ielems     = egroup_.getIndices ();
  const int   ielemCount = ielems.size         ();

  isActive_.resize ( ielemCount );

  isActive_ = 1;
 
}


GradientDamageModel::~GradientDamageModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void GradientDamageModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps  = props  .findProps ( myName_ );
  Properties  matProps = myProps.findProps ( MATERIAL_PROP );

  material_->configure ( matProps, globdat );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void GradientDamageModel::getConfig ( const Properties& conf,
                                      const Properties&  globdat )  const
{
  Properties  myConf  = conf  .makeProps ( myName_ );
  Properties  matConf = myConf.makeProps ( MATERIAL_PROP );

  material_->getConfig ( matConf, globdat );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool GradientDamageModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  // compute the tangent stiffness matrix and internal force vector

  if ( action == Actions::GET_MATRIX0 )
  {
    Ref<MatrixBuilder>  mbuilder;

    Vector  state;
    Vector  force;

    // Get the current displacements.

    StateVector::get ( state, dofs_, globdat );

    // Get the matrix builder and the internal force vector.

    params.get ( mbuilder, ActionParams::MATRIX0 );
    params.get ( force,    ActionParams::INT_VECTOR );

    getMatrix_ ( *mbuilder, force, state );

    return true;
  }

  // compute mass matrix 

  if ( action == Actions::GET_MATRIX2 )
  {
    Ref<MatrixBuilder> mbuilder;

    params.get ( mbuilder, ActionParams::MATRIX2 );
    
    getMatrix2_( *mbuilder );

    return true;
  }

  // initialisation: in this model, do constrain the eqv dofs
  // to make sure interpolation order of eqv = order of disp dofs -1

  if ( action == Actions::INIT )
  {
    // Get the constraints associated with the DOF space.

    Ref<Constraints>  cons = Constraints::get ( dofs_, globdat );

    setConstraints_ ( *cons );

    return true;
  }

  // things to be done when the solution converges, ie the end of 
  // a load step

  if ( action == Actions::COMMIT )
  {
    material_->commit ();

    return true;
  }

  // used with the StepModule
  // to advance to new load step or resolve with the same load vector

  if ( action == Actions::CHECK_COMMIT )
  {
    //checkCommit_ ( params );
    return true;
  }

  // post processing: compute some table for visualize stress, 
  // strain, damage ...

  if ( action == Actions::GET_TABLE )
  {
    return getTable_ ( params, globdat );
  }

  return false;
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void GradientDamageModel::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   state )

{
  IdxVector   ielems     = egroup_.getIndices ();

  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();
  const int   strCount   = STRAIN_COUNTS[rank_];

  // number of nonlocal equivalent strain dofs (scalar field)
  
  const int   eqvCount   = nodeCount;

  // number of displacement dofs
  
  const int   dispCount  = nodeCount * rank_;

  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );

  // current vector state:
  // non-local eqv strain and displacement
  
  Vector      eqv        ( eqvCount  );
  Vector      disp       ( dispCount );

  Matrix      stiff      ( strCount, strCount );
  Vector      stress     ( strCount );
  Vector      strain     ( strCount );

  // internal force vector
  // one for disp. dofs and one for equi. dofs

  Vector      elemForce1 ( dispCount );
  Vector      elemForce2 ( eqvCount  );

  // element stiffness matrices (four components)

  Matrix      elemMat1   ( dispCount, dispCount ); // disp-disp
  Matrix      elemMat2   ( dispCount, eqvCount  ); // disp-eqv 
  Matrix      elemMat3   ( eqvCount,  dispCount ); // eqv -disp
  Matrix      elemMat4   ( eqvCount,  eqvCount  ); // eqv -eqv
  
  // B matrices
  // one for disp. dofs and one for equi strain dofs
  
  Matrix      bd         ( strCount, dispCount );
  Matrix      bdt        = bd.transpose ();

  Matrix      be         ( rank_, eqvCount );
  Matrix      bet        = be.transpose ();
  
  IdxVector   inodes     ( nodeCount );
  IdxVector   eqvDofs    ( eqvCount  );
  IdxVector   dispDofs   ( dispCount );

  Vector      ipWeights  ( ipCount   );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;
  MChain4     mc4; 
  
  double      wip;
  double      equiStrain;

  double      c = material_->giveLengthScale ();

  // get the elastic stiffness matrix and the history
  // of the material for one time !!!
  // Will be used in case of damage grows.

  Matrix       stiff0 = material_->giveElasticMatrix ( );

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {

    if ( ! isActive_[ie] )
    {      
      continue;
    }

    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes,   ielem  );
    nodes_.getSomeCoords ( coords,   inodes );
    dofs_->getDofIndices ( eqvDofs,  inodes, eqvTypes_ );
    dofs_->getDofIndices ( dispDofs, inodes, dispTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Compute the shape functions and transpose

    Matrix N  = shape_->getShapeFunctions ();
    
    // Get the non-local equivalent strains and the displacements in
    // the element nodes.

    eqv  = select ( state, eqvDofs );
    disp = select ( state, dispDofs );

    // Initialization the internal forces
    // and the tangent stiffness matrix
    
    elemForce1 = 0.0;
    elemForce2 = 0.0;

    elemMat1   = 0.0;
    elemMat2   = 0.0;
    elemMat3   = 0.0;
    elemMat4   = 0.0;    

    // Loop on integration points 
    
    for ( int ip = 0; ip < ipCount; ip++ )
    {
      // get the corresponding material point ipoint
      // of element ielem integration point ip from the mapping

      int ipoint = ipMpMap_ ( ielem, ip );

      // Compute the B-matrix for this integration point.
      // it is the B-matrix of displacement dofs

      getShapeGrads_ ( bd, grads(ALL,ALL,ip) );

      // get B-matrix associated with eqv dofs
      
      be =  grads(ALL,ALL,ip);

      // Compute the strain for this integration point.

      matmul ( strain, bd, disp );

      // Store the regular strain components.

      strain_(slice(BEGIN,strCount),ipoint) = strain;

      // Compute the non-local equivalent strain at this GP
      // and store in the last position of strain_ 

      Vector Nip               = N(ALL,ip);
      
      strain_(strCount,ipoint) = dot ( Nip, eqv );

      // Update the material model.

      material_->update ( stress, stiff, strain_(ALL,ipoint), ipoint );

      // ----------------------------------------------
      //   Compute components of the stiffness matrix
      // ----------------------------------------------
      
      // the disp-disp stiffness, eqv-eqv andd eqv-disp always here

      wip         = ipWeights[ip];
      elemMat1   += wip * mc3.matmul ( bdt, stiff, bd );
      elemMat4   += wip * ( matmul ( Nip, Nip ) + c * mc2. matmul ( bet, be )   );
     
      Vector dEpsBardEps = material_->getdEpsBardEps ( strain );
      
      elemMat3   -=  wip * mc2.matmul ( matmul( Nip, dEpsBardEps ), bd );


      // disp-eqv stiffness matrix, only for loading case      

      if (  material_->isLoading ( ipoint ) )
      {
        double    kappa  = material_->giveHistory     ( ipoint );
	    double    dOmega = material_->getdOmegadKappa ( kappa  );
	  
	    elemMat2 -=  wip * dOmega *  mc3.matmul ( bdt, stiff0, 
                                                  matmul( strain, Nip ) );
      }
      
      // Compute internal forces

      equiStrain  = material_-> getEquiStrain ( strain );

      elemForce1 +=  wip * mc1.matmul ( bdt, stress );
      elemForce2 -=  wip * equiStrain * Nip;

    }  // end of loop on integration points

    elemForce2   +=  mc1.matmul ( elemMat4, eqv );

    // Assembly ...

    mbuilder.addBlock ( dispDofs, dispDofs, elemMat1 );
    mbuilder.addBlock ( dispDofs, eqvDofs , elemMat2 );
    mbuilder.addBlock ( eqvDofs , dispDofs, elemMat3 );
    mbuilder.addBlock ( eqvDofs , eqvDofs , elemMat4 );

    select ( force, dispDofs ) += elemForce1;
    select ( force, eqvDofs  ) += elemForce2;
  }
}


//-----------------------------------------------------------------------
//   getMatrix2_
//-----------------------------------------------------------------------

// compute the mass matrix 
// current implementation: consistent mass matrix

void GradientDamageModel::getMatrix2_

    ( MatrixBuilder&          mbuilder )
{
  IdxVector   ielems     = egroup_.getIndices  ();

  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();
  const int   dofCount   = rank_ * nodeCount;

  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );
  
  Matrix      elemMat    ( dofCount, dofCount );

  Matrix      R          ( rank_, rank_ );

  Matrix      sfuncs     = shape_->getShapeFunctions ();
  Matrix      N          ( rank_, rank_ * nodeCount );
  Matrix      Nt         = N.transpose ( ); 

  IdxVector   inodes     ( nodeCount );
  IdxVector   idofs      ( dofCount  );

  Vector      ipWeights  ( ipCount   );

  MChain3     mc3;

  double      rho = material_->giveRho ( );

  R = 0.0;
 
  for ( int i = 0; i < rank_ ; i++ )
  {
    R(i,i) = rho;
  }

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem    );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Assemble the element matrix and the internal force vector.

    elemMat   = 0.0;

    for ( int ip = 0; ip < ipCount; ip++ )
    {
      // compute matrix of shpae function N       

      getShapeFuncs_ ( N, sfuncs(ALL,ip) );

      // Add the contribution of this integration point.

      elemMat   += ipWeights[ip] * mc3.matmul ( Nt, R, N );
    }

    // Add the element secant matrix to the global stiffness matrix.

    mbuilder.addBlock ( idofs, idofs, elemMat );

  }
}

//-----------------------------------------------------------------------
//   setConstraints_
//-----------------------------------------------------------------------


void GradientDamageModel::setConstraints_ ( Constraints& cons )
{
  using jive::util::Printer;
  using jive::geom::Geometries;

  IdxVector     ielems     = egroup_.getIndices ();

  const int     ielemCount = ielems   .size     ();
  const int     nodeCount  = shape_->nodeCount  ();

  const String  shapeGeom  = shape_->getGeometry ();

  IdxVector     inodes     ( nodeCount );
  IdxVector     idofs      ( nodeCount );
  IdxVector     jdofs      ( 2 );
  Vector        coeffs     ( 2 );

  IntMatrix     cdofs;


  // only do this constraint for high order elements
  // not for linear elements

  if      ( shapeGeom == Geometries::LINE &&
	    nodeCount == 2 )
  {
    return;
  }
  else if ( shapeGeom == Geometries::TRIANGLE &&
	    nodeCount == 3 )
  {
    return;
  }
  else if ( shapeGeom == Geometries::SQUARE &&
	    nodeCount == 4 )
  {
    return;
  }
  else if ( shapeGeom == Geometries::TETRAHEDRON &&
	    nodeCount == 4 )
  {
    return;
  }

  System::info() << myName_
		 << " : initializing constraints ...\n";

  // Determine which equivalent strain DOFS are to be constrained. The
  // first row of the cdofs matrix contains the DOF indices to be
  // constrained, and the other two rows contain the DOF indices of
  // the two master nodes. Note that the DOF indices are local with
  // respect to an element. They will be translated to global DOF
  // indices below.

  if      ( shapeGeom == Geometries::LINE &&
	    nodeCount == 3 )
  {
    cdofs.resize ( 3, 1 );

    // Constrain the mid node.

    cdofs(0,0) = 1;
    cdofs(1,0) = 0;
    cdofs(2,0) = 2;
  }
  else if ( shapeGeom == Geometries::TRIANGLE &&
	    nodeCount == 6 )
  {
    cdofs.resize ( 3, 3 );

    // Constrain the mid nodes on each edge.

    cdofs(0,0) = 1;
    cdofs(1,0) = 0;
    cdofs(2,0) = 2;

    cdofs(0,1) = 3;
    cdofs(1,1) = 2;
    cdofs(2,1) = 4;

    cdofs(0,2) = 5;
    cdofs(1,2) = 4;
    cdofs(2,2) = 0;
  }
  else if ( shapeGeom == Geometries::SQUARE &&
	    nodeCount == 8 )
  {
    cdofs.resize ( 3, 4 );

    // Constrain the mid nodes on each edge.

    cdofs(0,0) = 1;
    cdofs(1,0) = 0;
    cdofs(2,0) = 2;

    cdofs(0,1) = 3;
    cdofs(1,1) = 2;
    cdofs(2,1) = 4;

    cdofs(0,2) = 5;
    cdofs(1,2) = 4;
    cdofs(2,2) = 6;

    cdofs(0,3) = 7;
    cdofs(1,3) = 6;
    cdofs(2,3) = 0;
  }
  else
  {
    throw IllegalInputException (
      getContext (),
      String::format (
	"unsupported shape geometry: %S%d",
	&shapeGeom,
	nodeCount
      )
    );
  } 

  // Determine the constraint coefficients.

  coeffs = 0.5;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element DOFs.

    elems_.getElemNodes  ( inodes, ielem );
    dofs_->getDofIndices ( idofs,  inodes, eqvTypes_ );

    // Add constraints.

    for ( int j = 0; j < cdofs.size(1); j++ )
    {
      // Get the global index of the DOF to be constrained.

      int  idof = idofs[cdofs(0,j)];

      // Check whether this DOF has already been constrained.

      if ( cons.isSlaveDof( idof ) )
      {
	continue;
      }

      // Get the global indices of the two master DOFs.

      jdofs[0] = idofs[cdofs(1,j)];
      jdofs[1] = idofs[cdofs(2,j)];

      // Add a linear constraint of the form:
      //
      //   u[idof] = u[jdofs[0]] * coeffs[0] +
      //             u[jdofs[1]] * coeffs[1] +

      cons.addConstraint ( idof, jdofs, coeffs );
    }
  }

  // This is for debugging only.

  //cons.printTo ( Printer::get() );

  //Printer::get() << "\n\n";

  //Printer::flush ();
}

//-----------------------------------------------------------------------
//   initializeIPMPMap_
//-----------------------------------------------------------------------


void  GradientDamageModel::initializeIPMPMap_ ( )

{
  jive::IdxVector   ielems     = egroup_.getIndices  ();

  const idx_t   ielemCount = ielems.size         ();
  const idx_t   ipCount    = shape_->ipointCount ();

        idx_t   ipoint     = 0;

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems[ie];

    // loop over integration points 

    for ( int ip = 0; ip < ipCount; ip++, ipoint++ )
    {
      ipMpMap_ ( ielem, ip ) = ipoint;
    }
  }
}


//-----------------------------------------------------------------------
//   getTable_
//-----------------------------------------------------------------------


bool GradientDamageModel::getTable_

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

  if ( name == "xoutTable" ) //&& table->getRowItems() == nodes_.getData() )
  {
    Vector  disp;
    String  contents;

    StateVector::get ( disp,     dofs_, globdat  );
    params.      get ( contents, "contentString" );

    
    getOutputData_ ( table, weights, contents, disp );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   getStrain_
//-----------------------------------------------------------------------


// void GradientDamageModel::getStrain_

//   ( XTable&        table,
//     const Vector&  weights )

// {
//   IdxVector  ielems     = egroup_.getIndices  ();

//   const int  ielemCount = ielems.size         ();
//   const int  nodeCount  = shape_->nodeCount   ();
//   const int  ipCount    = shape_->ipointCount ();
//   const int  strCount   = STRAIN_COUNTS[rank_];

//   Matrix     sfuncs     = shape_->getShapeFunctions ();

//   Matrix     ndStrain   ( nodeCount, strCount + 1 );
//   Vector     ndWeights  ( nodeCount );

//   IdxVector  inodes     ( nodeCount );
//   IdxVector  jcols      ( strCount + 1 );

//   int        ipoint;


//   // Add the columns for the normal strain components to the table.

//   switch ( strCount )
//   {
//   case 1:

//     jcols[0] = table.addColumn ( "e_xx" );

//     break;

//   case 3:

//     jcols[0] = table.addColumn ( "e_xx" );
//     jcols[1] = table.addColumn ( "e_yy" );
//     jcols[2] = table.addColumn ( "e_xy" );

//     break;

//   case 6:

//     jcols[0] = table.addColumn ( "e_xx" );
//     jcols[1] = table.addColumn ( "e_yy" );
//     jcols[2] = table.addColumn ( "e_zz" );
//     jcols[3] = table.addColumn ( "e_xy" );
//     jcols[4] = table.addColumn ( "e_yz" );
//     jcols[5] = table.addColumn ( "e_xz" );

//     break;

//   default:

//     throw Error (
//       JEM_FUNC,
//       "unexpected number of strain components: " +
//       String ( strCount )
//     );
//   }

//   // Add an extra column for the non-local equivalent strain.

//   jcols[strCount] = table.addColumn ( "e_eq" );

//   // Iterate over all elements assigned to this model.

//   ipoint = 0;

//   for ( int ie = 0; ie < ielemCount; ie++ )
//   {
//     // Get the global element index.

//     int  ielem = ielems[ie];

//     // Get the element nodes.

//     elems_.getElemNodes ( inodes, ielem );

//     // Iterate over the integration points.

//     ndStrain  = 0.0;
//     ndWeights = 0.0;

//     for ( int ip = 0; ip < ipCount; ip++, ipoint++ )
//     {
//       // Extrapolate the integration point strains to the nodes using
//       // the transposed shape functions.

//       ndStrain  += matmul ( sfuncs(ALL,ip), strain_(ALL,ipoint) );
//       ndWeights += sfuncs(ALL,ip);
//     }

//     // Increment the table weights. When the complete table has been
//     // filled, Jive will divide each row in the table by the
//     // corresponding table weight. In this way the strain components
//     // are automatically averaged over all elements that are attached
//     // to a node. The weight vector is initially zero.

//     select ( weights, inodes ) += ndWeights;

//     // Add the strains to the table.

//     table.addBlock ( inodes, jcols, ndStrain );
//   }
// }

//-----------------------------------------------------------------------
//   getHistory_
//-----------------------------------------------------------------------


// void GradientDamageModel::getHistory_

//   ( XTable&          table,
//     const Vector&    weights )
// {
//   IdxVector  ielems     = egroup_.getIndices  ();

//   const int  ielemCount = ielems.size         ();
//   const int  nodeCount  = shape_->nodeCount   ();
//   const int  ipCount    = shape_->ipointCount ();

//   Matrix     coords     ( rank_, nodeCount );
//   Vector     ipWeights  ( ipCount );

//   IdxVector  inodes     ( nodeCount );
//   IdxVector  jcols      ( 1 );

//   Matrix     ndHistory  ( nodeCount, 1 );
//   Vector     ndWeights  ( nodeCount    );

//   Matrix     sfuncs     = shape_->getShapeFunctions();

//   // Add a column for the history parameter to the table.

//   jcols[0] = table.addColumn ( "hist" );

//   // Get the history from the material

//   Vector  hist = material_->giveOmega(  );

//   // Iterate over all elements assigned to this model.

//   int ipoint = 0;

//   for ( int ie = 0; ie < ielemCount; ie++ )
//   {
//     // Get the global element index.

//     int  ielem = ielems[ie];

//     // Get the element coordinates.

//     elems_.getElemNodes  ( inodes, ielem  );

//     ndHistory = 0.0;
//     ndWeights = 0.0;

//     for ( int ip = 0; ip < ipCount; ip++, ipoint++ )
//     {
//       ndHistory(ALL,0) += abs ( sfuncs(ALL,ip) ) * hist[ipoint] ;
//       ndWeights +=  abs ( sfuncs(ALL,ip) );
//     }

//     select ( weights, inodes ) += ndWeights;
    
//     table.addBlock ( inodes, jcols, ndHistory );
//   }
// }

//-----------------------------------------------------------------------
//   getOutputData_
//-----------------------------------------------------------------------

// all data at Gauss points are calculated but only those specified in
// the vtk data block in the pro file will be written to vtu files.
// For example  vtk.data = "pstress | damage | stress_xx" -> write all 3 principal stresses,
// damage and sigma_xx.
// For history data, refer to historyNames_ in phase-field material classes to know which data
// are available.

void GradientDamageModel::getOutputData_

  ( Ref<XTable>    table,
    const Vector&  weights,
    const String&  contents,
    const Vector&  state )
          
{
  MChain2    mc2;

  // table filler related stuff

  TbFiller   tbFiller   ( rank_ );

  //StringVector hisNames =  materials_[0]->getHistoryNames ();

  Slice      iistrain   = tbFiller.announce ( "strain.tensor" ); // iistrain  = [0 1 2 3]
  Slice      iistress   = tbFiller.announce ( "stress.tensor" ); // iistrain  = [0 1 2 3]
  //Slice      iipstress  = tbFiller.announce ( "pstress.diag"  ); // iipstress = [4 5 6]
  Slice      iidamage   = tbFiller.announce ( "damage"        ); // [7]
  //Slice      iihistory  = tbFiller.announce ( hisNames        ); // material-dependent history variables

  Vector     ipValues   ( tbFiller.typeCount() ); // typeCount() = # of types = 8 in 2D

  Vector     strain     ( ipValues[iistrain]   );
  Vector     stress     ( ipValues[iistress]   );
  //Vector     pstress    ( ipValues[iipstress]  );
  Vector     dam        ( ipValues[iidamage]   );
  //Vector     history    ( ipValues[iihistory]  );


  // Let TbFiller find out which columns of ndValues to write to 
  // which columns of the table (based on filter in input file)

  IdxVector  i2table;
  IdxVector  jcols;

  tbFiller . setFilter   ( contents ); 
  tbFiller . prepareTable( i2table, jcols, table ); 
  // System::out() << strain <<"HAHHAHAHAH\n";
  // System::out() << stress <<"hehe\n";
  // System::out() << pstress <<"hehe\n";
  // 
  
  IdxVector  ielems     = egroup_.getIndices  ();

  const int  ielemCount  = ielems.size         ();
  const int  nodeCount   = shape_->nodeCount   ();
  const int  ipCount     = shape_->ipointCount ();
  const int  strCount    = STRAIN_COUNTS[rank_];
  const int  dispCount   = nodeCount * rank_;
   
  Matrix     ndValuesOut  ( nodeCount, i2table.size() );
  Matrix     ndValuesOut1 ( nodeCount, i2table.size() );
  Vector     ipValuesOut  ( i2table.size() );

  Matrix      stiff      ( strCount, strCount  );
  Matrix      bd         ( strCount, dispCount );
  Vector      ndWeights  ( nodeCount           ); 
  IdxVector   inodes     ( nodeCount           );


  Properties  params;

  // get all damage 

  Vector  hist = material_->giveOmega(  );

  idx_t       ipoint, igpoint = 0;

  Matrix      sfuncs     = shape_->getShapeFunctions ();

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    int  ielem = ielems[ie];

    elems_.getElemNodes ( inodes, ielem );

    ndValuesOut = 0.;
    ndWeights   = 0.;

    // Loop on integration points 
    
    for ( idx_t ip = 0; ip < ipCount; ip++, igpoint++ )
    {
      ipoint      = ipMpMap_( ielem, ip );
      
 
      strain      = strain_(slice(BEGIN,strCount),igpoint);
      dam[0]      = hist[igpoint]; // non-local dam

      material_->update ( stress, stiff, strain , ipoint );
 
      //funcs = abs ( funcs );

       // apply the filter now, only cols specified by i2table are 
       // written to the table 
      ipValuesOut  = ipValues[i2table]; 



      matmul (ndValuesOut1, sfuncs(ALL,ip), ipValuesOut );


      ndValuesOut += ndValuesOut1;
      ndWeights   += sfuncs(ALL,ip);
       //System::out() << " 2" << noda << "\n";
    }  // end of loop on integration points
    
    select ( weights, inodes ) += ndWeights;

    // Add the stresses to the table.

    table->addBlock ( inodes, jcols, ndValuesOut );
  }
}

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void GradientDamageModel::checkCommit_

  ( const Properties&  params )

{
  System::info() << myName_ << " : check commit ...\n";

  IdxVector   ielems     = egroup_.getIndices ();

  const int   ielemCount = ielems.size         ();
  const int   ipCount    = shape_->ipointCount ();

  bool        restart    = false; 

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // do not check element already disabled

    if ( ! isActive_[ie] )
    {
      continue;
    }

    int ielem = ielems[ie];

    for ( int ip = 0; ip < ipCount; ip++ )
    {
      int ipoint = ipMpMap_ (ielem, ip);

      if ( material_->isFullyDamaged ( ipoint ) )
      {
        isActive_[ie] = 0;
        restart       = true;

        break;
      }
    }
 
    // only allow one removed element per load step (can used with  
    // large load increment)

    // or allow many elements be removed per load step (then the load 
    // increment must be small): comment the next statement

    if ( restart ) break;
  }

  // if there is an element fully damaged, then
  // restart the same load step, ie do not accept the current
  // solution 

  if ( restart )
  {
    params.set ( "accept", false );
  }
 
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newGradientDamageModel
//-----------------------------------------------------------------------


static Ref<Model>     newGradientDamageModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<GradientDamageModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareGradientDamageModel
//-----------------------------------------------------------------------


void declareGradientDamageModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "GradientDamage", & newGradientDamageModel );
}
