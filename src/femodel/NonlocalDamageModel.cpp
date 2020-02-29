/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the nonlocal integral damage model
 *  using finite elements.
 *
 *  This models applies for element group. This design allows
 *  for example, nonlocal damage model and elastic model exist
 *  at the same time. The element group can contain subgroups
 *  for different materials. 
 *
 *  For efficiency, only mesh of one element type is allowed.
 *
 *  Remark: In the input file, must use
 *  matrix.type = "Sparse"; for non-standard node by node assembly
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
#include <jive/geom/GridBoxManager.h>
#include <jive/geom/boxes.h>
#include <jive/fem/ElementGroup.h>

#include "NonlocalDamageModel.h"
#include "FemModels.h"
#include "util/TbFiller.h"

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  NonlocalDamageModel::DOF_NAMES[3]     = { "dx", "dy", "dz" };

const char*  NonlocalDamageModel::SHAPE_PROP       = "shape";
const char*  NonlocalDamageModel::MATERIAL_PROP    = "materials";
const char*  NonlocalDamageModel::WEIGHT_FUNC_PROP = "weightFunc";
const char*  NonlocalDamageModel::RADIUS_PROP      = "radius";
const char*  NonlocalDamageModel::AVERAGING_PROP   = "averaging";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


NonlocalDamageModel::NonlocalDamageModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super ( name )

{
  using jive::util::joinNames;
  using jive::util::FuncUtils;
  using jive::geom::IShapeFactory;
  using jive::StringVector;

  Properties    myProps = props.getProps  ( myName_ );
  Properties    myConf  = conf .makeProps ( myName_ );

  const String  context = getContext ();

  // Get the element group assigned to this model.

  egroup_ = ElemGroup::get ( myConf, myProps, globdat, context );

  // then get elements and nodes from this group

  elems_  = egroup_.getElements ();
  nodes_  = elems_ .getNodes    ();

  numElem_= egroup_.size        ();

  ielems_.resize ( numElem_ );
  ielems_ = egroup_.getIndices  ();

  // compute some constant quantities

  rank_      = nodes_ .rank            ( );
  nodeCount_ = elems_.maxElemNodeCount ( );
  dofCount_  = nodeCount_ * rank_;
  strCount_  = STRAIN_COUNTS[rank_];

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

  // ----------------------------------------------------
  //     Create an internal shape object 
  // ----------------------------------------------------

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

  // ----------------------------------------------------
  //     create degree of freedom space
  // ----------------------------------------------------

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  dofTypes_.resize ( rank_ );

  for ( int i = 0; i < rank_; i++ )
  {
    dofTypes_[i] = dofs_->addType ( DOF_NAMES[i] );
  }
  
  System::out() << dofTypes_;

  dofs_->addDofs ( elems_.getUniqueNodesOf ( ielems_ ),
                   dofTypes_ );

  // Compute the total number of integration points.

  ipCount_     = shape_->ipointCount ();

  int  totalIpCount = numElem_ * ipCount_;

  // Allocate memory for the strains, including the non-local
  // equivalent strain.

  strain_.resize ( strCount_ + 1, totalIpCount );

  //strain_ = 0.0;

  // ------------------------------------------------------
  //       create material objects
  // ------------------------------------------------------

  StringVector mats;

  myProps.get ( mats ,MATERIAL_PROP );
  myConf. set ( MATERIAL_PROP, mats );

  int numMats = mats.size ( );

  materials_. resize  ( numMats  );
  //elemMatMap_.resize  ( numElem_ );

  for ( int iMat = 0; iMat < numMats; iMat++ )
  {
    String matName                = mats[iMat]; 

    Ref<Material>         tmpMat  = newMaterial ( matName, myConf, myProps, globdat );
    Ref<DamageMaterial>   damMat  = dynamicCast<DamageMaterial> ( tmpMat );

    if ( damMat == NIL )
    {
      throw IllegalInputException (
        context,
        "material type is not a damage material"
      );
    }

    Properties  matProps = myProps.findProps ( matName );
    Properties  matConf  = myConf.makeProps  ( matName );

    damMat->configure ( matProps, globdat );
    damMat->getConfig ( matConf , globdat );

    // get the elements of this material

    Properties    imatProps = myProps.getProps ( matName );
    Properties    imatConf  = myConf.makeProps ( matName );

    ElemGroup egroup    = ElemGroup::get ( imatConf, imatProps, globdat, context );

    IdxVector eIndices  = egroup.getIndices ( ); 

    int elemCount = eIndices.size ();

    for ( int ielem = 0; ielem < elemCount; ielem++ )
    {
      elemMatMap_[eIndices[ielem]] = iMat;
    }

    damMat->allocPoints  ( elemCount * ipCount_ );

    materials_[iMat] = damMat;
  }

  //System::out() << elemMatMap_ <<"\n"; 
  //System::out() << egroup_.size()<<"\n";

  // initalize the mapping between integration points
  // and the material points

  initializeIPMPMap_ ();

  // Initialize the weight function that is used for the non-local
  // integration. This function has a single argument name "r" that is
  // the distance between two integration points.

  weightFunc_ = FuncUtils::newFunc

    ( "r", WEIGHT_FUNC_PROP, myProps, globdat );

  FuncUtils::getConfig ( myConf, weightFunc_, WEIGHT_FUNC_PROP );

  myProps.find ( radius_   , RADIUS_PROP,  0.0, maxOf( radius_) );

  // Select the correct function for computing the B-matrix.

  getShapeGrads_ = getShapeGradsFunc ( rank_ );
  getShapeFuncs_ = getShapeFunc      ( rank_ );

  // This is to signal that the integration point pairs have not yet
  // been updated; this will be done when needed.

  updated_ = false;
  isStep_  = false; // default option: do not use step back 

  myProps.find ( isStep_, "isStep" );

  eqWeights_.resize   ( totalIpCount );
  eqWeights_ = 0.0;

  interactMap_.resize ( totalIpCount ); 

  isActive_.resize ( numElem_ );
  isActive_ = 1;
}


NonlocalDamageModel::~NonlocalDamageModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void NonlocalDamageModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void NonlocalDamageModel::getConfig ( const Properties& conf ) const
{
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool NonlocalDamageModel::takeAction

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

    Ref<MBuilder>  mbuilder;

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
    Ref<MBuilder> mbuilder;

    params.get ( mbuilder, ActionParams::MATRIX2 );

    return true;
  }

  // This action is called when solver gets convergence
  // Most often, it is time when material updates history variables

  if ( action == Actions::COMMIT )
  {
    for ( int iMat = 0; iMat < materials_.size ( ); iMat++ )
    { 
      materials_[iMat]->commit ();
    }

    return true;
  }

  // used with the StepModule
  // to advance to new load step or resolve with the same load vector

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


void NonlocalDamageModel::getMatrix_

  ( MBuilder&       mbuilder,
    const Vector&   force,
    const Vector&   disp )

{
  Cubix       grads      ( rank_, nodeCount_, ipCount_ );
  Matrix      coords     ( rank_, nodeCount_           );
  Cubix       grads2     ( rank_, nodeCount_, ipCount_ );
  Matrix      coords2    ( rank_, nodeCount_           );

  Matrix      elemMat    ( dofCount_, dofCount_ );
  Matrix      corrMat    ( dofCount_, dofCount_ );
  
  Vector      elemForce  ( dofCount_ );
  Vector      elemDisp   ( dofCount_ );

  Matrix      stiff      ( strCount_, strCount_ );

  Vector      stress     ( strCount_ );
  Vector      strain     ( strCount_ );

  Matrix      b          ( strCount_, dofCount_ );
  Matrix      bt         = b.transpose ();

  Matrix      b2         ( strCount_, dofCount_ );

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   inodes2    ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );
  IdxVector   idofs2     ( dofCount_  );

  Vector      ipWeights  ( ipCount_   );
  Vector      ipWeights2 ( ipCount_   );

  MChain1     mc1;
  MChain3     mc3;

  double      temp = 0.0;
  double      wip;

  int         iMat, ipoint, igpoint = 0;

  Ref<DamageMaterial>   eMat; 

  Flex<IPointPair>      ipointPairs;

  // Calculate the strain for this updated displacement:
  //
  //  - strain vector stored at the beginning of vector NonlocalDamageModel::strain_
  //  - local equivalent strain
  //  - nonlocal equivalent strain by weighted averaging: stored at the end 
  //                                        of vector NonlocalDamageModel::strain_

  calcStrain_ ( disp );

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < numElem_; ie++ )
  {
    int  ielem = ielems_[ie];

    // ignore disabled element (fully damaged element)

    if ( ! isActive_[ie] )
    {
      continue;
    }

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem             );
    nodes_.getSomeCoords ( coords, inodes            );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Get the displacements at the element nodes.

    elemDisp = select ( disp, idofs );

    // get correct material of current element

    iMat   = elemMatMap_[ielem]; 
    eMat   = materials_[iMat];

    Matrix      stiff0 = eMat->giveElasticMatrix ( ); 

    // The element matrix and the internal force vector.

    elemMat   = 0.0;
    elemForce = 0.0;

    // Loop on integration points 

    for ( int ip = 0; ip < ipCount_; ip++, igpoint++ )
    {
      // Compute the B-matrix for this integration point.

      getShapeGrads_ ( b, grads(ALL,ALL,ip) );

      // Update the material model and update the total integration
      // point index.

      ipoint = ipMpMap_ (ielem,ip); //System::out() << "ipoint: " << ipoint << "\n\n";

      eMat->update ( stress, stiff, strain_(ALL,igpoint), ipoint );

      // Add the contribution of this integration point.

      wip        = ipWeights[ip];

      elemMat   += wip * mc3.matmul ( bt, stiff, b );
      elemForce += wip * mc1.matmul ( bt, stress   );

      // Assembly the correction part of the consistent tangent stiffness
      // matrix (Jirasek et al.), only for loading case 
      // assembly point by point
	
      if ( eMat->isLoading ( ipoint ) )
      {
	    int     j, k;
        int     je, jelem;

        double  weight2; 

	    MChain2 mc2;

        double    kappa  = eMat->giveHistory     ( ipoint );
	    double    dOmega = eMat->getdOmegadKappa ( kappa  );
	    double    wd     = dOmega * wip;

        // normal strain vector at the current integration point ipoint

	    strain  = strain_(slice(BEGIN,strCount_), igpoint); //HERE !!!

	// Look for integration points in neighbouring of
	// the current integration point "igpoint"

        ipointPairs  = interactMap_[igpoint];

        int nipCount = ipointPairs.size ();

	// Then, loop on neighboring Gauss points of igpoint
	
	for ( int jp = 0; jp < nipCount; jp++ )
	{
	  const IPointPair& ipp = ipointPairs[jp];
	 
          // find the integration point interacted with ipoint
          // and the weight associated with it

          if ( igpoint == ipp.ipoint )
	  {
	    j       = ipp.jpoint; 
	    weight2 = ipp.weight2; 
	  }
	  else
	  {
	    j       = ipp.ipoint; 
	    weight2 = ipp.weight1; 
	  }

	  // element contains this GP j and its global index

	  je    =  j / ipCount_; jelem = ielems_[je];

          // if element je is fully damaged, then
          // ignore its integration point

          if ( ! isActive_[je] )
          {
            continue;
          }

	  // compute the derivative of local-equivalent strain w.r.t
          // strain tensor at the interacted integration point, i.e j

          // REMARK: assumed that the calculation of derivative of
          // equivalent strain does not depend on material constant !!!
          // Corrected this for multi material domains !!! 

          int                  jMat  = elemMatMap_[jelem]; 
          Ref<DamageMaterial>  ejMat = materials_ [jMat]; 

	  Vector dEpsBardEps = ejMat->getdEpsBardEps ( strain_(slice(BEGIN,strCount_),j) );

          temp  = ( weight2 / eqWeights_[igpoint] ) * wd;

	  // compute the local index of GP j in element jelem
	  // to get correct Bj matrix

	  if ( je == 0 )
	  {
	    k = j;
	  }
	  else
	  {
	    k = j % (ipCount_ * je);
	  }

	  // if je != ie, need to compute B matrix and
	  // dofs of je. Otherwise, reuse available stuffs

	  if ( jelem != ielem )
	  {
	    // get dofs, B-matrix of this element jelem
 
	    elems_.getElemNodes  ( inodes2, jelem   );
	    nodes_.getSomeCoords ( coords2, inodes2 );

	    shape_->getShapeGradients ( grads2, ipWeights2, coords2 );

	    // get the B matrix at the integration point j of element jelem

	    getShapeGrads_ ( b2, grads2 ( ALL, ALL, k ) );

	    dofs_->getDofIndices ( idofs2, inodes2, dofTypes_ );

	    // Compute the correction stiffness of this pair of
	    // integration points

	    corrMat = -temp * mc3.matmul ( bt, 
                                           mc2.matmul ( stiff0, matmul ( strain, dEpsBardEps ) ), 
                                           b2 );

	    // Then add this to the global stiffness matrix.

	    mbuilder.addBlock ( idofs, idofs2, corrMat );

	  }
	  else
	  {	    
	    getShapeGrads_ ( b2, grads ( ALL, ALL, k ) );

	    corrMat = -temp * mc3.matmul ( bt, 
                                           mc2.matmul ( stiff0, matmul ( strain, dEpsBardEps ) ), 
                                           b2 );

	    mbuilder.addBlock ( idofs, idofs, corrMat );
	  }
 	  
        }      // end of loop on neighboring integration points
      }
    }          // end of loop on integration points 

    // Add the element secant matrix to the global stiffness matrix.

    mbuilder.addBlock ( idofs, idofs, elemMat );

    // Add the element force vector to the global force vector.

    select ( force, idofs ) += elemForce;
  }
}


//-----------------------------------------------------------------------
//   calcStrain_
//-----------------------------------------------------------------------

void NonlocalDamageModel::calcStrain_ ( const Vector&  disp )
{
  // Gradient of shape functions for all nodes at all integration points

  Cubix       grads      ( rank_, nodeCount_, ipCount_ );
  Matrix      coords     ( rank_, nodeCount_           );

  // Nodal displacement of 1 element, strain at 1 GP and weights

  Vector      elemDisp   ( dofCount_ );
  Vector      ipStrain   ( strCount_ );
  Vector      ipWeights  ( ipCount_  );

  // Local equivalent strain at all integration points 

  Vector      eqStrain   ( ipCount_ * numElem_  );

  Matrix      b          ( strCount_, dofCount_ );

  // Connectivity list of 1 element, and position of its DOFs in the
  // global system 

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );

  // Initialize the integration point pairs if necessary.

  if ( ! updated_ )
  {
    initIPPairs_ ();
  }

  strain_     = 0.0;  // have to do this due to accumulation later
  eqStrain    = 0.0;

  int  ipoint = 0;

  // Iterate over all elements assigned to this model to calculate the
  // local equivalent strains.

  for ( int ie = 0; ie < numElem_; ie++ )
  {
    int  ielem = ielems_[ie]; 

    if ( ! isActive_[ie] )
    {
      continue;
    }

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem  );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Get the displacements at the element nodes.

    elemDisp = select ( disp, idofs );

    // get the correct material of current element

    int                   iMat = elemMatMap_[ielem];
    Ref<DamageMaterial>   eMat = materials_[iMat];

    // loop over integration points 

    for ( int ip = 0; ip < ipCount_; ip++, ipoint++ )
    {
      // Compute the B-matrix for this integration point.

      getShapeGrads_ ( b, grads(ALL,ALL,ip) );

      // Compute the strain for this integration point.

      matmul ( ipStrain, b, elemDisp );

      // Store the normal strain components.

      strain_(slice(BEGIN,strCount_),ipoint) = ipStrain;

      // Calculate the local equivalent strain in this integration
      // point.

      eqStrain[ipoint] = eMat->getEquiStrain ( ipStrain ) ;
    }

  }

  // System::out() << "local equivalent strain : \n"
  //              << eqStrain << "\n\n";

  // compute the nonlocal equivalent strain 


  getNonlocalStrainNSym_ ( eqStrain );


  // System::out() << "non-local equivalent strain : \n"
  //              << strain_(ieq,ALL) << "\n\n";
}

//-----------------------------------------------------------------------
//   getNonlocalStrainNSym_
//-----------------------------------------------------------------------

void NonlocalDamageModel::getNonlocalStrainNSym_

( const Vector& eqStrain )

{
  for ( int i = 0; i < ipPairs_.size(); i++ )
  {
    IPointPair&  ipp = ipPairs_[i];

    int     ip = ipp.ipoint;
    int     jp = ipp.jpoint;

    int     ie    =  ip / ipCount_;
    int     je    =  jp / ipCount_;

    if ( ! isActive_[ie] || ! isActive_[je] )
    {
      continue;
    }

    double  w1 = ipp.weight1;
    double  w2 = ipp.weight2;

    if ( ip != jp )
    {
      strain_(strCount_,ip) += w2 * eqStrain[jp];
      strain_(strCount_,jp) += w1 * eqStrain[ip];
    }
    else
    {
      strain_(strCount_,ip) += w2 * eqStrain[jp];
    }
  }

  strain_(strCount_,ALL) /= eqWeights_;

}

//-----------------------------------------------------------------------
//   getNonlocalStrainSym_
//-----------------------------------------------------------------------

void NonlocalDamageModel::getNonlocalStrainSym_

( const Vector& eqStrain )

{
  for ( int i = 0; i < ipPairs_.size(); i++ )
  {
    IPointPair&  ipp = ipPairs_[i];

    int     ip = ipp.ipoint;
    int     jp = ipp.jpoint;

    double  w1 = ipp.weight1;
    double  w2 = ipp.weight2;

    if ( ip != jp )
    {
      strain_(strCount_,ip) += w2 * eqStrain[jp];
      strain_(strCount_,jp) += w1 * eqStrain[ip];
    }
    else
    {
      strain_(strCount_,ip) +=  
          ( 1.0 - eqWeights_[ip] + ipp.weight1 ) * eqStrain[jp];
    }
  }

}

//-----------------------------------------------------------------------
//   initIPPairs_
//-----------------------------------------------------------------------


void NonlocalDamageModel::initIPPairs_ ()
{
  using jive::geom::makeBoxes;
  using jive::geom::getBoxCenters;
  using jive::geom::BoxManager;
  using jive::geom::GridBoxManager;

  Ref<BoxManager>  boxman = newInstance<GridBoxManager> ( rank_ );

  Matrix       boxes      ( rank_ * 2, numElem_ );
  Matrix       centers    ( rank_,     numElem_ );

  Matrix       ndCoords   ( rank_, nodeCount_ );
  Matrix       ipCoords1  ( rank_, ipCount_   );
  Matrix       ipCoords2  ( rank_, ipCount_   );

  Vector       ipWeights1 ( ipCount_ );
  Vector       ipWeights2 ( ipCount_ );

  IdxVector    inodes     ( nodeCount_ );
  IdxVector    jnodes     ( nodeCount_ );
  IdxVector    jelems     ( numElem_   );

  int          ipoint;
  int          jpoint;


  System::info() << myName_
		 << " : updating integration point pairs ...\n";

  ipPairs_.clear ();

  // Get the bounding boxes of the elements. Each box is represented
  // by a vector containing the minimum and maximum element  coordinates.

  elems_.getSomeElemBoxes ( boxes, ielems_ );

  // Get the box centers.

  getBoxCenters ( centers, boxes );

  // Make new boxes with the correct size.

  makeBoxes ( boxes, centers, ::sqrt( (double) rank_ ) * radius_ );

  // Add the boxes to the box manager.

  boxman->addBoxes ( boxes );

  // Iterate over the elements assigned to this model.

  for ( int ie = 0; ie < numElem_; ie++ )
  {
    int  ielem = ielems_[ie];

    elems_.getElemNodes  ( inodes,   ielem  );
    nodes_.getSomeCoords ( ndCoords, inodes );

    // Calculate the global integration point coordinates and weights.

    shape_->getGlobalIntegrationPoints ( ipCoords1,  ndCoords );
    shape_->getIntegrationWeights      ( ipWeights1, ndCoords );

    // Find the integration point pairs within the current element.
    // ipoint is the global number of the first integration point of 
    // an element

    ipoint = ie * ipCount_;

    findIPPairs_ ( radius_,
		   ipoint,     ipoint,
		   ipCoords1,  ipCoords1,
		   ipWeights1, ipWeights1 );

    // Get the neighbors of this element.

    int  n = boxman->findBoxNeighbors ( jelems, ie );

    for ( int j = 0; j < n; j++ )
    {
      // Get the local element index of this neighbor element.

      int  je = jelems[j];

      // This is to avoid finding the same pair of integration points
      // twice.

      if ( ie <= je )
      {
	continue;
      }

      // Calculate the global integration point coordinates and
      // weights.

      elems_.getElemNodes  ( jnodes,   je     );
      nodes_.getSomeCoords ( ndCoords, jnodes );

      shape_->getGlobalIntegrationPoints ( ipCoords2,  ndCoords );
      shape_->getIntegrationWeights      ( ipWeights2, ndCoords );

      // Find the integration point pairs in the elements ielem and
      // jelem.

      jpoint = je * ipCount_;

      findIPPairs_ ( radius_,
		     ipoint,     jpoint,
		     ipCoords1,  ipCoords2,
		     ipWeights1, ipWeights2 );
    }
  }

  // compute the weights 

  for ( int i = 0; i < ipPairs_.size(); i++ )
  {
    const IPointPair&  ipp = ipPairs_[i];

    double  w1 = ipp.weight1;
    double  w2 = ipp.weight2;

    if ( ipp.ipoint != ipp.jpoint )
    {
      eqWeights_[ipp.ipoint]   += w2;
      eqWeights_[ipp.jpoint]   += w1;
    }
    else
    {
      eqWeights_[ipp.ipoint]   += w2;
    }
  }


  System::info() << myName_
		 << " : number of integration point pairs = "
		 << ipPairs_.size() << '\n';

  updated_ = true;


  // debug the integration point pairs

  /*
  for ( int i = 0; i < ipPairs_.size(); i++ )
  {
    const IPointPair&  ipp = ipPairs_[i];

    System::out() << " integration point pair: "
                  << "(" << ipp.ipoint << ", " << ipp.jpoint << ")"
                  << "\n";
  }  

  System::out() << " \n\n";

  // debug the map of IP => integration point pairs 

  Flex<IPointPair> ipointPairs;

  for ( int ip = 0; ip < ielemCount * ipCount_; ip++ )
  {
    ipointPairs = interactMap_[ip];

    for ( int jp = 0; jp < ipointPairs.size(); jp++ )
    {
      const IPointPair& ipp = ipointPairs[jp];

      System::out() << " integration point pairs of point  " << ip << ": "
		    << "(" << ipp.ipoint << ", " << ipp.jpoint << ")"
		    << "\n";

      int i   =  ipp.ipoint; 
      int j   =  ipp.jpoint;   
      int k, m;
 
      j       = ( ip == i ) ? j : i; 

      m       = j / ipCount_;

      if ( m == 0 )
      {
	k = j;
      }
      else
      {
	k = j % (ipCount_ * m );
      }

      System::out() << " point  " << ip << " interacts with " << j
                    << " which is IP indexed " << k
                    << "  and who belongs to element numbered "    << m << "\n"; 
    }

    System::out() << " \n";
  }
  */
}


//-----------------------------------------------------------------------
//   findIPPairs_
//-----------------------------------------------------------------------


void NonlocalDamageModel::findIPPairs_

  ( double         radius,
    int            ipoint,
    int            jpoint,
    const Matrix&  coords1,
    const Matrix&  coords2,
    const Vector&  weights1,
    const Vector&  weights2 )

{
  using jem::numeric::norm2;

  const int   ipCount = coords1.size (1);
  const int   jpCount = coords2.size (1);
  const int   jpoint0 = jpoint;

  IPointPair        ipp;

  for ( int ip = 0; ip < ipCount; ip++, ipoint++ )
  {
    const Vector  u = coords1(ALL,ip);

    jpoint = jpoint0;

    for ( int jp = 0; jp < jpCount; jp++, jpoint++ )
    {
      // This is to avoid finding the same pair twice.

      if ( ipoint < jpoint )
      {
	continue;
      }

      Vector  v = coords2(ALL,jp);
      double  d = norm2 ( u - v );

      if ( d < radius_ )
      {
	double  w = weightFunc_->eval ( d );

	ipp.ipoint  = ipoint;
	ipp.jpoint  = jpoint;
	ipp.dist    = d;
	ipp.weight1 = w * weights1[ip];
	ipp.weight2 = w * weights2[jp];

        // add the newly found pair into the list

	ipPairs_.pushBack ( ipp );

        // add this to the map also, attention not to add twice for pairs (i,i)
	
        if ( ipoint != jpoint )
        { 
	  interactMap_[ipoint].pushBack ( ipp );
	  interactMap_[jpoint].pushBack ( ipp );
        }
        else
        {
          interactMap_[ipoint].pushBack ( ipp );
        }
      }
    }

  }
}

//-----------------------------------------------------------------------
//   initializeIPMPMap_
//-----------------------------------------------------------------------

void  NonlocalDamageModel::initializeIPMPMap_ ( )

{
        int   ipoint;

  IdxVector   matMap ( materials_.size ( ) );

  matMap = 0;

  for ( int ie = 0; ie < numElem_; ie++ )
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

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void NonlocalDamageModel::checkCommit_

  ( const Properties&  params )

{
  if ( !isStep_ )
  {
    params.set ( "accept", true );

    return;
  }

  const int             ipCount_    = shape_->ipointCount ();
  bool                  restart    = false; 
  Ref<DamageMaterial>   eMat;

  for ( int ie = 0; ie < numElem_; ie++ )
  {
    int  ielem = ielems_[ie];

    // do not check element already disabled

    if ( ! isActive_[ielem] )
    {
      continue;
    }

    // get correct material

    int matID = elemMatMap_[ielem];
    eMat      = materials_[matID];

    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      int ipoint = ipMpMap_ (ielem, ip);

      if ( eMat->isFullyDamaged ( ipoint ) )
      {
        isActive_[ielem] = 0;

        // found fully damaged element, then should recompute
        // for this load increment

        restart          = true;

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

//-----------------------------------------------------------------------
//   getTable_
//-----------------------------------------------------------------------


bool NonlocalDamageModel::getTable_

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


// //-----------------------------------------------------------------------
// //   getStrain_
// //-----------------------------------------------------------------------


// void NonlocalDamageModel::getStrain_

//   ( XTable&        table,
//     const Vector&  weights,
//     const Vector&  disp )

// {
//   Matrix     sfuncs     = shape_->getShapeFunctions ();

//   Matrix     ndStrain   ( nodeCount_, strCount_ + 1 );
//   Vector     ndWeights  ( nodeCount_ );

//   IdxVector  inodes     ( nodeCount_ );
//   IdxVector  jcols      ( strCount_ + 1 );

//   int        ipoint;


//   // Add the columns for the normal strain components to the table.

//   switch ( strCount_ )
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
//       String ( strCount_ )
//     );
//   }

//   // Add an extra column for the non-local equivalent strain.

//   jcols[strCount_] = table.addColumn ( "e_eq" );

//   // Update the strains.

//   calcStrain_ ( disp );

//   // Iterate over all elements assigned to this model.

//   ipoint = 0;

//   for ( int ie = 0; ie < numElem_; ie++ )
//   {
//     int  ielem = ielems_[ie];

//     // Get the element nodes.

//     elems_.getElemNodes ( inodes, ielem );

//     // Iterate over the integration points.

//     ndStrain  = 0.0;
//     ndWeights = 0.0;

//     for ( int ip = 0; ip < ipCount_; ip++, ipoint++ )
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

// //-----------------------------------------------------------------------
// //   getStress_
// //-----------------------------------------------------------------------


// void NonlocalDamageModel::getStress_

//   ( XTable&        table,
//     const Vector&  weights )

// {
//   Matrix     sfuncs     = shape_->getShapeFunctions ();

//   Matrix     ndStress   ( nodeCount_, strCount_ );
//   Vector     ndWeights  ( nodeCount_ );

//   Vector     stressIp   ( strCount_ );

//   IdxVector  inodes     ( nodeCount_ );
//   IdxVector  jcols      ( strCount_  );

//   MChain1    mc1;

//   int        ipoint, iMat, igpoint = 0;

//   Ref<DamageMaterial> eMat;


//   // Add the columns for the stress components to the table.

//   switch ( strCount_ )
//   {
//   case 1:

//     jcols[0] = table.addColumn ( "s_xx" );

//     break;

//   case 3:

//     jcols[0] = table.addColumn ( "s_xx" );
//     jcols[1] = table.addColumn ( "s_yy" );
//     jcols[2] = table.addColumn ( "s_xy" );

//     break;

//   case 6:

//     jcols[0] = table.addColumn ( "s_xx" );
//     jcols[1] = table.addColumn ( "s_yy" );
//     jcols[2] = table.addColumn ( "s_zz" );
//     jcols[3] = table.addColumn ( "s_xy" );
//     jcols[4] = table.addColumn ( "s_yz" );
//     jcols[5] = table.addColumn ( "s_xz" );

//     break;

//   default:

//     throw Error (
//       JEM_FUNC,
//       "unexpected number of stress components: " +
//       String ( strCount_ )
//     );
//   }

//   // Iterate over all elements assigned to this model.

//   for ( int ie = 0; ie < numElem_; ie++ )
//   {
//     int  ielem = ielems_[ie];

//     elems_.getElemNodes ( inodes, ielem );

//     // get correct material

//     iMat   = elemMatMap_[ielem];
//     eMat   = materials_[iMat];

//     Matrix stiff0  = eMat->giveElasticMatrix ();
//     Vector omega   = eMat->giveOmega         (); 

//     // Iterate over the integration points.

//     ndStress  = 0.0;
//     ndWeights = 0.0;

//     for ( int ip = 0; ip < ipCount_; ip++, igpoint++ )
//     {
//       // Extrapolate the integration point stresses to the nodes using
//       // the transposed shape functions.

//       stressIp   = mc1.matmul ( stiff0, strain_(slice(BEGIN,strCount_),igpoint) );

//       ipoint     = ipMpMap_ (ielem,ip); 

//       stressIp  *= 1.0 - eMat->giveOmega ( ipoint );

//       ndStress  += matmul ( sfuncs(ALL,ip), stressIp );
//       ndWeights += sfuncs(ALL,ip);
//     }

//     select ( weights, inodes ) += ndWeights;

//     // Add the stresses to the table.

//     table.addBlock ( inodes, jcols, ndStress );
//   }
// }


// //-----------------------------------------------------------------------
// //   getHistory_
// //-----------------------------------------------------------------------


// void NonlocalDamageModel::getHistory_

//   ( XTable&          table,
//     const Vector&    weights )
// {
//   Matrix     coords     ( rank_, nodeCount_ );
//   Vector     ipWeights  ( ipCount_ );

//   IdxVector  inodes     ( nodeCount_ );
//   IdxVector  jcols      ( 1 );

//   Matrix     ndHistory  ( nodeCount_, 1 );
//   Vector     ndWeights  ( nodeCount_    );

//   Ref<DamageMaterial> eMat;

//   int                 iMat;

//   Matrix     sfuncs     = shape_->getShapeFunctions();

//   // Add a column for the history parameter to the table.

//   jcols[0] = table.addColumn ( "hist" );

//   // Iterate over all elements assigned to this model.

//   for ( int ie = 0; ie < numElem_; ie++ )
//   {
//     int  ielem = ielems_[ie];

//     // Get the element coordinates.

//     elems_.getElemNodes  ( inodes, ielem  );

//     // get correct material

//     iMat   = elemMatMap_[ielem];
//     eMat   = materials_[iMat];

//     // Get the history from the material

//     //Vector  hist = eMat->giveOmega(  ); //System::out() << "omega " << hist <<"\n";

//     ndHistory = 0.0;
//     ndWeights = 0.0;

//     for ( int ip = 0; ip < ipCount_; ip++ )
//     {
//       int ipoint = ipMpMap_ (ielem,ip);

//       ndHistory(ALL,0) += abs ( sfuncs(ALL,ip) ) * eMat->giveOmega ( ipoint );
//       ndWeights        += abs ( sfuncs(ALL,ip) );
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

void NonlocalDamageModel::getOutputData_

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
  
   
  Matrix     ndValuesOut  ( nodeCount_, i2table.size() );
  Matrix     ndValuesOut1 ( nodeCount_, i2table.size() );
  Vector     ipValuesOut  ( i2table.size() );

  Matrix      stiff      ( strCount_, strCount_  );
//  Matrix      bd         ( strCount_, dispCount );
  Vector      ndWeights  ( nodeCount_           ); 
  IdxVector   inodes     ( nodeCount_           );


  Ref<DamageMaterial> eMat;

  Properties  params;

  idx_t       iMat, ipoint, igpoint = 0;

  Matrix      sfuncs     = shape_->getShapeFunctions ();

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < numElem_; ie++ )
  {
    int  ielem = ielems_[ie];

    elems_.getElemNodes ( inodes, ielem );

    // get correct material

    iMat   = elemMatMap_[ielem];
    eMat   = materials_[iMat];

    Matrix stiff0  = eMat->giveElasticMatrix ();
    Vector omega   = eMat->giveOmega         (); 

    ndValuesOut = 0.;
    ndWeights   = 0.;

    // Loop on integration points 
    
    for ( idx_t ip = 0; ip < ipCount_; ip++, igpoint++ )
    {
      ipoint      = ipMpMap_( ielem, ip );
      
 
      strain      = strain_(slice(BEGIN,strCount_),igpoint);
      dam[0]      = omega[igpoint]; // non-local dam

      stress      = matmul ( stiff0, strain );
      stress     *= 1.0 - eMat->giveOmega ( ipoint );
 
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
  return newInstance<NonlocalDamageModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareDamageModel
//-----------------------------------------------------------------------


void declareNonLocalDamageModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "NonLocalDamage", & newDamageModel );
}
