/*
 *
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *
 *  This class implemens the linear elasticity model  using finite elements.
 *  The number of elements and the integration scheme/element is fixed.
 *
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *
 *  Updates (what and who):
 *
 *  - October 2014: change number of stress/strain components of 2D problems to
 *  4 (before was 3). VPN
 *
 */

#include <jem/base/Array.h>
#include <jem/base/Error.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/System.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/util/Properties.h>
#include <jem/io/PrintWriter.h>
#include <jem/io/FileWriter.h>

#include <jive/util/utilities.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/XTable.h>
#include <jive/util/Assignable.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/geom/IShapeFactory.h>
#include <jive/geom/InternalShape.h>
#include <jive/geom/Geometries.h>
#include <jive/fem/ElementGroup.h>

#include "femodel/FemModels.h"
#include "util/utilities.h"
#include "util/TbFiller.h"
#include "material/Material.h"
#include "material/HookeMaterial.h"
#include "material/OrthotropicMaterial.h"

using namespace jem::io;

using namespace jem;

using jem::numeric::MatmulChain;
using jem::util::Properties;
using jem::io::PrintWriter;
using jem::io::FileWriter;


using jive::Vector;
using jive::IdxVector;
using jive::Cubix;
using jive::util::XDofSpace;
using jive::util::XTable;
using jive::util::Assignable;
using jive::algebra::MatrixBuilder;
using jive::model::Model;
using jive::geom::IShape;
using jive::fem::NodeSet;
using jive::fem::ElementSet;
using jive::fem::ElementGroup;

// some typedef to avoid typing

typedef ElementSet              ElemSet;
typedef ElementGroup            ElemGroup;
typedef MatmulChain<double,3>   MChain3;
typedef MatmulChain<double,2>   MChain2;
typedef MatmulChain<double,1>   MChain1;



class ElasticityModel : public Model
{
 public:

  typedef ElasticityModel    Self;
  typedef Model              Super;

  static const char*         DOF_NAMES[3];
  static const char*         SHAPE_PROP;
  static const char*         MATERIAL_PROP;

                       ElasticityModel

    ( const String&       name,
      const Properties&   conf,
      const Properties&   props,
      const Properties&   globdat );

  virtual void         configure

    ( const Properties&   props,
      const Properties&   globdat );

  virtual void         getConfig

    ( const Properties&   conf )      const;

  virtual bool         takeAction

    ( const String&       action,
      const Properties&   params,
      const Properties&   globdat );

 protected:

  virtual              ~ElasticityModel ();

 private:

  void                 getIntForce_

    ( const Vector&       force,
      const Vector&       disp )       const;

  void                 getMatrix_

    ( MatrixBuilder&      mbuilder,
      const Vector&       force,
      const Vector&       disp )       const;

  void                 getMatrix2_

    ( MatrixBuilder&      mbuilder );

  void                 computeVolume_ () const;

  bool                 getTable_

    ( const Properties&   params,
      const Properties&   globdat );

  void                 getHistory_

    ( XTable&             table,
      const Vector&       weights );


  void                 getOutputData_

    ( Ref<XTable>             table,
      const Vector&           weights,
      const String&           contents,
      const Vector&           disp  );    
  
  void                 getStress_

   ( XTable&        table,
     const Vector&  weights,
     const Vector&  disp );


 private:

  Assignable<ElemGroup>   egroup_;
  Assignable<ElemSet>     elems_;
  Assignable<NodeSet>     nodes_;

  int                     rank_;

  Ref<IShape>             shape_;

  Ref<XDofSpace>          dofs_;
  IdxVector               dofTypes_;

  Matrix                  stiff0_;

  ShapeGradsFunc          getShapeGrads_;
  ShapeFunc               getShapeFuncs_;

  IdxVector               ielems_;

  int                     ielemCount_;
  int                     nodeCount_;
  int                     ipCount_;
  int                     strCount_;
  int                     dofCount_;

  double                  thickness_;       // thickness for plane stress

  int                     outNode_;
  Ref<PrintWriter>        out_;

  Array<Ref<Material> >   materials_;
  IdxVector               elemMatMap_;

  double                  rho_;
  double                  g_;
};

// -------------------------------------------------------------
//    static constants
// -------------------------------------------------------------

const char* ElasticityModel::DOF_NAMES[3]  = {"dx","dy","dz"};
const char* ElasticityModel::SHAPE_PROP    = "shape";
const char* ElasticityModel::MATERIAL_PROP = "material";

// -------------------------------------------------------------
//    constructor
// -------------------------------------------------------------

ElasticityModel::ElasticityModel

   ( const String&       name,
     const Properties&   conf,
     const Properties&   props,
     const Properties&   globdat ) : Super(name)
{
  using jive::util::joinNames;
  using jive::geom::IShapeFactory;

  Properties  myProps = props.getProps ( myName_ );
  Properties  myConf  = conf.makeProps ( myName_ );

  const String context = getContext();

  egroup_ = ElemGroup::get ( myConf, myProps, globdat, context );

  elems_  = egroup_.getElements ( );
  nodes_  = elems_.getNodes     ( );
  rank_   = nodes_.rank         ( );

  ielemCount_ = egroup_.size    ( );

  ielems_.resize ( ielemCount_ );
  ielems_ = egroup_.getIndices ();

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

  shape_  = IShapeFactory::newInstance(
    joinNames (myName_, SHAPE_PROP ),
    conf,
    props );

  //System::out() << shape_->getIntegrationScheme();

  // compute some constants

  nodeCount_  = shape_->nodeCount    ();
  ipCount_    = shape_->ipointCount  ();
  strCount_   = STRAIN_COUNTS[rank_];
  dofCount_   = rank_ * nodeCount_;

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

  elems_.checkSomeElements ( context, ielems_, nodeCount_ );

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  dofTypes_.resize( rank_ );

  for ( int i = 0; i < rank_; i++ )
  {
    dofTypes_[i] = dofs_->addType ( DOF_NAMES[i]);
  }

  dofs_->addDofs ( elems_.getUniqueNodesOf ( ielems_ ), dofTypes_ );


  // ------------------------------------------------------
  //       create material objects
  // ------------------------------------------------------

  StringVector mats;

  myProps.get ( mats ,MATERIAL_PROP );
  myConf. set ( MATERIAL_PROP, mats );

  int numMats = mats.size ( );
  materials_.  resize  ( numMats     );

  //EJ: set size to the total number of elements (on this processor).
  elemMatMap_. resize  ( elems_.size() );

  elemMatMap_ = -1;

  for ( int iMat = 0; iMat < numMats; iMat++ )
  {
    String matName         = mats[iMat];
    Ref<Material> tmpMat   = newMaterial ( matName, myConf, myProps, globdat );
    Properties    matProps = myProps.findProps ( matName );
    Properties    matConf  = myConf.makeProps  ( matName );

    tmpMat->configure ( matProps, globdat );
    tmpMat->getConfig ( matConf,  globdat );

    // get the elements of this material

    ElementGroup  egroup    = ElementGroup::get ( matConf, matProps, globdat, context );
    IdxVector     eIndices  = egroup.getIndices ();

    int elemCount = eIndices.size ();

    for ( int ielem = 0; ielem < elemCount; ielem++ )
    {
      //EJ: use the element index.
      elemMatMap_[eIndices[ielem]] = iMat;
    }
    materials_[iMat] = tmpMat;

    // store the materials in globdat so that other models can access to
    // Other models: Discontinuous Galerkin interface model

    globdat.set ( String::format( "mat%d", iMat ), tmpMat);
  }
    
  //System::out() << "BulkModel::elemMatMap_ " << elemMatMap_ << "\n";

  // also store the elemMatMap in globdat

  globdat.set ( "elemMatMap" , elemMatMap_ );

  // correct function to compute B matrix

  getShapeGrads_  = getShapeGradsFunc        ( rank_ );
  getShapeFuncs_  = getShapeFunc             ( rank_ );

  thickness_ = 1;
  myProps.find ( thickness_ , "thickness" );
  myConf. set  ( "thickness", thickness_  );

  // compute the volume

  computeVolume_();

  outNode_ = -1;

  int node;

  if ( myProps.find ( node, "node" ) )
  {
    outNode_ = nodes_.findNode ( node );

    out_     = newInstance<PrintWriter> (
		  newInstance<FileWriter> ( "stress-disp.dat" ) );
  }

  stiff0_.resize(strCount_,strCount_);
  
  rho_ = 0.; 
  g_   = 0.;
  
  myProps.find ( rho_  , "density" );
  myProps.find ( g_  ,   "g" );

  myConf.set ( "density", rho_ );
  myConf.set ( "g", g_ );

}

ElasticityModel::~ElasticityModel()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void ElasticityModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void ElasticityModel::getConfig ( const Properties& conf ) const
{
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool ElasticityModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  // This is called by explicit time integration module

  if ( action == Actions::GET_INT_VECTOR )
  {
    Vector  disp;
    Vector  force;

    // Get the current displacements.

    StateVector::get ( disp, dofs_, globdat );

    // Get the matrix builder and the internal force vector.

    params.get ( force,    ActionParams::INT_VECTOR );

    getIntForce_ ( force, disp );

    return true;
  }


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

  if ( action == Actions::GET_MATRIX2 )
  {
    Ref<MatrixBuilder> mbuilder;

    params.get ( mbuilder, ActionParams::MATRIX2 );

    getMatrix2_( *mbuilder );

    return true;
  }

  if ( action == Actions::GET_TABLE )
  {
    return getTable_ ( params, globdat );
  }

  /*
   * Why calculating stress here?? [Phu]
  if ( action == Actions::COMMIT )
  {
    Vector  disp;

    StateVector::get ( disp, dofs_, globdat );

    getStress_ ( globdat, disp );

    return true;
  }
  */

  return false;
}

//-----------------------------------------------------------------------
//   getIntForce_
//-----------------------------------------------------------------------


void ElasticityModel::getIntForce_

  ( const Vector&   force,
    const Vector&   disp ) const

{
  using jem::numeric::matmul;

  Cubix       grads      ( rank_, nodeCount_, ipCount_ );

  Matrix      coords     ( rank_, nodeCount_ );

  Vector      elemForce  ( dofCount_ );
  Vector      elemDisp   ( dofCount_ );

  Vector      strain     ( strCount_ );
  Vector      stress     ( strCount_ );
  Matrix      stiff      ( strCount_, strCount_ );

  Matrix      b          ( strCount_, dofCount_  );
  Matrix      bt         = b.transpose ();

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );

  Vector      ipWeights  ( ipCount_   );

  MChain1     mc1;

  idx_t          iMat;
  Ref<Material>  eMat;

  double         wip;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount_; ie++ )
  {
    int  ielem = ielems_[ie];
    elems_.getElemNodes  ( inodes, ielem             );
    nodes_.getSomeCoords ( coords, inodes            );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );
    shape_->getShapeGradients ( grads, ipWeights, coords );
    elemDisp = select ( disp, idofs );

    ipWeights *= thickness_;
    // get the correct material

    //EJ: use the element index and not the group index.
    iMat    = elemMatMap_[ielem];
    eMat    = materials_[iMat];

    // Assemble the element matrix.

    elemForce = 0.0;

    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      getShapeGrads_ ( b, grads(ALL,ALL,ip) );
      matmul ( strain, b, elemDisp     );
      eMat->update ( stress, stiff, strain, 0  );
      wip = ipWeights[ip];
      elemForce += wip * mc1.matmul ( bt, stress   );
    }

    // Add the element force vector to the global force vector.

    select ( force, idofs ) += elemForce;
  }
}

//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void ElasticityModel::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   disp ) const

{
  using jem::numeric::matmul;

  Cubix       grads      ( rank_, nodeCount_, ipCount_ );

  Matrix      coords     ( rank_,     nodeCount_ );
  Matrix      elemMat    ( dofCount_, dofCount_  );
  Vector      elemForce  ( dofCount_ );
  Vector      BtSigma    ( dofCount_ );
  Vector      elemDisp   ( dofCount_ );

  Vector      strain     ( strCount_ );
  Vector      stress     ( strCount_ );
  Matrix      stiff      ( strCount_, strCount_ );

  Matrix      b          ( strCount_, dofCount_  );
  Matrix      bt         = b.transpose ();
  
  Matrix      sfuncs     = shape_->getShapeFunctions ();
  Matrix      N          ( rank_, rank_ * nodeCount_ );
  Matrix      Nt         = N.transpose ( );

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );

  Vector      ipWeights  ( ipCount_   );

  MChain1     mc1;
  MChain3     mc3;
  
  Vector      selfW ( rank_ );

  selfW         = 0.; 
  selfW[1]      = -rho_ * g_;

  idx_t          iMat;
  Ref<Material>  eMat;

  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    // Get the global element index.

    idx_t  ielem = ielems_[ie];

    //System::out() << ielem << "\n";

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem             );
    nodes_.getSomeCoords ( coords, inodes            );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    //System::out() << "solid nodes " << inodes << "\n";
    // Compute the spatial derivatives of the element shapoe
    // functions.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    ipWeights *= thickness_;

    // Get the displacements at the element nodes.

    elemDisp = select ( disp, idofs );

    // get the correct material

    //EJ: use the element index and not the group index.
    iMat    = elemMatMap_[ielem];
    
    //System::out() << "BulkModel: iMat of " << ielem << ": " << iMat << "\n";
    
    eMat    = materials_[iMat];

    // Assemble the element matrix.

    elemMat   = 0.0;
    elemForce = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Compute the B-matrix for this integration point.

      getShapeGrads_ ( b, grads(ALL,ALL,ip) );

      // Compute the strain vector of this integration point

      matmul ( strain, b, elemDisp     );

      eMat->update ( stress, stiff, strain, 0  );

      // Compute the element force vector.
      // Compute the stiffness matrix

      double wip = ipWeights[ip];

      matmul ( BtSigma, bt, stress   );

      elemForce += wip * BtSigma;
      elemMat   += wip * mc3.matmul ( bt, stiff, b );
      
      // self weight load is considered as internal force, with a minus sign

      if ( g_ > 0. ) 
      {
        getShapeFuncs_ ( N, sfuncs(ALL,ip) );
        elemForce  -= wip  *     matmul ( Nt, selfW      );
      }
    }

    // Add the element matrix to the global stiffness matrix.

    mbuilder.addBlock ( idofs, idofs, elemMat );

    // Add the element force vector to the global force vector.

    select ( force, idofs ) += elemForce;
  }
}

//-----------------------------------------------------------------------
//   getMatrix2_
//-----------------------------------------------------------------------

// compute the mass matrix
// consistent mass matrix
// If explicit time integration schemes are used, this consistent mass matrix
// will be automatically by jem/jive to a lumped mass.

void ElasticityModel::getMatrix2_

    ( MatrixBuilder&          mbuilder )
{
  Cubix       grads      ( rank_, nodeCount_, ipCount_ );
  Matrix      coords     ( rank_, nodeCount_ );

  Matrix      elemMat    ( dofCount_, dofCount_ );

  Matrix      R          ( rank_, rank_ );

  Matrix      sfuncs     = shape_->getShapeFunctions ();
  Matrix      N          ( rank_, rank_ * nodeCount_ );
  Matrix      Nt         = N.transpose ( );

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );

  Vector      ipWeights  ( ipCount_   );

  MChain2     mc2;
  MChain3     mc3;

  double      rho;
  //R = 0.0;

  //for ( int i = 0; i < rank_ ; i++ )
  //{
  //  R(i,i) = rho_;
  //}

  //System::out() << "Get mass matrix ...\n";
  idx_t          iMat;
  Ref<Material>  eMat;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount_; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems_[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem    );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    ipWeights *= thickness_;

    //EJ: use the element index and not the group index.
    iMat    = elemMatMap_[ielem];
    eMat    = materials_[iMat];

    rho = eMat->giveRho ();

    // Assemble the element mass matrix

    elemMat   = 0.0;

    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      // compute matrix of shpae function N

      getShapeFuncs_ ( N, sfuncs(ALL,ip) );

      // Add the contribution of this integration point.

      elemMat   += ipWeights[ip] * rho * mc2.matmul ( Nt, N );
      //elemMat   += ipWeights[ip] * mc3.matmul ( Nt, R, N );
    }

    // Add the element secant matrix to the global stiffness matrix.

    mbuilder.addBlock ( idofs, idofs, elemMat );

  }
  System::out() << "Get mass matrix ...done\n";
}

//-----------------------------------------------------------------------
//   computeEqvElemSize_
//-----------------------------------------------------------------------

void ElasticityModel::computeVolume_ () const
{
  Cubix       grads      ( rank_, nodeCount_, ipCount_ );
  Matrix      coords     ( rank_, nodeCount_ );

  IdxVector   inodes     ( nodeCount_ );
  Vector      ipWeights  ( ipCount_   );

  double      totalArea(0.);
  double      area;
  double      minArea(1e6);

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount_; ie++ )
  {
    int  ielem = ielems_[ie];

    elems_.getElemNodes  ( inodes, ielem    );
    nodes_.getSomeCoords ( coords, inodes );
    shape_->getShapeGradients ( grads, ipWeights, coords );

    area = sum ( ipWeights );
    totalArea += area;
    minArea = min (  minArea, area );
  }

  double he = ( rank_ == 2 ) ? sqrt(minArea) : cbrt(minArea);

  System::out() << "The area of          " << myName_ << " is: " << totalArea <<"\n";
  System::out() << "Minimum area         " << myName_ << " is: " << minArea   <<"\n";
  System::out() << "Minimum element size " << myName_ << " is: " << he        <<"\n";
}

//-----------------------------------------------------------------------
//   getTable_
//-----------------------------------------------------------------------


bool ElasticityModel::getTable_

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

  if ( name == "stress" &&
       table->getRowItems() == nodes_.getData() )
  {
    Vector  disp;

    StateVector::get ( disp, dofs_, globdat );
  
    getStress_ ( *table, weights, disp );

    return true;
  }

 
  
  if ( name == "damage" && table->getRowItems() == nodes_.getData() )
  {
    getHistory_ ( *table, weights );

    return true;
  }

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
//   getOutputData_
//-----------------------------------------------------------------------


void ElasticityModel::getOutputData_

  ( Ref<XTable>    table,
    const Vector&  weights,
    const String&  contents,
    const Vector&  disp )

{
  using jem::numeric::matmul;

  Cubix      grads      ( rank_, nodeCount_, ipCount_ );
  Matrix     coords     ( rank_, nodeCount_ );
  Vector     elemDisp   ( dofCount_ );
  Matrix     b          ( strCount_, dofCount_ );
  IdxVector  inodes     ( nodeCount_ );
  IdxVector  idofs      ( dofCount_ );
  Vector     ndWeights  ( nodeCount_ ); 
  Vector     ipWeights  ( ipCount_   );

  MChain2    mc2;

  Ref<Material>       eMat;
  idx_t               iMat;

  Matrix  sfuncs   = shape_->getShapeFunctions ();

  // table filler related stuff

  TbFiller   tbFiller   ( rank_ );

  Slice      iistrain   = tbFiller.announce ( "strain.tensor" ); // iistrain  = [0 1 2 3]
  Slice      iistress   = tbFiller.announce ( "stress.tensor" ); // iistrain  = [0 1 2 3]
  Slice      iipstress  = tbFiller.announce ( "pstress.diag"  ); // iipstress = [4 5 6]
  //Slice      iidamage   = tbFiller.announce ( "damage"        ); // [7]


  Vector     ipValues   ( tbFiller.typeCount() ); // typeCount() = # of types = 8 in 2D

  Vector     strain     ( ipValues[iistrain]   );
  Vector     stress     ( ipValues[iistress]   );
  Vector     pstress    ( ipValues[iipstress]  );
  //Vector     dam        ( ipValues[iidamage]   );


  // Let TbFiller find out which columns of ndValues to write to 
  // which columns of the table (based on filter in input file)

  IdxVector  i2table;
  IdxVector  jcols;

  tbFiller . setFilter   ( contents ); 
  tbFiller . prepareTable( i2table, jcols, table ); 
   
  Matrix     ndValuesOut ( nodeCount_, i2table.size() );
  Vector     ipValuesOut ( i2table.size() );

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount_; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems_[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem  );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Get the displacements in the nodes of this element.

    elemDisp = select ( disp, idofs );

    // get the correct material

    //EJ: use the element index and not the group index.
    iMat       = elemMatMap_[ielem];
    eMat       = materials_[iMat];

    ndValuesOut = 0.;
    ndWeights   = 0.;

    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      Vector Nip ( sfuncs(ALL,ip) );

      getShapeGrads_ ( b, grads(ALL,ALL,ip) );
      matmul ( strain, b, elemDisp );
      eMat->update ( stress, stiff0_, strain, 0  );
      computeEigenValues2D ( pstress, stress );

      //sfuncs(ALL,ip) = abs ( sfuncs(ALL,ip) );

      // apply the filter now, only cols specified by i2table are 
      // written to the table 
      ipValuesOut  = ipValues[i2table]; 
      ndValuesOut += matmul ( Nip, ipValuesOut );
      ndWeights   += Nip;
    }

    select ( weights, inodes ) += ndWeights;
    table->addBlock ( inodes, jcols, ndValuesOut );
  }
}



//-----------------------------------------------------------------------
//   getHistory
//-----------------------------------------------------------------------


void ElasticityModel::getHistory_

  ( XTable&        table,
    const Vector&  weights )

{
  Matrix     ndStrain   ( nodeCount_, 1 );

  ndStrain = 0.;

  IdxVector  inodes     ( nodeCount_ );
  IdxVector  jcols      ( 1          );

  // Add the columns for the strain components to the table.

 jcols[0] = table.addColumn ( "damage" );

  for ( int ie = 0; ie < ielemCount_; ie++ )
  {
    int  ielem = ielems_[ie];
    elems_.getElemNodes  ( inodes, ielem  );
    table.addBlock ( inodes, jcols, ndStrain );
  }
}


//-----------------------------------------------------------------------
//   getStress_
//-----------------------------------------------------------------------


void ElasticityModel::getStress_

  ( XTable&        table,
    const Vector&  weights,
    const Vector&  disp )

{
  using jem::numeric::matmul;

  Matrix     sfuncs     = shape_->getShapeFunctions ();

  Cubix      grads      ( rank_, nodeCount_, ipCount_ );


  Matrix     ndStress   ( nodeCount_, strCount_ );
  Vector     ndWeights  ( nodeCount_ );

  Matrix     coords     ( rank_, nodeCount_ );
  Matrix     gcoords    ( rank_, ipCount_   );

  Vector     elemDisp   ( dofCount_ );
  Vector     ipWeights  ( ipCount_   );
  

  Vector     stressIp   ( strCount_ );
  Vector     strainIp   ( strCount_ );

  Matrix     b          ( strCount_, dofCount_ );

  IdxVector  inodes     ( nodeCount_ );
  IdxVector  idofs      ( dofCount_ );
  IdxVector  jcols      ( strCount_  );

  MChain2    mc2;

  Ref<Material>       eMat;
  idx_t               iMat, ig(0);

  // Add the columns for the stress components to the table.

  switch ( strCount_ )
  {
  case 1:

    jcols[0] = table.addColumn ( "stress_xx" );

    break;

  case 4:

    jcols[0] = table.addColumn ( "stress_xx" );
    jcols[1] = table.addColumn ( "stress_yy" );
    jcols[2] = table.addColumn ( "stress_zz" );
    jcols[3] = table.addColumn ( "stress_xy" );

    break;

  case 6:

    jcols[0] = table.addColumn ( "stress_xx" );
    jcols[1] = table.addColumn ( "stress_yy" );
    jcols[2] = table.addColumn ( "stress_zz" );
    jcols[3] = table.addColumn ( "stress_xy" );
    jcols[4] = table.addColumn ( "stress_yz" );
    jcols[5] = table.addColumn ( "stress_xz" );

    break;

  default:

    throw Error (
      JEM_FUNC,
      "unexpected number of stress components: " +
      String ( strCount_ )
    );
  }


  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems_[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem  );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    shape_->getShapeGradients  ( grads, ipWeights, coords );

    // Get the displacements in the nodes of this element.

    elemDisp = select ( disp, idofs );

    // get the correct material

    //EJ: use the element index and not the group index.
    iMat    = elemMatMap_[ielem];
    eMat    = materials_[iMat];

    ndStress  = 0.0;
    ndWeights = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++, ig++ )
    {
      getShapeGrads_ ( b, grads(ALL,ALL,ip) );

      matmul         ( strainIp, b, elemDisp );

      eMat->update ( stressIp, stiff0_, strainIp, 0  );


      ndStress  += matmul ( sfuncs(ALL,ip), stressIp );
      ndWeights += sfuncs(ALL,ip);
    }

    select ( weights, inodes ) += ndWeights;

    // Add the stresses to the table.

    table.addBlock ( inodes, jcols, ndStress );

  }

}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newESolidModel
//-----------------------------------------------------------------------


static Ref<Model>     newElasticityModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<ElasticityModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareSolidModel
//-----------------------------------------------------------------------


void declareElasticityModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "Elasticity", & newElasticityModel );
}




