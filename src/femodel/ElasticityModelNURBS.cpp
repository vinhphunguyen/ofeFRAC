/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the linear elasticity model
 *  using finite elements. 
 *  The number of elements and the integration scheme/element
 *  is fixed.
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
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
#include <jive/geom/StdSquare.h>
#include <jive/geom/InternalShape.h>
#include <jive/geom/ParametricArea.h>
#include <jive/geom/Geometries.h>
#include <jive/fem/ElementGroup.h>

#include "femodel/FemModels.h"
#include "util/utilities.h"
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



class ElasticityModelNURBS : public Model
{
 public:

  typedef ElasticityModelNURBS    Self;
  typedef Model              Super;

  static const char*         DOF_NAMES[3];
  static const char*         SHAPE_PROP;
  static const char*         MATERIAL_PROP;

                       ElasticityModelNURBS

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

  virtual              ~ElasticityModelNURBS ();

 private:

  void                 getMatrix_

    ( MatrixBuilder&      mbuilder,
      const Vector&       force,
      const Vector&       disp )       const;

  bool                 getTable_

    ( const Properties&   params,
      const Properties&   globdat );

  void                 getHistory_

    ( XTable&             table,
      const Vector&       weights );

  void                 getStress_

    ( XTable&                 table,
      const Vector&           weights,
      const Vector&           disp  );

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

  int                     outNode_;
  Ref<PrintWriter>        out_;
};

// -------------------------------------------------------------
//    static constants
// -------------------------------------------------------------

const char* ElasticityModelNURBS::DOF_NAMES[3]  = {"dx","dy","dz"};
const char* ElasticityModelNURBS::SHAPE_PROP    = "shape";
const char* ElasticityModelNURBS::MATERIAL_PROP = "material";

// -------------------------------------------------------------
//    constructor
// -------------------------------------------------------------

ElasticityModelNURBS::ElasticityModelNURBS

   ( const String&       name,
     const Properties&   conf,
     const Properties&   props,
     const Properties&   globdat ) : Super(name)
{
  using jive::util::joinNames;
  using namespace jive::geom;

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
  
  Properties  shapeProps = myProps.getProps ( "shape" );
  String      shapeName;

  shapeProps.get( shapeName, "type" );

  if      ( shapeName == "9B" )
  {
    Ref<StdShape>  sshape = newInstance<StdSquare9B> ();

    shape_ = newInstance<ParametricArea> (
                   "quad",
                    StdSquare::getGaussScheme ( 3, 3 ),
                    sshape
                  );
  }
  //else if ( shapeName == "16B" )
  //{
  //  Ref<StdShape>  sshape = newInstance<StdSquare16B> ();

  //  shape_ = newInstance<ParametricArea> (
  //                 "quad",
  //                  StdSquare::getGaussScheme ( 4, 4 ),
  //                  sshape
  //                );
  //}
  //else if ( shapeName == "3x2B" )
  //{
  //  Ref<StdShape>  sshape = newInstance<StdSquare3x2B> ();

  //  shape_ = newInstance<ParametricArea> (
  //                 "quad",
  //                  StdSquare::getGaussScheme ( 3, 2 ),
  //                  sshape
  //                );
  //}
  //else if ( shapeName == "4x2B" )
  //{
  //  Ref<StdShape>  sshape = newInstance<StdSquare4x2B> ();

  //  shape_ = newInstance<ParametricArea> (
  //                 "quad",
  //                  StdSquare::getGaussScheme ( 4, 2 ),
  //                  sshape
  //                );
  //}
  //else if ( shapeName == "4x3B" )
  //{
  //  Ref<StdShape>  sshape = newInstance<StdSquare4x3B> ();

  //  shape_ = newInstance<ParametricArea> (
  //                 "quad",
  //                  StdSquare::getGaussScheme ( 4, 3 ),
  //                  sshape
  //                );
  //}
  //else if ( shapeName == "5x2B" )
  //{
  //  Ref<StdShape>  sshape = newInstance<StdSquare5x2B> ();

  //  shape_ = newInstance<ParametricArea> (
  //                 "quad",
  //                  StdSquare::getGaussScheme ( 5, 2 ),
  //                  sshape
  //                );
  //}
  //else if ( shapeName == "5x3B" )
  //{
  //  Ref<StdShape>  sshape = newInstance<StdSquare5x3B> ();

  //  shape_ = newInstance<ParametricArea> (
  //                 "quad",
  //                  StdSquare::getGaussScheme ( 5, 3 ),
  //                  sshape
  //                );
  //}
  else
  {
    throw IllegalInputException (
      context,
      String::format ( "unsupported continuum shape object")
    );
  }


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

  dofs_->addDofs (
    elems_.getUniqueNodesOf ( ielems_ ),
    dofTypes_
  );

  // Create a material model object.

  Ref<Material>     tmpMat = newMaterial ( MATERIAL_PROP, myConf,
                                            myProps, globdat );

  Ref<HookeMaterial>       isotrMa = dynamicCast<HookeMaterial>       ( tmpMat );
  Ref<OrthotropicMaterial> orthoMa = dynamicCast<OrthotropicMaterial> ( tmpMat );

  Properties  matProps = myProps.findProps ( MATERIAL_PROP );
  Properties  matConf  = myConf.makeProps  ( MATERIAL_PROP );


  stiff0_.resize ( strCount_, strCount_ );
  
 if ( isotrMa != NIL)
 {
  isotrMa->configure    ( matProps, globdat );
  isotrMa->getConfig    ( matConf , globdat );

  // get the elastic modulii matrix once

  stiff0_ = isotrMa->getStiffMat ( ); 
 }
 else if ( orthoMa != NIL)
 {
  orthoMa->configure    ( matProps, globdat );
  orthoMa->getConfig    ( matConf , globdat );

  stiff0_ = orthoMa->getStiffMat ( ); 
 }

  // correct function to compute B matrix

  getShapeGrads_  = getShapeGradsFunc        ( rank_ );

  outNode_ = -1;

  int node;

  if ( myProps.find ( node, "node" ) )
  {
    outNode_ = nodes_.findNode ( node );

    out_     = newInstance<PrintWriter> (
		  newInstance<FileWriter> ( "stress-disp.dat" ) );
  }

}

ElasticityModelNURBS::~ElasticityModelNURBS()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void ElasticityModelNURBS::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void ElasticityModelNURBS::getConfig ( const Properties& conf ) const
{
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool ElasticityModelNURBS::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;


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

    return true;
  }


  if ( action == Actions::GET_TABLE )
  {
    return getTable_ ( params, globdat );
  }

  return false;
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void ElasticityModelNURBS::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   disp ) const

{
  using jem::numeric::matmul;

  Cubix       grads      ( rank_, nodeCount_, ipCount_ );

  Matrix      coords     ( rank_, nodeCount_ );

  Matrix      elemMat    ( dofCount_, dofCount_  );
  Vector      elemForce  ( dofCount_ );
  Vector      elemDisp   ( dofCount_ );

  Vector      strain     ( strCount_ );
  Vector      stress     ( strCount_ );

  Matrix      b          ( strCount_, dofCount_  );
  Matrix      bt         = b.transpose ();

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );

  Vector      ipWeights  ( ipCount_   );

  MChain1     mc1;
  MChain3     mc3;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount_; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems_[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem             );
    nodes_.getSomeCoords ( coords, inodes            );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    //System::out() << "solid nodes " << inodes << "\n";
    // Compute the spatial derivatives of the element shapoe
    // functions.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Get the displacements at the element nodes.

    elemDisp = select ( disp, idofs );

    // Assemble the element matrix.

    elemMat   = 0.0;
    elemForce = 0.0;

    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      // Compute the B-matrix for this integration point.

      getShapeGrads_ ( b, grads(ALL,ALL,ip) ); 

      // Compute the strain vector of this integration point

      matmul ( strain, b, elemDisp     ); 

      matmul ( stress, stiff0_, strain );

      // Compute the element force vector.
      // Compute the stiffness matrix

      double wip = ipWeights[ip];

      elemForce += wip * mc1.matmul ( bt, stress     );
      elemMat   += wip * mc3.matmul ( bt, stiff0_, b );
    }

    // Add the element matrix to the global stiffness matrix.

    mbuilder.addBlock ( idofs, idofs, elemMat );

    // Add the element force vector to the global force vector.

    select ( force, idofs ) += elemForce;
  }
}

//-----------------------------------------------------------------------
//   getTable_
//-----------------------------------------------------------------------


bool ElasticityModelNURBS::getTable_

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

    
    if ( outNode_ != -1 )
    {

      int ndof  = dofs_->typeCount();

      Vector       state;
      IdxVector    jtypes ( ndof );
      int          dof;

      StateVector::get ( state, dofs_, globdat );

      for ( int i = 0; i < ndof; i++ ) jtypes[i] = i;

      for ( int i = 0; i < ndof; i++ ) 
      {
        dof = dofs_->findDofIndex ( outNode_, jtypes[i] );

        print ( *out_, state[dof] ); 
        print ( *out_, "   " );
      }

      double  sigmaXX = table->getValue ( outNode_, 0 );

      print ( *out_, sigmaXX );
      out_->printLine();
      out_->flush();
    }

    return true;
  }

  return false;
}

//-----------------------------------------------------------------------
//   getStress_
//-----------------------------------------------------------------------


void ElasticityModelNURBS::getStress_

  ( XTable&        table,
    const Vector&  weights,
    const Vector&  disp )

{
  using jem::numeric::matmul;

  Cubix      grads      ( rank_, nodeCount_, ipCount_ );

  Matrix     ndStress   ( nodeCount_, strCount_ );
  Vector     ndWeights  ( nodeCount_ );
  Vector     stressIp   ( strCount_ );
  Vector     ipWeights  ( ipCount_   );

  Matrix     coords     ( rank_, nodeCount_ );
  Vector     elemDisp   ( dofCount_ );

  Matrix     b          ( strCount_, dofCount_ );

  IdxVector  inodes     ( nodeCount_ );
  IdxVector  idofs      ( dofCount_ );
  IdxVector  jcols      ( strCount_ );

  MChain2    mc2;


  // Add the columns for the stress components to the table.

  switch ( strCount_ )
  {
  case 1:

    jcols[0] = table.addColumn ( "sxx" );

    break;

  case 3:

    jcols[0] = table.addColumn ( "sxx" );
    jcols[1] = table.addColumn ( "syy" );
    jcols[2] = table.addColumn ( "sxy" );

    break;

  case 6:

    jcols[0] = table.addColumn ( "sxx" );
    jcols[1] = table.addColumn ( "syy" );
    jcols[2] = table.addColumn ( "szz" );
    jcols[3] = table.addColumn ( "sxy" );
    jcols[4] = table.addColumn ( "syz" );
    jcols[5] = table.addColumn ( "sxz" );

    break;

  default:

    throw Error (
      JEM_FUNC,
      "unexpected number of stress components: " +
      String ( strCount_ )
    );
  }

  Matrix  sfuncs   = shape_->getShapeFunctions ();

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

    ndStress  = 0.0;
    ndWeights = 0.0;

    for ( int ip = 0; ip < ipCount_; ip++ )
    {
      getShapeGrads_ ( b, grads(ALL,ALL,ip) );

      stressIp   = mc2.matmul ( stiff0_, b, elemDisp );

      sfuncs(ALL,ip) = abs ( sfuncs(ALL,ip) );

      ndStress  += matmul ( sfuncs(ALL,ip), stressIp );
      ndWeights += sfuncs(ALL,ip) ;
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


static Ref<Model>     newElasticityModelNURBS

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<ElasticityModelNURBS> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareSolidModel
//-----------------------------------------------------------------------


void declareElasticityModelNURBS ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "ElasticityNURBS", & newElasticityModelNURBS );
}



