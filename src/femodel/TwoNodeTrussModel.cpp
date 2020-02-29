/*
 *
 *  Copyright (C) 2017 Monash University. All rights reserved.
 *
 *  This class implements a mechanical FE model of two-node truss elements. It is used
 *  with a 1D J2 plasticity material to model reinforcement steel bars.
 *  It will be used with existing continuum models such as VariationalDamMechModel and
 *  VariationalDamPhaseModel to model cracking of reinforced concrete structures.
 *
 *  Author: V.P. Nguyen, phu.nguyen@monash.edu
 *  Date: 13 November 2017
 *
 *  Updates (what and who):
 *
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
#include "material/OneDimJ2PlasticityMat.h"

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



class TwoNodeTrussModel : public Model
{
 public:

  typedef TwoNodeTrussModel    Self;
  typedef Model              Super;

  static const char*         DOF_NAMES[3];
  static const char*         MATERIAL_PROP;

                       TwoNodeTrussModel

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

  virtual              ~TwoNodeTrussModel ();

 private:

  void                 getMatrix_

    ( MatrixBuilder&      mbuilder,
      const Vector&       force,
      const Vector&       disp )       const;

  void                 computeInitialLengths_ ( );

 private:

  Assignable<ElemGroup>   egroup_;
  Assignable<ElemSet>     elems_;
  Assignable<NodeSet>     nodes_;

  Ref<XDofSpace>          dofs_;
  IdxVector               dofTypes_;

  //ShapeGradsFunc          getShapeGrads_;
  //ShapeFunc               getShapeFuncs_;

  IdxVector               ielems_;

  idx_t                   rank_;
  idx_t                   ielemCount_;
  idx_t                   nodeCount_;
  //idx_t                   ipCount_;
  //idx_t                   strCount_;
  idx_t                   dofCount_;

  double                  area_;       // area
  Vector                  lengths0_;   // initial lengths of all elements

  Ref<OneDimJ2PlasticityMat >         mat_;
};

// -------------------------------------------------------------
//    static constants
// -------------------------------------------------------------

const char*  TwoNodeTrussModel::DOF_NAMES[3] = { "dx", "dy", "dz" };
const char*  TwoNodeTrussModel::MATERIAL_PROP = "material";

// -------------------------------------------------------------
//    constructor
// -------------------------------------------------------------

TwoNodeTrussModel::TwoNodeTrussModel

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

  // compute some constants

  nodeCount_  = 2;
  dofCount_   = rank_ * nodeCount_;

  // Make sure that each element has the same number of nodes as the
  // shape object.

  elems_.checkSomeElements ( context, ielems_, nodeCount_ );

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  idx_t typeCount = dofs_->typeCount ( );

  System::out() << "dof type count " << typeCount << "\n";

  if ( typeCount == 0 )
  {
    dofTypes_.resize( rank_ );
    
    for ( idx_t i = 0; i < rank_; i++ )
    {
      dofTypes_[i] = dofs_->addType ( DOF_NAMES[i] );
    }
      
    dofs_->addDofs ( elems_.getUniqueNodesOf ( ielems_ ), dofTypes_ );
  }
  else
  {
    dofTypes_.resize( rank_ );
    for ( int i = 0; i < rank_; i++ )
    {
      dofTypes_[i] = i; 
    }
  }

  // ------------------------------------------------------
  //       create material objects
  // ------------------------------------------------------

  Properties    matProps = myProps.findProps ( MATERIAL_PROP );
  Properties    matConf  = myConf.makeProps  ( MATERIAL_PROP );

  Ref<OneDimJ2PlasticityMat> tmpMat   = newInstance <OneDimJ2PlasticityMat> ( matConf, matProps,  globdat );

  mat_ = tmpMat;

  mat_->allocPoints ( ielemCount_ );

  area_ = 1;
  myProps.find ( area_ , "area" );
  myConf. set  ( "area", area_  );

  computeInitialLengths_ ( );
}

TwoNodeTrussModel::~TwoNodeTrussModel()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void TwoNodeTrussModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void TwoNodeTrussModel::getConfig ( const Properties& conf ) const
{
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool TwoNodeTrussModel::takeAction

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

  if ( action == Actions::COMMIT )
  {
   //System::out() << myName_ << ": " << "computing stiffness matrix ... \n";

    mat_->commit ( );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void TwoNodeTrussModel::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   disp ) const

{
  using jem::numeric::matmul;

  Matrix      coords     ( rank_, nodeCount_ );

  Matrix      elemMat    ( dofCount_, dofCount_  );
  Vector      elemForce  ( dofCount_ );
  Vector      elemDisp   ( dofCount_ );

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount_  );

  double      l, l0, x1, x2, y1, y2, x21, y21;
  double      strain, stress, Cep, stiff, stre;
  double      c, s, c2, s2, cs;

  Properties  params;

  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems_[ie];

    //System::out() << ielem << "\n";

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem             );
    nodes_.getSomeCoords ( coords, inodes            );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );
    
    elemDisp = select ( disp, idofs );

    //System::out() << "solid nodes " << inodes << "\n";

    x1  = coords(0,0) ; y1  = coords(1,0) ; 
    x2  = coords(0,1) ; y2  = coords(1,1) ; 
    
    x21    = x2 + elemDisp[2] -  x1 - elemDisp[0]; 
    y21    = y2 + elemDisp[3] -  y1 - elemDisp[1]; 
    
    l0     = lengths0_[ie];
    l      = sqrt ( x21 * x21 + y21 * y21 );
    
    strain = ( l - l0 ) / l0;

    params.set ( "strain" , strain );
    params.set ( "i" ,      ie     );

    mat_->update ( params );
    
    params.get ( Cep ,     "tangent" );
    params.get ( stress ,  "stress"  );


    stiff  = Cep    * area_ / l0;
    stre   = stress * area_;

    c   = ( x2 - x1 ) / l0;
    s   = ( y2 - y1 ) / l0;
    c2  = c * c;
    s2  = s * s;
    cs  = c * s;

    elemMat(0,0) = c2;  elemMat(0,1) = cs;  elemMat(0,2) = -c2; elemMat(0,3) = -cs;
    elemMat(1,0) = cs;  elemMat(1,1) = s2;  elemMat(1,2) = -cs; elemMat(1,3) = -s2;
    elemMat(2,0) = -c2; elemMat(2,1) = -cs; elemMat(2,2) = c2;  elemMat(2,3) = cs;
    elemMat(3,0) = -cs; elemMat(3,1) = -s2; elemMat(3,2) = cs;  elemMat(3,3) = s2;

    elemMat *= stiff;

    //elemForce    = matmul ( elemMat, elemDisp ); 

    elemForce[0] = -c;
    elemForce[1] = -s;
    elemForce[2] =  c;
    elemForce[3] =  s;

    elemForce *= stre;

    //System::out() << "l0: " << elemMat << "\n"; 
    //System::out() << "l0: " << elemForce << "\n"; 
    //System::out() << "Cep: " << Cep << "\n"; 

    // Add the element matrix to the global stiffness matrix.
    // Add the element force vector to the global force vector.

    mbuilder.addBlock ( idofs, idofs, elemMat );

    select ( force, idofs ) += elemForce;
  }
}

//-----------------------------------------------------------------------
//   computeInitialLengths_
//-----------------------------------------------------------------------

void TwoNodeTrussModel::computeInitialLengths_ ( ) 
{
  Matrix      coords     ( rank_, nodeCount_ );

  IdxVector   inodes     ( nodeCount_ );

  double      x21, y21;

  lengths0_.resize ( ielemCount_ );

  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < ielemCount_; ie++ )
  {
    idx_t  ielem = ielems_[ie];

    elems_.getElemNodes  ( inodes, ielem    );
    nodes_.getSomeCoords ( coords, inodes   );

    x21  = coords(0,1) -  coords(0,0) ; 
    y21  = coords(1,1) -  coords(1,0) ; 

    lengths0_[ie] = sqrt ( x21 * x21 + y21 * y21 );
  }
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newTwoNodeTrussModel
//-----------------------------------------------------------------------


static Ref<Model>     newTwoNodeTrussModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<TwoNodeTrussModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareSolidModel
//-----------------------------------------------------------------------


void declareTwoNodeTrussModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "TwoNodeTruss", & newTwoNodeTrussModel );
}




