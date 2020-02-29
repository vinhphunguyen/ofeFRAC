
#include <jem/base/limits.h>
#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/base/ClassTemplate.h>
#include <jem/base/IllegalInputException.h>
#include <jem/io/ObjectInput.h>
#include <jem/io/ObjectOutput.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/util/Event.h>
#include <jive/util/error.h>
#include <jive/util/utilities.h>
#include <jive/util/Table.h>
#include <jive/util/TableUtils.h>
#include <jive/util/ItemSet.h>
#include <jive/util/DofSpace.h>
#include <jive/util/Constraints.h>
#include <jive/util/VectorManager.h>
#include <jive/algebra/VectorSpace.h>
#include <jive/model/Actions.h>
#include <jive/model/DummyModel.h>
#include <jive/model/StateVector.h>
#include <jive/model/ModelFactory.h>
#include <jive/implict/ArclenActions.h>

#include "UserModels.h"
#include "ArclenModel.h"

extern "C"
{
  #include <math.h>
}


using jive::IntVector;
using jive::util::VectorManager;
using jive::model::Actions;
using jive::model::ActionParams;
using jive::model::StateVector;
using jive::implict::ArclenActions;
using jive::implict::ArclenParams;


//=======================================================================
//   class ArclenModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  ArclenModel::TYPE_NAME       = "Arclen";

const char*  ArclenModel::MODEL_PROP      = "model";
const char*  ArclenModel::MAX_ITER_PROP   = "maxIter";
const char*  ArclenModel::MIN_INCR_PROP   = "minIncr";
const char*  ArclenModel::MAX_INCR_PROP   = "maxIncr";
const char*  ArclenModel::ARC_LENGTH_PROP = "arcLength";
const char*  ArclenModel::LOAD_SCALE_PROP = "loadScale";
const char*  ArclenModel::WGT_TABLE_PROP  = "weightTable";


const char*  ArclenModel::DELTA_STATE_    = "deltaState0";

const int    ArclenModel::U_LOAD_         = 1 << 0;
const int    ArclenModel::U_WEIGHTS_      = 1 << 1;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


ArclenModel::ArclenModel

  ( const String&      name,
    const Ref<Model>&  child ) :

    Super  ( name  ),
    child_ ( child )

{
  using jive::model::DummyModel;

  if ( child_ == NIL )
  {
    child_ = newInstance<DummyModel> ( name );
  }

  istep_     = 0;
  updated_   = 0;
                                // or IDC (indirect displacement control)
  maxIter_   = 4;
  minIncr_   = 1.0e-3;
  maxIncr_   = 1.0e+1;
  loadScale_ = 0.0;
  arcLength_ = 0.0;
}


ArclenModel::~ArclenModel ()
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool ArclenModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  if ( dofs_ == NIL )
  {
    if ( action == Actions::INIT )
    {
      child_->takeAction ( action, params, globdat );

      init_ ( globdat );

      return true;
    }
    else
    {
      return child_->takeAction ( action, params, globdat );
    }
  }

  if ( action == ArclenActions::GET_ARC_FUNC )
  {
    evalArcFunc_ ( params, globdat );

    return true;
  }

  if ( action == ArclenActions::GET_UNIT_LOAD )
  {
    getUnitLoad_ ( params, globdat );

    return true;
  }

  if ( action == Actions::COMMIT )
  {
    child_->takeAction ( action, params, globdat );

    commit_ ( params, globdat );

    return true;
  }

  return child_->takeAction ( action, params, globdat );
}


//-----------------------------------------------------------------------
//   findModel
//-----------------------------------------------------------------------


Model* ArclenModel::findModel ( const String& name ) const
{
  if ( name == myName_ )
  {
    return const_cast<Self*> ( this );
  }
  else
  {
    return child_->findModel ( name );
  }
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void ArclenModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  using jem::maxOf;

  if ( props.contains( myName_ ) )
  {
    Properties  myProps = props.findProps ( myName_ );

    myProps.find ( maxIter_,   MAX_ITER_PROP, 1,          1000 );
    myProps.find ( minIncr_,   MIN_INCR_PROP, 0.0,        maxOf( minIncr_ ) );
    myProps.find ( maxIncr_,   MAX_INCR_PROP, minIncr_,   maxOf( maxIncr_ ) );
    myProps.get  ( arcLength_, ARC_LENGTH_PROP );
    myProps.find ( loadScale_, LOAD_SCALE_PROP );

    if ( myProps.find( wtblName_, WGT_TABLE_PROP ) )
    {
      updated_ &= ~U_WEIGHTS_;
    }
  }

  child_->configure ( props, globdat );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void ArclenModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );

  myConf.set ( MAX_ITER_PROP,   maxIter_   );
  myConf.set ( MIN_INCR_PROP,   minIncr_   );
  myConf.set ( MAX_INCR_PROP,   maxIncr_   );
  myConf.set ( ARC_LENGTH_PROP, arcLength_ );
  myConf.set ( LOAD_SCALE_PROP, loadScale_ );
  myConf.set ( WGT_TABLE_PROP,  wtblName_  );

  child_->getConfig ( conf, globdat );
}

//-----------------------------------------------------------------------
//   setMaxIter
//-----------------------------------------------------------------------


void ArclenModel::setMaxIter ( int count )
{
  JEM_PRECHECK ( count > 0 );

  maxIter_ = count;
}


//-----------------------------------------------------------------------
//   setLoadIncr
//-----------------------------------------------------------------------


void ArclenModel::setLoadIncr ( double incr )
{
  //loadIncr_ = incr;

  if ( istep_ > 0 )
  {
    istep_ = 0;
  }
}


//-----------------------------------------------------------------------
//   setLoadScale
//-----------------------------------------------------------------------


void ArclenModel::setLoadScale ( double scale )
{
  loadScale_ = scale;

  if ( istep_ > 0 )
  {
    istep_ = 0;
  }
}


//-----------------------------------------------------------------------
//   setIncrRange
//-----------------------------------------------------------------------


void ArclenModel::setIncrRange

  ( double  minIncr,
    double  maxIncr )

{
  JEM_PRECHECK ( minIncr >= 0.0 && minIncr < maxIncr );

  minIncr_ = minIncr;
  maxIncr_ = maxIncr;
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Model> ArclenModel::makeNew

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  using jive::util::joinNames;
  using jive::model::newModel;

  Ref<Model>  child =

    newModel ( joinNames( name, MODEL_PROP ),
	       conf, props, globdat );

  return newInstance<Self> ( name, child );
}


//-----------------------------------------------------------------------
//   init_
//-----------------------------------------------------------------------


void ArclenModel::init_ ( const Properties& globdat )
{
  dofs_    = DofSpace    :: get ( globdat, getContext() );
  cons_    = Constraints :: get ( dofs_,   globdat );
  istep_   = 0;
  updated_ = 0;

  connect_ ();
}


//-----------------------------------------------------------------------
//   initLoad_
//-----------------------------------------------------------------------


void ArclenModel::initLoad_ ( const Properties& globdat )
{
  Properties  params;

  load_.resize ( dofs_->dofCount() );

  load_ = 0.0;

  params.set ( ActionParams::EXT_VECTOR, load_ );

  child_->takeAction  ( Actions::GET_EXT_VECTOR, params, globdat );
  dofs_ ->resetEvents ();

  updated_ |= U_LOAD_;
}


//-----------------------------------------------------------------------
//   initWeights_
//-----------------------------------------------------------------------


void ArclenModel::initWeights_ ( const Properties& globdat )
{
  using jem::select;
  using jive::util::zeroSlaveDofs;
  using jive::util::Table;
  using jive::util::TableUtils;


  weights_.resize ( dofs_->dofCount() );

  if ( wtblName_.size() == 0 )
  {
    weights_ = 1.0;
  }
  else
  {
    Ref<Table>  table  =

      Table::get ( wtblName_, dofs_->getItems(),
                   globdat,   getContext() );

    const int   nnz    = table->size ();

    Array<idx_t>   idofs  ( nnz );
    Vector         values ( nnz );


    TableUtils::toVector ( values, idofs, *table, *dofs_ );

    weights_ = 0.0;
    select ( weights_, idofs ) = values;
  }

  zeroSlaveDofs ( weights_, *cons_ );

  dofs_->resetEvents ();
  cons_->resetEvents ();

  updated_ |= U_WEIGHTS_;
}


//-----------------------------------------------------------------------
//   evalArcFunc_
//-----------------------------------------------------------------------


void ArclenModel::evalArcFunc_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jem::max;
  using jem::System;

  double  fvalue;
  double  lambda;
  double  jac11;

  //int     iiter;

  //params.get ( iiter, ArclenParams::IITER );

  if ( vspace_ == NIL )
  {
    vspace_ = VectorSpace::get ( dofs_, globdat );
  }

  if ( ! (updated_ & U_LOAD_) )
  {
    initLoad_    ( globdat );
  }

  if ( ! (updated_ & U_WEIGHTS_) )
  {
    initWeights_ ( globdat );
  }

  vtmp_.resize ( weights_.size() );

  Vector  jac10;
  Vector  u, u0;

  params      .get    ( jac10, ArclenParams::JACOBIAN10 );
  StateVector::get    ( u,  dofs_, globdat );
  StateVector::getOld ( u0, dofs_, globdat );

  vtmp_ = u - u0;

  lambda = vspace_->product ( weights_, vtmp_ );
  jac10  = weights_;

  if ( lambda < 0.0 )
  {
    lambda = -lambda;
    jac10 *= -1.0;
  }

  fvalue = lambda - arcLength_;
  jac11  = 0.0;

  params.set ( ArclenParams::JACOBIAN11, jac11  );
  params.set ( ArclenParams::ARC_FUNC,   fvalue );
}


//-----------------------------------------------------------------------
//   getUnitLoad_
//-----------------------------------------------------------------------


void ArclenModel::getUnitLoad_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::util::sizeError;

  Vector  f;

  if ( ! (updated_ & U_LOAD_) )
  {
    initLoad_ ( globdat );
  }

  params.get ( f, ArclenParams::UNIT_LOAD );

  if ( f.size() != load_.size() )
  {
    sizeError ( getContext(),
		"unit load vector", f.size(), load_.size() );
  }

  f = load_;
}


//-----------------------------------------------------------------------
//   commit_
//-----------------------------------------------------------------------


void ArclenModel::commit_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jem::maxOf;
  using jem::numeric::axpy;

  Vector  u, u0, du;

  double  lambda;
  int     iiter;

  params.get ( iiter,  ArclenParams::IITER      );
  params.get ( lambda, ArclenParams::LOAD_SCALE );

  //loadIncr_  = lambda - loadScale_;
  //loadScale_ = lambda;

  // Adjust the load increment according to the number of iterations.

  /*
  if ( iiter > 0 )
  {
    else
    {
      loadIncr_ *= (double) maxIter_ / (double) iiter;
    }
    double n    = ( iiter - optIter_ ) / 4.0;
    loadIncr_ *= ::pow ( 0.5, n );
  }
  */

  StateVector  ::get       ( u,                dofs_, globdat );
  StateVector  ::getOld    ( u0,               dofs_, globdat );
  VectorManager::getVector ( du, DELTA_STATE_, dofs_, globdat );

  axpy ( du, u, -1.0, u0 );

  if ( istep_ < maxOf( istep_ ) )
  {
    istep_++;
  }
}


//-----------------------------------------------------------------------
//   connect_
//-----------------------------------------------------------------------


void ArclenModel::connect_ ()
{
  using jem::util::connect;

  connect ( dofs_->newSizeEvent,   this, & Self::dofsChanged_ );
  connect ( dofs_->newOrderEvent,  this, & Self::dofsChanged_ );
  connect ( cons_->newStructEvent, this, & Self::consChanged_ );

  dofs_->resetEvents ();
  cons_->resetEvents ();
}


//-----------------------------------------------------------------------
//   dofsChanged_
//-----------------------------------------------------------------------


void ArclenModel::dofsChanged_ ()
{
  updated_ = 0;
}


//-----------------------------------------------------------------------
//   consChanged_
//-----------------------------------------------------------------------


void ArclenModel::consChanged_ ()
{
  updated_ &= ~U_WEIGHTS_;
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareArclenModel
//-----------------------------------------------------------------------


void declareArclenModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( ArclenModel::TYPE_NAME,
			  & ArclenModel::makeNew );
}
