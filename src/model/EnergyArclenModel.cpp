/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implements the energy release control method
 *  (Miguel Gutierrez 2004). This solver is best suited to
 *  problem involving snapback behavior. It starts with load
 *  control and then switch to, automatically, based on the
 *  sudden change of iterations, energy based arclength control.
 *  
 *  This solver when used in combination with the StepModule
 *  is able to capture the whole complicated behavior.
 *
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *
 *  Improved by Frans van der Meers
 *
 */

#include <jem/base/limits.h>
#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/base/IllegalOperationException.h>
#include <jem/base/ClassTemplate.h>
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
#include <jive/util/Globdat.h>
#include <jive/algebra/VectorSpace.h>
#include <jive/model/Actions.h>
#include <jive/model/DummyModel.h>
#include <jive/model/StateVector.h>
#include <jive/model/ModelFactory.h>
#include <jive/implict/ArclenActions.h>

#include "UserModels.h"
#include "EnergyArclenModel.h"

extern "C"
{
  #include <math.h>
}

using jem::io::endl;
using jem::numeric::axpy;
using jive::IntVector;
using jive::util::VectorManager;
using jive::model::Actions;
using jive::model::ActionParams;
using jive::model::StateVector;
using jive::implict::ArclenActions;
using jive::implict::ArclenParams;


//=======================================================================
//   class EnergyArclenModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  EnergyArclenModel::TYPE_NAME       = "EnergyArclen";

const char*  EnergyArclenModel::MODEL_PROP      = "model";
const char*  EnergyArclenModel::ARC_FUNC_PROP   = "arcFunc";
const char*  EnergyArclenModel::OPT_ITER_PROP   = "optIter";
const char*  EnergyArclenModel::SWT_ITER_PROP   = "swtIter";
const char*  EnergyArclenModel::MIN_INCR_PROP   = "minIncr";
const char*  EnergyArclenModel::MAX_INCR_PROP   = "maxIncr";
const char*  EnergyArclenModel::LOAD_INCR_PROP  = "loadIncr";
const char*  EnergyArclenModel::LOAD_SCALE_PROP = "loadScale";
const char*  EnergyArclenModel::EXIT_FRAC_PROP  = "exitFraction";
const char*  EnergyArclenModel::WGT_TABLE_PROP  = "weightTable";
const char*  EnergyArclenModel::STEP_ADJUST_PROP = "stepAdjustment";


const char*  EnergyArclenModel::DELTA_STATE_    = "deltaState0";

const int    EnergyArclenModel::U_LOAD_         = 1 << 0;
const int    EnergyArclenModel::U_WEIGHTS_      = 1 << 1;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


EnergyArclenModel::EnergyArclenModel

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
  arcFunc_   = ERC;
  optIter_   = 4;
  minIncr_   = 1.0e-3;
  maxIncr_   = 1.0e+1;
  loadIncr_  = 1.0;        // load increment at the first iteration of a load step
  loadScale_ = 0.0;
  oldScale_  = 0.0;
  arcLength_ = 0.0;

  maxLoadScale_ = 0.0;
  exitFraction_ = 0.0;          // terminate computation when load approaches zero

  isLoadControl_ = true;        // starts with load control

  stepAdjust_    = "lastValue"; // other is smooth

  nformat_.setScientific     (   );
  nformat_.setFractionDigits ( 4 );

  down_          = false;
}


EnergyArclenModel::~EnergyArclenModel ()
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool EnergyArclenModel::takeAction

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

  if ( action == "CHECK_COMMIT" )
  {
    checkCommit_ ( params, globdat );

    return true;
  }

  if ( action == Actions::COMMIT )
  {
    child_->takeAction ( action, params, globdat );

    commit_ ( params, globdat );

    return true;
  }

  // this REDUCE_STEP action is thrown by StepModule

  if ( action == "REDUCE_STEP" )
  {
    reduceStep_ ( params, globdat );

    return true;
  }

  return child_->takeAction ( action, params, globdat );
}


//-----------------------------------------------------------------------
//   findModel
//-----------------------------------------------------------------------


Model* EnergyArclenModel::findModel ( const String& name ) const
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


void EnergyArclenModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  using jem::maxOf;

  if ( props.contains( myName_ ) )
  {
    Properties  myProps = props.findProps ( myName_ );

    String      fnType;


    if ( myProps.find( fnType, ARC_FUNC_PROP ) )
    {
      if      ( fnType == "ERC" )
      {
        arcFunc_ = ERC;
      }
      else
      {
        myProps.propertyError (
          ARC_FUNC_PROP,
          String::format (
            "invalid arc-length function type: %S; "
            "valid types are "
            "`ERC\' (Energy Release Control) ",
            & fnType
          )
        );
      }
    }

    myProps.find ( optIter_,   OPT_ITER_PROP,
                   1,          1000 );

    myProps.find ( swtIter_,   SWT_ITER_PROP,
                   1,          1000 );

    myProps.find ( minIncr_,   MIN_INCR_PROP,
                   0.0,        maxOf( minIncr_ ) );

    myProps.find ( maxIncr_,   MAX_INCR_PROP,
                   minIncr_,   maxOf( maxIncr_ ) );

    myProps.find ( loadIncr_,  LOAD_INCR_PROP );

    myProps.find ( loadScale_, LOAD_SCALE_PROP );

    myProps.find ( exitFraction_, EXIT_FRAC_PROP,
                   0.0,        1.0 );

    myProps.find ( stepAdjust_, STEP_ADJUST_PROP );

    myProps.find ( switchBack_, "switchBack" );

    if ( myProps.find( wtblName_, WGT_TABLE_PROP ) )
    {
      updated_ &= ~U_WEIGHTS_;
    }
  }

  arcLength_ = loadIncr_;

  child_->configure ( props, globdat );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void EnergyArclenModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );

  String      fnType;

  if ( arcFunc_ == ERC )
  {
    fnType = "ERC";
  }
  else
  {
    fnType = "unknown";
  }

  myConf.set ( ARC_FUNC_PROP,    fnType        );
  myConf.set ( OPT_ITER_PROP,    optIter_      );
  myConf.set ( SWT_ITER_PROP,    swtIter_      );
  myConf.set ( MIN_INCR_PROP,    minIncr_      );
  myConf.set ( MAX_INCR_PROP,    maxIncr_      );
  myConf.set ( LOAD_INCR_PROP,   loadIncr_     );
  myConf.set ( LOAD_SCALE_PROP, loadScale_     );
  myConf.set ( EXIT_FRAC_PROP,   exitFraction_ );
  myConf.set ( WGT_TABLE_PROP,   wtblName_     );
  myConf.set ( STEP_ADJUST_PROP, stepAdjust_   );

  child_->getConfig ( conf, globdat );
}


//-----------------------------------------------------------------------
//   setArcFunc
//-----------------------------------------------------------------------


void EnergyArclenModel::setArcFunc ( ArcFunc func )
{
  arcFunc_ = func;
}


//-----------------------------------------------------------------------
//   setMaxIter
//-----------------------------------------------------------------------


void EnergyArclenModel::setMaxIter ( int count )
{
  JEM_PRECHECK ( count > 0 );

  optIter_ = count;
}


//-----------------------------------------------------------------------
//   setLoadIncr
//-----------------------------------------------------------------------


void EnergyArclenModel::setLoadIncr ( double incr )
{
  loadIncr_ = incr;

  if ( istep_ > 0 )
  {
    istep_ = 0;
  }
}


//-----------------------------------------------------------------------
//   setLoadScale
//-----------------------------------------------------------------------


void EnergyArclenModel::setLoadScale ( double scale )
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


void EnergyArclenModel::setIncrRange

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


Ref<Model> EnergyArclenModel::makeNew

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


void EnergyArclenModel::init_ ( const Properties& globdat )
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


void EnergyArclenModel::initLoad_ ( const Properties& globdat )
{  
  load_.resize ( dofs_->dofCount() );
  load_ = 0.0;

  Properties  params;
  params.set ( ActionParams::EXT_VECTOR, load_ );

  child_->takeAction  ( Actions::GET_EXT_VECTOR, params, globdat );
  dofs_ ->resetEvents ();

  updated_ |= U_LOAD_;
}


//-----------------------------------------------------------------------
//   initWeights_
//-----------------------------------------------------------------------


void EnergyArclenModel::initWeights_ ( const Properties& globdat )
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
    Vector      values ( nnz );


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


void EnergyArclenModel::evalArcFunc_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jem::max;
  using jem::isTiny;
  using jem::System;

  double  fvalue;
  double  lambda, lambda0;

  Vector  jac10;
  double  jac11;

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

  // for the elastic branch, use the standard arclen method

  if ( isLoadControl_ )
  {
    // load control with arclen function phi = delta lambda - delta l = 0

    params.get ( jac10,   ArclenParams::JACOBIAN10     );
    params.get ( lambda,  ArclenParams::LOAD_SCALE     );
    params.get ( lambda0, ArclenParams::OLD_LOAD_SCALE );

    jac10  = 0.0;
    jac11  = 1.0;

    fvalue = ( lambda - lambda0 )  - arcLength_ ;
  }

  // Energy Release Control Method

  else
  {
    Vector  u, u0, du, tmp;

    StateVector::get    ( u,  dofs_, globdat );
    StateVector::getOld ( u0, dofs_, globdat );

    const int uSize = u.size ();

    du. resize ( uSize );
    tmp.resize ( uSize );

    params.get ( jac10,   ArclenParams::JACOBIAN10     );
    params.get ( lambda,  ArclenParams::LOAD_SCALE     );
    params.get ( lambda0, ArclenParams::OLD_LOAD_SCALE );

    jac10  =  0.5 * lambda0 * load_;
    jac11  = -0.5 * vspace_->product ( u0, load_ );

    axpy ( du, u, -1.0, u0 );

    tmp = lambda0 * du - ( lambda - lambda0 ) * u0;

    double dEnergy   = 0.5 * vspace_->product ( tmp, load_ );

    // total stored energy is computed and stored 

    energy0_  = 0.5 * vspace_->product ( u0, load_ );
    energy0_ *= lambda0;

    fvalue = dEnergy - arcLength_ ;
  }

  params.set ( ArclenParams::JACOBIAN11, jac11  );
  params.set ( ArclenParams::ARC_FUNC,   fvalue );
}


//-----------------------------------------------------------------------
//   getUnitLoad_
//-----------------------------------------------------------------------


void EnergyArclenModel::getUnitLoad_

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


void EnergyArclenModel::commit_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jem::maxOf;

  int   iiter;

  params.get ( iiter,   ArclenParams::IITER );

  oldScale_ = loadScale_;

  //if ( loadIncr_ < 0 ) down_ = true;

  // arc-length control

  if ( !isLoadControl_ )
  {
    if ( loadIncr_ > 0. && arcLength_ < switchBack_ )
    {
      arcLength_     = 2 * loadIncr_;
      isLoadControl_ = true;

      System::out() << "Switch back to load control ...with load increment "
	            << arcLength_ << "\n";
    }
    else
    {
      // adjust the path following parameter Delta tau

      double n = ( iiter - optIter_ ) / 4.0;

      arcLength_ *= ::pow ( 0.5, n );

      // check if it falls within allowable interval (defined by user)

      if        ( arcLength_ < minIncr_ )
      {
        arcLength_ = minIncr_;
      }
      else if   ( arcLength_ > maxIncr_ )
      {
        arcLength_ = maxIncr_;
      }

      System::out() << "path following parameter: " 
                    << nformat_.print ( arcLength_ ) << "\n"
                    << ", energy0:                " 
                    << nformat_.print ( energy0_ )   << "\n\n";
      }

  }

  // increase the load step number

  if ( istep_ < maxOf( istep_ ) )
  {
    istep_++;
  }

  // store highest loadScale_ value

  if ( loadScale_ > maxLoadScale_ )
  { 
    maxLoadScale_ = loadScale_;
  }
}

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void EnergyArclenModel::checkCommit_

  ( const Properties&  params,
    const Properties&  globdat )

{
  Vector  u, u0, du, tmp;

  double  lambda;
  double  dlambda;
  double  lambda0;
  int     iiter;

  triedLarge_ = false;

  // need adapted StepModule for this! [fpm]

  if ( params.find ( iiter, "iterCount" ) )
  {
    params.get ( lambda,  "loadScale" );
    params.get ( dlambda, "loadIncr"  );

    lambda0    = lambda - dlambda;

    loadIncr_  = lambda - oldScale_;
    loadScale_ = lambda;

    StateVector  ::get       ( u,                dofs_, globdat );
    StateVector  ::getOld    ( u0,               dofs_, globdat );
    VectorManager::getVector ( du, DELTA_STATE_, dofs_, globdat );

    axpy ( du, u, -1.0, u0 );

    if ( iiter > 0 && isLoadControl_ )
    {
      // check the number of iterations
      // if bigger than the desired number (usually 4 or 5), then
      // switch to enery release control

      const int ls = load_.size ();

      tmp.resize ( du.size () );
      tmp = lambda0 * du - loadIncr_ * u0;
      tmp.reshape ( ls );
      u0. reshape ( ls );

      double dEnergy   = 0.5 * vspace_->product ( tmp, load_ );

      energy0_  = 0.5 * vspace_->product (  u0, load_ );
      energy0_ *= lambda0;

      System::out() << "EnergyArclen, dEnergy      = " 
                    << nformat_.print ( dEnergy ) << endl 
                    << "              Total energy = " 
                    << nformat_.print ( energy0_) << endl << endl;

      // switch to energy-based control
      
      if ( iiter >= swtIter_ )
      {
        System::out() << "Switching to arclength ...\n\n";

        isLoadControl_ = false;
        arcLength_     = dEnergy;

        if ( arcLength_ > maxIncr_ )
        {
          arcLength_ = maxIncr_;
          System::warn() << "Dissipated energy: " 
                         << nformat_.print ( arcLength_ ) 
                         << " > maxIncr_: " 
                         << nformat_.print ( maxIncr_ )
	                 << "  use reduced denergy. You may increase maxIncr_\n";
        }

        if ( arcLength_ < minIncr_ )
        {
          // energy MUST be significant

          System::warn() << "Not switching after all ... " 
                         << " since dissipated energy: " 
                         << nformat_.print ( arcLength_ ) 
                         << " < minIncr_: " 
                         << nformat_.print ( minIncr_ ) << "\n";

          isLoadControl_ = true;
          arcLength_     = loadIncr_; 
        }
      }
      else
      {
        lastArcL_ = dEnergy; // used when divergence occured in load control
	                     // then use lastArcL_ for arc-length control
      }
    }    // end of if isLoadControl
  }
}

//-----------------------------------------------------------------------
//   reduceStep_
//-----------------------------------------------------------------------

void EnergyArclenModel::reduceStep_

  ( const Properties&  params,
    const Properties&  globdat )

{
  if ( isLoadControl_ )
  {
    // load control, then switch to energy arc-length method

    isLoadControl_ = false;
    arcLength_     = lastArcL_; 

    if ( arcLength_ > maxIncr_ )
    {
      arcLength_ = maxIncr_;

      System::warn() << "Path following parameter reduced,"
                     << "  you may want to increase maxIncr.\n\n";
    }

    if ( arcLength_ < minIncr_ )
    {
      // energy MUST be significant

      System::warn() << "Path following parameter increased!\n\n";
      arcLength_ = minIncr_; 
    }

    System::out() << " Switching to arclength with step size "
                  << arcLength_ << "\n";
  }
  else
  {
    // arc-length control, then adapt step size

    System::out() << "Adapting step size...\n"
                  << "Current path following parameter: " << arcLength_ << " "
                  << "with load increment: " << loadIncr_ << "\n\n";

    if ( arcLength_ > 0.1 * minIncr_ ) 
    {
      // reduce arc-length

      arcLength_ *= .55;

      System::out() << "Reduced path following parameter: " 
                    << arcLength_ << "\n\n";
    }
    else
    {
      if ( ! triedLarge_ )
      {
        // try with large value 

        triedLarge_ = true;

        arcLength_ = 1.1 * maxIncr_ ;

        System::out() << "Increased path following parameter: " 
                      << arcLength_ << "\n\n";
      }
      else
      {
	// it's hopeless (whole range has been tried for this load step)

        //System::out() << myName_ << " Everything from large to small "
                                // << "arc-length has been tried...\n\n";

        //params.set( "RAMP", true );
      }
    }
  }

}


//-----------------------------------------------------------------------
//   connect_
//-----------------------------------------------------------------------


void EnergyArclenModel::connect_ ()
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


void EnergyArclenModel::dofsChanged_ ()
{
  updated_ = 0;
}


//-----------------------------------------------------------------------
//   consChanged_
//-----------------------------------------------------------------------


void EnergyArclenModel::consChanged_ ()
{
  updated_ &= ~U_WEIGHTS_;
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareEnergyArclenModel
//-----------------------------------------------------------------------


void declareEnergyArclenModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( EnergyArclenModel::TYPE_NAME,
                          & EnergyArclenModel::makeNew );
}
