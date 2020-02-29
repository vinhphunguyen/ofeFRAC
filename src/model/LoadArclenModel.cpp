/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implements a model that combines
 *  load and arclen control in one model. This model
 *  must be used with FlexArclenModule to have a flexible path following solver.
 *
 *  Basic ideas:
 *
 *    1. Load control: a so-called unit external force vector
 *       is computed once. The load scale lambda is updated
 *       by a constant amount (can be changed later). During
 *       this stage, in FlexArclenModule, the NonlinModule is being
 *       used. This is the case until divergence occurs (pass the peak)
 *
 *   2.  Arclen control (based on energy released). Now, in FlexArclenModule,
 *       the ArclenModule is active so that lambda is now an unknown.
 *       This is the case until hardening branch is detected, then
 *       switch back to (1), Load control
 *
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 January 2009
 *
 */

#include <jem/base/limits.h>
#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/base/Float.h>
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
#include <jive/implict/SolverInfo.h>


#include "UserModels.h"
#include "LoadArclenModel.h"
#include "module/SolverNames.h"


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
//   class LoadArclenModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  LoadArclenModel::TYPE_NAME       = "LoadArclen";

const char*  LoadArclenModel::MODEL_PROP      = "model";
const char*  LoadArclenModel::ARC_FUNC_PROP   = "arcFunc";
const char*  LoadArclenModel::OPT_ITER_PROP   = "optIter";
const char*  LoadArclenModel::SWT_ITER_PROP   = "swtIter";
const char*  LoadArclenModel::SWT_ENER_PROP   = "swtEner";
const char*  LoadArclenModel::MIN_INCR_PROP   = "minIncr";
const char*  LoadArclenModel::MAX_INCR_PROP   = "maxIncr";
const char*  LoadArclenModel::LOAD_INCR_PROP  = "loadIncr";
const char*  LoadArclenModel::LOAD_SCALE_PROP = "loadScale";
const char*  LoadArclenModel::EXIT_FRAC_PROP  = "exitFraction";
const char*  LoadArclenModel::REDUCTION_PROP  = "reduction";

const char*  LoadArclenModel::DELTA_STATE_    = "deltaState0";

const int    LoadArclenModel::U_LOAD_         = 1 << 0;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


LoadArclenModel::LoadArclenModel

  ( const String&      name,
    const Ref<Model>&  child ) :

    Super  ( name  ),
    child_ ( newInstance<LoadArclenBCs>() ),
    out_  ( System::info( name ) )

{
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
  swtEner_   = Float::MAX_VALUE;
  gC_        = Float::MAX_VALUE;

  maxLoadScale_ = 0.0;
  exitFraction_ = 0.0;          // terminate computation when load approaches zero
  reduction_    = 1.0;

  isLoadControl_ = true;        // starts with load control
  onceDown_      = false;

  nformat_.setScientific     (   );
  nformat_.setFractionDigits ( 4 );

}


LoadArclenModel::~LoadArclenModel ()
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool LoadArclenModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  if ( action == Actions::INIT )
  {
    init_ ( globdat );

    return true;
  }

  // compute the external force vector

  if ( action == Actions::GET_EXT_VECTOR )
  {
    System::out() << " Evaluating the external force ...\n\n";

    child_->computeExternalForce ( params, globdat );
    System::out() << " Evaluating the external force done.\n";

    return  true;
  }

  // compute arclength function

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

 
  if ( action == SolverNames::CHECK_COMMIT )
  {
    checkCommit_ ( params, globdat );
    checkSwitch_ ( params, globdat );

    return true;
  }


  if ( action == Actions::COMMIT )
  {
    commit_ ( params, globdat );

    return true;
  }

  // this REDUCE_STEP action is thrown by FlexArclenModule

  if ( action == SolverNames::REDUCE_STEP )
  {
    reduceStep_ ( params, globdat );

    return true;
  }

  if ( action == SolverNames::GET_DENERGY )
  {
    arcLength_     = getReleasedEnergy_ ( params, globdat );
    isLoadControl_ = false;

    if        ( arcLength_ < minIncr_ )
    {
      arcLength_ = minIncr_;
    }
    else if   ( arcLength_ > maxIncr_ )
    {
      arcLength_ = maxIncr_;
    }

    return true;
  }

  // this action is called in FlexArclenModule,
  // when switch from load control to arclen control

  if ( action == SolverNames::TO_ARCL )
  {
    toArclControl_ ( params, globdat );

    child_->toArclControl ( params, globdat );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   findModel
//-----------------------------------------------------------------------


Model* LoadArclenModel::findModel ( const String& name ) const
{
 return const_cast<Self*> ( this );
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void LoadArclenModel::configure

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
    myProps.find ( swtEner_,   SWT_ENER_PROP );

    myProps.find ( exitFraction_, EXIT_FRAC_PROP,
                   0.0,        1.0 );

    myProps.find ( reduction_, REDUCTION_PROP,
                   0.0,        1.0 );

    myProps.find ( gC_,  "fractureEnergy" );

    // configuring the child object of class DispArclenBCs

    Properties childProps = myProps.findProps ( "loadObj" );
    child_->configure ( childProps, globdat );
  }

  arcLength_ = loadIncr_;
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void LoadArclenModel::getConfig

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
  myConf.set ( LOAD_SCALE_PROP,  loadScale_    );
  myConf.set ( EXIT_FRAC_PROP,   exitFraction_ );

  Properties childConf = myConf.makeProps ( "loadObj" );
  child_->getConfig ( childConf, globdat );
}


//-----------------------------------------------------------------------
//   setArcFunc
//-----------------------------------------------------------------------


void LoadArclenModel::setArcFunc ( ArcFunc func )
{
  arcFunc_ = func;
}


//-----------------------------------------------------------------------
//   setMaxIter
//-----------------------------------------------------------------------


void LoadArclenModel::setMaxIter ( int count )
{
  JEM_PRECHECK ( count > 0 );

  optIter_ = count;
}


//-----------------------------------------------------------------------
//   setLoadIncr
//-----------------------------------------------------------------------


void LoadArclenModel::setLoadIncr ( double incr )
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


void LoadArclenModel::setLoadScale ( double scale )
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


void LoadArclenModel::setIncrRange

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


Ref<Model> LoadArclenModel::makeNew

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


void LoadArclenModel::init_ ( const Properties& globdat )
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


void LoadArclenModel::initLoad_ ( const Properties& globdat )
{

  System::out() << "Initialize unit load ... \n";

  Properties  params;

  load_.resize ( dofs_->dofCount() );

  load_ = 0.0;

  params.set ( ActionParams::EXT_VECTOR, load_ );

  child_->getUnitLoad ( params, globdat );

  dofs_ ->resetEvents ();

  updated_ |= U_LOAD_;
}


//-----------------------------------------------------------------------
//   evalArcFunc_
//-----------------------------------------------------------------------


void LoadArclenModel::evalArcFunc_

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

  // Energy Release Control Method

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

/*
  using jem::Float;

  if ( Float::isNaN( fvalue ) )
  {
    System::out() << "fvalue" << fvalue << ", arcLength_" << arcLength_ << "\n";
  }
*/

  params.set ( ArclenParams::JACOBIAN11, jac11  );
  params.set ( ArclenParams::ARC_FUNC,   fvalue );
}


//-----------------------------------------------------------------------
//   getUnitLoad_
//-----------------------------------------------------------------------


void LoadArclenModel::getUnitLoad_

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

// This is called when NonlinModule got convergence
// or ArclenModule did converge

void LoadArclenModel::commit_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jem::maxOf;
  using jive::model::ActionParams;

  //using jive::util::Globdat;

  oldScale_ = loadScale_;

  if ( isLoadControl_ )
  {
    // increase the load scale 

    loadScale_ += loadIncr_;

    // set to globdat so that the child model
    // compute correct fext = scale * unit force

    globdat.set ( ActionParams::SCALE_FACTOR, loadScale_ );
    globdat.set ( "OldLoadScale", oldScale_              );

    out_    << " old load scale: " << oldScale_ 
            << ", current load scale: " << loadScale_ 
            << ", load increment: " << loadIncr_ <<"\n"
            << ", energy0: " << energy0_ <<"\n";

    out_.flush();

    child_->commit ( isLoadControl_, globdat );
  }

  if ( !isLoadControl_ )
  {
    int   iiter;

    params.get ( iiter,   ArclenParams::IITER );

    double n    = ( iiter - optIter_ ) / 4.0;
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

    out_  << "path following parameter: " 
          << nformat_.print ( arcLength_ ) << "\n"
          << ", energy0:                " 
          << nformat_.print ( energy0_ )   << "\n\n";

    out_.flush();

    // loads change sign, a sign of hardening, then
    // switch back to load control

    double lambda0; 

    params      .get    ( loadScale_,  ArclenParams::LOAD_SCALE     );
    params      .get    ( lambda0,     ArclenParams::OLD_LOAD_SCALE );

    double dlambda = loadScale_ - lambda0;

    if ( dlambda < 0 ) onceDown_ = true;

    // This switching is problem dependent !!!

    if ( ( dlambda > 30. ) & ( onceDown_ ) )
    {
      System::out() << "Loads change sign, going back to load control...\n";

      System::out() << "load scale increment: " << dlambda << "\n";

      loadIncr_ = 30.0;

      toLoadControl_ ( params, globdat );

      onceDown_ = false;
    }
  }

  // increase the increment (load) step number

  if ( istep_ < maxOf( istep_ ) )
  {
    istep_++;
  }

  globdat.set ( SolverNames::IS_NONLIN, isLoadControl_ );
}

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void LoadArclenModel::checkCommit_

  ( const Properties&  params,
    const Properties&  globdat )

{
}

//-----------------------------------------------------------------------
//   reduceStep_
//-----------------------------------------------------------------------

// For now, only adjust the load increment for Arclen control
// by either (1) reduce the amount of released energy
// or (2) switch back to load control since the
// load is going up again

void LoadArclenModel::reduceStep_

  ( const Properties&  params,
    const Properties&  globdat )

{
  if ( isLoadControl_ )
  {
  }
  else
  {
    // arclen control, then adapt step size

    if ( arcLength_ > reduction_ * minIncr_ ) 
    {
      arcLength_ *= reduction_; // reduce arc-length

      System::out() << "Reduced path following parameter: "
                    << arcLength_ << "\n\n";

      params.set ( SolverNames::IS_NONLIN, false );
    }

  }
}

//-----------------------------------------------------------------------
//   toLoadControl_
//-----------------------------------------------------------------------

void LoadArclenModel::toLoadControl_

  ( const Properties&  params,
    const Properties&  globdat )

{
  // set constant load increment for coming
  // load control regime
  // how to define this value wisely???
  // it has unit of force, while the path following parameter
  // arcLength_ has unit of energy !!!

  //loadIncr_      = 5.;

  isLoadControl_ = true;

  // set IS_NONLIN to true so that FlexArclenModule
  // knows it should now use NonlinModule to solve the system

  globdat .set ( SolverNames::IS_NONLIN, true );

  // set to globdat so that the child model
  // compute correct fext = scale * unit force

  globdat.set  ( ActionParams::SCALE_FACTOR, loadScale_ + loadIncr_ );

  System::out() << myName_ << " Switching back to loadControl "
                << "with incr: " << loadIncr_ << "\n\n";
}


//-----------------------------------------------------------------------
//   connect_
//-----------------------------------------------------------------------


void LoadArclenModel::connect_ ()
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


void LoadArclenModel::dofsChanged_ ()
{
  updated_ = 0;
}


//-----------------------------------------------------------------------
//   consChanged_
//-----------------------------------------------------------------------


void LoadArclenModel::consChanged_ ()
{
}

//-----------------------------------------------------------------------
//   getReleasedEnergy_
//-----------------------------------------------------------------------

// Attention:
// Calling this when NonlinModule diverged,
// to compute the energy released gives ZERO since
// u = u0

double LoadArclenModel::getReleasedEnergy_

 ( const Properties&  params,
   const Properties&  globdat )

{
  using jem::maxOf;
  using jive::implict::SolverInfo;
  
  // so far, using NonlinModule so load_ is not yet initialized!!!

  if ( ! (updated_ & U_LOAD_) )
  {
    initLoad_ ( globdat );
  }

  Vector  u, u0, du;

  double  lambda;
  double  dlambda;
  double  lambda0;
  int     iiter;

  lambda0 = child_->getLoadScale0();
  dlambda = child_->getLoadScale () - lambda0;

  //loadIncr_  = lambda - oldScale_;
  //loadScale_ = lambda;

  System::out() << "load scale: " << lambda
                << ", load incr : " << loadIncr_
                << ", old load scale: " << lambda0 << "\n";

  StateVector  ::get       ( u,                dofs_, globdat );
  StateVector  ::getOld    ( u0,               dofs_, globdat ); // u = u0 !!!
  VectorManager::getVector ( du, DELTA_STATE_, dofs_, globdat );

  axpy ( du, u, -1.0, u0 );

  energy0_  = 0.5 * dot ( u0, load_ );
  energy0_ *= lambda0;

  //const int ls = load_.size ();

  Vector tmp ( du.size () );

  tmp = lambda0 * du - dlambda * u0; 

  //tmp.reshape ( ls );
  //u0. reshape ( ls );

  if ( vspace_ == NIL )
  {
    vspace_ = VectorSpace::get ( dofs_, globdat );
  }

  double dEnergy = 0.5 * vspace_->product ( tmp, load_ );

  System::out() << "released energy: " << dEnergy << "\n";

  return dEnergy;
}


//-----------------------------------------------------------------------
//   toArclControl_
//-----------------------------------------------------------------------

void LoadArclenModel::toArclControl_

  ( const Properties&  params,
    const Properties&  globdat )

{
  isLoadControl_ = false;
  arcLength_     = lastArcl_;

  if        ( arcLength_ < minIncr_ )
  {
    arcLength_ = minIncr_;
  }
  else if   ( arcLength_ > maxIncr_ )
  {
    arcLength_ = maxIncr_;
  }

  out_ << "Switching to arclength control with dEnergy: " 
       << arcLength_ << endl;
}

//-----------------------------------------------------------------------
//   checkSwitch_
//-----------------------------------------------------------------------

void LoadArclenModel::checkSwitch_

  ( const Properties&  params,
    const Properties&  globdat )

{
  int    iiter;

  if ( isLoadControl_ )
  {
    lastArcl_ = getReleasedEnergy_( params, globdat );

    params.get ( iiter, "iterCount" );

    out_ << "LoadControl: dEnergy " << lastArcl_ << 
            ", loadScale " << child_->getLoadScale() << endl;

    // store value anyway

    if ( lastArcl_ > swtEner_ ||
         ( iiter > swtIter_ && lastArcl_ > minIncr_ ) ||
         ( energy0_ > gC_ )
       ) 
    {
      params.set ( SolverNames::DO_SWITCH, "please" );
    }

    // check if it's going down 

    double lambda0 = ::fabs ( child_->getLoadScale0() );
    double dlambda = ::fabs ( child_->getLoadScale () ) - lambda0;

    onceDown_ |= dlambda < 0.;

    if (onceDown_) out_ << "Going down !!! \n";
  }
  else
  {
    lastArcl_ = arcLength_;
  }
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareLoadArclenModel
//-----------------------------------------------------------------------


void declareLoadArclenModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( LoadArclenModel::TYPE_NAME,
                          & LoadArclenModel::makeNew );
}
