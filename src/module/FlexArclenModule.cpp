/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implements a flexible path following
 *  method using both displacement control and energy based 
 *  arclength control. This is very helpful in tracing
 *  complex equilibrium paths.
 *
 *  Basic idea:
 *
 *    Start the simulation with load control
 *    until first snapback is detected (divergence encountered)
 *    then resolve the same load step with energy-based arclength
 *    control with amount of released energy computed. Proceeding
 *    with this arclength control till hardening branch is 
 *    detected, then switch back to displacement control.
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 28 January 2009
 *
 *  Status:
 *
 *   30 January 2009: works with transition from load to arclen
 *   26 February 2009:
 *     - implemented CHECK_COMMIT for crack propagation
 *     - switch back to disp due to high residual (with XArclenModule)
 *     - TmpNonlin: do the first solve with arclen,
 *       and switch to nonlin for remainder of step if not accepted
 *     - reduceStep for both strategies
 *     - when minIncr reached, switch solver
 *
 *
 */

#include <jem/base/limits.h>
#include <jem/base/System.h>
#include <jem/base/Exception.h>
#include <jem/base/IllegalOperationException.h>
#include <jem/util/Properties.h>
#include <jem/util/StringUtils.h>
#include <jive/util/Globdat.h>
#include <jive/util/utilities.h>
#include <jive/implict/SolverInfo.h>
#include <jive/implict/NonlinModule.h>
#include <jive/implict/ArclenModule.h>
#include <jive/app/ModuleFactory.h>

#include "FlexArclenModule.h"
#include "SolverNames.h"


using jem::IllegalOperationException;
using jem::System;
using jem::newInstance;
using jem::io::endl;
using jem::util::StringUtils;
using jive::implict::SolverInfo;
using jive::util::Globdat;
using jive::util::joinNames;


//=======================================================================
//   class FlexArclenModule
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* FlexArclenModule::TYPE_NAME       = "FlexArclen";
const char* FlexArclenModule::TMP_NONLIN_PROP = "doTmpNonlin";
const char* FlexArclenModule::NONLIN          = "nonLin";
const char* FlexArclenModule::ARCLEN          = "arcLen";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

FlexArclenModule::FlexArclenModule 

  ( const String&  name,
    Ref<NonlinModule>  solver1,
    Ref<ArclenModule>  solver2 ) :

      Super         ( name    ),
      nonlinSolver_ ( solver1 ),
      arclenSolver_ ( solver2 ),
      currentSolver_( solver1 ),
      istep_        ( 0       ), 
      istep0_       ( 0       ),
      doTmpNonlin_  ( false   ),
      isTmpNonlin_  ( false   ),
      doneTmpNonlin_( false   )

{
}

FlexArclenModule::~FlexArclenModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status FlexArclenModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  model_ = Model::get ( globdat, getContext() );
  istep_ = istep0_ = 0;

  globdat.set ( Globdat::TIME_STEP, istep_ );

  nonlinSolver_->init ( conf, props, globdat );
  arclenSolver_->init ( conf, props, globdat );

  return OK;
}


//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------


Module::Status FlexArclenModule::run ( const Properties& globdat )
{
  using jem::maxOf;
  using jem::Exception;
  using jem::Ref;
  using jive::model::Model;

  Properties  info      = SolverInfo::get ( globdat );
  Properties  params;

  bool        accept    = true;
  String      doSwitch;

  info.clear ();

  currentSolver_->advance ( globdat );
  
  // Advance the time/load step number

  istep_ = istep0_;

  if ( istep_ < maxOf( istep_ ) )
  {
    istep_++;
  }

  globdat.set ( Globdat::TIME_STEP, istep_ );

  // try to find a solution

  try
  {
    currentSolver_->solve ( info, globdat );

    converged_ = true;
  }
  catch ( const jem::Exception ex )
  {
    System::out() << "FlexArclenModule caught a solver exception: \n" 
                  << "   " << ex.name() << "\n"
                  << "   " << ex.what() << "\n\n";

    converged_ = false;
  }

  int         iterCount;
  info.find ( iterCount , SolverInfo::ITER_COUNT );

  if ( converged_ )
  {
    // Ask the model whether to accept converged solution
    // (and check for switch)

    params.set( "iterCount", iterCount );

    model_->takeAction ( SolverNames::CHECK_COMMIT, params, globdat );

    params.find ( accept, SolverNames::ACCEPT );

    if ( accept )
    {
      // The final solution for this time step has been obtained
      // (possibly switch from tmpNonlin)

      commit_ ( globdat );
    }
    else
    {
      // We have an equilibrium solution but it is not accepted:
      // continue with this time step
      // (possibly switch with tmpNonlin)

      continue_ ( globdat );
    }

    // in CHECK_COMMIT, a check for switch has also been done:
    // do switch if the model has asked for it
    // (except when in tmpNonlin)

    params.find ( doSwitch, SolverNames::DO_SWITCH );

    if ( doSwitch == "please" && ! isTmpNonlin_ )
    {
      switchStrategy_ ( globdat );
    }

    if ( converged_ && accept ) isTmpNonlin_ = false;
  }
  else   // not converged
  {
    // 1. discard solution

    cancel_ ( globdat );

    // 2. find proper strategy to try again
   
    // ****
    // if the XArclenModule has quit because high residual scale factor
    // --> switch to nonlin

    if ( iterCount == 1 && switchToNonlin_ ( globdat ) ) 
    {
      System::out() << "Returning to the same load step with "
        << currentSolver_->getName() << "\n";
    }
    
    // ****
    // if tmpNonlin was tried, try again without

    else if ( isTmpNonlin_ )
    {
      doneTmpNonlin_ = true;
      isTmpNonlin_   = false;

      switchToArclen_ ( globdat );

      System::out() << "Returning to the same load step "
                    << "with tmpNonlin option turned off\n";
    }

    // ****
    // try to reduce the step size

    else if ( reduceStep_ ( globdat ) )
    {
      System::out() << "Returning to the same load step "
                    << "with reduced step\n";
    }

    // ****
    // switch solver

    else if ( switchStrategy_ ( globdat ) )
    {
      System::out() << "Returning to the same load step, now with "
        << currentSolver_->getName() << "\n";
    }

    // ****
    // give it up

    else
    {
      System::out() << "\n*** FlexArclenModule is out of inspiration,\n"
        "    can't find a winning strategy for this step. Sorry.\n\n";

      return EXIT;
    }
  }
  // set flag in globdat for output indicating solution quality

  globdat.set ( "var.accepted", converged_ && accept );

  return OK;
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void FlexArclenModule::shutdown ( const Properties& globdat )
{
  nonlinSolver_->shutdown ( globdat );
  arclenSolver_->shutdown ( globdat );
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void FlexArclenModule::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties myProps = props.findProps ( myName_ );

  myProps.find ( doTmpNonlin_, TMP_NONLIN_PROP );

  nonlinSolver_->configure ( props, globdat );
  arclenSolver_->configure ( props, globdat );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void FlexArclenModule::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf  = conf.makeProps ( myName_ );

  myConf.set( "type" , TYPE_NAME );

  myConf.set( TMP_NONLIN_PROP, doTmpNonlin_ );

  nonlinSolver_->getConfig ( conf, globdat );
  arclenSolver_->getConfig ( conf, globdat );
}

//-----------------------------------------------------------------------
//   commit_
//-----------------------------------------------------------------------

void FlexArclenModule::commit_

  ( const Properties&       globdat )

{
  // store this solution and proceed with next time step
  // model->takeAction(COMMIT) is called by solver!

  currentSolver_->commit ( globdat );

  istep0_ = istep_;

  // reset variables

  nCancels_       = 0;
  nContinues_     = 0;
  nonlinTriedAll_ = false;
  arclenTriedAll_ = false;
  doneTmpNonlin_  = false;

  // switch back to Arlcen if tmpNonlin option was used

  if ( isTmpNonlin_ )
  {
    switchToArclen_ ( globdat );
  }

  // isTmpNonlin_ flag is reset later!
}

//-----------------------------------------------------------------------
//   continue_
//-----------------------------------------------------------------------

void FlexArclenModule::continue_

  ( const Properties&       globdat )

{
  // the solver is not canceled, because the current solution which
  // (converged but not accepted) is a good starting point
  
  System::out() << "Continuing with this load step\n\n";

  istep_ = istep0_;

  globdat.set ( Globdat::TIME_STEP, istep_ );

  ++nContinues_;

  if ( doTmpNonlin_ 
       && ! ( isTmpNonlin_ || doneTmpNonlin_ ) 
       && currentSolver_ == arclenSolver_ )
  {
    isTmpNonlin_ = switchToNonlin_ ( globdat );
  }
}

//-----------------------------------------------------------------------
//   cancel_
//-----------------------------------------------------------------------

void FlexArclenModule::cancel_

  ( const Properties&       globdat )

{
  // discard solution, go back to the beginning of the time step
  // model->takeAction(CANCEL) is called by solver!

  ++nCancels_;

  nContinues_ = 0;

  currentSolver_->cancel ( globdat );

  istep_ = istep0_;

  globdat.set ( Globdat::TIME_STEP, istep_ );
}

//-----------------------------------------------------------------------
//   switchStrategy_
//-----------------------------------------------------------------------

bool FlexArclenModule::switchStrategy_

  ( const Properties&       globdat )

{
  return ( switchToNonlin_(globdat) || switchToArclen_(globdat) );
}

//-----------------------------------------------------------------------
//   switchToNonlin_
//-----------------------------------------------------------------------

bool FlexArclenModule::switchToNonlin_

  ( const Properties&       globdat )

{
  if ( nonlinTriedAll_ || currentSolver_ == nonlinSolver_ ) return false;

  Properties params;

  params.set ( SolverNames::CONVERGED, converged_ );

  model_->takeAction ( SolverNames::TO_DISP, params, globdat );
    
  currentSolver_ = nonlinSolver_;

  return true;
}

//-----------------------------------------------------------------------
//   switchToArclen_
//-----------------------------------------------------------------------

bool FlexArclenModule::switchToArclen_

  ( const Properties&       globdat )

{
  if ( arclenTriedAll_ || currentSolver_ == arclenSolver_ ) return false;

  Properties  params;

  double      loadScale;

  model_->takeAction ( SolverNames::TO_ARCL, params, globdat );

  // give loadScale from model to solver 

  params.get( loadScale, "OldLoadScale" );

  arclenSolver_->setLoadScale ( loadScale );

  currentSolver_ = arclenSolver_;

  return true;
}

//-----------------------------------------------------------------------
//   reduceStep_
//-----------------------------------------------------------------------

bool FlexArclenModule::reduceStep_

  ( const Properties&       globdat )

{
  Properties  params;

  bool        reduced = false;

  model_->takeAction ( SolverNames::REDUCE_STEP, params, globdat );

  params.find ( reduced, SolverNames::DONE );

  if ( ! reduced )
  {
    // time step can't be reduced further for this solver, 
    // don't switch to it again this time step

    currentSolver_ == nonlinSolver_ ?  nonlinTriedAll_ = true : arclenTriedAll_ = true ;
  }

  return reduced;
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------

Ref<Module>  FlexArclenModule::makeNew

  ( const String&           name,
    const Properties&       conf,
    const Properties&       globdat,
    const Properties&       props )

{
  Properties  myConf  = conf .makeProps ( name );
  Properties  myProps = props.findProps ( name );

  Ref<NonlinModule>  solver1 = newInstance<NonlinModule> 

                                  ( joinNames ( name, NONLIN ) );

  Ref<ArclenModule>  solver2 = newInstance<ArclenModule> 

                                  ( joinNames ( name, ARCLEN ) );


  return newInstance<Self> ( name, solver1, solver2 );
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareFlexArclenModule
//-----------------------------------------------------------------------

void declareFlexArclenModule ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( FlexArclenModule::TYPE_NAME,
                         & FlexArclenModule::makeNew );
}

