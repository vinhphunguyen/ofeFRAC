/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implements the step module. This module
 *  check if a converged solution is accepted or not. In case not,
 *  it omits the converged solution and solves the current load step
 *  again with some changes such as reduced step size or recently
 *  introduced crack.
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
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
#include <jive/model/Model.h>
#include <jive/implict/SolverInfo.h>
#include <jive/implict/NonlinModule.h>
#include <jive/implict/ArclenModule.h>
#include <jive/app/ModuleFactory.h>

#include "StepModule.h"


using jem::IllegalOperationException;
using jem::System;
using jem::newInstance;
using jem::io::endl;
using jem::util::StringUtils;
using jive::util::Globdat;
using jive::util::joinNames;
using jive::implict::NonlinModule;
using jive::implict::ArclenModule;



//=======================================================================
//   class StepModule
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* StepModule::TYPE_NAME = "Step";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

StepModule::StepModule 

  ( const String&  name,
    Ref<SolverModule> solver ) :

      Super   ( name   ),
      solver_ ( solver )

{
  istep_  = 0;
  istep0_ = 0; 
}

StepModule::~StepModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status StepModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  istep_ = istep0_ = 0;

  globdat.set ( Globdat::TIME_STEP, istep_ );

  return solver_->init ( conf, props, globdat );
}


//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------


Module::Status StepModule::run ( const Properties& globdat )
{
  using jem::maxOf;
  using jem::Ref;
  using jem::System;
  using jem::io::endl;
  using jive::model::Model;
  using jive::implict::SolverInfo;

  Ref<Model>  model     = Model::get ( globdat, getContext() );

  Properties  info      = SolverInfo::get ( globdat );
  Properties  params;

  bool        accept    = true;
  bool        converged = true;

  info.clear ();

  solver_->advance ( globdat );

  // Advance the time/load step number in a robust way (avoiding overflow).

  istep_ = istep0_;

  if ( istep_ < maxOf( istep_ ) )
  {
    istep_++;
  }

  globdat.set ( Globdat::TIME_STEP, istep_ );

  // try to find a solution

  try
  {
    solver_->solve ( info, globdat );
  }
  catch ( const jem::Exception& ex )
  {
    System::out() << "StepModule, solver exception: " << ex.name() << "\n"
                  << "                              " << ex.what() << "\n\n";

    // solver has not converged

    converged = false;

    // try to change the step size 

    if ( model->takeAction ( "REDUCE_STEP", params, globdat ) )
    {
      // stop computation if step size change is hopeless

      bool ramp;

      if ( params.find( ramp, "RAMP" ) )
      {
        if ( ramp )
        {
            System::out() << "StepModule: solver did not converge. No\n";
            return EXIT;
        }
      }
    }
    else
    {
      // cancel previous solution and throw exception from solver

      System::out() << "No time step reduction implemented.\n";

      solver_->cancel ( globdat );

      // Restore the time/load step.

      istep_ = istep0_;

      globdat.set ( Globdat::TIME_STEP, istep_ );

      throw;
    }
  }

  if ( converged )
  {
    // Ask the model whether to accept converged solution

    int         iterCount; 
    double      loadScale;
    double      loadIncr ;

    info.find ( iterCount, SolverInfo::ITER_COUNT );
    info.find ( loadScale, SolverInfo::LOAD_SCALE );
    info.find ( loadIncr , SolverInfo::LOAD_INCR  );

    params.set ( "iterCount", iterCount );
    params.set ( "loadScale", loadScale );
    params.set ( "loadIncr",  loadIncr  );

    model->takeAction ( "CHECKCOMMIT", params, globdat );

    params.find ( accept, "accept" );

    if ( accept )
    {
      // Done with this time step

      solver_->commit ( globdat );

      istep0_ = istep_;
    }
    else
    {
      // Continue with this solution

      // Not canceling here may be dangerous, but in some (most/all?) cases
      // it saves some iterations [fpm]

      System::out() << "Continuing with the same load step\n\n";

      istep_ = istep0_;

      globdat.set ( Globdat::TIME_STEP, istep_ );
    }
  }

  else  // not converged
  {
    // Discard this solution.

    System::out() << "Returning to the same load step\n\n";

    solver_->cancel ( globdat );

    istep_ = istep0_;

    globdat.set ( Globdat::TIME_STEP, istep_ );
  }

  return OK;
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void StepModule::shutdown ( const Properties& globdat )
{
  solver_->shutdown ( globdat );
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void StepModule::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  solver_->configure ( props, globdat );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void StepModule::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  solver_->getConfig ( conf, globdat );
}

//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------

Ref<Module>  StepModule::makeNew

  ( const String&           name,
    const Properties&       conf,
    const Properties&       props,
    const Properties&       globdat )
{

  Properties  myConf  = conf .makeProps ( name );
  Properties  myProps = props.findProps ( name );

  Ref<SolverModule>  solver;

  Properties childSolProps = myProps.findProps ( "childSolver"  );

  String solverName;

  childSolProps.find ( solverName, "type" );

  if       ( solverName == "Nonlin" )
  {
    solver = newInstance<NonlinModule> ( joinNames ( name, "childSolver" ) );
  }
  else if  ( solverName == "Arclen" )
  {
    solver = newInstance<ArclenModule> ( joinNames ( name, "childSolver" ) );
  }
  else
  {
    solver = NIL;
  }

  return newInstance<Self> ( name, solver );
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareStepModule
//-----------------------------------------------------------------------

void declareStepModule ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( StepModule::TYPE_NAME,  
                         & StepModule::makeNew );

}

