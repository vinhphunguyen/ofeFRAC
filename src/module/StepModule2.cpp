//#define FRANK_DEBUG
#include <jem/base/limits.h>
#include <jem/base/System.h>
#include <jem/base/Exception.h>
#include <jem/base/IllegalOperationException.h>
#include <jem/util/Properties.h>
#include <jive/util/Globdat.h>
#include <jem/util/StringUtils.h>
#include <jive/util/utilities.h>
#include <jive/model/Model.h>
#include <jive/model/Actions.h> 
#include <jive/implict/SolverInfo.h>
#include <jive/implict/NonlinModule.h>
#include <jive/implict/ArclenModule.h>
#include <jive/app/ModuleFactory.h>

#include "StepModule2.h"

using namespace jem;

using jem::util::Properties; 
using jive::util::Globdat;
using jem::util::Properties;
using jive::util::joinNames;
using jive::implict::NonlinModule;
using jive::implict::ArclenModule;

//=======================================================================
//   class StepModule2
//=======================================================================
//in case solver does not converge - output is written
//with gmsh and fd plot
//but the last written step is not OK!!!!
//this is not a converged but the not 
//converged result!!!!!!!


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* StepModule2::TYPE_NAME = "Step2";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


StepModule2::StepModule2 

( const String&  name,
  Ref<SolverModule> solver )

 :  Super   ( name   ),
    solver_ (solver)
{
  using jem::NIL;

  //JEM_PRECHECK ( solver != NIL );

  myName_ = solver->getName ();
  istep_  = 0;
  istep0_ = 0;
}




StepModule2::~StepModule2 ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status StepModule2::init

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


Module::Status StepModule2::run ( const Properties& globdat )
{
  using jem::maxOf;
  using jem::Ref;
  using jem::System;
  using jem::io::endl;
  using jive::model::Model;
  using jive::implict::SolverInfo;
  using jive::model::Actions;

  Properties  info = SolverInfo::get ( globdat );


  info.clear ();

  solver_->advance ( globdat );

  // Advance the time/load step number in a robust way (avoiding
  // overflow).

  istep_ = istep0_;

  if ( istep_ < maxOf( istep_ ) )
  {
    istep_++;
  }

  globdat.set ( Globdat::TIME_STEP, istep_ );

  //refernece to models that should take action

  Ref<Model>  model = Model::get ( globdat, getContext() );

  //make a properties object to go with the STORE_RESIDUAL option

  Properties params_residual;
  double dum = 0;
  params_residual.set("residual", dum);

  try
  {
    solver_->solve ( info, globdat );

    //get residual and store it in params)residual

    Properties info = Globdat::getVariables("solverInfo", globdat);

    info.find(dum,SolverInfo::RESIDUAL);

    params_residual.set("residual", dum);

    //call a model to store the residual

    model->takeAction("STORE_RESIDUAL", params_residual, globdat);

    solver_->commit ( globdat );

    //check if max displacement is reached
    //if this check is implemented in a model

    Properties params_check_max_disp;
    model->takeAction("CHECK_MAX_DISP", params_check_max_disp, globdat);

    bool max_disp_reached = false;
    params_check_max_disp.find(max_disp_reached, "max_disp_reached");

    if(max_disp_reached)
      {
	System::warn() << "Maximum displacement is reached\n";
	return EXIT;
      }
    else
      {
	istep0_ = istep_; 
      }
  }
  catch ( const jem::Exception& ex )
  {
    //make parameter object to transfer data back

    Properties  params_step;

    //call a model to do something now (change stepsize)

    model->takeAction("REDUCE_STEP", params_step, globdat);

   //get residual and store it in params)residual

    Properties info = Globdat::getVariables("solverInfo", globdat);

    info.find(dum,SolverInfo::RESIDUAL);

    params_residual.set("residual", dum);

    //call a model to store the residual

    model->takeAction("STORE_RESIDUAL", params_residual, globdat);


    bool ramp = false;
    if(params_step.find(ramp, "RAMP"))
      {
	if(ramp)
	  {
	    //call a model to do something now (change stepsize)
	    model->takeAction("WRITE_GMSH_FLASH", params_step, globdat);

	    throw IllegalOperationException("AdaptLoadScaleModel","Solver did not converge. No ");     
	  
	  }
      }

    //stop solver

    solver_->cancel ( globdat );

    // Restore the time/load step.

    istep_ = istep0_;
    

    globdat.set ( Globdat::TIME_STEP, istep_ );

    //throw;
  }

  // Ask the model whether to accept this new solution.



  Properties  params;
  bool        accept = true;
  //ra removed - because Commit eliminates higly damaged elements
  model->takeAction ( "CHECK_COMMIT", params, globdat );

  params.find ( accept, "accept" );


  if ( ! accept )
  {
    System::err() << "In StepModule2 there is no answer on CHECK_COMMIT implemented!";
    //throw IllegalOperationException("StepModule2::run","In StepModule2 there is no answer on CHECK_COMMIT implemented!");     

  }
  else
  {
    System::err() << "In StepModule2 there is no answer on CHECK_COMMIT implemented!";
    //throw IllegalOperationException("StepModule2::run","In StepModule2 there is no answer on CHECK_COMMIT implemented!");     

  }



  return OK;
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void StepModule2::shutdown ( const Properties& globdat )
{
  solver_->shutdown ( globdat );
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void StepModule2::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  solver_->configure ( props, globdat );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void StepModule2::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  solver_->getConfig ( conf, globdat );
}

//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------

Ref<Module>  StepModule2::makeNew

  ( const String&           name,
    const Properties&       conf,
    const Properties&       globdat,
    const Properties&       props )
{

  Properties  myConf  = conf .makeProps ( name );
  Properties  myProps = props.findProps ( name );

  Ref<SolverModule>  solver;

  Properties childSolProps = myProps.findProps ( "childSolver"  );


  String solverName;

  childSolProps.find ( solverName, "type" );
  System::out() << solverName;

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

void declareStepModule2 ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( StepModule2::TYPE_NAME,
                         & StepModule2::makeNew );
}

