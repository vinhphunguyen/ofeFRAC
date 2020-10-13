
#include <jive/app/ChainModule.h>
#include <jive/app/ControlModule.h>
#include <jive/app/OutputModule.h>
#include <jive/app/ReportModule.h>
#include <jive/app/Application.h>
#include <jive/app/InfoModule.h>
#include <jive/app/SampleModule.h>
#include <jive/app/UserconfModule.h>
#include <jive/geom/declare.h>
#include <jive/model/declare.h>
#include <jive/femodel/declare.h>
#include <jive/fem/declare.h>
#include <jive/fem/InputModule.h>
#include <jive/fem/InitModule.h>
#include <jive/fem/MPInputModule.h>
#include <jive/fem/PartitionModule.h>
#include <jive/fem/ShapeModule.h>
#include <jive/implict/LinsolveModule.h>
#include <jive/implict/NonlinModule.h>
#include <jive/implict/ArclenModule.h>
#include <jive/implict/Park3Module.h>
#include <jive/algebra/declare.h>
#include <jive/implict/declare.h>
#include <jive/app/declare.h>
#include <jive/gl/declare.h>
#include <jive/gl/DisplayModule.h>
#include <jive/gl/FemViewModule.h>
#include <jive/gl/GraphModule.h>

#include "model/UserModels.h"
#include "module/UserModules.h"
#include "femodel/FemModels.h"
#include "module/StepModule.h"
#include "module/StepModule2.h"
#include "module/ParaviewModule.h"
#include "declareShapes.h"

using namespace jem;

using jive::app::Application;
using jive::app::Module;
using jive::app::ChainModule;
using jive::app::OutputModule;
using jive::app::InfoModule;
using jive::app::ControlModule;
using jive::app::ReportModule;
using jive::app::SampleModule;
using jive::app::UserconfModule;
using jive::fem::InputModule;
using jive::fem::InitModule;
using jive::fem::MPInputModule;
using jive::fem::PartitionModule;
using jive::fem::ShapeModule;
using jive::implict::ArclenModule;
using jive::implict::NonlinModule;
using jive::implict::LinsolveModule;
using jive::implict::Park3Module;
using jive::gl::FemViewModule;
using jive::gl::GraphModule;
using jive::gl::DisplayModule;

//-----------------------------------------------------------------------
//   mainModule
//-----------------------------------------------------------------------


Ref<Module> mainModule ()
{
  Ref<ChainModule>    chain = newInstance<ChainModule> ();
  Ref<ControlModule>  ctrl;

  // Declare internal shapes, models and matrix builders. These
  // functions essentially store pointers to construction functions
  // that are called when Jive needs to create a shape, model or
  // matrix builder of a particular type.

  jive::geom:: declareIShapes     ();
  jive::geom::declareShapes       ();
  jive::model::declareModels      ();
  jive::fem::declareMBuilders     ();
  jive::implict::declareModels    ();
  jive::algebra::declareMBuilders ();
  jive::femodel::declareModels    ();

  // declare models/modules created by users

  declareShapes                   ();
  declareUserModels               ();
  declareUserModules              ();
  declareFemModels                ();

  // Declare all modules that are added dynamically
  // specified in *.pro file, not known at compile time

  jive::app     ::declareModules ();
  jive::implict ::declareModules ();
  jive::gl      ::declareModules ();

  // --------------------------------------------------------------
  //  add default modules

  chain->pushBack ( newInstance<MPInputModule>   ( ) );  // parallel computing
  chain->pushBack ( newInstance<PartitionModule> ( ) );  // must have !!!
  chain->pushBack ( newInstance<InputModule>     ( ) );  // read the input

  //EJ: This module is not needed when using the PartitionModule.
  //chain->pushBack ( newInstance<MeshpartModule>  ( ) );  // partition mesh into blocks

  chain->pushBack ( newInstance<ShapeModule>     ( ) );
  chain->pushBack ( newInstance<InitModule>      ( ) );
  chain->pushBack ( newInstance<InfoModule>      ( ) );

  // add user defined modules. This allows flexible input since the above is
  // fixed but the solver, the graph module change from problems to problems.

  chain->pushBack ( newInstance<UserconfModule> ( "extraModules" ) );


  // The ControlModule controls when the program stops. The
  // FemViewModule handles visualization.


   ctrl = newInstance<ControlModule> (
      "control",
        NIL,
        ControlModule::FG_MODE
        );

  ctrl ->runWhile ( "i > 0" );
  chain->pushBack ( ctrl );

  // Finally, the chain is wrapped in a ReportModule that prints some
  // overall information about the current calculation.

  //return newInstance<ReportModule> ( "report", chain );

  // Mac OS

 return newInstance<DisplayModule> ( newInstance<ReportModule> ( "report", chain ));
}


//-----------------------------------------------------------------------
//   main
//-----------------------------------------------------------------------


int main ( int argc, char** argv )
{
  return Application::pexec ( argc, argv, & mainModule );
}
