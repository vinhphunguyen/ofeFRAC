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
 */

#ifndef FLEX_ARCLEN_MODULE_H
#define FLEX_ARCLEN_MODULE_H

#include <jive/app/Module.h>
#include <jive/model/Model.h>


namespace jive
{
  namespace implict
  {
    class ArclenModule;
    class NonlinModule;
  }
}


using jem::Ref;
using jem::String;
using jem::NIL;
using jem::util::Properties;
using jive::app::Module;
using jive::model::Model;
using jive::implict::NonlinModule;
using jive::implict::SolverModule;
using jive::implict::ArclenModule;



//-----------------------------------------------------------------------
//   class FlexArclenModule
//-----------------------------------------------------------------------


class FlexArclenModule : public Module
{
 public:

  typedef FlexArclenModule  Self;
  typedef Module            Super;

  static const char*        TYPE_NAME;
  static const char*        NONLIN;
  static const char*        ARCLEN;
  static const char*        TMP_NONLIN_PROP;

  explicit                  FlexArclenModule

    ( const String&           name    = "FlexArclen",
      Ref<NonlinModule>       solver1 = NIL,
      Ref<ArclenModule>       solver2 = NIL );

  virtual Status            init

    ( const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual Status            run

    ( const Properties&       globdat );

  virtual void              shutdown

    ( const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       props,
      const Properties&       globdat )    const;

  static Ref<Module>        makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       globdat,
      const Properties&       props );

 protected:

  virtual                  ~FlexArclenModule  ();

 private:

  virtual bool              switchStrategy_

    ( const Properties&       globdat );

  virtual bool              switchToNonlin_

    ( const Properties&       globdat );

  virtual bool              switchToArclen_

    ( const Properties&       globdat );

  virtual bool              reduceStep_

    ( const Properties&       globdat );

  void                      commit_ 
    
    ( const Properties&       globdat );

  void                      continue_ 
    
    ( const Properties&       globdat );

  void                      cancel_
    
    ( const Properties&       globdat );

 private:

  Ref<NonlinModule>         nonlinSolver_;   // used for standard disp control
  Ref<ArclenModule>         arclenSolver_;   // used for arclen control

  Ref<SolverModule>         currentSolver_;

  Ref<Model>                model_;

  int                       istep_;
  int                       istep0_;

  int                       nCancels_;
  int                       nContinues_;

  bool                      nonlinTriedAll_;
  bool                      arclenTriedAll_;
  bool                      converged_;

  // some flags for trying to switch to nonlin after crack propagation

  bool                      doTmpNonlin_;    // general flag from input

  bool                      isTmpNonlin_;    // current state

  bool                      doneTmpNonlin_;  // do only once 
};


#endif
