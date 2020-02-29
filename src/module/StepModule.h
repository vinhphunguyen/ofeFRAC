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

#ifndef STEP_MODULE_H
#define STEP_MODULE_H

#include <jive/implict/SolverModule.h>


using jem::Ref;
using jem::String;
using jem::NIL;
using jem::util::Properties;
using jive::app::Module;
using jive::implict::SolverModule;


//-----------------------------------------------------------------------
//   class StepModule
//-----------------------------------------------------------------------


class StepModule : public Module
{
 public:

  typedef StepModule        Self;
  typedef Module            Super;

  static const char*        TYPE_NAME;

  explicit                  StepModule

    ( const String&           name   = "step",
      Ref<SolverModule>       solver = NIL );

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

  virtual                  ~StepModule  ();


 private:

  Ref<SolverModule>         solver_;

  int                       istep_;
  int                       istep0_;

};


#endif
