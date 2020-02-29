#ifndef STEP_MODULE2_H
#define STEP_MODULE2_H

#include <jive/implict/SolverModule.h>


using jem::Ref;
using jem::String;
using jem::util::Properties;
using jive::app::Module;
using jive::implict::SolverModule;


//-----------------------------------------------------------------------
//   class StepModule
//-----------------------------------------------------------------------


class StepModule2 : public Module
{
 public:

  static const char*        STEP_REDUCTION;
  static const char*        TYPE_NAME;

  typedef StepModule2       Self;
  typedef Module            Super;


  explicit                  StepModule2

    ( const String&           name,
      Ref<SolverModule>       solver );

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

  virtual                  ~StepModule2  ();


 private:

  Ref<SolverModule>         solver_;

  int                       istep_;
  int                       istep0_;

};


#endif
