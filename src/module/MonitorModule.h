#ifndef MONITOR_MODULE_H
#define MONITOR_MODULE_H

#include <jive/app/Module.h>


using jem::String;
using jem::util::Properties;
using jive::app::Module;


//-----------------------------------------------------------------------
//   class MonitorModule
//-----------------------------------------------------------------------


class MonitorModule : public Module
{
 public:

  typedef MonitorModule     Self;
  typedef Module            Super;

  static const char*        NODE_PROP;
  static const char*        FORCE_NAMES[3];


  explicit                  MonitorModule

    ( const String&           name = "monitor" );

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
      const Properties&       globdat )        const;


 protected:

  virtual                  ~MonitorModule   ();


 private:

  int                       nodeID_;

};


#endif
