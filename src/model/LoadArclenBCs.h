/*
 *  Copyright (C) 2012, Ton Duc Thang University.
 *
 * This class implements a model to define imposed external forces
 * on certain nodes in a computation with a combined load control
 * and arclength control. It is used as a child of class 
 * LoadArclenModel and is able to switch between the two different
 * strategies.
 *
 * Author(s):
 *
 *   Vinh Phu Nguyen, nvinhphu@gmail.com
 *   Saigon, Vietnam October 2012
 *
 */

#ifndef LOAD_ARCLEN_BCS_H
#define LOAD_ARCLEN_BCS_H


#include <jem/base/array/utilities.h>
#include <jem/base/array/select.h>
#include <jem/base/System.h>
#include <jem/util/Properties.h>
#include <jem/base/Error.h>

#include <jive/fem/NodeGroup.h>
#include <jive/model/Actions.h>
#include <jive/util/Assignable.h>

namespace jive
{
  namespace util
  {
    class Constraints;
    class XDofSpace;
  }
}

using namespace jem;

using jem::Error;
using jem::util::Properties;
using jem::io::endl;
using jive::Vector;
using jive::IntVector;
using jive::StringVector;
using jive::fem::NodeGroup;
using jive::fem::NodeSet;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::util::XDofSpace;

// ======================================================
//     class LoadArclenBCs
// ======================================================

class LoadArclenBCs : public Object
{
 public:

  typedef LoadArclenBCs      Self;

  static const char*         TYPE_NAME;
  static const char*         NODES_PROP;
  static const char*         DOF_PROP;
  static const char*         LOAD_PROP;
  static const char*         INCR_PROP;

                            LoadArclenBCs ();

  virtual void         configure

    ( const Properties&   myProps,
      const Properties&   globdat );

  virtual void         getConfig

    ( const Properties&   myConf,
      const Properties&   globdat )      const;

  void                 toLoadControl

    ( const Properties&  params,
      const Properties&  globdat );

  void                 toArclControl

    ( const Properties&  params,
      const Properties&  globdat );

  void                 getUnitLoad

    ( const Properties&  params,
      const Properties&  globdat );

  void                 computeExternalForce

    ( const Properties&  params,
      const Properties&  globdat );

  void                 commit

    ( const bool         isDispControl,
      const Properties&  globdat );

  int                  getLoadedDofs

    ( IntVector&          idofs )        const;

  inline  void         setLoadIncr

    ( double val );

  inline double        getLoadIncr  ()  const;
  inline double        getLoadScale ()  const;
  inline double        getLoadScale0()  const;

 private:

  void                 initializeRefLoad_ ( const Properties& globdat );

 protected:

  virtual              ~LoadArclenBCs ();

 private:

  Assignable<NodeSet>     nodeSet_;

  Ref<XDofSpace>          dofs_;

  double                  loadIncr_;
  double                  loadScale_; 
  double                  loadScale0_;

  int                     nnodes_;

  IntVector               nodes_;
  IntVector               index_;
  Vector                  refLoad_; //reference load vector g0
  StringVector            dofTypes_;
};

// ======================================================
//     implementation of inline functions
// ======================================================


inline void LoadArclenBCs::setLoadIncr

 ( double val )
{
  loadIncr_ = val;
}

inline double LoadArclenBCs::getLoadScale  () const 
{ 
  return loadScale_; 
}

inline double LoadArclenBCs::getLoadScale0 () const 
{ 
  return loadScale0_; 
}


#endif



