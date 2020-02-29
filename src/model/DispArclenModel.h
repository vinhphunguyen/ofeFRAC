/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implements a model that combines
 *  displacement and arclen control in one model. This model
 *  must be used with FlexArclenModule to have a
 *  flexible path following solver.
 *
 *  Basic ideas:
 *
 *    1. Disp control: The simulation starts with displacement control
 *       where some nodes are constrained. During
 *       this stage, in FlexArclenModule, the NonlinModule is being
 *       used. This is the case until divergence occurs (snapback).
 *       Then, one node is defined as master and others as its slave.
 *       A force is then applied on this master node so that the
 *       unit external force required by ArclenModule can be
 *       constructed. 
 *
 *   2.  Arclen control (based on energy released). Now, in FlexArclenModule,
 *       the ArclenModule is active so that lambda is now an unknown.
 *       This is the case until hardening branch is detected, then
 *       switch back to (1), Load control
 *
 *   Constrained dofs are managed by a child object which is an instance
 *   of class DispArclenBCs
 *
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 31 January 2009
 *
 */


#ifndef DISP_ARCLEN_MODEL_H
#define DISP_ARCLEN_MODEL_H

#include <jem/util/Flex.h>
#include <jem/io/NumberFormat.h>
#include <jem/io/Writer.h>
#include <jive/Array.h>
#include <jive/model/Model.h>
#include <jive/fem/NodeGroup.h>
#include <jive/model/Actions.h>
#include <jive/model/ModelFactory.h>

#include "DispArclenBCs.h"

namespace jive
{
  namespace util
  {
    class Constraints;
    class XDofSpace;
  }
}

namespace jive
{
  namespace model
  {
    class Model;
  }
}

namespace jive
{
  namespace algebra
  {
    class VectorSpace;
  }
}

using namespace jem;

using jem::util::Properties;
using jem::util::Flex;
using jem::io::NumberFormat;
using jem::io::Writer;
using jive::Vector;
using jive::IntVector;
using jive::algebra::VectorSpace;
using jive::model::Model;
using jive::StringVector;
using jive::fem::NodeGroup;
using jive::fem::NodeSet;
using jive::model::Model;
using jive::util::XDofSpace;
using jive::util::DofSpace;
using jive::util::Constraints;

//-----------------------------------------------------------------------
//   class DispArclenModel
//-----------------------------------------------------------------------


class DispArclenModel : public Model
{
 public:

  typedef DispArclenModel   Self;
  typedef Model             Super;

  static const char*        TYPE_NAME;

  static const char*        CONSTRAINT_PROP;
  static const char*        OPT_ITER_PROP;
  static const char*        SWT_ITER_PROP;
  static const char*        SWT_ENER_PROP;
  static const char*        MIN_INCR_PROP;
  static const char*        MAX_INCR_PROP;
  static const char*        DISP_INCR_PROP;
  static const char*        MIN_DISP_PROP;
  static const char*        LOAD_SCALE_PROP;
  static const char*        REDUCTION_PROP;
  static const char*        PREFER_DC_PROP;

  explicit                  DispArclenModel

    ( const String&           name  = "arclen",
      const Ref<Model>&       child = NIL );

  virtual Model*            findModel

    ( const String&           name )         const;

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       conf,
      const Properties&       globdat )      const;

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );

  void                      setMaxIter

    ( int                     count );

  inline int                getMaxIter      () const;

  void                      setLoadIncr

    ( double                  incr );

  inline double             getLoadIncr     () const;

  void                      setIncrRange

    ( double                  minIncr,
      double                  maxIncr );

  inline double             getMinIncr      () const;
  inline double             getMaxIncr      () const;


  static Ref<Model>         makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );


 protected:

  virtual                  ~DispArclenModel  ();


 private:

  void                      init_

    ( const Properties&       globdat );

  void                      initLoad_

    ( const Properties&       globdat );

  void                      evalArcFunc_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      getUnitLoad_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      commit_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      checkCommit_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      checkSwitch_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      reduceStep_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      connect_        ();

  void                      dofsChanged_    ();
  void                      consChanged_    ();

  double                    getReleasedEnergy_

    ( const Properties&  params,
      const Properties&  globdat );

  void                      toArclControl_

    ( const Properties&  params,
      const Properties&  globdat );

  void                      toDispControl_

    ( const Properties&  params,
      const Properties&  globdat );

 private:

  static const int          U_LOAD_;
  static const char*        DELTA_STATE_;


  Ref<DispArclenBCs>        child_;
  Ref<DofSpace>             dofs_;
  Ref<VectorSpace>          vspace_;

  int                       updated_;
  Vector                    load_;

  // stored quantities

  double                    arcLength_;
  double                    lastArcl_;
  double                    energy0_;
  double                    totalDiss_;

  // flags

  bool                      isDispControl_;
  bool                      triedLarge_;
  bool                      onceDown_;

  // fancy info to the scrren

  NumberFormat              nformat_;
  Writer&                   out_;

  /* the following members are constant input variables 
   *
   * optIter:     optimal number of iterations (for adaptive stepping)
   * swtIter:     number of iterations for which to switch to arc-length
   * swtEnergy:   amount of energy for which to switch to arc-length
   * reduction:   factor for reduction of increments
   * minIncr:     minimum energy increment
   * maxIncr:     maximum energy increment
   * dispIncr:    initial displacement increment
   * minDispIncr: minimum (absolute) displacement increment
   * preferDisp:  switch to displacement control when hardening
   *
   */

  int                       optIter_;
  int                       swtIter_;
  double                    swtEner_;
  double                    reduction_;
  double                    minIncr_;
  double                    maxIncr_;
  double                    dispIncr0_; 
  double                    minDispIncr_; 
  double                    gC_; 
  bool                      preferDisp_;
};




//#######################################################################
//   Implementation
//#######################################################################


//-----------------------------------------------------------------------
//   getMaxIter
//-----------------------------------------------------------------------


inline int DispArclenModel::getMaxIter () const
{
  return optIter_;
}


//-----------------------------------------------------------------------
//   get(Min|Max)Incr
//-----------------------------------------------------------------------


inline double DispArclenModel::getMinIncr () const
{
  return minIncr_;
}


inline double DispArclenModel::getMaxIncr () const
{
  return maxIncr_;
}



#endif
