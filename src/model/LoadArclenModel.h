/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implements a model that combines
 *  load and arclen control in one model. This model
 *  must be used with FlexArclenModule to have a
 *  flexible path following solver.
 *
 *  Basic ideas:
 *
 *    1. Load control: a so-called unit external force vector
 *       is computed once. The load scale lambda is updated
 *       by a constant amount (can be changed later). During
 *       this stage, in FlexArclenModule, the NonlinModule is being
 *       used. This is the case until divergence occurs (pass the peak)
 *
 *   2.  Arclen control (based on energy released). Now, in FlexArclenModule,
 *       the ArclenModule is active so that lambda is now an unknown.
 *       This is the case until hardening branch is detected, then
 *       switch back to (1), Load control
 *
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 January 2009
 *
 */


#ifndef LOAD_ARCLEN_MODEL_H
#define LOAD_ARCLEN_MODEL_H

#include <jem/util/Flex.h>
#include <jem/io/NumberFormat.h>
#include <jem/io/Writer.h>
#include <jive/Array.h>

#include "LoadArclenBCs.h"

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
using jive::algebra::VectorSpace;
using jive::util::DofSpace;
using jive::util::Constraints;
using jive::model::Model;


//-----------------------------------------------------------------------
//   class LoadArclenModel
//-----------------------------------------------------------------------


class LoadArclenModel : public Model
{
 public:

  typedef LoadArclenModel   Self;
  typedef Model             Super;

  static const char*        TYPE_NAME;

  static const char*        MODEL_PROP;
  static const char*        ARC_FUNC_PROP;
  static const char*        OPT_ITER_PROP;
  static const char*        SWT_ITER_PROP;
  static const char*        SWT_ENER_PROP;
  static const char*        MIN_INCR_PROP;
  static const char*        MAX_INCR_PROP;
  static const char*        LOAD_INCR_PROP;
  static const char*        LOAD_SCALE_PROP;
  static const char*        EXIT_FRAC_PROP;
  static const char*        REDUCTION_PROP;

  enum                      ArcFunc
  {
                               ERC     // Energy Release Control
  };


  explicit                  LoadArclenModel

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


  void                      setArcFunc

    ( ArcFunc                 func );

  inline ArcFunc            getArcFunc      () const;

  void                      setMaxIter

    ( int                     count );

  inline int                getMaxIter      () const;

  void                      setLoadIncr

    ( double                  incr );

  inline double             getLoadIncr     () const;

  void                      setLoadScale

    ( double                  scale );

  inline double             getLoadScale    () const;

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

  virtual                  ~LoadArclenModel  ();


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

  void                      toLoadControl_

    ( const Properties&  params,
      const Properties&  globdat );

  void                      toArclControl_

    ( const Properties&  params,
      const Properties&  globdat );

 private:

  static const int          U_LOAD_;

  static const char*        DELTA_STATE_;


  Ref<LoadArclenBCs>        child_;
  Ref<DofSpace>             dofs_;
  Ref<Constraints>          cons_;
  Ref<VectorSpace>          vspace_;

  int                       istep_;
  int                       updated_;

  // name of the constrained function, currently, there is only one

  ArcFunc                   arcFunc_;

  // the name of table contains the dofs where the constraints are applied.
  // Most of the time, all dofs are used.

  String                    wtblName_;

  String                    stepAdjust_;

  Vector                    load_;

  Flex<double>              stepSizes_;

  // optimal number of iterations

  int                       optIter_;

  // number of iterations for which to switch to arc-length
  // and switch energy

  int                       swtIter_;
  double                    swtEner_;

  // stop computation when ( load << maximum of load over time )
  // varies between [0,1]

  double                    exitFraction_; 
  double                    reduction_;

  // arc-length bounds

  double                    minIncr_;
  double                    maxIncr_;

  // stored quantities

  double                    loadIncr_;
  double                    loadScale_;
  double                    oldScale_;
  double                    maxLoadScale_; // stored highest loadScale_ value
  double                    arcLength_;
  double                    lastArcl_;
  double                    energy0_;
  double                    gC_;

  // flags

  bool                      isLoadControl_;
  bool                      triedLarge_;
  bool                      onceDown_;

  NumberFormat              nformat_;
  Writer&                   out_;
};




//#######################################################################
//   Implementation
//#######################################################################

//-----------------------------------------------------------------------
//   getArcFunc
//-----------------------------------------------------------------------


inline LoadArclenModel::ArcFunc

  LoadArclenModel::getArcFunc () const

{
  return arcFunc_;
}


//-----------------------------------------------------------------------
//   getMaxIter
//-----------------------------------------------------------------------


inline int LoadArclenModel::getMaxIter () const
{
  return optIter_;
}


//-----------------------------------------------------------------------
//   getLoadIncr
//-----------------------------------------------------------------------


inline double LoadArclenModel::getLoadIncr () const
{
  return loadIncr_;
}


//-----------------------------------------------------------------------
//   getLoadScale
//-----------------------------------------------------------------------


inline double LoadArclenModel::getLoadScale () const
{
  return loadScale_;
}


//-----------------------------------------------------------------------
//   get(Min|Max)Incr
//-----------------------------------------------------------------------


inline double LoadArclenModel::getMinIncr () const
{
  return minIncr_;
}


inline double LoadArclenModel::getMaxIncr () const
{
  return maxIncr_;
}



#endif
