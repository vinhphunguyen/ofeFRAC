/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implements the energy release control method
 *  (Miguel Gutierrez 2004). This solver is best suited to
 *  problem involving snapback behavior. It starts with load
 *  control and then switch to, automatically, based on the
 *  sudden change of iterations, energy based arclength.
 *  
 *  This solver when used in combination with the StepModule
 *  is able to capture the whole complicated behavior.
 *
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *
 *  Improved by Frans van der Meers
 *
 */


#ifndef ENERGY_ARCLEN_MODEL_H
#define ENERGY_ARCLEN_MODEL_H

#include <jem/util/Flex.h>
#include <jem/io/NumberFormat.h>
#include <jive/Array.h>

namespace jive
{
  namespace util
  {
    class DofSpace;
    class Constraints;
  }
}

namespace jive
{
  namespace algebra
  {
    class VectorSpace;
  }
}

namespace jive
{
  namespace model
  {
    class Model;
  }
}


using namespace jem;

using jem::util::Properties;
using jem::util::Flex;
using jem::io::NumberFormat;
using jive::Vector;
using jive::algebra::VectorSpace;
using jive::util::DofSpace;
using jive::util::Constraints;
using jive::model::Model;


//-----------------------------------------------------------------------
//   class EnergyArclenModel
//-----------------------------------------------------------------------


class EnergyArclenModel : public Model
{
 public:

  typedef EnergyArclenModel Self;
  typedef Model             Super;

  static const char*        TYPE_NAME;

  static const char*        MODEL_PROP;
  static const char*        ARC_FUNC_PROP;
  static const char*        OPT_ITER_PROP;
  static const char*        SWT_ITER_PROP;
  static const char*        MIN_INCR_PROP;
  static const char*        MAX_INCR_PROP;
  static const char*        LOAD_INCR_PROP;
  static const char*        LOAD_SCALE_PROP;
  static const char*        EXIT_FRAC_PROP;
  static const char*        WGT_TABLE_PROP;
  static const char*        STEP_ADJUST_PROP;

  enum                      ArcFunc
  {
                               ERC     // Energy Release Control
  };


  explicit                  EnergyArclenModel

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

  virtual                  ~EnergyArclenModel  ();


 private:

  void                      init_

    ( const Properties&       globdat );

  void                      initLoad_

    ( const Properties&       globdat );

  void                      initWeights_

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

  void                      reduceStep_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      connect_        ();

  void                      dofsChanged_    ();
  void                      consChanged_    ();


 private:

  static const int          U_LOAD_;
  static const int          U_WEIGHTS_;

  static const char*        DELTA_STATE_;


  Ref<Model>                child_;
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
  Vector                    weights_;

  Flex<double>              stepSizes_;

  // optimal number of iterations

  int                       optIter_;

  // number of iterations for which to switch to arc-length

  int                       swtIter_;

  // stop computation when ( load << maximum of load over time )
  // varies between [0,1]

  double                    exitFraction_; 

  // arc-length bounds

  double                    minIncr_;
  double                    maxIncr_;

  // stored quantities

  double                    loadIncr_;
  double                    loadScale_;
  double                    oldScale_;
  double                    maxLoadScale_; // stored highest loadScale_ value
  double                    arcLength_;
  double                    lastArcL_;
  double                    energy0_;
  double                    switchBack_; // val to define when to switch back to load control 

  // flags

  bool                      isLoadControl_;
  bool                      triedLarge_;
  bool                      down_;

  NumberFormat              nformat_;
};




//#######################################################################
//   Implementation
//#######################################################################

//-----------------------------------------------------------------------
//   getArcFunc
//-----------------------------------------------------------------------


inline EnergyArclenModel::ArcFunc

  EnergyArclenModel::getArcFunc () const

{
  return arcFunc_;
}


//-----------------------------------------------------------------------
//   getMaxIter
//-----------------------------------------------------------------------


inline int EnergyArclenModel::getMaxIter () const
{
  return optIter_;
}


//-----------------------------------------------------------------------
//   getLoadIncr
//-----------------------------------------------------------------------


inline double EnergyArclenModel::getLoadIncr () const
{
  return loadIncr_;
}


//-----------------------------------------------------------------------
//   getLoadScale
//-----------------------------------------------------------------------


inline double EnergyArclenModel::getLoadScale () const
{
  return loadScale_;
}


//-----------------------------------------------------------------------
//   get(Min|Max)Incr
//-----------------------------------------------------------------------


inline double EnergyArclenModel::getMinIncr () const
{
  return minIncr_;
}


inline double EnergyArclenModel::getMaxIncr () const
{
  return maxIncr_;
}



#endif
