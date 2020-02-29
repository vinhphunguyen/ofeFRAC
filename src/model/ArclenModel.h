#ifndef ARCLEN_MODEL_H
#define ARCLEN_MODEL_H

#include <jive/Array.h>
#include <jive/algebra/VectorSpace.h>
#include <jive/util/DofSpace.h>
#include <jive/util/Constraints.h>
#include <jive/model/Model.h>


using namespace jem;

using jem::util::Properties;
using jive::Vector;
using jive::algebra::VectorSpace;
using jive::util::DofSpace;
using jive::util::Constraints;
using jive::model::Model;


//-----------------------------------------------------------------------
//   class ArclenModel
//-----------------------------------------------------------------------


class ArclenModel : public Model
{
 public:

  typedef ArclenModel       Self;
  typedef Model             Super;

  static const char*        TYPE_NAME;

  static const char*        MODEL_PROP;
  static const char*        MAX_ITER_PROP;
  static const char*        MIN_INCR_PROP;
  static const char*        MAX_INCR_PROP;
  static const char*        ARC_LENGTH_PROP;
  static const char*        LOAD_SCALE_PROP;
  static const char*        WGT_TABLE_PROP;


  explicit                  ArclenModel

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

  virtual                  ~ArclenModel  ();


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

  int                       maxIter_;

  double                    minIncr_;
  double                    maxIncr_;
  double                    loadScale_;
  double                    arcLength_;

  String                    wtblName_;

  Vector                    vtmp_;
  Vector                    load_;
  Vector                    weights_;

};




//#######################################################################
//   Implementation
//#######################################################################

//-----------------------------------------------------------------------
//   getMaxIter
//-----------------------------------------------------------------------


inline int ArclenModel::getMaxIter () const
{
  return maxIter_;
}


//-----------------------------------------------------------------------
//   getLoadIncr
//-----------------------------------------------------------------------


inline double ArclenModel::getLoadIncr () const
{
  return arcLength_;
}


//-----------------------------------------------------------------------
//   getLoadScale
//-----------------------------------------------------------------------


inline double ArclenModel::getLoadScale () const
{
  return loadScale_;
}


//-----------------------------------------------------------------------
//   get(Min|Max)Incr
//-----------------------------------------------------------------------


inline double ArclenModel::getMinIncr () const
{
  return minIncr_;
}


inline double ArclenModel::getMaxIncr () const
{
  return maxIncr_;
}



#endif
