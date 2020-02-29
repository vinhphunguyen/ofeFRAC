/*
 *  Copyright (C) 2008 TU Delft. All rights reserved.
 *  
 *  This class implements a mixed mode traction-seperation law
 *  for cohesive crack models with dummy stiffness (Turon et al. 2006). 
 *  For use in InterfaceModel.
 *  
 *  Author: F.P. van der Meer
 *  Date: June 2008
 *
 */


#ifndef TURON_COHESIVE_MATERIAL_H 
#define TURON_COHESIVE_MATERIAL_H

#include <jem/base/String.h>
#include <jem/util/Flex.h>

#include "CohesiveMaterial.h"
//#include "TractionLaw.h"

using jem::String;
using jem::util::Flex;


// =======================================================
//  class TuronCohesiveMaterial
// =======================================================


class TuronCohesiveMaterial : public virtual CohesiveMaterial
{
 public:

  static const char*      F2T_PROP;
  static const char*      F6_PROP;
  static const char*      G_I_PROP;
  static const char*      G_II_PROP;
  static const char*      G_III_PROP;
  static const char*      DUMMY_PROP;
  static const char*      ETA_PROP;
  static const double     MIN_SECANT;
  static const double     OMEGA_MAX;
  static const double     EPS;

  /*
   *  constructor
   */

  explicit                TuronCohesiveMaterial

    ( const int             rank,
      const Properties&     globdat );

  /*
   *  configure and getConfig (from input file)
   */

  virtual void            configure

    ( const Properties&     props,
      const Properties&     globdat );

  virtual void            getConfig

    ( const Properties&     conf,
      const Properties&     globdat )         const;

  /*
   *  compute the traction (t) and cohesive tangent stiffness matrix
   *  (stiff) at material point mpoint given the displacement jump (jump)
   *   jump[0] = crack opening displacement
   *   jump[1] = crack sliding displacement
   */

  virtual void            update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump,
      int                   mpoint );
  
  virtual void            update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump0,
      const Vector&         djump,
      int                   mpoint );

  virtual void            elasticUpdate

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump );

  /*
   * update for moonen material
   */

  virtual void            update

    ( Vector&               traction,
      Matrix&               stiff,
      const Matrix&         jumpStiff,
      const Vector&         tEff,
      int                   mpoint );

  /*
   *  Called when the Newton Raphson converged
   *  Swap newHist_ to oldHist_ 
   */

  virtual void            commit           ();
  
  virtual bool            justTractionFree

    ( const int             i ) const;


  /**
   *  Allocate for history variables if it is not yet 
   *  initialized. Otherwise, extend it by appending at the end.
   *  Initial value of damage is optional
   */

  void                    allocPoints

    ( int                   count,
      double                dam = 0. );

  virtual void            deallocPoints

    ( int                   count );

  /*
   *  Return history variables at material point ipoint
   */

  inline double           giveDissipation   ( int ipoint ) const;
  
  inline virtual double   giveTensileStrength  (  ) const;

  inline double           giveHistory       ( int ipoint ) const;

  inline int              isLoading         ( int ipoint ) const;

  virtual int             wasLoading        ( int ipoint ) const;

  /*
   * Return number of integration points
   */

  virtual int             ipointCount  (  ) const
  {
    return preHist_.loading.size ( ) ;
  }

  double                  getDummy () const {return dummy_;}

 protected:

  virtual                ~TuronCohesiveMaterial   ();

  double                  getDelta0_ ( const double beta ) const;

  void                    allocPoints_

     ( const int             count,
       const double          dam,
       const int             loading );


 protected:

  double                  f2t_;
  double                  f6_;
  double                  gI_;
  double                  gII_;
  double                  gIII_;
  double                  eta_;

  mutable double          f2t2_;
  mutable double          f62_;
  double                  deltaN02_;
  double                  deltaS02_;
  double                  deltaN0F_;
  double                  deltaS0F_;

  double                  dummy_;

  // history variable (equivalent crack opening), involve in time

  struct                  hist_
  {
    Flex<double>            damage;     
    Flex<int>               loading;     
    Flex<double>            dissipation;
  };

  hist_                   preHist_;      // history of the previous load step
  hist_                   newHist_;      // history of the current iteration
  hist_*                  latestHist_;   // points to latest hist
};

// -------------------------------------------------------------------
//  giveHistory
// -------------------------------------------------------------------

inline double  TuronCohesiveMaterial::giveHistory ( int ipoint ) const
{
  return latestHist_->damage[ipoint];
}

// -------------------------------------------------------------------
//  giveTensileStrength
// -------------------------------------------------------------------

inline double  TuronCohesiveMaterial::giveTensileStrength () const
{
  return f2t_;;
}

// -------------------------------------------------------------------
//  isLoading
// -------------------------------------------------------------------

inline int     TuronCohesiveMaterial::isLoading ( int ipoint ) const
{
  return latestHist_->loading[ipoint];
}

// -------------------------------------------------------------------
//  wasLoading
// -------------------------------------------------------------------

inline int     TuronCohesiveMaterial::wasLoading ( int ipoint )  const
{
  return preHist_.loading[ipoint];
}

// -------------------------------------------------------------------
//  giveDissipation
// -------------------------------------------------------------------

inline double  TuronCohesiveMaterial::giveDissipation ( int ipoint ) const
{
  return latestHist_->dissipation[ipoint];
}


#endif 
