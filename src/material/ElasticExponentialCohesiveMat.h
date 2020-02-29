/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  This class implements a 2D elastic-exponential cohesive law presented in
 *  "Modelling of cohesive crack growth in concrete structures with the extended
 *  finite element method" of Unger et al, CMAME, 2007.
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 21 December 2007
 *
 */


#ifndef ELASTIC_EXPONENTIAL_COHESIVE_MATERIAL_H
#define ELASTIC_EXPONENTIAL_COHESIVE_MATERIAL_H

#include <jem/base/String.h>
#include <jem/util/Flex.h>

#include "XCohesiveMat.h"


using jem::String;
using jem::util::Flex;


// =======================================================
//  class ElasticExponentialCohesiveMat
// =======================================================


class ElasticExponentialCohesiveMat : public XCohesiveMat
{
 public:

  static const char*      BETA_PROP;
  static const char*      TENSILE_PROP;
  static const char*      ENERGY_PROP;
  static const char*      DUMMY_PROP;
  static const char*      LIMIT_SHIFT_PROP;
  
  /*
   *  constructor
   */

  explicit                ElasticExponentialCohesiveMat

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
      const Properties&     globdat  )         const;

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

    ( Vector&             traction,
      Matrix&             stiff,
      const Vector&       jump0,
      const Vector&       djump,
      int                 mpoint );

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
      const Vector&         tEff,
      const Matrix&         jumpStiff,
      int                   mpoint );

  /*
   * This is used in explicit dynamics  
   */

  virtual void            update

    ( Vector&               traction,
      const Vector&         jump,
      int                   mpoint );

  /*
   *  Called when the Newton Raphson converged
   *  Swap newHist_ to oldHist_ 
   */

  virtual void            commit           ();


  /*
   *  Allocate for history variables if it is not yet 
   *  initialized. Otherwise, extend it by appending at the
   *  end 
   */

  void                    allocPoints

    ( int                   count,
      double                dam = 0. );


  /*
   *  Return history variables at material point ipoint
   */

  inline double           giveHistory       ( int ipoint ) const;

  /*
   *  Return the tensile stress of the material, ft
   */

  inline double           giveTensileStress ( )            const;
  
  void                    initShift 

    ( const int             ipoint,
      const Vector&         traction );

  inline double           giveTensileStrength ( ) const;
  
  virtual double          evalFailure
    
    ( const Vector&         sigmaN,
      int                   mpoint ) const;

  /*
   * Return number of integration points
   */

  virtual int             ipointCount  (  ) const
  {
    return preHist_.loading.size ( ) ;
  }

  virtual bool            justTractionFree

    ( const int             i ) const;
  
  void                    scaleDissipationIncrement

    ( const double     factor,
      const int        ipoint );

 protected:

  virtual                ~ElasticExponentialCohesiveMat   ();

 private:

 private:

  int                     rank_;

  double                  beta_;  // weight of t he  tangential jump
  double                  dummy_; // Kp
  double                  ft_;    // tensile strength
  double                  g1c_;
  double                  delta0_;
  double                  beta2_;
  double                  mftOverGf_;

  // history variable (equivalent crack opening), involve in time

  struct                  hist_
  {
    Flex<double>            eqvOpenMax ;
    Flex<int>               loading;
  };

  hist_                   preHist_;      // history of the previous load step
  hist_                   newHist_;      // history of the current iteration

  Flex<Vector>            shift_;     
  Vector                  damage_;
};

// -------------------------------------------------------------------
//  giveHistory
// -------------------------------------------------------------------

inline double  ElasticExponentialCohesiveMat::giveHistory ( int ipoint ) const
{
  return damage_[ipoint];
}

// -------------------------------------------------------------------
//  giveTensileStress
// -------------------------------------------------------------------

inline double  ElasticExponentialCohesiveMat::giveTensileStress (  ) const
{
  return ft_;
}

inline double  ElasticExponentialCohesiveMat::giveTensileStrength (  ) const
{
  return ft_;
}

#endif 
