/*
 *  Copyright (C) 2009 TU Delft. All rights reserved.
 *  
 *  This class implements a mixed mode traction-seperation law
 *  for cohesive crack models with dummy stiffness (Turon et al. 2006)
 *  with shifted origin to mimick rigid behavior (Hille et al. 2009). 
 *  For use in XFEMModel.
 *  
 *  Author: F.P. van der Meer
 *  Date: April 2009
 *
 */


#ifndef SHIFTED_TURON_C_MATERIAL_WEIBUL_H 
#define SHIFTED_TURON_C_MATERIAL_WEIBUL_H

#include "TuronXCohesiveMat.h"

// =======================================================
//  class TuronXCohesiveWeibulMat
// =======================================================


class TuronXCohesiveWeibulMat : public TuronXCohesiveMat
{
 public:

  typedef TuronXCohesiveMat  Super;

  static const char*      WEIBUL_MOD_PROP;
  static const char*      SIGMA_MIN_PROP;

  /*
   *  constructor 
   */

  explicit                TuronXCohesiveWeibulMat

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

  // ------- overloaded functions ------------------------

  /*
   * Check failure in crack tip element
   */

  double                  evalFailure

    ( const Vector&       sigmaN,
      int                 mpoint )           const;
  
  /*
   * This is used in explicit dynamics  
   */

  virtual void            update

    ( Vector&               traction,
      const Vector&         jump,
      int                   mpoint );

  /*
   * allocate memory to store shift value for each ip, 
   * and call Super::allocate 
   */

  void                    allocPoints

    ( int                   count,
      double                dam = 0. );

  // ------- newly defined functions --------------------

  /**
   * initialize translation
   */

 protected:

  virtual                ~TuronXCohesiveWeibulMat ();


 private:

  double                  m_; // Weibul modulus
  double                  sigmaMin_; 

  Flex<double>            f2ts_;     
  Flex<double>            f6ts_;     
};

#endif 
