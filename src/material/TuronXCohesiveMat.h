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


#ifndef SHIFTED_TURON_C_MATERIAL_H 
#define SHIFTED_TURON_C_MATERIAL_H

#include "XCohesiveMat.h"
#include "TuronCohesiveMaterial.h"

// =======================================================
//  class TuronXCohesiveMat
// =======================================================


class TuronXCohesiveMat : public TuronCohesiveMaterial,
                          public XCohesiveMat
{
 public:

  typedef TuronCohesiveMaterial  Super;

  static const char*      LIMIT_SHIFT_PROP;

  /*
   *  constructor 
   */

  explicit                TuronXCohesiveMat

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

  virtual double          evalFailure

    ( const Vector&     sigmaN,
      int               mpoint )           const;
  
  /*
   * scale the latest dissipation increment for this ipoint with factor
   */

  void                    scaleDissipationIncrement

    ( const double     factor,
      const int        ipoint );

  /*
   * give translated jump vector to Super::update
   */

  virtual void            update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump,
      int                   mpoint );
  /*
   * This is used in explicit dynamics  
   */

  virtual void            update

    ( Vector&               traction,
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
   * allocate memory to store shift value for each ip, 
   * and call Super::allocate 
   */

  void                    allocPoints

    ( int                   count,
      double                dam = 0. );

  void                    deallocPoints

    ( int                   count );

    
  // ------- newly defined functions --------------------

  /**
   * initialize translation
   */

  void                    initShift 

    ( const int             ipoint,
      const Vector&         traction );

  /**
   * get shift
   */

  inline Vector           getShift

    ( int                   ipoint ) const;

  inline double           giveTensileStrength ( ) const;

 protected:

  virtual                ~TuronXCohesiveMat ();


 protected:

  bool                    limitShift_;
  Flex<Vector>            shift_;     
};

// -------------------------------------------------------------------
//  getShift
// -------------------------------------------------------------------

inline Vector  TuronXCohesiveMat::getShift ( int ipoint ) const
{
  return shift_[ipoint];
}

inline double  TuronXCohesiveMat::giveTensileStrength (  ) const
{
  return f2t_;
}

#endif 
