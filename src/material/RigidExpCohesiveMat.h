/*
 *
 */


#ifndef RIGID_EXP_COH_MATERIAL_H 
#define RIGID_EXP_COH_MATERIAL_H

#include <jem/base/String.h>
#include <jem/util/Flex.h>

#include "XCohesiveMat.h"

using jem::util::Flex;

// =======================================================
//  class RigidExpCohesiveMat
// =======================================================


class RigidExpCohesiveMat : public XCohesiveMat
{
 public:

  static const char*      TENSILE_STRENTH_PROP;
  static const char*      FRACTURE_ENERGY_PROP;
  static const char*      SHEAR_STIFFNESS_PROP;

  explicit                RigidExpCohesiveMat

    ( const int             rank,
      const Properties&     globdat );

  virtual void            configure

    ( const Properties&     props,
      const Properties&     globdat );

  virtual void            getConfig

    ( const Properties&     conf,
      const Properties&     globdat )      const;

  virtual double          evalFailure

    ( const Vector&     sigmaN,
      int               mpoint )           const;
  
  void                    scaleDissipationIncrement

    ( const double     factor,
      const int        ipoint );

  virtual void            update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump,
      int                   mpoint );
  
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


  void                    allocPoints

    ( int                   count,
      double                dam = 0. );

  void                    deallocPoints

    ( int                   count );

  virtual bool            justTractionFree

    ( const int             i ) const;

  void                    initShift 

    ( const int             i, 
      const Vector&         traction );
    
  inline double           giveTensileStrength ( ) const;
  
  inline double           giveHistory         
      
    ( const int             i ) const;

 protected:

  virtual                ~RigidExpCohesiveMat ();

 protected:

  Flex<double>           ft_;
  double                 gf_;
  double                 dint_;
  double                 ft0_;

  Flex<double>           damage_;       // for visualising cracks in a continuum 
                                        // formulation.
  
  struct                  hist_
  {
    Flex<double>            opening;     
    Flex<int>               loading;     
    Flex<double>            dissipation;
  };

  hist_                   preHist_;      // history of the previous load step
  hist_                   newHist_;      // history of the current iteration
  hist_*                  latestHist_;   // points to latest hist

};

inline double  RigidExpCohesiveMat::giveHistory ( const int ip ) const
{
  return damage_[ip];
}

inline double  RigidExpCohesiveMat::giveTensileStrength (  ) const
{
  return ft0_;
}

#endif 
