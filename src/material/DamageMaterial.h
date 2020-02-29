/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the isotropic elastic-damage material
 *  This represents the material at a point in space.
 *  It is implemented in such a way that can be used for any
 *  discretisation methods, say finite elements, EFG and so on.
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *
 */

#ifndef DAMAGEMATERIAL_H
#define DAMAGEMATERIAL_H

#include <jive/Array.h>
#include <jem/base/String.h>
#include <jem/util/Flex.h>
#include <vector>

#include "HookeMaterial.h"
#include "VonMisesStrain.h"


using jem::String;
using jem::util::Flex;
using jive::Vector;
using jive::IntVector;
using std::vector;


enum SofteningLaw {
  Perfect,
  Linear,
  Exponential3Params,
  Exponential2Params,
  Exponential2ParamsReg,
  Hyperbolic,
  Power,
  ExpoEnergy
};

enum EquiStrainDef {
  Mazars,
  VonMises,
  Rankine,
  StrainEnergy
};

//-----------------------------------------------------------------------
//   class ElastoDamageMaterial
//-----------------------------------------------------------------------

/*  The IsoDamageMaterial implements an isotropic linear elastic based damage
 *  material: sigma = (1-omega):D:epsilon
 *  
 * Start: 18 August 2007:
 * 
 * 1. Softening laws implemented:
 *  + Perfect softening (16/10/2007)
 *  + Linear softening
 *  + Modified power softening
 *  + Exponential softening
 *  + Exponential softening
 *  + Hyperbolic softening
 * 2. Equivalent strain definition
 *  + Positive strain or Mazars
 *  + Modified von Mises
 *  + Rankine-type 
 * 3. Energy-based regularised softening law for exponential law done
 *  + Milan Jirasek's version with Gf
 *  + de Vree's version with crack band with (solution depends on this param!!!)
 */

class DamageMaterial : public Material
{
 public:

  typedef  DamageMaterial Self;

  static const char*      SOFTENING_PROP;
  static const char*      EQUISTRAIN_PROP;

  static const char*      KAPPAI_PROP;
  static const char*      KAPPAC_PROP;
  static const char*      ALPHA_PROP;
  static const char*      BETA_PROP;
  static const char*      ETA_PROP;
  static const char*      B_PROP;
  static const char*      LENGTH_PROP;

  static const char*      TENSILE_PROP;
  static const char*      FRACTURE_ENERGY_PROP;

  static const char*      CRACK_WIDTH_PROP;
  static const char*      REMOVE_DAMAGE_PROP;

  static const char*      MAZARS_EQUI_STRAIN;
  static const char*      MISES_EQUI_STRAIN;
  static const char*      RANKINE_EQUI_STRAIN;
  static const char*      ENERGY_STRAIN;

  static const char*      PERFECT_SOFTENING;
  static const char*      LINEAR_SOFTENING;
  static const char*      EXPONENT1_SOFTENING; 
  static const char*      EXPONENT2_SOFTENING;
  static const char*      EXPONENT3_SOFTENING;
  static const char*      POWER_SOFTENING;
  static const char*      HYPERBOLIC_SOFTENING;
  static const char*      EXPOENERGY_SOFTENING;

  static const double     CRITICAL_DAMAGE;

  static vector<String>   equiStrainDefs;
  static vector<String>   softeningLawDefs;

  explicit                DamageMaterial

    ( int                   rank,
      const Properties&     globdat );

  virtual void            configure

    ( const Properties&     props,
      const Properties&     globdat );

  virtual void            getConfig

    ( const Properties&     conf,
      const Properties&     globdat )         const;

  // update the stress, the reduced elastic modulii
  // given the updated strain at GP ipoint

  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      int                   ipoint );

  // the same as above but overloaded for mesh adjusted damage model
  // he is the element characteristic element length 

  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      int                   ipoint,
      double                he );
  
  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain0,
      const Vector&         dstrain,
      int                   ipoint,
      double                he );


  /*
   *  Called when the Newton Raphson converged
   *  Swap newHist_ to oldHist_ 
   */

  virtual void            commit           ();

  virtual void            allocPoints

    ( int                   ipCount );

  /*
   * Damage class overwrites this function (defined as normal virtual
   * function in Material class). Other material classes that do not deal with
   * localisation do not have to worry about this function.
   **/

  double                  checkLocalisation

    (       Vector&         normal,
      const Vector&         stress,
      const Matrix&         tangent,
            int             ip ) const;

  // ----------------------------------------------------------
  //   return the internal variable: kappa variable
  // ----------------------------------------------------------
  
  inline double           giveHistory      ( int ipoint ) const;
  inline Flex<double>     giveHistory      (            ) const;

  // return omega, used for visualization the damage evolution

  inline Vector           giveOmega        ()             const;
  inline double           giveOmega        ( int ipoint ) const;

  inline double           giveLengthScale  ()             const;
  inline double           giveRho          ()             const;

  // compute the equivalent strain given a strain vector

  double                  getEquiStrain

    ( const Vector& strain ) ;

  // compute derivatives of equivalent strain w.r.t strain
  // vector

  Vector                  getdEpsBardEps

    ( const Vector& strain ) ;

  Matrix                  giveElasticMatrix ()const;

  int                     isLoading         ( int ip )    const;

  // compute the derivative of Omega w.r.t kappa

  double                  getdOmegadKappa

    ( double kappa )                                      const;

  double                  getdOmegadKappa

    ( double kappa, double he )                           const;

  /*
   * return true if the material point mpoint is completely damaged
   */
 
  bool                    isFullyDamaged     ( int mpoint ) const;
 
 protected:

  virtual                ~DamageMaterial   ();

 private:

  // compute the damage variable Omega from kappa

  double                  damageEvolution_

    ( double kappa )                 const;

  // overloaded version for mesh adjusted damage model

  double                  damageEvolution_

    ( double kappa,
      double he )                    const;

  // ---------------------------------------------------------------
  //  compute the derivatives of equivalent strain
  //  w.r.t the strain vector
  //  Mazars formulas
  // ---------------------------------------------------------------

  Vector                  getdEquiStraindEpsilon1_

    ( const Vector&   strain );


  // ---------------------------------------------------------------
  //  compute the derivatives of equivalent strain
  //  w.r.t the strain vector
  //  Rankine formulas
  // ---------------------------------------------------------------

  Vector                  getdEquiStraindEpsilon3_

   ( const Vector&   stress,
     const double    young,
     const Matrix&   De );

  
 private:

  // parameters for softening and equivalent strain def

  SofteningLaw            softening_;
  EquiStrainDef           equiStrain_;

  double                  kappaI_;    // damage threshold
  double                  kappaC_;    // softening modulus
  double                  alpha_;
  double                  beta_;
  double                  b_;
  double                  s_;

  double                  eta_;
  
  double                  threshold_;

  double                  c_;

  double                  ft_;        // tensile strength
  double                  gf_;        // fracture energy

  double                  lambda_;    // crack band width 

  // history variables

  struct                  hist_
  {
    Flex<double>            eqveps ;     // non-local equivalent strain history variable
    Flex<int>               loading;     // dk/deqveps: loading/unloading index 
  };

  hist_                   preHist_;      // history of the previous load step
  hist_                   newHist_;      // history of the current iteration

  Ref<HookeMaterial>      elasticMat_;
  Matrix                  elasticMod_;

  VonMisesStrain          vonMisesEqv_;

}; // end of class IsoDamageMaterial


// -------------------------------------------------------------------
//  giveOmega
// -------------------------------------------------------------------

inline Vector  DamageMaterial::giveOmega()  const
{
  const int s =  preHist_.eqveps.size ( );

  Vector omega ( s );

  for ( int i = 0; i < s; i++ )
  {
    omega[i] =  damageEvolution_ ( newHist_.eqveps[i] ); 
  }

  return omega;
}

// -------------------------------------------------------------------
//  giveOmega(int ip)
// -------------------------------------------------------------------

inline double  DamageMaterial::giveOmega ( int ipoint )  const
{
  return damageEvolution_ ( newHist_.eqveps[ipoint] ); 
}

// -------------------------------------------------------------------
//  giveHistory
// -------------------------------------------------------------------

inline double  DamageMaterial::giveHistory ( int ipoint ) const
{
  return newHist_.eqveps[ipoint];
}

inline Flex<double>  DamageMaterial::giveHistory ( ) const
{
  return newHist_.eqveps;
}

// -------------------------------------------------------------------
//  isLoading
// -------------------------------------------------------------------

inline int    DamageMaterial::isLoading ( int ipoint ) const
{
  return newHist_.loading[ipoint];
}

// -------------------------------------------------------------------
//  giveLengthScale
// -------------------------------------------------------------------

inline double  DamageMaterial::giveLengthScale   ( ) const
{
  return c_;
}

// -------------------------------------------------------------------
//  giveRho
// -------------------------------------------------------------------

inline double  DamageMaterial::giveRho           ( ) const
{
  return elasticMat_->giveRho ();
}

#endif

// *********************************************
//    
//    definition of some helper functions
//    
// ********************************************

// ----------------------------------------------------------------------
//  perfect softening 
// ----------------------------------------------------------------------

double                  perfectSoftening

  ( double Ki,
    double K );

double                  getDerivPerfectSoftening

    ( double Ki,
      double K );

// ----------------------------------------------------------------------
//  linear softening 
// ----------------------------------------------------------------------

double                  linearSoftening

  ( double Ki,
    double Kc,
    double K );

double                  getDerivLinearSoftening

    ( double Ki,
      double Kc,
      double K );

// ----------------------------------------------------------------------
//  exponential softening1 ( Ki, alpha, beta, K) 
// ----------------------------------------------------------------------

double                  exponentialSoftening1

  ( double Ki,
    double alpha,
    double beta,
    const double K );

double                  getDerivExponentialSoftening1

  ( double Ki,
    double alpha,
    double beta,
    double K );

// overloaded for mesh adjusted softenting 

double                  exponentialSoftening1

  ( double Ki,
    double alpha,
    double gf, double ft,
    double K, 
    double he );

double                  getDerivExponentialSoftening1

  ( double Ki,
    double alpha,
    double gf, double ft,
    double K,
    double he );

// ----------------------------------------------------------------------
//  exponential softening2 
// ----------------------------------------------------------------------2

double                  exponentialSoftening2

  ( double Ki,
    double Kc,
    double K );

double                  getDerivExponentialSoftening2

    ( double Ki,
      double Kc,
      double K );

// overloaded for mesh adjusted softening modulus

double                  exponentialSoftening2

  ( double Ki, double Kc,
    double K, double lamda,
    double he );

double                  getDerivExponentialSoftening2

    ( double Ki,double Kc,
      double K, double lamda,
      double he );

// Jirasek regularised exponential softening law

double                  exponentialSoftening3

  ( double Ki, double K,
    double gf, double ft, 
    double he );

double                  getDerivExponentialSoftening3

   ( double Ki, double K,
     double gf, double ft, 
     double he );


// ----------------------------------------------------------------------
//  power softening 
// ----------------------------------------------------------------------

double                  powerSoftening

  ( double Ki,
    double Kc,
    double alpha,
    double beta,
    double K );

double                  getDerivPowerSoftening

  ( double Ki,
    double Kc,
    double alpha,
    double beta,
    double K );

// ----------------------------------------------------------------------
//  hyperbolic softening 
// ----------------------------------------------------------------------

double                  hyperbolicSoftening

  ( double Ki,
    double b,
    double K );

double                  getDerivHyperbolicSoftening

  ( double Ki,
    double b,
    double K );


// ----------------------------------------------------------------------
//  exponential strain energy softening 
// ----------------------------------------------------------------------

double                     expoEnergySoftening

 ( double kappa, double kappaI,
   double kappaC, double s );

double                     getDerivExpoEnergySoftening

 ( double kappa, double kappaI,
   double kappaC, double s );

// ----------------------------------------------------------------------
//  equivalent strain definition of Rankine
// ----------------------------------------------------------------------

double                  equiStrainRankine

  ( const double  young,
    const Vector& princ );




// ----------------------------------------------------------------
//  compute the principal values of a symmetric tensor
// ----------------------------------------------------------------

Vector                  getPrincipalValues

  ( const Vector&       v );

//-----------------------------------------------------------------------
//   getPrincipalValues + return invariants
//-----------------------------------------------------------------------

void                    getPrincipalValues

  ( Vector&             I,
    Vector&             vI,
    const Vector&       v );
  

// ----------------------------------------------------------------
//  compute the derivatives of three principals w.r.t a vector
// ----------------------------------------------------------------

Matrix                  getDerivOfPrincipalValues

  ( const Vector&      v );

// Compute the ramp function <x> = 1/2(x+x)

double                  evalMcAuley

  ( const double );

double                  evalHeaviside

  ( const double );
