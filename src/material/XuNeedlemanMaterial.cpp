/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the Xu-Needleman traction-seperation law
 *  for two dimensions.
 *
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 21 December 2007
 *
 */

#include <jem/base/System.h>
#include <jem/base/Float.h>
#include <jem/base/Error.h>
#include <jem/base/array/utilities.h>
#include <jem/base/array/operators.h>
#include <jem/base/array/tensor.h>
#include <jem/base/array/intrinsics.h>
#include <jem/numeric/utilities.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>

extern "C"
{
  #include  <math.h>
}

#include "HookeMaterial.h"
#include "XuNeedlemanMaterial.h"

using namespace jem;
using namespace jem::io;

using jem::numeric::MatmulChain;
using jem::util::Properties;
using jive::Vector;

const double       e = 2.718281828459;

//=======================================================================
//   class XuNeedlemanMaterial
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  XuNeedlemanMaterial::BETA_PROP       = "beta";
const char*  XuNeedlemanMaterial::TENSILE_PROP    = "sigmaMax";
const char*  XuNeedlemanMaterial::SHEAR_PROP      = "tauMax";
const char*  XuNeedlemanMaterial::OPEN_DISP_PROP  = "critOpen";
const char*  XuNeedlemanMaterial::SHEAR_DISP_PROP = "critSlid";
const char*  XuNeedlemanMaterial::R_PROP          = "r";
const char*  XuNeedlemanMaterial::Q_PROP          = "q";
const char*  XuNeedlemanMaterial::ENERGY_PROP     = "energy";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------

XuNeedlemanMaterial::XuNeedlemanMaterial

    ( const int             rank,
      const Properties&     globdat )

 : CohesiveMaterial ( rank, globdat )
{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );
  rank_    = rank;

  critOpen_       = 0.;
  critSlid_       = 0.;
  beta_           = 0.0;
  r_              = 0.0;
  phin_           = 0.;
  tensileStress_  = 0.;
  shearStress_    = 0.;
}


XuNeedlemanMaterial::~XuNeedlemanMaterial ()
{} 

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void XuNeedlemanMaterial::configure 

 ( const Properties& props,
   const Properties& globdat  )
{
  using jem::maxOf;

  bool energy  = props.find  ( phin_,           ENERGY_PROP,    0.0, maxOf( phin_           ) );
  bool tensile = props.find  ( tensileStress_,  TENSILE_PROP,   0.0, maxOf( tensileStress_ ) );
  bool shear   = props.find  ( shearStress_,    SHEAR_PROP,     0.0, maxOf( shearStress_   ) );
  bool open    = props.find  ( critOpen_,       OPEN_DISP_PROP, 0.0, maxOf( critOpen_      ) );
  bool slid    = props.find  ( critSlid_,       SHEAR_DISP_PROP,0.0, maxOf( critSlid_      ) );

  props.get ( beta_,           BETA_PROP,      0.0, maxOf( beta_           ) );
  props.get ( r_,              R_PROP,         0.0, maxOf( r_              ) );
  props.get ( q_,              Q_PROP,         0.0, maxOf( q_             ) );

  if      ( open && slid && energy )
  {
  }
  else if ( energy && shear && tensile )
  {
    critOpen_ = phin_      / ( tensileStress_ * e                );
    critSlid_ = q_ * phin_ / ( shearStress_   * sqrt ( 0.5 * e ) );
  }
  else if ( tensile && open  && shear )
  {
    phin_     = e * tensileStress_ * critOpen_;
    critSlid_ = q_ * phin_ / ( shearStress_   * sqrt ( 0.5 * e ) );
  }
  else
  {
    System::out() << "XuNeedleman cohesive law: incorrect parameter set\n";
    JEM_ASSERT(1==0);
  }

  //  phit_ = sqrt(0.5 * e) * shearStress_ * critSlid_;
  //  q_    = phit_ / phin_;

  critOpenInv_ = 1. / critOpen_;
  critSlidInv_ = 1. / critSlid_;
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void XuNeedlemanMaterial::getConfig 

  ( const Properties& conf,
    const Properties& globdat ) const
{
  conf.set ( BETA_PROP,       beta_           );
  conf.set ( R_PROP,          r_              );
  conf.set ( Q_PROP,          q_              );
  conf.set ( SHEAR_PROP,      shearStress_    );
  conf.set ( TENSILE_PROP,    tensileStress_  );
  conf.set ( OPEN_DISP_PROP,  critOpen_       );
  conf.set ( SHEAR_DISP_PROP, critSlid_       );  
  conf.set ( ENERGY_PROP,     phin_           );  
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void  XuNeedlemanMaterial::update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump,
      int                   mpoint )
{
  double    kappa;
  int       load;

  // the normal and tangential crack opening

  double    crOpen     = jump[0];
  double    crSlid     = jump[1];
 
  // the effective crack opening 

  double    eqvOpen    = ::sqrt ( crOpen * crOpen + 
                                  beta_ * beta_ * crSlid * crSlid );

  // the loading function 

  double    eqvOpenMax = preHist_.eqvOpenMax[mpoint];
  double    f          = eqvOpen - eqvOpenMax;

  // f == 0 (exactly equals to zero), at the beginning of every 
  // loading step, this condition holds for loading integration point

  if ( ::abs (f) < 1e-14 )
  {
    kappa = eqvOpenMax;
    load  = preHist_.loading[mpoint];
  }
  else 
  {
    if ( f > 0. )
    {
      kappa = eqvOpen;
      load  = 1;
    }
    else 
    {
      kappa = eqvOpenMax;
      load  = 0;
    }
  }

  //if ( crOpen < 0. )
  //{
  //  System::out() << crOpen << "interpenetration!!!\n";
  //}

  // compute the traction vector and the tangent stiffness matrix

  getTraction_      ( traction, jump, eqvOpen, eqvOpenMax, load );
  getTangentMatrix_ ( stiff,    jump, eqvOpen, eqvOpenMax, load );

  // update history variables

  newHist_.eqvOpenMax[mpoint]  = kappa;
  newHist_.loading[mpoint]     = load;

  // debug

//  System::out () << "jump       : " << jump       << "\n" 
//                 << "eqvOpen    : " << eqvOpen    << "\n"
//                << "eqvOpenMax : " << eqvOpenMax << "\n"
//                << "loading    : " << load       << "\n\n";
// if ( load == 0 )
// {
//    System::out () << "jump       : " << jump       << "\n" 
//                 << "eqvOpen    : " << eqvOpen    << "\n"
//                << "eqvOpenMax : " << eqvOpenMax << "\n";
// }

  if ( jem::Float::isNaN( sum(    stiff ) ) ||
       jem::Float::isNaN( sum( traction ) ) )
  { 
    System::out() << "Invalid matrix element" <<
      "\n traction: " << traction <<
      "\n stiff   :\n" << stiff << endl;
  }
}

//-----------------------------------------------------------------------
//   update (moonen material)
//-----------------------------------------------------------------------

void  XuNeedlemanMaterial::update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         tEff,
      const Matrix&         jumpStiff,
      int                   mpoint )
{
  throw Error ( JEM_FUNC,
      "Moonen material not defined with Rankine law" );
}


//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void  XuNeedlemanMaterial::elasticUpdate

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump
    )
{
}

// --------------------------------------------------------------------
//  commit
// --------------------------------------------------------------------

void  XuNeedlemanMaterial::commit()
{
  newHist_.eqvOpenMax. swap ( preHist_.eqvOpenMax );
  newHist_.loading.    swap ( preHist_.loading    );
}

// --------------------------------------------------------------------
//  allocPoints
// --------------------------------------------------------------------

void  XuNeedlemanMaterial::allocPoints 

  ( int count, double dam )
{
  // not yet allocated

  if ( preHist_.eqvOpenMax.size ( ) == 0 ) 
  {
    preHist_.eqvOpenMax.resize ( count );  
    preHist_.loading.   resize ( count );

    newHist_.eqvOpenMax.resize ( count );
    newHist_.loading.   resize ( count );

    preHist_.eqvOpenMax = 0.0;
    preHist_.loading    = 1  ;
    newHist_.eqvOpenMax = 0.0;
    newHist_.loading    = 0  ;
  }

  // or appending at the end of the vector

  else
  {
    for ( int i = 0; i < count; i++ )
    {
      preHist_.eqvOpenMax.pushBack ( 0.0 );
      preHist_.loading.   pushBack ( 1   );

      newHist_.eqvOpenMax.pushBack ( 0.0 );  
      newHist_.loading.   pushBack ( 0   );  
    }
  }
}

//-----------------------------------------------------------------------
//   getTraction
//-----------------------------------------------------------------------

void XuNeedlemanMaterial::getTraction_

 ( const Vector& traction,
   const Vector& jump,
   double eqvOpen, double eqvOpenMax,
   int    loading ) const

{
  double crOpen = jump[0];
  double crSlid = jump[1];

  traction = 0.;

  if ( loading )
  {
    double normalFrac = crOpen * critOpenInv_;
    double tangenFrac = crSlid * critSlidInv_;

    double normalExp  = ::exp ( - normalFrac );
    double tangenExp  = ::exp ( - tangenFrac * tangenFrac );

    double rminus1Inv = 1. / ( r_ - 1. );

    traction[0] = -phin_ * critOpenInv_ * normalExp * ( 
	           (-r_ + normalFrac) * (1. - q_) * rminus1Inv - 
		   (q_ + (r_ - q_) * rminus1Inv * (normalFrac - 1.) ) * tangenExp
	          );

    traction[1] = 2. * phin_ * tangenFrac * critSlidInv_ * ( q_ + ( r_ - q_ ) * rminus1Inv  
	          * normalFrac ) * normalExp * tangenExp;
  }
  else
  {
    double frac = eqvOpen / eqvOpenMax;

    double sCrOpen = crOpen / frac;
    double sCrSlid = crSlid / frac;

    double normalFrac = sCrOpen * critOpenInv_;
    double tangenFrac = sCrSlid * critSlidInv_;

    double normalExp  = exp ( - normalFrac );
    double tangenExp  = exp ( - tangenFrac * tangenFrac );

    double rminus1Inv = 1. / ( r_ - 1. );

    traction[0] = -phin_ * critOpenInv_ * normalExp * ( 
	           (-r_ + normalFrac) * (1. - q_) * rminus1Inv - 
		   (-q_ + (r_-q_) * rminus1Inv * (normalFrac - 1.) ) * tangenExp
	          ) * frac ;

    traction[1] = 2. * phin_ * tangenFrac * critSlidInv_ * ( q_ + ( r_ - q_ ) * rminus1Inv  
	          * normalFrac ) * normalExp * tangenExp * frac;
  }
}

//-----------------------------------------------------------------------
//   getTangentMatrix
//-----------------------------------------------------------------------

void XuNeedlemanMaterial::getTangentMatrix_

 ( const Matrix& tangent,
   const Vector& jump,
   double eqvOpen, double eqvOpenMax,
   int    loading ) const

{ 
  double crOpen = jump[0];
  double crSlid = jump[1];

  double T11, T12, T21, T22;

  tangent = 0.;

  if ( eqvOpen < 1e-8 )
  {
    loading = 1;
  }

  if ( loading )
  {
    double normalFrac = crOpen * critOpenInv_;
    double tangenFrac = crSlid * critSlidInv_;

    double normalExp  = exp ( - normalFrac );
    double tangenExp  = exp ( - tangenFrac * tangenFrac );

    double exp1exp2   = normalExp * tangenExp;

    double rminus1Inv = 1. / ( r_ - 1. );

    T11 = phin_ * critOpenInv_ * critOpenInv_ * normalExp * ( 
	           (1. - q_) * rminus1Inv * (-r_ + normalFrac - 1. ) -
		   (q_ + (r_-q_) * rminus1Inv * (normalFrac -2.) ) * tangenExp
		 );

    T12 = 2. * phin_ * critOpenInv_ * (
	          -q_ + (r_ - q_) * rminus1Inv * ( 1. - normalFrac )
	          ) * tangenFrac * critSlidInv_ * exp1exp2;

    T21 = 2. * phin_ * critSlidInv_ * critOpenInv_ * tangenFrac * exp1exp2 
		     * ( (r_ - q_) * rminus1Inv * ( 1. - normalFrac ) - q_ );

    T22 = 2. * phin_ * critSlidInv_* critSlidInv_ * 
                 (q_ + (r_ - q_) * rminus1Inv * normalFrac ) * 
		 (1. - 2. * tangenFrac * tangenFrac) * exp1exp2;
  }
  else
  {
    double frac = eqvOpen / eqvOpenMax;

    double sCrOpen = crOpen / frac;
    double sCrSlid = crSlid / frac;

    double normalFrac = sCrOpen * critOpenInv_;
    double tangenFrac = sCrSlid * critSlidInv_;

    double normalExp  = exp ( - normalFrac );
    double tangenExp  = exp ( - tangenFrac * tangenFrac );

    double rminus1Inv = 1. / ( r_ - 1. );

    double tn = -phin_ / critOpen_ * normalExp * ( 
	           (-r_ + normalFrac) * (1. - q_) * rminus1Inv - 
		   (q_ + (r_-q_) * rminus1Inv * (normalFrac - 1.) ) * tangenExp
	          );

    double tt = 2. * phin_ * tangenFrac * critSlidInv_ * ( q_ + ( r_ - q_ ) * rminus1Inv  
	           * normalFrac ) * normalExp * tangenExp;

    double den         = 1. / ( eqvOpen * eqvOpenMax );
    double mul         = crOpen * crSlid;
    double eqvOpen2Inv = 1. / ( eqvOpen * eqvOpen ); 
    double beta2       = beta_ * beta_ ;

    double t11 = phin_ * critOpenInv_ * normalExp * ( 
	           (1. - q_) * rminus1Inv * (-r_ + normalFrac - 1. ) -
		   (q_ + (r_-q_) * rminus1Inv * (normalFrac -2.) ) * tangenExp
		 ) * frac ;

    double t12 = 2. * phin_ * critOpenInv_ * (
	          -q_ + (r_ - q_) * rminus1Inv * ( 1. - normalFrac )
	          ) * crSlid * normalExp * tangenExp / critSlidInv_ * frac;

    double t21 = 2. * phin_ * critSlidInv_ * critOpenInv_ * tangenFrac * 
                    normalExp * tangenExp * 
		    ( (r_-q_) * rminus1Inv * ( 1. - normalFrac ) - q_ ) * frac;

    double t22 = 2. * phin_ * critSlidInv_* critSlidInv_ * 
                 (q_ + (r_ - q_) * rminus1Inv * normalFrac ) * 
		 (1. - 2. * tangenFrac * tangenFrac) * 
		 normalExp * tangenExp * frac;

    //System::out() << tn << endl 
    //              << tt << endl
    //		    << t11 << " " << t12 << " " << t21 << " " << t22 << endl;
                   
    T11 = crOpen * den * tn + 
          beta2  * eqvOpen2Inv * crSlid * crSlid * t11 - 
	  mul    * eqvOpen2Inv * t12;

    T12 = beta_  * beta_  * crSlid * den * tn + 
          crOpen * crOpen * eqvOpen2Inv * t12 -
	  beta2  * mul    * eqvOpen2Inv * t11;

    T21 = crOpen * den * tt + 
          beta2  * crSlid * crSlid * eqvOpen2Inv * t21 -
	  mul    * eqvOpen2Inv * t22;

    T22 = beta_  * beta_  * crSlid * den * tt +
          crOpen * crOpen * eqvOpen2Inv * t22 -
	  beta2  * mul    * eqvOpen2Inv * t21;
  }

  tangent(0,0) = T11; tangent(0,1) = T12;
  tangent(1,0) = T21; tangent(1,1) = T22;
}

//-----------------------------------------------------------------------
//   justTractionFree_ 
//-----------------------------------------------------------------------

bool  XuNeedlemanMaterial::justTractionFree

  ( const int   i ) const

{
  return false;
}

//-----------------------------------------------------------------------
//   evalFailure 
//-----------------------------------------------------------------------

//double   XuNeedlemanMaterial::evalFailure
//
//    ( const Vector&       stress,
//      const Ref<Material> bulkMat )  const
//
//{
//  Vector     priStress ( rank_ );
//
//  return  priStress[0] / tensileStress_;
//}
//
//bool XuNeedlemanMaterial::isFailure
//
//    ( const Vector&       stress,
//      const Ref<Material> bulkMat )  const
//{
//  return  evalFailure( stress, bulkMat ) > 1.;
//}


//bool XuNeedlemanMaterial::evalFailure 
//
//     ( const Vector& traction ) const
//{
//  double tn = traction[0];
//  double tt = traction[1];
//
//         tn = tn > 0 ? tn : 0.;
//
//  double bet2 = beta_ * beta_;
//
//  double teff = sqrt ( tt * tt / bet2 + tn * tn );
//
//  return teff > tensileStress_;
//}
