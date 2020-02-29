/*
 *  Copyright (C) 2008 TU Delft. All rights reserved.
 *
 *  This class implements a mixed mode traction-seperation law
 *  for cohesive crack models with dummy stiffness (Turon et al. 2006). 
 *  For use in InterfaceModel (2D or 3D).
 *
 *  Author: F.P. van der Meer
 *  Date: June 2008
 *
 */

#include <jem/base/System.h>
#include <jem/base/array/utilities.h>
#include <jem/base/array/select.h>
#include <jem/base/array/operators.h>
#include <jem/base/array/tensor.h>
#include <jem/base/array/intrinsics.h>
#include <jem/io/PrintWriter.h>
#include <jem/util/Properties.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>

#include "OrthotropicMaterial.h"
#include "TuronCohesiveMaterial.h"

extern "C"
{
  #include <math.h>
}

using namespace jem;
using namespace jem::io;

using jem::numeric::matmul;
using jem::numeric::norm2;
using jem::numeric::MatmulChain;
using jem::util::Properties;
using jive::Vector;

//=======================================================================
//   class TuronCohesiveMaterial
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  TuronCohesiveMaterial::F2T_PROP      = "f2t";
const char*  TuronCohesiveMaterial::F6_PROP       = "f6";
const char*  TuronCohesiveMaterial::G_I_PROP      = "gI";
const char*  TuronCohesiveMaterial::G_II_PROP     = "gII";
const char*  TuronCohesiveMaterial::G_III_PROP    = "gIII";
const char*  TuronCohesiveMaterial::ETA_PROP      = "eta";
const char*  TuronCohesiveMaterial::DUMMY_PROP    = "dummy";

const double TuronCohesiveMaterial::OMEGA_MAX     = 1. - 1.e-12;
const double TuronCohesiveMaterial::EPS           = 1.e-10;

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------

TuronCohesiveMaterial::TuronCohesiveMaterial

  ( const int          rank,
    const Properties&  globdat )

  : CohesiveMaterial ( rank, globdat )

{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );

  eta_  = 0.;

  // some default values for optional input

  dummy_       = 1.e6;
}


TuronCohesiveMaterial::~TuronCohesiveMaterial ()
{} 

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void TuronCohesiveMaterial::configure

  ( const Properties& props,
    const Properties& globdat )

{
  using jem::maxOf;

  props.get  ( f2t_  , F2T_PROP,  0., maxOf( f2t_  ) );
  props.get  ( f6_   , F6_PROP,   0., maxOf( f6_   ) );
  props.get  ( gI_   , G_I_PROP,  0., maxOf( gI_   ) );
  props.get  ( gII_  , G_II_PROP, 0., maxOf( gII_  ) );
  props.find ( dummy_, DUMMY_PROP,0., maxOf( dummy_) );

  props.find ( eta_  ,   ETA_PROP,  .5, maxOf( eta_  ) );

  // normally g3c = g2c

  gIII_ = gII_;

  props.find  ( gIII_  , G_III_PROP, 0., maxOf( gIII_  ) );

  // onset and final displacement in normal and tangential
  // directions, 0 for onset, F for final.
  // N for normal and S for shear

  double deltaNF = 2. *  gI_ / f2t_;
  double deltaSF = 2. * gII_ / f6_;
  double deltaN0 = f2t_ / dummy_;
  double deltaS0 = f6_  / dummy_;

  deltaN02_ = deltaN0 * deltaN0;
  deltaS02_ = deltaS0 * deltaS0;
  deltaN0F_ = deltaN0 * deltaNF;
  deltaS0F_ = deltaS0 * deltaSF;

  f2t2_ = f2t_ * f2t_;
  f62_  = f6_  * f6_;
}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void TuronCohesiveMaterial::getConfig

  ( const Properties& conf , 
    const Properties& globdat ) const

{
  conf.set (   F2T_PROP,  f2t_   );
  conf.set (    F6_PROP,  f6_    );
  conf.set (   G_I_PROP,  gI_    );
  conf.set (  G_II_PROP,  gII_   );
  conf.set (  G_III_PROP, gIII_  );
  conf.set (   ETA_PROP,  eta_   );
  conf.set ( DUMMY_PROP,  dummy_ );
  conf.set ( "f2t2",      f2t2_ );
}

//-----------------------------------------------------------------------
//   update (regular)
//-----------------------------------------------------------------------

void  TuronCohesiveMaterial::update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump,
      int                   mpoint )
{
  int    loading;
  double damage;
  double beta;

  bool   tension = jump[0] > 0;

  // compute equivalent displacement jump

  Vector jumpCorr( rank_ ) ;  

  jumpCorr     = jump;
  jumpCorr[0] *= tension;      // jumpCorr does not count for negative normal jump


  double delta  = norm2 ( jumpCorr );
  double deltaS = norm2 ( jump[slice(1,END)] );

  // update the damage

  if ( deltaS > 1e-16 ) // original code: if (deltaS > 0.)
  {
    beta = deltaS / ( deltaS + jumpCorr[0] );
  }
  else
  {
    beta = 0.;
  }

  double B     = beta * beta;
         B    /= 1. + 2.*B - 2.*beta;
  double Beta  = ::pow( B , eta_ );

  double delta0  = deltaN02_ + ( deltaS02_ - deltaN02_ ) * Beta;
         delta0  = sqrt( delta0 );
  double deltaF  = deltaN0F_ + ( deltaS0F_ - deltaN0F_ ) * Beta;
         deltaF /= delta0;

  double hisDam = ( delta - delta0 ) / ( deltaF - delta0 );
  double damMax = preHist_.damage[mpoint];

  if      ( jem::numeric::abs ( hisDam - damMax ) < 1.e-10 )
  {
    loading = preHist_.loading[mpoint];
  }
  else if ( hisDam > damMax )
  {
    loading = 1;
  }
  else
  {
    loading = 0;
    hisDam  = damMax;
  }

  if ( hisDam > 0. )
  {
    if ( loading )
    {
      damage = hisDam * deltaF / delta;
    }
    else
    {
      damage  = hisDam * deltaF;
      damage /= delta0 + hisDam * ( deltaF - delta0 );
    }
  }
  else
  {
    hisDam = damage =  0.;
  }

  // do not allow damage=1.0, use a threshold for this
  // OMEGA_MAX = 0.99 for example.

  if ( hisDam >= OMEGA_MAX )
  {
    hisDam = damage = OMEGA_MAX;
    loading = 0;
  }

  // compute secant stiffness and traction

  stiff = 0;

  stiff(0,0) = ( 1. - damage * tension ) * dummy_;

  for ( int i = 1; i < rank_; ++i )
  {
    stiff(i,i) = ( 1. - damage ) * dummy_;
  }

  traction = matmul( stiff, jump );

  // compute the tangent stiffness matrix

  if ( loading )
  {
    double H  = deltaF * delta0 ;
           H /= ( deltaF - delta0 ) * delta * delta * delta;

    TensorIndex i,j;

    stiff(i,j) -= dummy_ * H * jumpCorr[i] * jumpCorr[j];

    // apocryphal term

    Vector dtdd   ( - dummy_ * jumpCorr );

    double   tmp  = deltaF - delta0;
             tmp *= tmp * delta;
    double dddd0  = - deltaF * ( deltaF - delta ) / tmp;
    double ddddF  = - delta0 * ( delta - delta0 ) / tmp;

    double dBedB  = eta_ * ::pow( B, eta_-1. );

    double dd0dB  = deltaS02_ - deltaN02_;
           dd0dB /= 2. * delta0; 
    double ddFdB  = deltaS0F_ - deltaN0F_ - deltaF * dd0dB;
           ddFdB /= delta0;

    Vector  dBdd  ( rank_ );

    if ( deltaS > 1.e-10 )
    {
             tmp  = jumpCorr[0] / deltaS / deltaS;
         dBdd[0]  = - 2. * tmp * B * B;

      if (   rank_ == 2 )
      {
         dBdd[1]  = - dBdd[0] * jumpCorr[0] / jumpCorr[1];
      }
      else
      {
         dBdd[1]  = 2 * tmp * tmp * B * B * jump[1];
         dBdd[2]  = 2 * tmp * tmp * B * B * jump[2];
      }
    }
    else
    {
      dBdd[0] = 0.;
      dBdd[slice(1,END)] = 1. / jumpCorr[0] / jumpCorr[0];
    }

    double  dddB  = dddd0 * dd0dB + ddddF * ddFdB;
            dddB *= dBedB;
    Vector  dddd  ( dddB * dBdd );

    stiff(i,j)   += dtdd[i] * dddd[j];
  }

  // compute dissipation per unit surface (incremental, Euler backward)

  double dG  = .5 * ( hisDam - preHist_.damage[mpoint] );
         dG *= deltaF * dummy_ * delta0;
  double  G  = preHist_.dissipation[mpoint] + dG;

  // update history variables

  newHist_.damage     [mpoint] = hisDam;
  newHist_.loading    [mpoint] = loading;
  newHist_.dissipation[mpoint] = G;

  latestHist_ = &newHist_;
  
 /* 
  System::out() << "jum: " << jump << "\n";
  System::out() << "damMax: " << damMax << "\n";
  System::out() << "HisDam: " << hisDam << "\n";
  System::out() << "dam: " << damage << "\n";
  System::out() << "tra: " << traction << "\n";
  System::out() << "sti: " << stiff << "\n";
  */
 
 
} 


void TuronCohesiveMaterial::elasticUpdate

( Vector& traction, Matrix& stiff, const Vector& jump )
{
  TensorIndex i,j;

  traction = dummy_ * jump;

  stiff(i,j) = where ( i == j, dummy_, 0. );
}

// end: update (regular)

//-----------------------------------------------------------------------
//   update (moonen material)
//-----------------------------------------------------------------------

void  TuronCohesiveMaterial::update

    ( Vector&               traction,
      Matrix&               stiff,
      const Matrix&         jumpStiff,
      const Vector&         tEff,
      int                   mpoint )

{ }

//-----------------------------------------------------------------------
//   update (regular)
//-----------------------------------------------------------------------

void  TuronCohesiveMaterial::update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump0,
      const Vector&         djump,
      int                   mpoint )
{
    Vector jump ( rank_ ); jump = jump0 + djump;
    update ( traction, stiff, jump, mpoint );
}

// --------------------------------------------------------------------
//  commit
// --------------------------------------------------------------------

void  TuronCohesiveMaterial::commit()

{
  newHist_.damage     . swap ( preHist_.damage      );
  newHist_.loading    . swap ( preHist_.loading     );
  newHist_.dissipation. swap ( preHist_.dissipation );

  latestHist_ = &preHist_;
}


// --------------------------------------------------------------------
//  allocPoints
// --------------------------------------------------------------------

void  TuronCohesiveMaterial::allocPoints 

    ( int     count, 
      double  dam   )

{
  allocPoints_ ( count, dam, 0 );
}

// --------------------------------------------------------------------
//  allocPoints_
// --------------------------------------------------------------------

void  TuronCohesiveMaterial::allocPoints_ 

    ( const int     count, 
      const double  dam,
      const int     loading)

{
  // not yet allocated

  if ( preHist_.damage.size ( ) == 0 ) 
  {
    preHist_.damage      .resize ( count );  
    preHist_.loading     .resize ( count );
    preHist_.dissipation .resize ( count );

    newHist_.damage      .resize ( count );
    newHist_.loading     .resize ( count );
    newHist_.dissipation .resize ( count );

    preHist_.damage      = dam;
    preHist_.loading     = loading  ;
    preHist_.dissipation = 0. ;

    newHist_.damage      = dam;
    newHist_.loading     = loading  ;
    newHist_.dissipation = 0. ;
  }

  // appending at the end of the vector

  else
  {
    for ( int i = 0; i < count; i++ )
    {
      preHist_.damage     .pushBack ( dam );  
      preHist_.loading    .pushBack (  loading  );  
      preHist_.dissipation.pushBack (  0. );

      newHist_.damage     .pushBack ( dam );  
      newHist_.loading    .pushBack (  loading  );  
      newHist_.dissipation.pushBack (  0. );
    }
  }
}

// --------------------------------------------------------------------
//  getDelta0
// --------------------------------------------------------------------

double  TuronCohesiveMaterial::getDelta0_ ( const double beta ) const

{
    double B       = beta * beta;
    B             /= 1. + 2.*B - 2.*beta;
    double Beta    = ::pow( B , eta_ );
    double delta02 = deltaN02_ + ( deltaS02_ - deltaN02_ ) * Beta;

    return sqrt( delta02 );
}

// --------------------------------------------------------------------
//  deallocPoints
// --------------------------------------------------------------------

void TuronCohesiveMaterial::deallocPoints ( int count )
{
  preHist_.damage     .popBack ( count );
  preHist_.loading    .popBack ( count );
  preHist_.dissipation.popBack ( count );

  newHist_.damage     .popBack ( count );
  newHist_.loading    .popBack ( count );
  preHist_.dissipation.popBack ( count );
}

//-----------------------------------------------------------------------
//   justTractionFree_ 
//-----------------------------------------------------------------------

bool  TuronCohesiveMaterial::justTractionFree

  ( const int   i ) const

{
  return (  latestHist_->damage[i] >= OMEGA_MAX );
}
