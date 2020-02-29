/*
 *  Copyright (C) 2009 TU Delft. All rights reserved.
 *  
 *  This class implements a mixed mode traction-seperation law
 *  for cohesive crack models with dummy stiffness (Turon et al. 2006)
 *  with shifted origin to mimick rigid behavior (Hille et al. 2009). 
 *  For use in (Dam)XFEMModel.
 *  
 *  Author: F.P. van der Meer
 *  Date: April 2009
 *
 */

#include <jem/base/Array.h>
#include <jem/base/array/operators.h>
#include <jem/base/array/tensor.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/algebra/matmul.h>

#include "TuronXCohesiveMat.h"

using namespace jem;

using jem::numeric::norm2;
using jem::numeric::matmul;

//=======================================================================
//   class TuronXCohesiveMat
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* TuronXCohesiveMat::LIMIT_SHIFT_PROP = "limitShift";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

TuronXCohesiveMat::TuronXCohesiveMat

  ( const int          rank,
    const Properties&  globdat )

  : CohesiveMaterial ( rank, globdat ), 
    TuronCohesiveMaterial ( rank, globdat ),
    XCohesiveMat ( rank, globdat ), 
    limitShift_ ( true )

{}

TuronXCohesiveMat::~TuronXCohesiveMat ()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void  TuronXCohesiveMat::configure

  ( const Properties&     props,
    const Properties&     globdat )

{
  Super::configure ( props, globdat );
  
  props.find ( limitShift_, LIMIT_SHIFT_PROP );
}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------

void  TuronXCohesiveMat::getConfig

    ( const Properties&     conf,
      const Properties&     globdat ) const

{
  Super::getConfig ( conf, globdat );

  conf.set ( LIMIT_SHIFT_PROP, limitShift_ );
}


//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void  TuronXCohesiveMat::update

  ( Vector&               traction,
    Matrix&               stiff,
    const Vector&         jump,
    int                   mpoint )

{
  Vector shiftedJump ( jump + shift_[mpoint] );

  Super::update ( traction, stiff, shiftedJump, mpoint );
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------
// In DG interface model, jump already subtracted from jump0

void  TuronXCohesiveMat::update

  ( Vector&               traction,
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

  Matrix stiff(rank_,rank_);

  stiff = 0;

  stiff(0,0) = ( 1. - damage * tension ) * dummy_;

  for ( int i = 1; i < rank_; ++i )
  {
    stiff(i,i) = ( 1. - damage ) * dummy_;
  }

  traction = matmul( stiff, jump );

  // compute dissipation per unit surface (incremental, Euler backward)

  double dG  = .5 * ( hisDam - preHist_.damage[mpoint] );
         dG *= deltaF * dummy_ * delta0;
  double  G  = preHist_.dissipation[mpoint] + dG;

  // update history variables

  newHist_.damage     [mpoint] = hisDam;
  newHist_.loading    [mpoint] = loading;
  newHist_.dissipation[mpoint] = G;

  latestHist_ = &newHist_;
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void  TuronXCohesiveMat::update

  ( Vector&               traction,
    Matrix&               stiff,
    const Vector&         jump0,
    const Vector&         djump,
    int                   mpoint )

{
  Vector shiftedJump ( jump0 + djump + shift_[mpoint] );

  //System::out() << "shift jump : " << shift_[mpoint] << "\n";
  //System::out() << " djump : " << djump << "\n";
  //System::out() << "jump0 : " << jump0 << "\n";
  //System::out() << "tot jump : " << shiftedJump << "\n";

  Super::update ( traction, stiff, shiftedJump, mpoint );
}


//-----------------------------------------------------------------------
//   elasticUpdate
//-----------------------------------------------------------------------

void  TuronXCohesiveMat::elasticUpdate

  ( Vector&               traction,
    Matrix&               stiff,
    const Vector&         jump )

{
  // NB: shift cannot be included because mpoint is unknown! 

  Super::elasticUpdate ( traction, stiff, jump );
}

//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------

void  TuronXCohesiveMat::allocPoints

  ( int      count,
    double   dam )

{
  for ( int i = 0; i < count; ++i )
  {
    shift_.pushBack ( Vector ( rank_ ) );
    
    shift_[i] = 0.;
  }

  Super::allocPoints_ ( count, dam, 1 );
}

// --------------------------------------------------------------------
//  scaleDissipationIncrement
// --------------------------------------------------------------------

void  TuronXCohesiveMat::scaleDissipationIncrement

  ( const double     factor,
    const int        ipoint )

{
  double dG  = newHist_.dissipation[ipoint] - preHist_.dissipation[ipoint];
         dG *= factor;

  if ( factor > 1. ) 
  {
    System::warn() << "scaleDissipationIncrement is scaling up" << endl;
  }

  newHist_.dissipation[ipoint] = preHist_.dissipation[ipoint] + dG ;
}

//-----------------------------------------------------------------------
//   evalFailure 
//-----------------------------------------------------------------------

double  TuronXCohesiveMat::evalFailure

    ( const Vector&       sigmaN,
      int                 mpoint )  const

{
  // System::out() << "tip stress in material frame " << sigmaN << endl;

  double tauN  = ( sigmaN[0] + abs(sigmaN[0]) ) / 2.;
  double tauS2 = sigmaN[1] * sigmaN[1];

  if ( rank_ == 3 )
  {
    tauS2 += sigmaN[2] * sigmaN[2];
  }

  double tauS = sqrt ( tauS2 );
  
  //double criterion;
  //if ( sigmaN[0] < 0. ) criterion = 0.;
  //if ( sigmaN[0] > 0. ) criterion = tauN / f2t_;

  //System::out() << "B-K failure criterion: " << criterion << endl;


  if ( tauS + tauN < EPS ) return 0.;

  double beta  = tauS / ( tauS + tauN );
  double B     = beta * beta;
         B    /= 1. + 2.*B - 2.*beta;
  double Beta  = ::pow( B , eta_ );

  double Fb2   = f2t2_ + ( f62_ - f2t2_ ) * Beta;
  double tauN2 = tauN * tauN;
  double tauB2 = tauN2 + tauS2;

  double criterion = tauB2 / Fb2;

  return  criterion;
}

//-----------------------------------------------------------------------
//   deallocPoints
//-----------------------------------------------------------------------

void  TuronXCohesiveMat::deallocPoints

  ( int      count )

{
  shift_.popBack ( count );

  Super::deallocPoints ( count );
}

//-----------------------------------------------------------------------
//   initShift
//-----------------------------------------------------------------------

void  TuronXCohesiveMat::initShift 

  ( const int             ipoint,
    const Vector&         traction )

{
  TensorIndex i;
  Vector shift ( shift_[ipoint] ); // shallow copy

  shift = traction / getDummy();

  bool tension   = shift[0] > 0;
       shift[0] *= tension;


  if ( limitShift_ )
  {
    double delta  = norm2 ( shift );

    if ( delta > 0. )
    {
      double deltaS = norm2 ( shift[slice(1,END)] );
      double beta   = deltaS / delta;

      // make sure the shift doesn't translate the origin beyond the peak

      double delta0 = getDelta0_ ( beta );

      if ( delta > delta0 )
      {
        shift    *= delta0 / delta;

        shift[i] += where ( shift[i] >= 0., EPS, -EPS );
      }
    }
  }
  
  //System::out () << " shift: " << shift << "\n";
}


