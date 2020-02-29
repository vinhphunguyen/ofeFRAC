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
#include <jem/io/PrintWriter.h>
#include <jem/io/FileWriter.h>
#include <jive/util/Printer.h>

//#include <random>

#include "TuronXCohesiveWeibulMat.h"

using namespace jem;

using jem::io::PrintWriter;
using jem::io::FileWriter;
using jem::numeric::norm2;
using jem::numeric::matmul;
using jive::util::Printer;

//=======================================================================
//   class TuronXCohesiveMat
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* TuronXCohesiveWeibulMat::WEIBUL_MOD_PROP = "weibulModulus";
const char* TuronXCohesiveWeibulMat::SIGMA_MIN_PROP  = "sigmaMin";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

TuronXCohesiveWeibulMat::TuronXCohesiveWeibulMat

  ( const int          rank,
    const Properties&  globdat )

  : CohesiveMaterial ( rank, globdat ),
    TuronXCohesiveMat ( rank, globdat )

{}

TuronXCohesiveWeibulMat::~TuronXCohesiveWeibulMat ()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void  TuronXCohesiveWeibulMat::configure

  ( const Properties&     props,
    const Properties&     globdat )

{
  Super::configure ( props, globdat );
  
  props.get ( m_,        WEIBUL_MOD_PROP );
  props.get ( sigmaMin_, SIGMA_MIN_PROP  );
}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------

void  TuronXCohesiveWeibulMat::getConfig

    ( const Properties&     conf,
      const Properties&     globdat ) const

{
  Super::getConfig ( conf, globdat );

  conf.set ( WEIBUL_MOD_PROP, m_ );
  conf.set ( SIGMA_MIN_PROP, sigmaMin_ );
}


//-----------------------------------------------------------------------

void  TuronXCohesiveWeibulMat::update

  ( Vector&               traction,
    const Vector&         jump,
    int                   mpoint )

{
  Vector shiftedJump ( jump + shift_[mpoint] );
  
  double deltaNF = 2. *  gI_ / f2ts_[mpoint];
  double deltaSF = 2. *  gII_/ f6ts_[mpoint];
  double deltaN0 = f2ts_[mpoint] / dummy_;
  double deltaS0 = f6ts_[mpoint] / dummy_;

  deltaN02_ = deltaN0 * deltaN0;
  deltaS02_ = deltaS0 * deltaS0;
  deltaN0F_ = deltaN0 * deltaNF;
  deltaS0F_ = deltaS0 * deltaSF;

  int    loading;
  double damage;
  double beta;

  bool   tension = shiftedJump[0] > 0;

  // compute equivalent displacement jump

  Vector jumpCorr( rank_ ) ;  

  jumpCorr     = shiftedJump;
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

  if      ( abs ( hisDam - damMax ) < 1.e-10 )
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
//   allocPoints
//-----------------------------------------------------------------------

void  TuronXCohesiveWeibulMat::allocPoints

  ( int      count,
    double   dam )

{
  Super::allocPoints ( count, dam );

  // random cohesive strength generation

/*
     
  double xrand;
  double value, f2t, f6t;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0, 1);

  int ipCount   = 2;
  int elemCount = count / ipCount;

  for ( int i = 0; i < elemCount; i++ )
  {
    xrand = dis(gen);
    value = pow(-log(xrand),1/m_) ;
    f2t   = f2t_ * value + sigmaMin_;
    f6t   = f6_  * value + sigmaMin_;
    
    f2ts_.pushBack ( f2t, ipCount );
    f6ts_.pushBack ( f6t, ipCount );
    //System::out() << f2t_ * pow(-log(xrand),1/m_) + sigmaMin_ << "\n";
  }


  Ref<PrintWriter> pdfFile   = newInstance<PrintWriter> (newInstance<FileWriter> ( "pdf1.dat" ));

  double sigmaC, f;

  for (int i = 0; i < f2ts_.size()/2; i++ )
  {
    sigmaC = f2ts_[2*i];
    f      = m_/f2t_ * pow ( (sigmaC-sigmaMin_)/f2t_, m_-1) * exp(- pow ( (sigmaC-sigmaMin_)/f2t_,m_) );
    *pdfFile << sigmaC << " " << f << "\n";
  }
  */
}

//-----------------------------------------------------------------------
//   evalFailure 
//-----------------------------------------------------------------------

double  TuronXCohesiveWeibulMat::evalFailure

    ( const Vector&       sigmaN,
      int                 mpoint )  const

{
  // System::out() << "tip stress in material frame " << sigmaN << endl;

  f2t2_ = f2ts_[mpoint]*f2ts_[mpoint]; 
  f62_  = f6ts_[mpoint]*f6ts_[mpoint]; 

  return  Super::evalFailure ( sigmaN, mpoint );
}

