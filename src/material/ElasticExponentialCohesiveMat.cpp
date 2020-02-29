/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens an exponential traction-seperation law
 *  for mode I problems. The model consists of a linear elastic
 *  branch (with dummy stiffness) and an exponential softening
 *  branch. This cohesive law was presented in
 *  "Modelling of cohesive crack growth in concrete structures with the extended
 *  finite element method" of Unger et al, CMAME, 2007.
 *
 *  History variable: maximum equivalent jump = sqrt{<un>^2 + beta^2*us^2}
 *  
 *  Author: V.P. Nguyen, phu.nguyen@adelaide.edu.au
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
#include <jem/numeric/algebra/utilities.h>

#include "HookeMaterial.h"
#include "ElasticExponentialCohesiveMat.h"

using namespace jem;
using namespace jem::io;

using jem::util::Properties;
using jem::numeric::norm2;
using jive::Vector;

//=======================================================================
//   class ElasticExponentialCohesiveMat
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  ElasticExponentialCohesiveMat::BETA_PROP     = "beta";
const char*  ElasticExponentialCohesiveMat::TENSILE_PROP  = "tensile";
const char*  ElasticExponentialCohesiveMat::ENERGY_PROP   = "energy";
const char*  ElasticExponentialCohesiveMat::DUMMY_PROP    = "dummy";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------

ElasticExponentialCohesiveMat::ElasticExponentialCohesiveMat

    ( const int             rank,
      const Properties&     globdat )

 : CohesiveMaterial ( rank, globdat ),
   XCohesiveMat     ( rank, globdat )
{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );
  rank_    = rank;

  beta_  = 0.0;
  ft_    = 0.0;
  g1c_   = 0.0;
  dummy_ = 1e5;
}


ElasticExponentialCohesiveMat::~ElasticExponentialCohesiveMat ()
{} 

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void ElasticExponentialCohesiveMat::configure 

 ( const Properties& props,
   const Properties& globdat  )
{
  using jem::maxOf;

  props.get   ( ft_,    TENSILE_PROP, 0.0, maxOf( ft_   ) );
  props.get   ( g1c_,   ENERGY_PROP,  0.0, maxOf( g1c_  ) );
  props.find  ( dummy_, DUMMY_PROP,   0.0, maxOf( dummy_) );
  props.find  ( beta_,  BETA_PROP,    0.0, maxOf( beta_ ) );

  delta0_     = ft_ / dummy_;
  beta2_      = beta_ * beta_;
  mftOverGf_  = -ft_/g1c_;
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void ElasticExponentialCohesiveMat::getConfig 

  ( const Properties& conf,
    const Properties& globdat ) const
{
  conf.set ( BETA_PROP,    beta_  );
  conf.set ( TENSILE_PROP, ft_    );
  conf.set ( ENERGY_PROP,  g1c_   );
  conf.set ( DUMMY_PROP,   dummy_ );  
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void  ElasticExponentialCohesiveMat::update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump,
      int                   mpoint )
{
  Vector    corJump ( jump.size() );
  double    crOpen  = jump[0];                       // the normal crack opening
  double    crSlid  = jump[1];
  corJump[0]        = crOpen;
  corJump[1]        = beta_*crSlid;
  bool      tension = crOpen > 0;
  corJump[0]       *= tension;                       // negative normal jump does not contribute
  double  eqvOpen   = norm2 ( corJump );             // the effective (equivalent) crack opening 

  // Get the history variable of the previous converged load step
  // and compute the loading function 

  double    eqvOpenMax = preHist_.eqvOpenMax[mpoint];
  double    f          = eqvOpen - eqvOpenMax;

  double    kappa = eqvOpenMax;
  int       load  = 0;

  // f == 0 (exactly equals to zero), at the beginning of every 
  // loading step, this condition holds for loading integration point

  if (  jem::numeric::abs (f) < 1e-16 )
  {
    load  = preHist_.loading[mpoint];
  }
  else 
  {
    if ( f > 0.0 )
    {
      kappa = eqvOpen;
      load  = 1;
    }
    else 
    {
      //System::out() << "Isn't it strange to unload in this case!!!\n";
    }
  }
  
  // update history variables

  newHist_.eqvOpenMax[mpoint]  = kappa;
  newHist_.loading[mpoint]     = load;

  // compute the traction vector and the tangent stiffness matrix

  traction = 0.;
  stiff    = 0.;

  double sigma, dsigma;

  if ( load )
  {
    if ( eqvOpen < delta0_ ) 
    { 
      sigma   = dummy_ * eqvOpen;
      dsigma  = dummy_;
      
      sigma /= eqvOpen;
    
      stiff(0,0) = dummy_;
      stiff(1,1) = dummy_ * beta2_;
    }
    else
    {
      sigma  = ft_ * exp ( mftOverGf_ * ( eqvOpen - delta0_ ) );
      dsigma = sigma * mftOverGf_ ;
    
      double tem = ( dsigma * eqvOpen - sigma ) / pow( eqvOpen, 3. );

      // diagonal terms

      stiff(0,0) = sigma / eqvOpen          + crOpen * crOpen * tem;
      stiff(1,1) = sigma / eqvOpen * beta2_ + beta2_ * beta2_ * crSlid * crSlid * tem;

      // coupling terms

      stiff(0,1) = crOpen * crSlid * beta2_ * tem;
      stiff(1,0) = stiff(0,1);

      sigma /= eqvOpen;
    }
    
    traction[0] = sigma * crOpen;
    traction[1] = sigma * crSlid * beta2_;
  }
  else     // secant unloading
  {
    if ( eqvOpenMax < delta0_ ) 
    { 
      sigma = dummy_ * eqvOpenMax;
    }
    else
    {
      sigma = ft_ * exp ( mftOverGf_ * ( eqvOpenMax - delta0_ )  );
    }
      
    sigma /= eqvOpenMax;
  
    traction[0] = sigma * crOpen;
    traction[1] = sigma * crSlid * beta2_;

    stiff(0,0) = sigma;
    stiff(1,1) = sigma * beta2_;
  }

  // special case of penetration : negative jump

  if (!tension)
  {
      traction[0] = dummy_ * crOpen;
      stiff(0,1)  = stiff(1,0) = 0.;
      stiff(0,0)  = dummy_;
  }

  // debuging stuff

  /*System::out () << "eqv crack opening: " << eqvOpen     << "\n"
                 << "open             : " << jump        << "\n" 
                 << "eqv crack max    : " << eqvOpenMax  << "\n" 
                 << "traction         : " << traction    << "\n"
                 << "local tangent    : " << stiff       << "\n";
   */
  // compute damage variable for visualizing the cracks in a continuum
  // formulation. Only update damage for loading points.
  // Be careful with negative traction.
  
  if  ( ( eqvOpen > delta0_ ) && ( load ) && (tension) )
  {
      double d = 1. - traction[0]/ft_;
      damage_[mpoint] = d;
  }
}


//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void  ElasticExponentialCohesiveMat::elasticUpdate

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump 
    )
{
}

//-----------------------------------------------------------------------
//   update (moonen material)
//-----------------------------------------------------------------------

void  ElasticExponentialCohesiveMat::update

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

void  ElasticExponentialCohesiveMat::update

  ( Vector&               traction,
    Matrix&               stiff,
    const Vector&         jump0,
    const Vector&         djump,
    int                   mpoint )

{
  Vector shiftedJump ( jump0 + djump + shift_[mpoint] );

  update ( traction, stiff, shiftedJump, mpoint );
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void  ElasticExponentialCohesiveMat::update

  ( Vector&               traction,
    const Vector&         jump,
    int                   mpoint )
{
}

// --------------------------------------------------------------------
//  commit
// --------------------------------------------------------------------

void  ElasticExponentialCohesiveMat::commit()
{
  newHist_.eqvOpenMax. swap ( preHist_.eqvOpenMax );
  newHist_.loading.    swap ( preHist_.loading    );
}

// --------------------------------------------------------------------
//  allocPoints
// --------------------------------------------------------------------

void  ElasticExponentialCohesiveMat::allocPoints 

  ( int count, double dam )
{
  // not yet allocated

  if ( preHist_.eqvOpenMax.size ( ) == 0 ) 
  {
    preHist_.eqvOpenMax.resize ( count );  
    preHist_.loading.resize ( count );

    newHist_.eqvOpenMax.resize ( count );
    newHist_.loading.resize ( count );

    preHist_.eqvOpenMax = 0.0;
    preHist_.loading = 0  ;
    newHist_.eqvOpenMax = 0.0;
    newHist_.loading = 0  ;
  }

  // appending at the end of the vector

  else
  {
    for ( int i = 0; i < count; i++ )
    {
      preHist_.eqvOpenMax.pushBack ( 0.0 );
      preHist_.loading.pushBack ( 0   );

      newHist_.eqvOpenMax.pushBack ( 0.0 );  
      newHist_.loading.pushBack ( 0   );  
    }
  }
  
  for ( int i = 0; i < count; ++i )
  {
    shift_.pushBack ( Vector ( rank_ ) );
    
    shift_[i] = 0.;
  }

  damage_.resize ( count );
  damage_ = 0.;
}

//-----------------------------------------------------------------------
//   initShift
//-----------------------------------------------------------------------

void  ElasticExponentialCohesiveMat::initShift 

  ( const int             ipoint,
    const Vector&         traction )

{
  TensorIndex i;
  Vector shift ( shift_[ipoint] ); // shallow copy

  shift = traction / dummy_;

  bool tension   = shift[0] > 0;
       shift[0] *= tension;
       
  if (beta_)     shift[1] /= beta2_;

  System::out () << traction << " " << shift << "\n";
}

//-----------------------------------------------------------------------
//   justTractionFree_ 
//-----------------------------------------------------------------------

bool  ElasticExponentialCohesiveMat::justTractionFree

  ( const int   i ) const

{
  return false;
}

//-----------------------------------------------------------------------
//   evalFailure 
//-----------------------------------------------------------------------

double  ElasticExponentialCohesiveMat::evalFailure

    ( const Vector&       sigmaN,
      int                 mpoint )  const

{
    return 0.;
}

// --------------------------------------------------------------------
//  scaleDissipationIncrement
// --------------------------------------------------------------------

void  ElasticExponentialCohesiveMat::scaleDissipationIncrement

  ( const double     factor,
    const int        ipoint )

{
}
