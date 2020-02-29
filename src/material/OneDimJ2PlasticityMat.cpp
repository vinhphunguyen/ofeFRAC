
#include <cmath>
#include <jem/base/System.h>
#include <jem/base/Error.h>
#include <jem/base/limits.h>
#include <jem/base/Float.h>
#include <jem/base/Error.h>
#include <jem/util/Properties.h>

#include "OneDimJ2PlasticityMat.h"

using namespace jem;
using jem::util::Properties;

// ------------------------------------------------
//  Constructor
// ------------------------------------------------

OneDimJ2PlasticityMat::OneDimJ2PlasticityMat 

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  props.get ( k_,      "plasticModulus" );
  props.get ( young_,   "young"         );
  props.get ( sigmaY_, "yieldStress"    );
  
  conf .set ( "plasticModulus", k_      );
  conf. set ( "young",          young_  );
  conf. set ( "yieldStress",    sigmaY_ );
}

// ------------------------------------------------
//  Destructor
// ------------------------------------------------

OneDimJ2PlasticityMat::~OneDimJ2PlasticityMat()
{}


// ------------------------------------------------
//  commit
// ------------------------------------------------

void   OneDimJ2PlasticityMat::update

    ( const Properties&      params )
{
  idx_t i;

  double strain, stress, eps, Cep;

  params.get ( i,         "i"        );
  params.get ( eps,        "strain"  );

  double  epsP0  = epsilonP0_[i]; 
  double  alpha0 = alpha0_   [i]; 

  double sigmaTrial = young_ * ( eps - epsP0 );
  double yield      = std::fabs ( sigmaTrial ) - ( sigmaY_ + k_ * alpha0 );

  if ( yield <= 0. )
  {
    stress = sigmaTrial;
    Cep    = young_;

    epsilonP_[i]  = epsP0;
    alpha_[i]     = alpha0;
  }
  else
  {
    double sigOfSigma = (sigmaTrial > 0) - (sigmaTrial < 0);
    double deltaGamma = yield / ( young_ + k_ );
    stress            = ( 1. - deltaGamma * young_ / std::fabs ( sigmaTrial ) ) * sigmaTrial;
    epsilonP_[i]      = epsP0  + deltaGamma * sigOfSigma;
    alpha_[i]         = alpha0 + deltaGamma;
    Cep               = young_ * k_ / ( young_ + k_ );
  }

  params.set ( "tangent", Cep    );
  params.set ( "stress",  stress );
}

// ------------------------------------------------
//  commit
// ------------------------------------------------

void   OneDimJ2PlasticityMat::commit()
{
  epsilonP0_ .swap ( epsilonP_  );
  alpha0_    .swap ( alpha_     );
}


void OneDimJ2PlasticityMat::allocPoints ( jem::idx_t count )
{
  epsilonP_ .resize ( count );       // plastic strain
  epsilonP0_.resize ( count );      // converged plastic strain
  alpha_    .resize ( count );          // hardening variable alpha
  alpha0_   .resize ( count );         //  converged value of alpha

  epsilonP_  = 0.;
  epsilonP0_ = 0.;
  alpha_     = 0.;
  alpha0_    = 0.;
}

