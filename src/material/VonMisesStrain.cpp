#include <jem/base/assert.h>
#include <jem/base/RuntimeException.h>
#include <jem/base/Array.h>
#include <jem/numeric/utilities.h>
#include <jem/base/System.h>

#include "VonMisesStrain.h"


// -------------------------------------------------------------------
//  constructor
// -------------------------------------------------------------------

VonMisesStrain::VonMisesStrain

( double k, double nu )

{
}

// -------------------------------------------------------------------
//  init
// -------------------------------------------------------------------

void VonMisesStrain::init

  ( double  k,
    double  nu )

{
  a_ = (k - 1.0) / (1.0 - 2.0 * nu);
  b_ =  0.5 / k;
  c_ = 12.0 * k / (1.0 + 2.0 * nu + nu * nu);
}

// ---------------------------------------------------------------------
//  double  operator() (I1, J2)
// ---------------------------------------------------------------------

double VonMisesStrain::operator ()

  ( const double  I1,
    const double  J2 )             const

{
  const double  d = a_ * a_ * I1 * I1 + c_ * J2;

  if ( d < 0 )
  {
    using namespace jem;
    throw RuntimeException (
      JEM_FUNC,
      " invalid parameters for von Mises equivalent strain !!! "
    );
  }

  return b_ * (a_ * I1 + ::sqrt( d ));
}

// ---------------------------------------------------------------------
//  Vector  operator() (dI1, dJ2)
// ---------------------------------------------------------------------

Vector  VonMisesStrain::operator () 

    ( const Vector& dI1dStrain, 
      const Vector& dJ2dStrain, 
      const double  I1,
      const double  J2 )            const

{
  const int s  = dI1dStrain.size();

  Vector    ret   ( s );

  double  d = a_ * a_ * I1 * I1 + c_ * J2;

  if ( d < 0 )
  {
    using namespace jem;
    throw RuntimeException (
      JEM_FUNC,
      " invalid parameters for von Mises equivalent strain !!! "
    );
  }

  d = ::sqrt(d);
    
  if ( d == 0.0 )
  {
    d = 1.0; 
  }
  
  const double  e      = b_ / d ;
  
  const double  deqdI1 = a_ * b_  + e * a_ * a_ * I1; 
  const double  deqdJ2 = e * c_ * 0.5;

  ret = deqdI1 * dI1dStrain + deqdJ2 * dJ2dStrain;

  return ret;
}
