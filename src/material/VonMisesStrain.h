#ifndef VON_MISES_STRAIN_H
#define VON_MISES_STRAIN_H


#include <jive/Array.h>

using jive::Vector;

//=======================================================================
//   class VonMisesStrain
//=======================================================================

class VonMisesStrain
{

public:

            VonMisesStrain 

   ( double k = 0.0, double nu = 0.0 );


  double    operator () 

    ( const double I1, 
      const double J2 )           const;

  Vector    operator () 

    ( const Vector& dI1dStrain, 
      const Vector& dJ2dStrain,
      const double  I1,
      const double  J2 )          const;

  void      init

    ( double          k,
      double          nu );


private:

  double    a_; 
  double    b_;
  double    c_; 

};

#endif
