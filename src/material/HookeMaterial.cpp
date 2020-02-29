/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the isotropic elastic material
 *  This represents the material at a point in space.
 *  It is implemented in such a way that can be used for any
 *  discretisation methods, say finite elements, EFG and so on.
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *
 *  7 October 2014: make De 4x4 matrix for 2D problems (VPN).
 *
 */


#include <jem/base/System.h>
#include <jem/base/limits.h>
#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/EigenUtils.h>

#include "util/utilities.h"
#include "HookeMaterial.h"

using namespace jem;
using namespace jem::io;
using jem::numeric::matmul;

const double one_third = 0.33333334;

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  HookeMaterial::YOUNG_PROP   = "young";
const char*  HookeMaterial::POISSON_PROP = "poisson";
const char*  HookeMaterial::RHO_PROP     = "rho";
const char*  HookeMaterial::STATE_PROP   = "state";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


HookeMaterial::HookeMaterial 

  ( int rank, const Properties& globdat )
    : Material ( rank, globdat )
{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );

  young_   = 1.0;
  poisson_ = 1.0;
  rho_     = 1.0;

  stiffMat_ .resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  stiffMat_ = 0.0;
}


HookeMaterial::~HookeMaterial ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void HookeMaterial::configure ( const Properties& props,
                                const Properties& globdat )
{
  using jem::maxOf;

  props.find ( young_,   YOUNG_PROP,   0.0, maxOf( young_ ) );
  props.find ( poisson_, POISSON_PROP, 0.0, 0.5 );
  
  double K, G;

  if ( ( props.find ( K, "bulk_modulus"  ) ) && 
       ( props.find ( G, "shear_modulus" ) ) )
  {
    young_   = 9*K*G / ( 3*K + G );
    poisson_ = ( 3*K - 2*G ) / 2. / ( 3*K + G );
  }

  //System::out() << props << "\n";

  if ( rank_ == 1 )
  {
    props.get ( area_, "area" );
  }

  props.find ( rho_, RHO_PROP, 0.0, maxOf( rho_ ) );

  // read problem type, plane stress ...

  String state;

  if ( rank_ == 2  )
  {
    props.get( state, STATE_PROP );

    if      ( state == "PLANE_STRAIN" )
    {
      state_ = PlaneStrain;
    }
    else if ( state == "PLANE_STRESS" )
    {
      state_ = PlaneStress;
    }
    else if ( state == "AXISYMMETRIC" )
    {
      state_ = AxiSymmetric;
    }
  }

  // compute the elastic moduli, only once time

  computeStiffMat_ ();
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void HookeMaterial::getConfig ( const Properties& conf, 
                                const Properties& globdat ) const
{
  conf.set ( YOUNG_PROP,   young_   );
  conf.set ( POISSON_PROP, poisson_ );
  conf.set ( STATE_PROP,   state_   );
  conf.set ( RHO_PROP,     rho_     );
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void HookeMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      int                   ipoint )

{
  // compute the elastic moduli

  stiff = stiffMat_;

  // compute the stress vector

  matmul ( stress, stiffMat_, strain );
  
}

//-----------------------------------------------------------------------
//   giveHistory
//-----------------------------------------------------------------------

double HookeMaterial::giveHistory(int ip) const
{
  return 0.0;
}

//-----------------------------------------------------------------------
//   getStiffMat
//-----------------------------------------------------------------------


Matrix HookeMaterial::getStiffMat() const
{
  return stiffMat_;
}
  
//-----------------------------------------------------------------------
//   computeStiffMat_
//-----------------------------------------------------------------------


void   HookeMaterial::computeStiffMat_ () 
{
  const int     n  = STRAIN_COUNTS[rank_];

  const double  e  = young_;
  const double  nu = poisson_;


  if      ( rank_ == 1 )
  {
    stiffMat_(0,0) = e * area_;
  }
  else if ( rank_ == 3 )
  {
    const double  a = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double  b = 0.5 * (1.0 - 2.0 * nu);
    const double  c = 1.0 - nu;

    stiffMat_(0,0) = a * c;
    stiffMat_(0,1) = stiffMat_(0,2) = a * nu;
    stiffMat_(1,1) = stiffMat_(0,0);
    stiffMat_(1,2) = stiffMat_(0,1);
    stiffMat_(2,2) = stiffMat_(0,0);
    stiffMat_(3,3) = a * b;
    stiffMat_(4,4) = stiffMat_(3,3);
    stiffMat_(5,5) = stiffMat_(3,3);

    // Copy lower triangle of the stress-strain matrix.

    for ( int i = 0; i < n; i++ )
    {
      for ( int j = 0; j < i; j++ )
      {
	    stiffMat_(i,j) = stiffMat_(j,i);
      }
    }
  }
  else if ( state_ == PlaneStrain )
  {
    const double  a = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double  b = 0.5 * (1.0 - 2.0 * nu);
    const double  c = 1.0 - nu;

    stiffMat_(0,0) = a * c;
    stiffMat_(0,1) = stiffMat_(1,0) = a * nu;
    stiffMat_(0,2) = a * nu;
    stiffMat_(1,1) = a * c;
    stiffMat_(1,2) = a * nu;
    stiffMat_(3,3) = a * b;
    stiffMat_(2,0) = a * nu; stiffMat_(2,1) = a * nu; stiffMat_(2,2) = a * c;
  }
  else if ( state_ == PlaneStress )
  {
    const double  a = e / (1.0 - nu * nu);

    stiffMat_(0,0) = a;
    stiffMat_(0,1) = stiffMat_(1,0) = a * nu;
    stiffMat_(1,1) = a;
    stiffMat_(3,3) = a * 0.5 * (1.0 - nu);
  }
  else
  {
    throw Error ( JEM_FUNC, "unexpected rank: " + String ( rank_ ) );
  }

}



//-----------------------------------------------------------------------
//   giveStrainVector (given FEM strain vector)
//-----------------------------------------------------------------------

Vector HookeMaterial:: giveStrainVector

    ( const Vector& strain ) const
   
{
  const int strCount = strain.size ();

  JEM_PRECHECK ( strCount == 1 ||
                 strCount == 3 ||
                 strCount == 6 );

  if      ( strCount == 1 )
  {
    return strain;
  }
  else if ( strCount == 6 )
  {
    Vector ret ( 6 );

    ret[0] = strain[0];
    ret[1] = strain[1];
    ret[2] = strain[2];
    
    ret[3] = 0.5 * strain[3];
    ret[4] = 0.5 * strain[4];
    ret[5] = 0.5 * strain[5];

    return ret;
  }
  else if ( state_ == PlaneStrain )
  {
    Vector ret ( 3 );

    ret[0] = strain[0];
    ret[1] = strain[1];
    ret[2] = 0.5 * strain[2];
    
    return ret;
  }

  // Plane stress, epsilon_zz put in the last position
  
  else 
  {
    Vector ret ( 4 );

    double epsXX = strain[0];
    double epsYY = strain[1];
    
    ret[0] = epsXX;
    ret[1] = epsYY;
    ret[2] = 0.5 * strain[2];
    ret[3] = - poisson_ / ( 1.0 - poisson_ ) * ( epsXX + epsYY );

    return ret;
  }

}


//-----------------------------------------------------------------------
//   giveStressVector
//-----------------------------------------------------------------------

Vector HookeMaterial:: giveStressVector

    ( const Vector& stress ) const

{
  // Plane strain, sigma_zz put in last position

  if ( state_ == PlaneStrain )
  {
    Vector ret ( 4 );

    double sigmaXX = stress[0];
    double sigmaYY = stress[1];

    ret[0]         = sigmaXX;
    ret[1]         = sigmaYY;
    ret[2]         = stress[2];
    ret[3]         = poisson_ * ( sigmaXX + sigmaYY );

    return ret;
  }
  else
  {
    return stress;
  }
}

//-----------------------------------------------------------------------
//   getFirstInvariant
//-----------------------------------------------------------------------

// strain = { e_xx, e_yy, e_zz, 2*e_xy, 2*e_yz, 2*e_zx}^T

double   HookeMaterial:: getFirstInvariant 

    ( const Vector&         v )          const 

{

  double ret;

  if      ( state_ ==  PlaneStrain )
  {
    ret = v[0] + v[1];
  }

  else if ( state_ ==  PlaneStress )
  {
    ret = (1.0 - 2.0 * poisson_) / (1.0 - poisson_)  * ( v[0] + v[1] );
  }

  else 
  {
    ret =  v[0] + v[1] + v[2];
  }

  return ret;
}

//-----------------------------------------------------------------------
//   getSecondInvariant
//-----------------------------------------------------------------------

 double   HookeMaterial::getSecondInvariant 

    ( const Vector&         v )          const
{

  double I2;
 
  if      ( state_ == PlaneStrain )
  {
    I2 = v[0] * v[1] - 0.25 * v[2] * v[2] ;
  }
  else if ( state_ == PlaneStress )
  {
    I2 = v[0] * v[1] +  (1.0 - 2.0 * poisson_) / (1.0 - poisson_) * (v[0] + v[1]) 
         - 0.25 * v[2] * v[2];
  }
  else
  {
    I2 = v[0] * v[1] + v[1] * v[2] + v[2] * v[0] - 0.25 * (
         v[3] * v[3] - v[4] * v[4] - v[5] * v[5] );
  }

  return I2;
}

//-----------------------------------------------------------------------
//   getThirdInvariant
//-----------------------------------------------------------------------

double  HookeMaterial::getThirdInvariant 

    ( const Vector&         v )          const
{
  const int strCount = v.size ();

  double I3;

  if      ( strCount == 4 )
  {
    I3 = v[0] * v[1] * v[3] - v[3] * v[2] * v[2];
  }  
  else if ( strCount == 6 )
  {
    I3 = v[0] * v[1] * v[2] + 2.0 * v[3] * v[4] * v[5] -
         v[0] * v[4] * v[4] - v[1]* v[5] * v[5] - v[2] * v[3] * v[3];
  }
  else
  {
    I3 = 0.0;
  }

  return I3;
}

//-----------------------------------------------------------------------
//   get2ndInvDevStrain 
//-----------------------------------------------------------------------

// 2d strain = { e_xx, e_yy, 2*e_xy}^T
// 3d strain = { e_xx, e_yy, e_zz, 2*e_xy, 2*e_yz, 2*e_zx}^T

double    HookeMaterial::get2ndInvDevStrain 

    ( const Vector&         v )          const
{
  double J2;

  if ( state_ ==  PlaneStrain )
  {
    J2 = one_third * ( v[0] * v[0] + v[1] * v[1] - v[0] * v[1])
                + 0.25 * v[2] * v[2];
  }

  if ( state_ ==  PlaneStress )
  {
    double fac = poisson_ / ( 1.0 - poisson_ );

    J2  = one_third * ( v[0] * v[0] + v[1] * v[1] + 
			       ( fac * fac + fac ) * ( v[0] + v[1] ) * ( v[0] + v[1] )
                             - v[0] * v[1] ) + 0.25 * v[2] * v[2];
  }

  else
  {
    J2 = one_third * ( v[0] * v[0] + v[1] * v[1] + v[2] * v[2]   -
	   	       v[0] * v[1] - v[1] * v[2] - v[2] * v[0] ) +
                       v[3] * v[3] + v[4] * v[4] + v[5] * v[5];
  }

  return J2;

}
  
// ---------------------------------------------------------------------
//  compute  the derivatives of first  invariant w.r.t a vector
// ---------------------------------------------------------------------

// 2d strain = { e_xx, e_yy, 2*e_xy}^T
// 3d strain = { e_xx, e_yy, e_zz, 2*e_xy, 2*e_yz, 2*e_zx}^T

Vector  HookeMaterial:: getDI1DStrain

  ( const Vector& strain )                      const
  
{
  const int strCount = strain.size();

  Vector I ( strCount );

  if ( state_ ==  PlaneStrain )
  {
    I[0] = 1.0;
    I[1] = 1.0;
    I[2] = 0.0;
  }

  if ( state_ ==  PlaneStress )
  {
    double f = (1.0 - 2.0 * poisson_) / (1.0 - poisson_);

    I[0] = f;
    I[1] = f;
    I[2] = 0.0;    
  }

  if ( strCount == 6)  
  {
    I[0] = 1.0;
    I[1] = 1.0;
    I[2] = 1.0;
    
    I[3] = 0.0;
    I[4] = 0.0;
    I[5] = 0.0;    
  }

  return I;
}

// ---------------------------------------------------------------------
//  compute  the derivatives of second  invariant w.r.t a vector
// ---------------------------------------------------------------------

Vector  HookeMaterial:: getDI2DStrain

  ( const Vector&      v )                 const
  
{
  const int strCount = v.size();

  Vector    ret ( strCount );

  if      ( state_ ==  PlaneStrain )
  {
    ret[0] =   v[1];
    ret[1] =   v[0];
    
    ret[2] = - v[2];  
  }
  
  else if ( state_ ==  PlaneStress )
  {
    double f = ( 2.0 * poisson_ ) / ( 1.0 - poisson_ );

    ret[0]   =  (1.0 - f) * v[1] - f * v[0];
    ret[1]   =  (1.0 - f) * v[0] - f * v[1];

    ret[2]   = - v[2];
  }
  else
  {
    ret[0] = v[1] + v[2];
    ret[1] = v[2] + v[0];
    ret[2] = v[0] + v[1];

    ret[3] = - v[3];
    ret[4] = - v[4];
    ret[5] = - v[5];
  }

  return ret;
}

// ---------------------------------------------------------------------
//  compute  the derivatives of third invariant w.r.t a vector
// ---------------------------------------------------------------------

Vector  HookeMaterial::getDI3DStrain

  ( const Vector&      v )               const
  
{
  const int strCount = v.size();

  Vector    ret ( strCount );

  if      ( state_ ==  PlaneStrain )
  {
    ret = 0.0;
  }
  else if ( state_ ==  PlaneStress )
  {
    double f = poisson_ / ( 1.0 - poisson_ );

    ret[0] = -f * v[1] * 2.0 * v[0] - f * v[1] * v[1] + f * 0.25 * v[2];
    ret[1] = -f * v[0] * 2.0 * v[1] - f * v[0] * v[0] + f * 0.25 * v[2];

    ret[0] =  f * ( v[0] + v[1] ) * v[2];
  }
  else
  {
    ret[0] = v[1] * v[2] - v[4] * v[4];
    ret[1] = v[2] * v[0] - v[5] * v[5];
    ret[2] = v[0] * v[1] - v[3] * v[3];

    ret[3] = 2.0 * ( v[4] * v[5] - v[2] * v[3] );
    ret[4] = 2.0 * ( v[5] * v[3] - v[0] * v[4] );
    ret[5] = 2.0 * ( v[3] * v[4] - v[1] * v[5] );
  }

  return ret;
}


// ---------------------------------------------------------------
//   getIandDIDStrain
// ---------------------------------------------------------------

// compute invariants, principal strains and their derivatives w.r.t
// the strains. 
// Implemented only for 3D for efficiency

void  HookeMaterial::getIandDIDStrain

  ( const Vector& v,
          Vector& I,   Vector& princ,
          Vector& dI1, Vector& dI2, Vector& dI3 )  const

{
  const int    strCount  = v.size();

  JEM_PRECHECK ( strCount == 6 );

  const double exx = v[0];
  const double eyy = v[1];
  const double ezz = v[2];

  const double exy = 0.5 * v[3];
  const double eyz = 0.5 * v[4];
  const double ezx = 0.5 * v[5];

  // compute invarianrs

  double I1 = exx + eyy + ezz;
  double I2 = exx * eyy + eyy * ezz + ezz * exx - 
              exy * exy - eyz * eyz - ezx * ezx;

  double I3 = exx * eyy * ezz + 2.0 * exy * eyz * ezx -
              exx * eyz * eyz - eyy * ezx * ezx - ezz * exy * exy;

  // store I1 and I2 in vector I

  I[0]   = I1;
  I[1]   = I2;

  // compute principal strains

  Vector temp;
  
  solveCubicEqua ( temp, 1.0, -I1, I2, -I3 );

  princ = 0.0;

  const int no = temp.size();

  if ( no > 0 )
  {
    princ[0] = temp[0];
  }
  
  if ( no > 1 )
  {
    princ[1] = temp[1];
  }
  
  if ( no > 2 )
  {
    princ[2] = temp[2];
  }
  
  // derivatives of I1 w.r.t strains

  dI1[0] = 1.0;
  dI1[1] = 1.0;
  dI1[2] = 1.0;
    
  dI1[3] = 0.0;
  dI1[4] = 0.0;
  dI1[5] = 0.0;    

  // derivatives of I2 w.r.t strains

  dI2[0] = eyy + ezz;
  dI2[1] = ezz + exx;
  dI2[2] = exx + eyy;

  dI2[3] = - 2.0 * exy;
  dI2[4] = - 2.0 * eyz;
  dI2[5] = - 2.0 * ezx;

  // derivatives of I3 w.r.t strains

  dI3[0] = eyy * ezz - eyz * eyz;
  dI3[1] = ezz * exx - ezx * ezx;
  dI3[2] = exx * eyy - exy * exy;

  dI3[3] = 2.0 * ( eyz * ezx - ezz * exy );
  dI3[4] = 2.0 * ( ezx * exy - exx * eyz );    
  dI3[5] = 2.0 * ( exy * eyz - eyy * ezx );    

}

// ------------------------------------------------------------------
//  getI1J2andGrads   
// ------------------------------------------------------------------

// 2d strain = { e_xx, e_yy, e_zz, 2*e_xy}^T
// 3d strain = { e_xx, e_yy, e_zz, 2*e_xy, 2*e_yz, 2*e_zx}^T

void   HookeMaterial::getI1J2andGrads   

  ( const Vector&       v,
    double& I1,
    double& J2, 
    Vector& dI1, Vector& dJ2 )  const

{
  const int    strCount  = v.size();

  if      ( state_ ==  PlaneStrain ) // e_zz = 0
  {

    const double exx = v[0];
    const double eyy = v[1];
    const double exy = v[3];

    I1 = exx + eyy;

    J2 = one_third * ( exx * exx + eyy * eyy - exx * eyy )
                       + 0.25 * exy * exy;

    dI1[0] = 1.0;
    dI1[1] = 1.0;
    dI1[2] = 0.0;
    dI1[3] = 0.0;

    dJ2[0] = one_third * ( 2.0 * exx - eyy );
    dJ2[1] = one_third * ( 2.0 * eyy - exx );
    dJ2[2] = one_third * ( -eyy - exx );
    dJ2[3] = 0.5 * exy ;
  }

  else if ( state_ ==  PlaneStress )
  {

    const double exx = v[0];
    const double eyy = v[1];
    //const double ezz = v[2]; // ???
    const double exy = v[3];


    double a = poisson_ / ( 1.0 - poisson_ );
    double c = 1.0 - a;;
    double b = a * a + a;

    I1 = c  * ( exx + eyy );
    J2 = one_third * ( exx * exx + eyy * eyy + b * ( exx + eyy ) * ( exx + eyy ) - exx * eyy ) + 0.25 * exy * exy;

    dI1[0] = c;
    dI1[1] = c;
    dI1[2] = 0.0;
    dI1[3] = 0.0;

    dJ2[0] = one_third * ( 2.0 * ( exx + b * ( exx + eyy ) ) - eyy );
    dJ2[1] = one_third * ( 2.0 * ( eyy + b * ( eyy + exx ) ) - exx );
    dJ2[2] = one_third * ( 2.0 * ( eyy + b * ( eyy + exx ) ) - exx );
    dJ2[3] = 0.5 * exy ;
  }

  else if ( strCount == 6 )
  {

    const double exx = v[0];
    const double eyy = v[1];
    const double ezz = v[2];
    const double exy = v[3];
    const double eyz = v[4];
    const double ezx = v[5];

    I1 = exx + eyy + ezz;

    J2 = one_third * ( exx * exx + eyy * eyy + ezz * ezz   -
	   	       exx * eyy - eyy * ezz - ezz * exx ) +
                       0.25 * ( exy * exy + eyz * eyz + ezx * ezx );

    dI1[0] = 1.0;
    dI1[1] = 1.0;
    dI1[2] = 1.0;

    dI1[3] = 0.0;
    dI1[4] = 0.0;
    dI1[5] = 0.0;

    dJ2[0] = one_third * ( 2.0 * exx - eyy - ezz );
    dJ2[1] = one_third * ( 2.0 * eyy - ezz - exx );
    dJ2[2] = one_third * ( 2.0 * ezz - exx - eyy );

    dJ2[3] = 0.5 * exy ;
    dJ2[4] = 0.5 * eyz ;
    dJ2[5] = 0.5 * ezx ;

  }

}

// ------------------------------------------------------------------
//  getI1andJ2
// ------------------------------------------------------------------

void   HookeMaterial::getI1andJ2

  ( const Vector&       v,
    double& I1,
    double& J2 )  const

{
  const int    strCount  = v.size();
  
  //System::out() << state_ << "\n";

  if      ( state_ ==  PlaneStrain )
  {

    const double exx = v[0];
    const double eyy = v[1];
    const double exy = v[2];

    I1 = exx + eyy;

    J2 = one_third * ( exx * exx + eyy * eyy - exx * eyy )
                       + 0.25 * exy * exy;

  }

  else if ( state_ ==  PlaneStress )
  {

    const double exx = v[0];
    const double eyy = v[1];
    const double exy = v[2];


    double a = poisson_ / ( 1.0 - poisson_ );
    double c = 1.0 - a;;
    double b = a * a + a;

    I1 = c  * ( exx + eyy );
    J2 = one_third * ( exx * exx + eyy * eyy + 
                       b * ( exx + eyy ) * ( exx + eyy ) -
                       exx * eyy ) + 0.25 * exy * exy;
  }

  else if ( strCount == 6 )
  {

    const double exx = v[0];
    const double eyy = v[1];
    const double ezz = v[2];
    const double exy = v[3];
    const double eyz = v[4];
    const double ezx = v[5];

    I1 =  exx + eyy + ezz;

    J2 = one_third * ( exx * exx + eyy * eyy + ezz * ezz   -
        	   	       exx * eyy - eyy * ezz - ezz * exx ) +
                       exy * exy + eyz * eyz + ezx * ezx;

  }
}

// ------------------------------------------------------------------
//  getI1andJ2Stress
// ------------------------------------------------------------------

void   HookeMaterial::getI1andJ2ForStress

  ( const Vector&       s,
    double& I1,
    double& J2 )  const

{
   const int    strCount  = s.size();

  if      ( state_ ==  PlaneStress )
  {
    const double sxx = s[0];
    const double syy = s[1];
    const double sxy = s[3];

    I1 = sxx + syy;

    J2 = one_third * ( sxx * sxx + syy * syy - sxx * syy ) + sxy * sxy;

  }
  else if ( state_ ==  PlaneStrain )
  {
  }
  else if ( strCount == 6 )
  {
    const double sxx = s[0];
    const double syy = s[1];
    const double szz = s[2];
    const double sxy = s[3];
    const double syz = s[4];
    const double szx = s[5];

    I1 =  sxx + syy + szz;

    J2 = one_third * ( sxx * sxx + syy * syy + szz * szz   -
        	   	       sxx * syy - syy * szz - szz * sxx ) +
                       sxy * sxy + syz * syz + szx * szx;
  }

}

// ------------------------------------------------------------------
//  getI1andJ2Stress
// ------------------------------------------------------------------

void   HookeMaterial::getI1andJ2ForStress

  ( const Vector&       v,
    const Vector&       S,
    double& I1,
    double& J2 )  const

{
   const int    strCount  = v.size();

  if      ( state_ ==  PlaneStress )
  {

    const double exx = v[0];
    const double eyy = v[1];
    const double exy = v[2];

    I1 = exx + eyy;

    J2 = one_third * ( exx * exx + eyy * eyy - exx * eyy )
                       + 0.25 * exy * exy;

    S[0] = exx - one_third * I1;
    S[1] = eyy - one_third * I1;
    S[2] = exy;

  }

  else if ( state_ ==  PlaneStrain )
  {

    const double exx = v[0];
    const double eyy = v[1];
    const double ezz = v[2];
    const double exy = v[3];


    I1 = exx + eyy + ezz;
    J2 = one_third * ( exx * exx + eyy * eyy + ezz * ezz   -
        	   	       exx * eyy - eyy * ezz - ezz * exx ) +
                       exy * exy ;
    
    S[0] = exx - one_third * I1;
    S[1] = eyy - one_third * I1;
    S[2] = ezz - one_third * I1;
    S[3] = exy;
  }

  else if ( strCount == 6 )
  {
    const double exx = v[0];
    const double eyy = v[1];
    const double ezz = v[2];
    const double exy = v[3];
    const double eyz = v[4];
    const double ezx = v[5];

    I1 =  exx + eyy + ezz;

    J2 = one_third * ( exx * exx + eyy * eyy + ezz * ezz   -
        	   	       exx * eyy - eyy * ezz - ezz * exx ) +
                       exy * exy + eyz * eyz + ezx * ezx;

    S = v;

    S[0] = exx - one_third*I1;
    S[1] = eyy - one_third*I1;
    S[2] = ezz - one_third*I1;

  }

}

//-----------------------------------------------------------------------
//   getPrincipalStrains
//-----------------------------------------------------------------------

Vector  HookeMaterial::getPrincipalStrains

  ( const Vector&       v )  const
  
{

  const int strCount = v.size();

  Vector  prinstr(3);          
  

  // PLANE PROBLEM

  if ( strCount == 4 ) 
  {

    double exx = v[0];
    double eyy = v[1];
    double exy = 0.5 * v[3];

    double poi = poisson_;

    // principal strains are root of a quadratic equation with det = d

    double d   = ( exx - eyy ) * ( exx - eyy ) + 4.0 * exy * exy;

    d = ::sqrt ( d );

    prinstr[0] = 0.5 * ( exx + eyy + d );
    prinstr[1] = 0.5 * ( exx + eyy - d );

    if ( state_ == PlaneStress )
    {
      prinstr[2] = poi / ( poi - 1.0 ) * ( exx + eyy );
    }
    else
    {
      prinstr[2] = 0.0;
    }
  }

  // THREE DIMENSIONS
  
  else 
  {
    
    double exx = v[0];
    double eyy = v[1];
    double ezz = v[2];

    double exy = 0.5 * v[3];
    double eyz = 0.5 * v[4];
    double ezx = 0.5 * v[5];

    double I1  = exx + eyy + ezz;
    double I2  = exx * eyy + eyy * ezz + ezz * exx - 
                 exy * exy - eyz * eyz - ezx * ezx; 

    double I3  = exx * eyy * ezz + 2.0 * exy * eyz * ezx -
                 exx * eyz * eyz - eyy * ezx * ezx - ezz * exy * exy;

    // solve the cubic equation:
    //  epsilon^3 -I1 * epsilon^2 + I2 * epsilon - I3 = 0
    // make sure there are always solutions for symmetric 
    // second order tensor. This is adopted from OOFEM.

    Vector temp;
    
    solveCubicEqua ( temp, 1.0, -I1, I2, -I3 );

    prinstr = 0.0;

    const int no = temp.size();

    if ( no > 0 )
    {
      prinstr[0] = temp[0];
    }

    if ( no > 1 )
    {
      prinstr[1] = temp[1];
    }

    if ( no > 2 )
    {
      prinstr[2] = temp[2];
    }

  }

  return prinstr;
}

//-----------------------------------------------------------------------
//   getPrincipalValues
//-----------------------------------------------------------------------

Vector  HookeMaterial::getPrincipalValues

  ( const Vector&       v )  const
  
{
  const int strCount = v.size();

  JEM_PRECHECK ( strCount == 1 ||
                 strCount == 3 ||
		             strCount == 4 ||
                 strCount == 6 );

  Vector ans;

  // Plane problem

  if      ( strCount == 3 )
  {
    double I1 = v[0] + v[1];
    double I2 = v[0] * v[1] - v[2] * v[2];

    double D  = I1 * I1 - 4.0 * I2;

    if ( D < 0 )
    {
      using namespace jem;
      throw Error (
      JEM_FUNC,
      " imaginary principal values !!! "
      );  
    }

    ans.resize(2);

    ans[0] = 0.5 * ( I1 + sqrt( D ) );
    ans[1] = 0.5 * ( I1 - sqrt( D ) );

    //System::out() << ans[0] << " " << ans[1] <<"\n\n";

  }

  // Plane stress : epsilon_zz non zero
  // Plane strain : sigma_zz non zero
  // In both cases, there are three principal values
  
  else if ( strCount == 4 )
  {
    double I1 = v[0] + v[1];
    double I2 = v[0] * v[1] - v[3] * v[3];

    double D  = I1 * I1 - 4.0 * I2;

    /*

    if ( D < 0 )
    {
      using namespace jem;
      throw Error (
      JEM_FUNC,
      " imaginary principal stress !!! "
      );
    }*/

    D = sqrt ( D );

    ans.resize(3);

    ans[0] = 0.5 * ( I1 + D );
    ans[1] = 0.5 * ( I1 - D );
    ans[2] = v[3];
  }
  
  // 3D problem
  
  else 
  {
    double I1 = getFirstInvariant  ( v );
    double I2 = getSecondInvariant ( v );
    double I3 = getThirdInvariant  ( v );

    Vector temp;
    
    solveCubicEqua ( temp, 1.0, -I1, I2, -I3 );

    ans.resize( 3 );
    ans = 0.0;

    const int no = temp.size();

    if ( no > 0 )
    {
      ans[0] = temp[0];
    }

    if ( no > 1 )
    {
      ans[1] = temp[1];
    }

    if ( no > 2 )
    {
      ans[2] = temp[2];
    }
  }
  
  return ans;
}


//-----------------------------------------------------------------------
//   getPrincipalValues + return invariants
//-----------------------------------------------------------------------

void    HookeMaterial:: getPrincipalValues

  ( Vector&             I,
    Vector&             vI,
    const Vector&       v ) const
  
{
  const int strCount = v.size();

  JEM_PRECHECK ( strCount == 1 ||
                 strCount == 3 ||
		             strCount == 4 ||
                 strCount == 6 );

  // Plane problem
  
  if      ( strCount == 3 )
  {
    double I1 = v[0] + v[1];
    double I2 = v[0] * v[1] - v[2] * v[2];

    double D  = I1 * I1 - 4.0 * I2;

    if ( D < 0 )
    {
      using namespace jem;
      throw Error (
      JEM_FUNC,
      " imaginary principal strains !!! "
      );  
    }

    vI.resize ( 2 );
    
    vI[0] = 0.5 * ( I1 + sqrt( D ) );
    vI[1] = 0.5 * ( I1 - sqrt( D ) );

    I.resize ( 2 );

    I[0]  = I1;
    I[1]  = I2;
  }

  // Plane stress : epsilon_zz non zero
  // Plane strain : sigma_zz non zero
  // In both cases, there are three principal values
  
  else if ( strCount == 4 )
  {
    double I1 = v[0] + v[1];
    double I2 = v[0] * v[1] - v[2] * v[2];

    double D  = I1 * I1 - 4.0 * I2;

    if ( D < 0 )
    {
      using namespace jem;
      throw Error (
      JEM_FUNC,
      " imaginary principal strains !!! "
      );  
    }

    vI.resize ( 3 );

    vI[0] = 0.5 * ( I1 + sqrt( D ) );
    vI[1] = 0.5 * ( I1 - sqrt( D ) );
    vI[2] = v[3];

    I.resize ( 3 );

    I[0]  = I1;
    I[1]  = I2;
    I[2]  = v[0] * v[1] * v[3] - v[3] * v[2] * v[2];
  }
  
  // 3D problem
  
  else 
  {
    const double I1 = getFirstInvariant  ( v );
    const double I2 = getSecondInvariant ( v );
    const double I3 = getThirdInvariant  ( v );
    
    Vector temp;

    // call the cubic solver works only for case
    // having three real roots (true for symmetric tensor like strain
    // and stress
    
    solveCubicEqua ( temp, 1.0, -I1, I2, -I3 );

    vI.resize( 3 );
    vI = 0.0;

    const int no = temp.size();

    if ( no > 0 )
    {
      vI[0] = temp[0];
    }

    if ( no > 1 )
    {
      vI[1] = temp[1];
    }

    if ( no > 2 )
    {
      vI[2] = temp[2];
    }

    I.resize ( 3 );

    I[0]  = I1;
    I[1]  = I2;
    I[2]  = I3;
  }
}

//-----------------------------------------------------------------------
//   getPrincipalStressAndDirections
//-----------------------------------------------------------------------

// each COLUMN of matrix dir = direction corresponding to each principal
// stress

void    HookeMaterial::getPrincipalStressAndDirections

    ( const Vector&      stress,
      const Vector&      sigmaI,
      const Matrix&      dir )               const
{
  const int dim = stress.size ();

  double sigmaX, sigmaY, sigmaXY;

  double sigmaZ  = 0.0;
  double sigmaYZ = 0.0;
  double sigmaZX = 0.0;

  if    ( dim == 4 )
  {
     sigmaX  = stress[0];
     sigmaY  = stress[1];
     sigmaXY = stress[3];

     if  ( state_ == PlaneStrain )
     {
       sigmaZ = poisson_ * ( sigmaX + sigmaY );
     } 
  }
  else
  {
     sigmaX  = stress[0];
     sigmaY  = stress[1];
     sigmaZ  = stress[2];

     sigmaXY = stress[3];
     sigmaYZ = stress[4];
     sigmaZX = stress[5];
  }

  Matrix sigma(3,3);

  sigma(0,0) = sigmaX;
  sigma(0,1) = sigmaXY;
  sigma(0,2) = sigmaZX;

  sigma(1,1) = sigmaY;
  sigma(1,2) = sigmaYZ;

  sigma(2,2) = sigmaZ;

  //System::out() << stress << "\n";

  using jem::numeric::EigenUtils;

  EigenUtils::symSolve ( sigmaI, dir, sigma);
}

//-----------------------------------------------------------------------
//   checkLocalisation
//-----------------------------------------------------------------------

double    HookeMaterial::checkLocalisation 

    ( Vector&              normal,
      const Vector&        stress,
      const Matrix&        tangent,
            int            ipoint )          const
{
     Vector sigmaI ( 3 );
     Matrix dir    ( 3, 3 );
     getPrincipalStressAndDirections ( stress, sigmaI, dir );

     int    i   = 0;
     double max = sigmaI[0];

     if ( max < sigmaI[1] ){
         max = sigmaI[1]; i = 1;
     }
     
     if ( max < sigmaI[2] ){
         max = sigmaI[2]; i = 2;
     }

     normal = dir(slice(BEGIN,rank_),i);

     //normal[0]=1.; normal[1]=0.;

     return max;
}



