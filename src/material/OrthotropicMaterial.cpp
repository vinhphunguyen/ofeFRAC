#include <jem/base/limits.h>
#include <jem/base/Error.h>
#include <jem/base/System.h>
#include <jem/base/Array.h>
#include <jem/base/PrecheckException.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/algebra/EigenUtils.h>
#include <jem/numeric/algebra/LUSolver.h>
#include <jem/numeric/utilities.h>


#include "util/utilities.h"
#include "OrthotropicMaterial.h"

using namespace jem;
using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jem::numeric::invert;
using jem::io::endl;

const double one_third = 0.33333333;

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  OrthotropicMaterial::YOUNG_1_PROP    = "young1";
const char*  OrthotropicMaterial::YOUNG_2_PROP    = "young2";
const char*  OrthotropicMaterial::YOUNG_3_PROP    = "young3";
const char*  OrthotropicMaterial::POISSON_12_PROP = "poisson12";
const char*  OrthotropicMaterial::POISSON_23_PROP = "poisson23";
const char*  OrthotropicMaterial::POISSON_31_PROP = "poisson31";
const char*  OrthotropicMaterial::SHEAR_12_PROP   = "shear12";
const char*  OrthotropicMaterial::SHEAR_23_PROP   = "shear23";
const char*  OrthotropicMaterial::SHEAR_31_PROP   = "shear31";
const char*  OrthotropicMaterial::RHO_PROP        = "rho";
const char*  OrthotropicMaterial::STATE_PROP      = "state";
const char*  OrthotropicMaterial::THETA_PROP      = "theta";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


OrthotropicMaterial::OrthotropicMaterial 

  ( const int          rank,
    const Properties&  globdat )

  : Material ( rank, globdat )

{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );

  pCount_  = 0;

  young1_    = 1.0;
  young2_    = 1.0;
  young3_    = 1.0;
  poisson12_ =  .0;
  poisson23_ =  .0;
  poisson31_ =  .0;
  shear12_   = 0.5;
  shear23_   = 0.5;
  shear31_   = 0.5;
  rho_       = 1.0;
  theta_     = 0.0;

  transformMat_     . resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  transformMat_     = 0.0;
  transformMatInv_  . resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  transformMatInv_  = 0.0;
  tt_               . ref    ( transformMatInv_.transpose() );
  stiffMat_         . resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  stiffMat_         = 0.0;
  materialCompMat_  . resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  materialCompMat_  = 0.0;
  materialStiffMat_ . resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  materialStiffMat_ = 0.0;
}


OrthotropicMaterial::~OrthotropicMaterial ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void OrthotropicMaterial::configure 

  ( const Properties& props,
    const Properties& globdat )
{
  using jem::maxOf;

  JEM_PRECHECK ( rank_ > 1 );

  props.get ( young1_   , YOUNG_1_PROP   , 0.0, maxOf( young1_ ) );
  props.get ( young2_   , YOUNG_2_PROP   , 0.0, maxOf( young2_ ) );
  props.get ( poisson12_, POISSON_12_PROP, 0.0, 0.5 );
  props.get ( shear12_  , SHEAR_12_PROP  , 0.0, maxOf( shear12_) );

  props.find ( rho_, RHO_PROP, 0.0, maxOf( rho_ ) );

  props.find ( theta_, THETA_PROP, -90. , 90. );

  if ( rank_ == 2  )
  {
    props.get ( state_, STATE_PROP);
  }
  else
  {
    state_ = "NOT_PLANE";
  }

  if ( rank_ == 3 || state_ == "PLANE_STRAIN" )
  {
    props.get ( poisson23_, POISSON_23_PROP, 0.0, 0.5 );

    if ( props.find ( young3_, YOUNG_3_PROP, 0.0, maxOf( young3_ ) ) )
    {
      // completely orthotropic material

      props.get ( shear23_  , SHEAR_23_PROP  , 0.0, maxOf( shear23_) );
      props.get ( poisson31_, POISSON_31_PROP, 0.0, 0.5 );
      props.get ( shear31_  , SHEAR_31_PROP  , 0.0, maxOf( shear31_) );
    }
    else
    {
      // transversely isotropic material

      young3_ = young2_;
      poisson31_ = poisson12_;
      shear31_ = shear12_;
      shear23_ = young2_ / ( 2. + 2. * poisson23_ );
    }
  }

  // compute the elastic stiffness matrix, only one time

  computeTransformMats_ ();
  computeStiffMat_ ();
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void OrthotropicMaterial::getConfig

  ( const Properties& conf , 
    const Properties& globdat ) const

{
  if ( rank_ == 2 )
  {
    conf.set ( STATE_PROP     , state_     );
  }

  conf.set ( YOUNG_1_PROP   , young1_    );
  conf.set ( YOUNG_2_PROP   , young2_    );
  conf.set ( POISSON_12_PROP, poisson12_ );
  conf.set ( SHEAR_12_PROP  , shear12_   );
  conf.set ( THETA_PROP     , theta_     );

  if ( rank_ == 3 || state_ == "PLANE_STRAIN" )
  {
    conf.set ( YOUNG_3_PROP   , young3_    );
    conf.set ( POISSON_23_PROP, poisson23_ );
    conf.set ( POISSON_31_PROP, poisson31_ );
    conf.set ( SHEAR_23_PROP  , shear23_   );
    conf.set ( SHEAR_31_PROP  , shear31_   );
  }
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void OrthotropicMaterial::update

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

double OrthotropicMaterial::giveHistory(int ip) const
{
  return 0.0;
}

//-----------------------------------------------------------------------
//   giveStressAtPoint 
//-----------------------------------------------------------------------

Vector OrthotropicMaterial::giveStressAtPoint 

 ( const int     ipoint, 
   const Vector& strain ) const

{
  MatmulChain<double,1>    mc1;

  Vector stress = mc1.matmul ( stiffMat_, strain );

  return stress;
}


//-----------------------------------------------------------------------
//   getStiffMat
//-----------------------------------------------------------------------


Matrix OrthotropicMaterial::getStiffMat() const
{
  return stiffMat_;
}
  
//-----------------------------------------------------------------------
//   computeTransformMats_
//-----------------------------------------------------------------------


void   OrthotropicMaterial::computeTransformMats_ () 
{
  const double  pi = 3.14159265;
  const double  c  = cos( theta_ * pi / 180.0 );
  const double  s  = sin( theta_ * pi / 180.0 );
  const double  sc = s*c;
  const double  c2 = c*c;
  const double  s2 = s*s;
  
  if ( rank_ == 3 ) 
  {
    transformMat_(0,0) = c2;
    transformMat_(0,1) = s2;
    transformMat_(0,3) = 2.0 * sc;

    transformMat_(1,0) = s2;
    transformMat_(1,1) = c2;
    transformMat_(1,3) = - 2.0 * sc;

    transformMat_(2,2) = 1.0;

    transformMat_(3,0) = - sc;
    transformMat_(3,1) = sc;
    transformMat_(3,3) = c2 - s2;

    transformMat_(4,4) = c;
    transformMat_(4,5) = - s;

    transformMat_(5,4) = s;
    transformMat_(5,5) = c;

    transformMatInv_(0,0) = c2;
    transformMatInv_(0,1) = s2;
    transformMatInv_(0,3) = - 2.0 * sc;

    transformMatInv_(1,0) = s2;
    transformMatInv_(1,1) = c2;
    transformMatInv_(1,3) = 2.0 * sc;

    transformMatInv_(2,2) = 1.0;

    transformMatInv_(3,0) = sc;
    transformMatInv_(3,1) = - sc;
    transformMatInv_(3,3) = c2 - s2;

    transformMatInv_(4,4) = c;
    transformMatInv_(4,5) = s;

    transformMatInv_(5,4) = - s;
    transformMatInv_(5,5) = c;
  }
  else if ( rank_ == 2 )
  {
    transformMat_(0,0) = c2;
    transformMat_(0,1) = s2;
    transformMat_(0,2) = 2.0 * sc;

    transformMat_(1,0) = s2;
    transformMat_(1,1) = c2;
    transformMat_(1,2) = - 2.0 * sc;

    transformMat_(2,0) = - sc;
    transformMat_(2,1) = sc;
    transformMat_(2,2) = c2 - s2;

    transformMatInv_(0,0) = c2;
    transformMatInv_(0,1) = s2;
    transformMatInv_(0,2) = - 2.0 * sc;

    transformMatInv_(1,0) = s2;
    transformMatInv_(1,1) = c2;
    transformMatInv_(1,2) = 2.0 * sc;

    transformMatInv_(2,0) = sc;
    transformMatInv_(2,1) = - sc;
    transformMatInv_(2,2) = c2 - s2;
  }
  else
  {
    throw Error ( JEM_FUNC, "unexpected rank: " + String ( rank_ ) );
  }
}


//-----------------------------------------------------------------------
//   computeStiffMat_
//-----------------------------------------------------------------------


void   OrthotropicMaterial::computeStiffMat_ () 
{
  const double  e1   = young1_;
  const double  e2   = young2_;
  const double  e3   = young3_;
  const double  nu12 = poisson12_;
  const double  nu23 = poisson23_;
  const double  nu31 = poisson31_;
  const double  g12  = shear12_;
  const double  g23  = shear23_;
  const double  g31  = shear31_;

  MatmulChain<double,3>    mChain;

  if ( rank_ == 3 ) 
  {
    materialCompMat_(0,0) = 1.0 / e1;
    materialCompMat_(1,1) = 1.0 / e2;
    materialCompMat_(2,2) = 1.0 / e3;

    materialCompMat_(0,1) = materialCompMat_(1,0) = -nu12 / e1;
    materialCompMat_(0,2) = materialCompMat_(2,0) = -nu31 / e1;
    materialCompMat_(1,2) = materialCompMat_(2,1) = -nu23 / e2;

    materialCompMat_(3,3) = 1.0 / g12;
    materialCompMat_(4,4) = 1.0 / g23;
    materialCompMat_(5,5) = 1.0 / g31;

    materialStiffMat_ = materialCompMat_;

    invert( materialStiffMat_ );
  }
  else if ( state_ == "PLANE_STRAIN" )
  {
    // NB: the plane strain materialCompMat_ contains the 
    // slice(0,3)-part of the full 3D compliance matrix
    // This cannot be used to compute strain from stress but it can be
    // for update and inversion to compute damaged stiffness matrix

    materialCompMat_(0,0) = 1.0 / e1;
    materialCompMat_(1,1) = 1.0 / e2;
    materialCompMat_(2,2) = 1.0 / e3;

    materialCompMat_(0,1) = materialCompMat_(1,0) = -nu12 / e1;
    materialCompMat_(0,2) = materialCompMat_(2,0) = -nu31 / e1;
    materialCompMat_(1,2) = materialCompMat_(2,1) = -nu23 / e2;

    materialStiffMat_ = materialCompMat_;

    invert( materialStiffMat_ );

    materialStiffMat_(ALL,  2) = 0.;
    materialStiffMat_(  2,ALL) = 0.;
    materialStiffMat_(  2,  2)   = g12;
  }
  else if ( state_ == "PLANE_STRESS" )
  {
    materialCompMat_(0,0) = 1.0 / e1;
    materialCompMat_(1,1) = 1.0 / e2;
    materialCompMat_(2,2) = 1.0 / g12;
    materialCompMat_(0,1) = materialCompMat_(1,0) = -nu12 / e1;
    materialCompMat_(0,2) = materialCompMat_(2,0) = 0.0;
    materialCompMat_(1,2) = materialCompMat_(2,1) = 0.0;

    materialStiffMat_ = materialCompMat_;

    invert( materialStiffMat_ );
  }
  else
  {
    throw Error ( JEM_FUNC, "unexpected rank: " + String ( rank_ ) );
  }
 
  // System::out() << "computed orthotropic stiffness matrix: " << endl <<
    // materialStiffMat_ << endl;
  
  stiffMat_ = mChain.matmul ( transformMatInv_, materialStiffMat_, tt_);

  // System::out() << "transformed orthotropic stiffness matrix: " << endl <<
    // stiffMat_ << endl;

}



//-----------------------------------------------------------------------
//   giveStrainVector
//-----------------------------------------------------------------------

Vector OrthotropicMaterial:: giveStrainVector

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
  else if ( state_ == "PLANE_STRAIN" )
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
    ret[3] = - poisson12_ / ( 1.0 - poisson12_ ) * ( epsXX + epsYY );

    return ret;
  }

}


//-----------------------------------------------------------------------
//   giveStressVector
//-----------------------------------------------------------------------

Vector OrthotropicMaterial:: giveStressVector

    ( const Vector& stress ) const
   
{
  // Plane strain, sigma_zz put in last position

  if ( state_ == "PLANE_STRAIN" )
  {
    Vector ret ( 4 );

    double sigmaXX = stress[0];
    double sigmaYY = stress[1];

    ret[0]         = sigmaXX;
    ret[1]         = sigmaYY;
    ret[2]         = stress[2];
    ret[3]         = poisson12_ * ( sigmaXX + sigmaYY );

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

double   OrthotropicMaterial:: getFirstInvariant 

    ( const Vector&         v )          const 

{
 
  double ret;

  if      ( state_ ==  "PLANE_STRAIN" )
  {
    ret = v[0] + v[1];
  }

  else if ( state_ ==  "PLANE_STRESS" )
  {
    ret = (1.0 - 2.0 * poisson12_) / (1.0 - poisson12_)  * ( v[0] + v[1] );
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

 double   OrthotropicMaterial::getSecondInvariant 

    ( const Vector&         v )          const
{

  double I2;
 
  if      ( state_ == "PLANE_STRAIN" )
  {
    I2 = v[0] * v[1] - 0.25 * v[2] * v[2] ;
  }
  else if ( state_ == "PLANE_STRESS" )
  {
    I2 = v[0] * v[1] +  (1.0 - 2.0 * poisson12_) / (1.0 - poisson12_) * (v[0] + v[1]) 
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

double  OrthotropicMaterial::getThirdInvariant 

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

double    OrthotropicMaterial::get2ndInvDevStrain 

    ( const Vector&         v )          const
{
  double J2;

  if ( state_ ==  "PLANE_STRAIN" )
  {
    J2 = one_third * ( v[0] * v[0] + v[1] * v[1] - v[0] * v[1])
                + 0.25 * v[2] * v[2];
  }

  if ( state_ ==  "PLANE_STRESS" )
  {
    double fac = poisson12_ / ( 1.0 - poisson12_ );

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

Vector  OrthotropicMaterial:: getDI1DStrain

  ( const Vector& strain )                      const
  
{
  const int strCount = strain.size();

  Vector I ( strCount );

  if ( state_ ==  "PLANE_STRAIN" )
  {
    I[0] = 1.0;
    I[1] = 1.0;
    I[2] = 0.0;    
  }

  if ( state_ ==  "PLANE_STRESS" )
  {
    double f = (1.0 - 2.0 * poisson12_) / (1.0 - poisson12_);

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

Vector  OrthotropicMaterial:: getDI2DStrain

  ( const Vector&      v )                 const
  
{
  const int strCount = v.size();

  Vector    ret ( strCount );

  if      ( state_ ==  "PLANE_STRAIN" )
  {
    ret[0] =   v[1];
    ret[1] =   v[0];
    
    ret[2] = - v[2];  
  }
  
  else if ( state_ ==  "PLANE_STRESS" )
  {
    double f = ( 2.0 * poisson12_ ) / ( 1.0 - poisson12_ );

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

Vector  OrthotropicMaterial::getDI3DStrain

  ( const Vector&      v )               const
  
{
  const int strCount = v.size();

  Vector    ret ( strCount );

  if      ( state_ ==  "PLANE_STRAIN" )
  {
    ret = 0.0;
  }
  else if ( state_ ==  "PLANE_STRESS" )
  {
    double f = poisson12_ / ( 1.0 - poisson12_ );

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

void  OrthotropicMaterial::getIandDIDStrain

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

// 2d strain = { e_xx, e_yy, 2*e_xy}^T
// 3d strain = { e_xx, e_yy, e_zz, 2*e_xy, 2*e_yz, 2*e_zx}^T

void   OrthotropicMaterial::getI1J2andGrads   

  ( const Vector&       v,
    double& I1,
    double& J2, 
    Vector& dI1, Vector& dJ2 )  const

{
  const int    strCount  = v.size();

  if      ( state_ ==  "PLANE_STRAIN" )
  {

    const double exx = v[0];
    const double eyy = v[1];
    const double exy = v[2];

    I1 = exx + eyy;

    J2 = one_third * ( exx * exx + eyy * eyy - exx * eyy )
                       + 0.25 * exy * exy;

    dI1[0] = 1.0;
    dI1[1] = 1.0;
    dI1[2] = 0.0;    

    dJ2[0] = one_third * ( 2.0 * exx - eyy );
    dJ2[1] = one_third * ( 2.0 * eyy - exx );
    dJ2[2] = 0.5 * exy ;

  }

  else if ( state_ ==  "PLANE_STRESS" )
  {

    const double exx = v[0];
    const double eyy = v[1];
    const double exy = v[2];


    double a = poisson12_ / ( 1.0 - poisson12_ );
    double c = 1.0 - a;;
    double b = a * a + a;

    I1 = c  * ( exx + eyy );
    J2 = one_third * ( exx * exx + eyy * eyy + 
                       b   * ( exx + eyy ) * ( exx + eyy ) - 
                       exx * eyy ) + 0.25 * exy * exy;

    dI1[0] = c;
    dI1[1] = c;
    dI1[2] = 0.0;   



    dJ2[0] = one_third * ( 2.0 * ( exx + b * ( exx + eyy ) ) - eyy );
    dJ2[1] = one_third * ( 2.0 * ( eyy + b * ( eyy + exx ) ) - exx );
    dJ2[2] = 0.5 * exy ;
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
                       exy * exy + eyz * eyz + ezx * ezx;

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

void   OrthotropicMaterial::getI1andJ2

  ( const Vector&       v,
    double& I1,
    double& J2 )  const

{
   const int    strCount  = v.size();

   const double one_third = 0.33333334;

  if      ( state_ ==  "PLANE_STRAIN" )
  {

    const double exx = v[0];
    const double eyy = v[1];
    const double exy = v[2];

    I1 = exx + eyy;

    J2 = one_third * ( exx * exx + eyy * eyy - exx * eyy )
                       + 0.25 * exy * exy;

  }

  else if ( state_ ==  "PLANE_STRESS" )
  {

    const double exx = v[0];
    const double eyy = v[1];
    const double exy = v[2];


    double a = poisson12_ / ( 1.0 - poisson12_ );
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

//-----------------------------------------------------------------------
//   getPrincipalStrains
//-----------------------------------------------------------------------

Vector  OrthotropicMaterial::getPrincipalStrains

  ( const Vector&       v )  const
  
{

  const int strCount = v.size();

  Vector  prinstr(3);          
  

  // PLANE PROBLEM

  if ( strCount == 3 ) 
  {

    double exx = v[0];
    double eyy = v[1];
    double exy = 0.5 * v[2];

    double poi = poisson12_;

    // principal strains are root of a quadratic equation with det = d

    double d   = ( exx - eyy ) * ( exx - eyy ) + 4.0 * exy * exy;

    if ( d < 0 )
    {
      using namespace jem;
      throw Error (
                   JEM_FUNC,
                   " imaginary principal strains !!! "
                   );  
    }

    d = ::sqrt ( d );

    prinstr[0] = 0.5 * ( exx + eyy + d );
    prinstr[1] = 0.5 * ( exx + eyy - d );

    if ( state_ == "PLANE_STRESS" )
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

Vector  OrthotropicMaterial::getPrincipalValues

  ( const Vector&       v )  const
  
{
  const int strCount = v.size();

  JEM_PRECHECK ( strCount == 1 ||
                 strCount == 3 ||
                 strCount == 4 ||
                 strCount == 6 );
  
  Vector ans;

  /*  int nonZeroFlag = 0;

  // Check if vector v is so small, just return zero

  for ( int i = 0; i < strCount; i++)
  {
    if ( ::fabs ( v[i] ) > 1.e-20 )
    {
      nonZeroFlag = 1;
    }
  }

  if ( nonZeroFlag == 0 )
  {
    const int ss = ( strCount == 3 ) ? 2 : 3; 
       
    ans.resize(ss);
    ans = 0.0;

    return ans;
    } */

  // else compute ...

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

    ans.resize(2);
    
    ans[0] = 0.5 * ( I1 + sqrt( D ) );
    ans[1] = 0.5 * ( I1 - sqrt( D ) );

    System::out() << ans[0] << " " << ans[1] <<"\n\n";

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

    ans.resize(3);

    ans[0] = 0.5 * ( I1 + sqrt( D ) );
    ans[1] = 0.5 * ( I1 - sqrt( D ) );
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

void    OrthotropicMaterial:: getPrincipalValues

  ( Vector&             I,
    Vector&             vI,
    const Vector&       v ) const
  
{
  const int strCount = v.size();

  JEM_PRECHECK ( strCount == 1 ||
                 strCount == 3 ||
                 strCount == 4 ||
                 strCount == 6 );

  /*int nonZeroFlag = 0;

  // Check if vector v is so small, just return zero

  for ( int i = 0; i < strCount; i++)
  {
    if ( ::fabs ( v[i] ) > 1.e-20 )
    {
      nonZeroFlag = 1;
    }
  }

  if ( nonZeroFlag == 0 )
  {
    const int ss = ( strCount == 3 ) ? 2 : 3; 
       
    vI.resize(ss);
    I. resize(ss);

    vI = 0.0;
    I  = 0.0;

    return;
  }
  */

  // else compute ...

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

void    OrthotropicMaterial::getPrincipalStressAndDirections

    ( const Vector&      stress,
      const Vector&      sigmaI,
      const Matrix&      dir )               const
{
  const int dim = stress.size ();

  double sigmaX, sigmaY, sigmaXY;

  double sigmaZ  = 0.0;
  double sigmaYZ = 0.0;
  double sigmaZX = 0.0;

  if    ( dim == 3 )
  {
     sigmaX  = stress[0];
     sigmaY  = stress[1];
     sigmaXY = stress[2];

     if  ( state_ == "PLANE_STRAIN" )
     {
       sigmaZ = poisson12_ * ( sigmaX + sigmaY );
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

     Matrix sigma(3,3);

     sigma(0,0) = sigmaX;
     sigma(0,1) = sigmaXY;
     sigma(0,2) = sigmaZX;

     sigma(1,1) = sigmaY;
     sigma(1,2) = sigmaYZ;

     sigma(2,2) = sigmaZ;

     using jem::numeric::EigenUtils;

     EigenUtils::symSolve ( sigmaI, dir, sigma);
  }

}


//-----------------------------------------------------------------------
//   getMaximumPrincipalStress
//-----------------------------------------------------------------------

double    OrthotropicMaterial::getMaximumPrincipalStress

    ( const Vector&      stress ) const
{
  const int dim = stress.size ();

  double sigmaX, sigmaY, sigmaXY;

  double sigmaZ  = 0.0;
  double sigmaYZ = 0.0;
  double sigmaZX = 0.0;

  if    ( dim == 3 )
  {
     sigmaX  = stress[0];
     sigmaY  = stress[1];
     sigmaXY = stress[2];

     if  ( state_ == "PLANE_STRAIN" )
     {
       sigmaZ = poisson12_ * ( sigmaX + sigmaY );
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

  using jem::numeric::EigenUtils;

  Vector sigmaI (3);
  EigenUtils::symSolve ( sigmaI, sigma );

  //System::out() << "local principal stresses :" << sigmaI <<"\n\n";

  // find the maximum principal stress

  double max = sigmaI[0];

  max        = max > sigmaI[1] ? max : sigmaI[1];
  max        = max > sigmaI[2] ? max : sigmaI[2];

  return max;
}

//-----------------------------------------------------------------------
//   getPrincipalStress
//-----------------------------------------------------------------------

Vector    OrthotropicMaterial::getPrincipalStress

    ( const Vector&      stress ) const
{
  const int dim = stress.size ();

  Vector sigmaI ( dim == 3 ? 2 : 3 );

  // two dimensions

  if    ( dim == 3 )
  {
     double sigmaX  = stress[0];
     double sigmaY  = stress[1];
     double sigmaXY = stress[2];

     double s1      = 0.5 * ( sigmaX + sigmaY );
     double s2      = 0.5 * ( sigmaX - sigmaY );
     double d0      = ::sqrt( s2 * s2 + sigmaXY * sigmaXY );

     sigmaI[0]      =  s1 + d0;
     sigmaI[1]      =  s1 - d0; 
  }

  // three dimensions

  else
  {
     double sigmaX  = stress[0];
     double sigmaY  = stress[1];
     double sigmaZ  = stress[2];

     double sigmaXY = stress[3];
     double sigmaYZ = stress[4];
     double sigmaZX = stress[5];

     Matrix sigma(3,3);

     sigma(0,0) = sigmaX;
     sigma(0,1) = sigmaXY;
     sigma(0,2) = sigmaZX;

     sigma(1,1) = sigmaY;
     sigma(1,2) = sigmaYZ;

     sigma(2,2) = sigmaZ;

     using jem::numeric::EigenUtils;
  
     EigenUtils::symSolve ( sigmaI, sigma );

     sort ( sigmaI ); 

     double temp = sigmaI[0];

     sigmaI[0]   = sigmaI[2];
     sigmaI[2]   = temp;

  }

  return sigmaI;
}

//-----------------------------------------------------------------------
//   getMaterialStress
//-----------------------------------------------------------------------

Vector    OrthotropicMaterial::getMaterialStress

    ( const Vector&      stress ) const
{

  const int dim = stress.size ();

  Vector   matStress ( dim );

  matmul ( matStress, transformMat_ , stress );

  return matStress;
}

//-----------------------------------------------------------------------
//   getMaterialStrain
//-----------------------------------------------------------------------

Vector    OrthotropicMaterial::getMaterialStrain

    ( const Vector&      strain ) const

{
  const int dim = strain.size ();

  Vector    matStrain ( dim );

  matmul ( matStrain, tt_, strain );

  return matStrain;
}

//-----------------------------------------------------------------------
//   getMaterialStiffMat 
//-----------------------------------------------------------------------

Matrix    OrthotropicMaterial::getMaterialStiffMat () const

{
  return materialStiffMat_;
}

//-----------------------------------------------------------------------
//   allocPoints 
//-----------------------------------------------------------------------

void     OrthotropicMaterial::allocPoints ( const int count ) 
{
  pCount_ += count;
}

void     OrthotropicMaterial::allocPoints

  ( const int        count,
    const Matrix&    transfer,
    const IntVector& oldPoints )

{
  allocPoints ( count );
}

