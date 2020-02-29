#ifndef ORTHOTROPICMATERIAL_H
#define ORTHOTROPICMATERIAL_H

#include <jem/base/String.h>
#include <jem/base/System.h>

#include "Material.h"


using jem::String;
using jem::System;
using jem::io::endl;
using jive::IntVector;


// =======================================================
//  class OrthotropicMaterial
// =======================================================

// This class implements an orthotropic or  transversely isotropic elastic material
// Frans van der Meer, February 2008
//
// for transversely isotropic material, input:
//
//   - young1
//   - young2
//   - poisson12
//   - poisson23
//   - shear12
//   - theta
//
// for plane stress, input:
//
//   - young1
//   - young2
//   - poisson12
//   - shear12
//   - theta

class OrthotropicMaterial : public Material
{
 public:

  static const char*      YOUNG_1_PROP;
  static const char*      YOUNG_2_PROP;
  static const char*      YOUNG_3_PROP;
  static const char*      POISSON_12_PROP;
  static const char*      POISSON_23_PROP;
  static const char*      POISSON_31_PROP;
  static const char*      SHEAR_12_PROP;
  static const char*      SHEAR_23_PROP;
  static const char*      SHEAR_31_PROP;
  static const char*      RHO_PROP;
  static const char*      STATE_PROP;
  static const char*      THETA_PROP;

  explicit                OrthotropicMaterial

    ( const int             rank,
      const Properties&     globdat );

                         ~OrthotropicMaterial ();

  virtual void            configure

    ( const Properties&     props,
      const Properties&     globdat );

  virtual void            getConfig

    ( const Properties&     conf,
      const Properties&     globdat )         const;

  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      int                   ipoint );

  virtual double          giveHistory       ( int ipoint ) const;

  Vector                  giveStressAtPoint 

    ( const int             ipoint,
      const Vector&         strain ) const;

  virtual void            allocPoints     ( const int count );

  virtual void            allocPoints

    ( const int        count,
      const Matrix&    transfer,
      const IntVector& oldPoints );

  inline virtual int      pointCount      () const;

  Matrix                  getStiffMat     () const;

  inline double           giveYoung1      () const;
  inline double           giveYoung2      () const;
  inline double           giveYoung3      () const;
  inline double           givePoisson12   () const;
  inline double           givePoisson23   () const;
  inline double           givePoisson31   () const;
  inline double           giveShear12     () const;
  inline double           giveShear23     () const;
  inline double           giveShear31     () const;
  inline double           giveRho         () const;
  inline double           giveTheta       () const;
  inline int              giveDim         () const;
  inline String           giveState       () const;

  // Given the engineering strain vector (compute from FE)
  // compute the normal strain vector (divide by two shear
  // strains) including add eps_zz
  // in case of plane stress problem
  
  Vector                  giveStrainVector

    ( const Vector&         strain )          const;


  // Given the stress vector computed from FE, build the
  // new one, need in case of plane strain problem (sigma_zz != 0)
  
  Vector                  giveStressVector

    ( const Vector&         strain )          const;


  //-----------------------------------------------------------------------
  //   compute the invariants of the strain tensor I1 I2 I3
  //-----------------------------------------------------------------------

  double                  getFirstInvariant 

    ( const Vector&         strain )          const;

  double                  getSecondInvariant 

    ( const Vector&         strain )          const;

  double                  getThirdInvariant 

    ( const Vector&         strain )          const;

  //-----------------------------------------------------------------------
  //   compute the second invariant of the deviatoric strain tensor J2
  //-----------------------------------------------------------------------

  double                  get2ndInvDevStrain 

    ( const Vector&         strain )          const;

  //-----------------------------------------------------------------------
  //   compute the derivatives w.r.t strain tensor  of
  //   invariants of the strain tensor
  //-----------------------------------------------------------------------

  Vector                  getDI1DStrain

    ( const Vector&         strain )          const;

  Vector                  getDI2DStrain

    ( const Vector&         strain )          const;

  Vector                  getDI3DStrain

    ( const Vector&         strain )          const;

  //-----------------------------------------------------------------------
  //   compute the derivatives w.r.t strain tensor of the 
  //   second invariant of the deviatoric strain tensor
  //-----------------------------------------------------------------------

  Vector                  getDJ2DStrain

    ( const Vector&         strain )          const;

  //-----------------------------------------------------------------------
  //   compute the I1, J2 and its derivatives, at the same time, of 
  //   strain tensor
  //   Usage: in the modified Von Mises equivalent strain
  //-----------------------------------------------------------------------

  void                    getI1J2andGrads   

  ( const Vector&       v,
    double& I1,
    double& J2, 
    Vector& dI1, Vector& dJ2 )                const;

  //-----------------------------------------------------------------------
  //   compute the I1 and J2 at the same time of strain tensor
  //   used in the modified Von Mises equivalent strain
  //-----------------------------------------------------------------------

  void                    getI1andJ2  

  ( const Vector&       v,
    double& I1,
    double& J2 )                              const;

  //-----------------------------------------------------------------------
  //   compute the I1, I2, dI1, dI2, dI3 and principal strains
  //-----------------------------------------------------------------------

  void                    getIandDIDStrain

  ( const Vector& v,
          Vector& I,   Vector& princ,
          Vector& dI1, Vector& dI2, Vector& dI3 )  const;


  //-----------------------------------------------------------------------
  //   getPrincipalStrains
  //-----------------------------------------------------------------------

  Vector                  getPrincipalStrains

    ( const Vector&       v )                 const;

  //-----------------------------------------------------------------------
  //   getPrincipalValues
  //-----------------------------------------------------------------------

  Vector                  getPrincipalValues

    ( const Vector&       v )                 const;


  //-----------------------------------------------------------------------
  //   getPrincipalValues + return invariants
  //-----------------------------------------------------------------------

  void                    getPrincipalValues

    ( Vector&             I,
      Vector&             vI,
      const Vector&       v )                const;

  //-----------------------------------------------------------------------
  //   compute principal stresses and directions
  //-----------------------------------------------------------------------

  void                    getPrincipalStressAndDirections

    ( const Vector&      stress,
      const Vector&      sigmaI,
      const Matrix&      dir )               const;

  //-----------------------------------------------------------------------
  //   compute maximum principal stress
  //-----------------------------------------------------------------------

  double                  getMaximumPrincipalStress

    ( const Vector&      stress )            const;

  //-----------------------------------------------------------------------
  //   compute principal stress and return the sorted vector
  //-----------------------------------------------------------------------

  Vector                  getPrincipalStress

    ( const Vector&      stress )            const;

  //-----------------------------------------------------------------------
  //   compute stress in material axes
  //-----------------------------------------------------------------------

  Vector                  getMaterialStress

    ( const Vector&      stress )            const;

  //-----------------------------------------------------------------------
  //   compute strain in material axes
  //-----------------------------------------------------------------------

  Vector                  getMaterialStrain

    ( const Vector&      strain )            const;

  //-----------------------------------------------------------------------
  //   get stiffness matrix in material frame
  //-----------------------------------------------------------------------
  
  Matrix                  getMaterialStiffMat     () const;

 protected:

  double                  young1_;
  double                  young2_;
  double                  young3_;
  double                  poisson12_;
  double                  poisson23_;
  double                  poisson31_;
  double                  shear12_;
  double                  shear23_;
  double                  shear31_;
  double                  rho_;
  double                  theta_;
  Matrix                  stiffMat_;
  Matrix                  materialStiffMat_;
  Matrix                  materialCompMat_;
  Matrix                  transformMat_;
  Matrix                  transformMatInv_;
  Matrix                  tt_;
  String                  state_;

  double                  area_;
  
 private:

  void                    computeStiffMat_();
  void                    computeTransformMats_();
  
 private:

  int                     pCount_;

};

//-----------------------------------------------------------------------
//   pointCount
//-----------------------------------------------------------------------

inline int    OrthotropicMaterial::pointCount  () const

{
  return pCount_;
}

//-----------------------------------------------------------------------
//   giveYoung1
//-----------------------------------------------------------------------

inline double OrthotropicMaterial::giveYoung1  () const
{
  return young1_;
}

//-----------------------------------------------------------------------
//   giveYoung2
//-----------------------------------------------------------------------

inline double OrthotropicMaterial::giveYoung2  () const
{
  return young2_;
}

//-----------------------------------------------------------------------
//   giveYoung3
//-----------------------------------------------------------------------

inline double OrthotropicMaterial::giveYoung3  () const
{
  return young3_;
}

//-----------------------------------------------------------------------
//   givePoisson12
//-----------------------------------------------------------------------

inline double OrthotropicMaterial::givePoisson12() const
{
  return poisson12_;
}

//-----------------------------------------------------------------------
//   givePoisson23
//-----------------------------------------------------------------------

inline double OrthotropicMaterial::givePoisson23() const
{
  return poisson23_;
}

//-----------------------------------------------------------------------
//   givePoisson31
//-----------------------------------------------------------------------

inline double OrthotropicMaterial::givePoisson31() const
{
  return poisson31_;
}

//-----------------------------------------------------------------------
//   giveShear12
//-----------------------------------------------------------------------

inline double OrthotropicMaterial::giveShear12() const
{
  return shear12_;
}

//-----------------------------------------------------------------------
//   giveShear23
//-----------------------------------------------------------------------

inline double OrthotropicMaterial::giveShear23() const
{
  return shear23_;
}

//-----------------------------------------------------------------------
//   giveShear31
//-----------------------------------------------------------------------

inline double OrthotropicMaterial::giveShear31() const
{
  return shear31_;
}

//-----------------------------------------------------------------------
//   giveTheta
//-----------------------------------------------------------------------

inline double OrthotropicMaterial::giveTheta() const
{
  return theta_;
}

//-----------------------------------------------------------------------
//   giveRho
//-----------------------------------------------------------------------

inline double OrthotropicMaterial::giveRho() const
{
  return rho_;
}

//-----------------------------------------------------------------------
//   giveDim
//-----------------------------------------------------------------------

inline int    OrthotropicMaterial::giveDim() const
{
  return rank_;
}

//-----------------------------------------------------------------------
//   giveState
//-----------------------------------------------------------------------

inline String    OrthotropicMaterial::giveState() const
{
  return state_;
}

#endif 
