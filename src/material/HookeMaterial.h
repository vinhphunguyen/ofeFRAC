/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the isotropic linear elastic material
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *
 */

#ifndef HOOKEMATERIAL_H
#define HOOKEMATERIAL_H

#include <jem/base/String.h>

#include "Material.h"
#include "util/utilities.h"


using jem::String;

/*
enum ProblemType {
  PlaneStrain,
  PlaneStress,
  AxiSymmetric
};
*/

// =======================================================
//  class HookeMaterial
// =======================================================

// This class implements an isotropic elastic material

class HookeMaterial : public Material
{
 public:

  static const char*      YOUNG_PROP;
  static const char*      POISSON_PROP;
  static const char*      RHO_PROP;
  static const char*      STATE_PROP;

  explicit                HookeMaterial

    ( int                   rank,
      const Properties&     globdat );

  virtual void            configure

    ( const Properties&     props,
      const Properties&     globdat );

  virtual void            getConfig

    ( const Properties&     conf,
      const Properties&     globdat ) const;

  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      int                   ipoint );

  virtual double          giveHistory     ( int ipoint ) const;

  Matrix                  getStiffMat     () const;

  inline double           giveYoung       () const;
  inline double           givePoisson     () const;
  inline double           giveRho         () const;
  inline int              giveDim         () const;
  inline ProblemType      giveState       () const;

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
  //   compute the I1 and J2 of stress vector
  //-----------------------------------------------------------------------

  void                    getI1andJ2ForStress  

  ( const Vector&       s,
    double& I1,
    double& J2 )                              const;

  
  //-----------------------------------------------------------------------
  //   compute the I1 and J2 at the same time of strain tensor
  //   used in the modified Von Mises equivalent strain
  //-----------------------------------------------------------------------

  void                    getI1andJ2ForStress  

  ( const Vector&       v,
    const Vector&       S,
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

  void                  getPrincipalValues

    ( Vector&             I,
      Vector&             vI,
      const Vector&       v )                const;


  void                  getPrincipalStressAndDirections

    ( const Vector&      stress,
      const Vector&      sigmaI,
      const Matrix&      dir )               const;
  
  double                checkLocalisation 

    ( Vector&              normal,
      const Vector&        stress,
      const Matrix&        tangent,
            int            ipoint )          const;

 protected:

  virtual                ~HookeMaterial   ();

 protected:

  double                  young_;
  double                  poisson_;
  Matrix                  stiffMat_;

  double                  area_;

 private:

  void                    computeStiffMat_();

 private:

  ProblemType             state_;
};

//-----------------------------------------------------------------------
//   giveYoung
//-----------------------------------------------------------------------

inline double HookeMaterial::giveYoung  () const
{
  return young_;
}

//-----------------------------------------------------------------------
//   givePoisson
//-----------------------------------------------------------------------

inline double HookeMaterial::givePoisson() const
{
  return poisson_;
}

//-----------------------------------------------------------------------
//   giveRho
//-----------------------------------------------------------------------

inline double HookeMaterial::giveRho() const
{
  return Material::rho_;
}

//-----------------------------------------------------------------------
//   giveDim
//-----------------------------------------------------------------------

inline int    HookeMaterial::giveDim() const
{
  return rank_;
}

//-----------------------------------------------------------------------
//   giveState
//-----------------------------------------------------------------------

inline ProblemType   HookeMaterial::giveState() const
{
  return state_;
}

#endif 
