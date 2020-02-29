#ifndef LARGEDISPUTILS_H 
#define LARGEDISPUTILS_H

#include <jive/Array.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>

#include "utilities.h"
#include "XtraUtilities.h"

using jive::Cubix ;
using jive::Vector;
using jive::Matrix;
using jem::numeric::matmul;
using jem::numeric::MatmulChain;

typedef MatmulChain<double,1>  MChain1;
typedef MatmulChain<double,2>  MChain2;
typedef MatmulChain<double,3>  MChain3;
typedef MatmulChain<double,4>  MChain4;

/*
 * Some utility functions for large displacement formulations.
 *
 * Authors:
 *
 * Frans van der Meer, TU Delft, The Netherlands
 * Vinh Phu Nguyen, Cardiff University, Wales, UK
 * TKM: 
 * 30-Aug-2019: modified B-matrix B0 = B * F computed for Total Lagrangian
 * 12-Sep-2019: "evalDeformationGradient" corrected
 */

//-----------------------------------------------------------------------
//   typedefs
//-----------------------------------------------------------------------

typedef void        (*BMatrixLinFunc)

  (       Matrix&       b,
    const Matrix&       f,
    const Matrix&       g );

typedef void        (*BMatrixNonlinFunc)

  (       Matrix&       b,
    const Matrix&       g );

// A pointer to a function that computes the spatial derivatives of
// the interpolation matrix. This is the so-called B-matrix.
// it points to the corresponding function for 1D, 2D and 3D case.

typedef void        (*ShapeGradsTLFunc)

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f );

//-----------------------------------------------------------------------
//   public functions
//-----------------------------------------------------------------------


BMatrixLinFunc        getBMatrixLinFunc

  ( int                 rank );

BMatrixNonlinFunc     getBMatrixNonlinFunc

  ( int                 rank );

//-----------------------------------------------------------------------
// Compute deformation gradient f 
//   from nodal displacements u and shape function gradients g
//-----------------------------------------------------------------------

void   evalDeformationGradient

  (       Matrix&    f,
    const Vector&    u,
    const Matrix&    g );

//-----------------------------------------------------------------------
// Compute Green-LagrangeStrain vector from deformation gradient
//-----------------------------------------------------------------------

void   getGreenLagrangeStrain
  
  (       Vector&    eps,
    const Matrix&    f );


//-----------------------------------------------------------------------
// Rotate normal vector n with deformation gradient F
//-----------------------------------------------------------------------

Vector rotateNormalVector

  ( const Vector&    n0,
    const Matrix&    f );

//-----------------------------------------------------------------------
// Update normal matrix 
//-----------------------------------------------------------------------

void   updateNormalMatrix

  (       Matrix&    nMat,
    const Vector&    nVec );

//-----------------------------------------------------------------------
// Update rotation matrix 
//-----------------------------------------------------------------------

void   updateRotationMatrix

  (       Matrix&    Q,
    const Vector&    nVec );

//-----------------------------------------------------------------------
// Rotate normal vector n0 with deformatin gradient f to compute
//   updated rotation matrix Q and derivative of angle phi w.r.t. f
//-----------------------------------------------------------------------

void  updateMatrices

  (       Matrix&  Q,
          Matrix&  dPhidF,
    const Vector&  n0,
    const Matrix&  f );

//-----------------------------------------------------------------------
// Compute 2D B matrix for linear part of K
//   from deformation gradient f and shape function gradients g
//-----------------------------------------------------------------------

void  getBMatrixLin2D

  (       Matrix&  b,
    const Matrix&  f,
    const Matrix&  g );

//-----------------------------------------------------------------------
// Compute 3D B matrix for linear part of K
//   from deformation gradient f and shape function gradients g
//-----------------------------------------------------------------------

void  getBMatrixLin3D

  (       Matrix&  b,
    const Matrix&  f,
    const Matrix&  g );

//-----------------------------------------------------------------------
// Assemble B matrix for nonlinear part of K
//   from shape function gradients g
// (more efficient to skip this and use addElemMatLargeD instead for 
//   direct assembly of additional term in K)
//-----------------------------------------------------------------------

void  getBMatrixNonlin2D

  (       Matrix&  b,
    const Matrix&  g );


void  getBMatrixNonlin3D

  (       Matrix&  b,
    const Matrix&  g );

//-----------------------------------------------------------------------
// Add part nonlinear in the displacement field to stiffness matrix k
//   from second piola kirchhoff stress tau and shape func gradients g
//-----------------------------------------------------------------------

void  addElemMatLargeD

  (       Matrix&  k,
    const Vector&  tau,
    const Matrix&  g,
    const double   w );

//-----------------------------------------------------------------------
// Compute nominal stress: stre=F*pk2
//-----------------------------------------------------------------------

void   getNominalStress

  (       Matrix& stre,
    const Vector& pk2,
    const Matrix& f );

//----------------------------------------------------------------------
// Premultiply stiffness term DB with F, while going from Voigt notation
//   to tensor notation (because F*pk2 is asymmetric)
//----------------------------------------------------------------------

void   getNominalStiff

  (       Cubix & FDB,
    const Matrix& DB,
    const Matrix& f );

//-----------------------------------------------------------------------
// Convert stress/strain from/to Voigt notation from/to tensor
//-----------------------------------------------------------------------

Matrix voigtToTensorStress ( const Vector& vec );

void   tensorToVoigtStrain 

  (       Vector& eps, 
    const Matrix& tens );


void   tensorToVoigtStress

(       Vector& s, 
    const Matrix& tens );


//-----------------------------------------------------------------------
// Convert any rectangular matrix to Column vector
//-----------------------------------------------------------------------

void getMatrixtoVector 
	( 	Vector& fvec,
	  const Matrix& tens );


//-----------------------------------------------------------------------
// Functions to get modified gradient matrix B0 = B * F
//-----------------------------------------------------------------------

// These functions compute the B-matrix given the matrix of shape function 1st 
// gradients (in global coordinates). For 1D, 2D (plane stress/strain) and 3D problems.

void                  get1DShapeGradsTL

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f );

void                  get2DShapeGradsTL

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f );

void                  get3DShapeGradsTL

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f );

// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.

ShapeGradsTLFunc        getShapeGradsTLFunc

  ( int                 rank );


//-----------------------------------------------------------------------
// evaluate geometric stiffness matrix
//-----------------------------------------------------------------------

void evalGeometricStiffness 
  (   Matrix& tens,
    const Matrix& grads,
    const Vector& stress );

//-----------------------------------------------------------------------
// get Cauchy-Green tensor and its inverse tensors (in Voigt notation), invariances from given Lagrange strain
//-----------------------------------------------------------------------

void getCauchyTensor 
  ( const Vector& strain,
   		Vector&  C,
   		Vector&  invC,
   		Matrix&  invCtensor,
   		double& I1,
   		double& I2,
   		double& I3 );

//-----------------------------------------------------------------------
// evaluate C_ijkl from given inv(C)
//-----------------------------------------------------------------------

void getInvCtensor 
  (  const Vector&   invC,
           Matrix&   invCtensor );

//-----------------------------------------------------------------------
// get inverse Cauchy-Green tensor (in Voigt notation) from given eigenvalues/vector of "C"
//-----------------------------------------------------------------------

void getCauchyTensorEigen 

  ( const Vector&  evals,
    const Matrix&  evecs,
    	Vector&  invCp,
    	Vector&  invCn,
      	Matrix&  invCtensorP,
      	Matrix&  invCtensorN );

//-----------------------------------------------------------------------
// get Beta tensors:beta_n_ij and beta_n_ijkl (in Voigt notation) from given eigenvector of "C"
//-----------------------------------------------------------------------

void getCauchyTensorEigen 

  ( const 	Matrix&  evecs,
      		Matrix&  Beta,
      		Cubix&  BetaTensor );

//-----------------------------------------------------------------------
// compute Cauchy stress from given pk2 stress
//-----------------------------------------------------------------------

void getCauchyStress 
  (       Vector&   cstress,
    const Vector&   pk2,
    const Matrix&   f );

//-----------------------------------------------------------------------
// compute determinant of a Matrix
//-----------------------------------------------------------------------

void getDeterminant 
  (       double&  J,
    const Matrix&  f);


//-----------------------------------------------------------------------
// update stiff matrix after local iteration for plane stress state
//-----------------------------------------------------------------------

void updateStiffPlaneStress 
  (       Matrix&   stiffS );

#endif
