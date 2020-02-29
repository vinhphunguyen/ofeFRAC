#include <jem/base/System.h>
#include <jem/base/Array.h>
#include <jem/base/array/tensor.h>
#include <jem/numeric/algebra/utilities.h>

#include "LargeDispUtilities.h"

// using namespace jem;

using jem::System;
using jem::ALL;
using jem::END;
using jem::TensorIndex;
using jem::numeric::norm2;
using jem::io::endl;

//-----------------------------------------------------------------------
//    evalDeformationGradient (2D and 3D)
//-----------------------------------------------------------------------

void  evalDeformationGradient

  (       Matrix&  f,
    const Vector&  u,
    const Matrix&  g )
{
  int rank = g.size(0);

  for ( int i = 0; i < rank; ++i )
  {
    for ( int j = 0; j < rank; ++j )
    {
      f(i,j) = dot ( g(j,ALL) , u[slice(i,END,rank)] );
    }
    f(i,i) += 1.;
  }
}
  

//-----------------------------------------------------------------------
//    getGreenLagrangeStrain (2D and 3D)
//-----------------------------------------------------------------------

void  getGreenLagrangeStrain

  (       Vector&  eps,
    const Matrix&  f )

{
  int rank = f.size(0); eps = 0.0;

  // Compute Green-Lagrange strain tensor: 1/2 * ( F^T*F - I )

  TensorIndex i,j,k;
  Matrix      tens( rank, rank );

  tens(i,j) = 0.5 * ( dot( f(k,i), f(k,j), k ) - where(i==j,1.,0.) );

  // Convert to vector (Voigt notation)
  
  tensorToVoigtStrain( eps, tens );
}

//-----------------------------------------------------------------------
//    rotateNormalVector
//-----------------------------------------------------------------------

Vector  rotateNormalVector

  ( const Vector&  n0,
    const Matrix&  f )

{
  // int rank = f.size(0);

  // Only implemented for 2D!!

  JEM_ASSERT ( f.size(0) == 2 );

  // Evaluate det(F)(F^T)^(-1)n

  Vector nUpdated(2);

  nUpdated[0] =   f(1,1)*n0[0] - f(1,0)*n0[1];
  nUpdated[1] = - f(0,1)*n0[0] + f(0,0)*n0[1] ;

  // Return normalized vector
  
  double l = norm2( nUpdated );
  
  nUpdated /= l;

  return nUpdated ;
};

//-----------------------------------------------------------------------
//    updateNormalMatrix
//-----------------------------------------------------------------------

void  updateNormalMatrix

  (       Matrix&  nMat,
    const Vector&  nVec )

{
  JEM_ASSERT ( nVec.size() == 2 );

  nMat.resize (2,3);

  nMat(0,0) = nVec[0];
  nMat(0,1) = 0.;
  nMat(0,2) = nVec[1];
  nMat(1,0) = 0.;
  nMat(1,1) = nVec[1];
  nMat(1,2) = nVec[0];
}

//-----------------------------------------------------------------------
//    updateRotationMatrix
//-----------------------------------------------------------------------

void  updateRotationMatrix

  (       Matrix&  Q,
    const Vector&  nVec )

{
  Q(0,0) =   nVec[0];
  Q(0,1) =   nVec[1];
  Q(1,0) =   nVec[1];
  Q(1,1) = - nVec[0];
}

//-----------------------------------------------------------------------
//    updateMatrices
//-----------------------------------------------------------------------

void  updateMatrices

  (       Matrix&  Q,
          Matrix&  dPhidF,
    const Vector&  n0,
    const Matrix&  f )

{
  // int rank = f.size(0);

  // Only implemented for 2D!!

  JEM_ASSERT ( f.size(0) == 2 );

  // Evaluate det(F)(F^T)^(-1)n

  double nBar0 =   f(1,1)*n0[0] - f(1,0)*n0[1];
  double nBar1 = - f(0,1)*n0[0] + f(0,0)*n0[1] ;

  // Normalize vector
  
  double l2 = nBar0 * nBar0 + nBar1 * nBar1;

  double l = sqrt( l2 );
  
  Vector n ( 2 );

  n[0] = nBar0 / l;
  n[1] = nBar1 / l;

  // Fill matrices

  // Vector n = rotateNormalVector ( n0, f );
  updateRotationMatrix ( Q, n );

  dPhidF(0,0) =   nBar0 * n0[1] / l2;
  dPhidF(0,1) = - nBar0 * n0[0] / l2;
  dPhidF(1,0) =   nBar1 * n0[1] / l2;
  dPhidF(1,1) = - nBar1 * n0[0] / l2;
}


//-----------------------------------------------------------------------
//    getBMatrixLin2D
//-----------------------------------------------------------------------

void  getBMatrixLin2D

  (       Matrix&  b,
    const Matrix&  f,
    const Matrix&  g )

{
  JEM_ASSERT   ( b.size(0) == 3 &&
                 g.size(0) == 2 &&
                 f.size(0) == 2 &&
                 f.size(1) == 2 &&
                 b.size(1) == 2 * g.size(1) );

  const int  nodeCount = g.size (1);

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    int  i0 = 2 * inode;
    int  i1 = i0 + 1;

    b(0,i0) = f(0,0) * g(0,inode);
    b(0,i1) = f(1,0) * g(0,inode);
    b(1,i0) = f(0,1) * g(1,inode);
    b(1,i1) = f(1,1) * g(1,inode);

    b(2,i0) = f(0,0) * g(1,inode) + f(0,1) * g(0,inode);
    b(2,i1) = f(1,1) * g(0,inode) + f(1,0) * g(1,inode);
  }
}

//-----------------------------------------------------------------------
//    getBMatrixLin3D
//-----------------------------------------------------------------------

void  getBMatrixLin3D

  (       Matrix&  b,
    const Matrix&  f,
    const Matrix&  g )

{
  JEM_ASSERT   ( b.size(0) == 6 &&
                 g.size(0) == 3 &&
                 f.size(0) == 3 &&
                 f.size(1) == 3 &&
                 b.size(1) == 3 * g.size(1) );

  const int  nodeCount = g.size (1);

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    int  i0 = 3 * inode;
    int  i1 = i0 + 1;
    int  i2 = i0 + 2;

    b(0,i0) = f(0,0) * g(0,inode);
    b(0,i1) = f(1,0) * g(0,inode);
    b(0,i2) = f(2,0) * g(0,inode);

    b(1,i0) = f(0,1) * g(1,inode);
    b(1,i1) = f(1,1) * g(1,inode);
    b(1,i2) = f(2,1) * g(1,inode);

    b(2,i0) = f(0,2) * g(2,inode);
    b(2,i1) = f(1,2) * g(2,inode);
    b(2,i2) = f(2,2) * g(2,inode);

    b(3,i0) = f(0,0) * g(1,inode) + f(0,1) * g(0,inode);
    b(3,i1) = f(1,0) * g(1,inode) + f(1,1) * g(0,inode);
    b(3,i2) = f(2,0) * g(1,inode) + f(2,1) * g(0,inode);

    b(4,i0) = f(0,1) * g(2,inode) + f(0,2) * g(1,inode);
    b(4,i1) = f(1,1) * g(2,inode) + f(1,2) * g(1,inode);
    b(4,i2) = f(2,1) * g(2,inode) + f(2,2) * g(1,inode);

    b(5,i0) = f(0,2) * g(0,inode) + f(0,0) * g(2,inode);
    b(5,i1) = f(1,2) * g(0,inode) + f(1,0) * g(2,inode);
    b(5,i2) = f(2,2) * g(0,inode) + f(2,0) * g(2,inode);
  }
}
 


//-----------------------------------------------------------------------
//    getBMatrixNonlin2D
//-----------------------------------------------------------------------

void  getBMatrixNonlin2D

  (       Matrix&  b,
    const Matrix&  g )

{
  JEM_ASSERT   ( b.size(0) == 4 &&
                 g.size(0) == 2 &&
                 b.size(1) == 2 * g.size(1) );

  const int  nodeCount = g.size (1);

  b = 0.0;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    int  i0 = 2 * inode;
    int  i1 = i0 + 1;

    b(0,i0) = g(0,inode);
    b(1,i0) = g(1,inode);

    b(2,i1) = g(0,inode);
    b(3,i1) = g(1,inode);
  }
}


//-----------------------------------------------------------------------
//    getBMatrixNonlin3D
//-----------------------------------------------------------------------

void  getBMatrixNonlin3D

  (       Matrix&  b,
    const Matrix&  g )

{
  JEM_ASSERT   ( b.size(0) == 9 &&
                 g.size(0) == 3 &&
                 b.size(1) == 3 * g.size(1) );

  const int  nodeCount = g.size (1);

  b = 0.0;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    int  i0 = 3 * inode;
    int  i1 = i0 + 1;
    int  i2 = i0 + 2;

    b(0,i0) = g(0,inode);
    b(1,i0) = g(1,inode);
    b(2,i0) = g(2,inode);

    b(3,i1) = g(0,inode);
    b(4,i1) = g(1,inode);
    b(5,i1) = g(2,inode);
    
    b(6,i2) = g(0,inode);
    b(7,i2) = g(1,inode);
    b(8,i2) = g(2,inode);
  }
}


//-----------------------------------------------------------------------
//   getBMatrixNonlinFunc
//-----------------------------------------------------------------------


BMatrixNonlinFunc getBMatrixNonlinFunc ( int rank )
{
  JEM_PRECHECK ( rank > 1 && rank <= 3 );

  if ( rank == 2 )
  {
    return & getBMatrixNonlin2D;
  }
  else
  {
    return & getBMatrixNonlin3D;
  }
}

//-----------------------------------------------------------------------
//   getBMatrixLinFunc
//-----------------------------------------------------------------------

BMatrixLinFunc getBMatrixLinFunc ( int rank )
{
  JEM_PRECHECK ( rank > 1 && rank <= 3 );

  if ( rank == 2 )
  {
    return & getBMatrixLin2D;
  }
  else
  {
    return & getBMatrixLin3D;
  }
}


//-----------------------------------------------------------------------
//    addElemMatLargeD
//-----------------------------------------------------------------------

void addElemMatLargeD

  (       Matrix& k,
    const Vector& tau,
    const Matrix& g,
    const double  w )

{
  int rank = g.size(0);
  int nn   = g.size(1);

  JEM_ASSERT    ( k.size(0) == rank*nn &&
                  k.size(1) == rank*nn );

  Matrix   t = voigtToTensorStress( tau );

  JEM_ASSERT    ( t.size(0) == rank &&
                  t.size(1) == rank );

  Matrix   btb ( nn, nn );

  TensorIndex in,jn,it,jt;

  btb(in,jn) = dot( g(jt,in), dot( t(it,jt), g(it,jn), it ), jt );

  for ( int in = 0; in < nn; ++in )
  {
    for ( int jn = 0; jn < nn; ++jn )
    {
      for ( int ix = 0; ix < rank; ++ix )
      {
        k( in*rank+ix, jn*rank+ix ) += w*btb(in,jn);
      }
    }
  }
}

//-----------------------------------------------------------------------
//    getNominalStress:: P = S*F' ==>? or S*F
//-----------------------------------------------------------------------

void   getNominalStress

  (       Matrix& stre,
    const Vector& pk2,
    const Matrix& f )

{
  JEM_ASSERT( f.size(1)    == 2 && 
              pk2.size()   == 3 &&
              stre.size(0) == 2 &&
              stre.size(1) == 2 );

  stre(0,0) = f(0,0)*pk2[0] + f(0,1)*pk2[2];
  stre(0,1) = f(0,0)*pk2[2] + f(0,1)*pk2[1];
  stre(1,0) = f(1,0)*pk2[0] + f(1,1)*pk2[2];
  stre(1,1) = f(1,0)*pk2[2] + f(1,1)*pk2[1];
}


//-----------------------------------------------------------------------
//    getNominalStiff
//-----------------------------------------------------------------------

void   getNominalStiff

  (       Cubix & nK,
    const Matrix& K,
    const Matrix& f )

{
  int nDof = nK.size(2);

  for ( int i = 0; i < nDof; ++i ) 
  {
    Matrix  thisNK = nK(ALL,ALL,i);
    getNominalStress ( thisNK, K(ALL,i), f );
  }
}



//-----------------------------------------------------------------------
//    voigtToTensorStress
//-----------------------------------------------------------------------

Matrix voigtToTensorStress ( const Vector& vec )

{
  Matrix tens;
  int nc = vec.size();

  tens = 0.0;

  if ( nc == 1 )
  {
    tens.resize( 1, 1 );
    tens(0,0) = vec[0];
  }
  else if ( nc == 3 )
  {
    tens.resize( 2, 2 );
    tens(0,0) = vec[0];
    tens(1,1) = vec[1];
    tens(0,1) = tens(1,0) = vec[2];
  }
  else if ( nc == 4 )
  {
    tens.resize( 3, 3 ); tens = 0.0;
    tens(0,0) = vec[0];
    tens(1,1) = vec[1];
    tens(2,2) = vec[2];
    tens(0,1) = tens(1,0) = vec[3];
  }
  else if ( nc == 6 )
  {
    tens.resize( 3, 3 );
    tens(0,0) = vec[0];
    tens(1,1) = vec[1];
    tens(2,2) = vec[2];
    tens(0,1) = tens(1,0) = vec[3];
    tens(1,2) = tens(2,1) = vec[4];
    tens(2,0) = tens(0,2) = vec[5];
  }
  return tens;
}

//-----------------------------------------------------------------------
//    tensorToVoigtStrain
//-----------------------------------------------------------------------

void tensorToVoigtStrain 

  (       Vector& eps,
    const Matrix& tens )

{
  int rank = tens.size(0);

  eps = 0.0;

  for ( int i = 0; i < rank; ++i )
  {
    eps[i] = tens(i,i);
  }

  if ( rank == 2 )
  {
    eps[3] = 2.*tens(0,1); //shear strain Exy in last position
  }
  else if ( rank == 3 )
  {
    eps[3] = 2.*tens(0,1);
    eps[4] = 2.*tens(1,2);
    eps[5] = 2.*tens(0,2);
  }
}


//-----------------------------------------------------------------------
//    tensorToVoigtStress
//-----------------------------------------------------------------------

void tensorToVoigtStress 

  (       Vector& s,
    const Matrix& tens )

{
  int rank = tens.size(0);

  s = 0.0;

  for ( int i = 0; i < rank; ++i )
  {
    s[i] = tens(i,i);
  }

  if ( rank == 2 )
  {
    s[3] = tens(0,1); //shear in last position
  }
  else if ( rank == 3 )
  {
    s[3] = tens(0,1);
    s[4] = tens(1,2);
    s[5] = tens(0,2);
  }
}


//-----------------------------------------------------------------------
// Convert any rectangular matrix to Column vector
//-----------------------------------------------------------------------

void getMatrixtoVector 
	( 	Vector& fvec,
	  const Matrix& tens )
{
  int n_i = tens.size(0);   int n_j = tens.size(1);  int nn = n_i * n_j;

  JEM_ASSERT( fvec.size() == nn );

  int ii = 0;

  for ( int i = 0; i < n_i; ++i )
  {
    for ( int j = 0; j < n_j; ++j )
    {
      fvec[ii] = tens(i,j); ii += 1;
    }
  }
}


//-----------------------------------------------------------------------
//   get1DShapeGradsTL
//-----------------------------------------------------------------------


void              get1DShapeGradsTL

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f )

{
  JEM_PRECHECK ( b.size(0) == 1 &&
                 g.size(0) == 1 &&
                 f.size(0) == 1 &&
                 b.size(1) == g.size(1) );

  b = g * f;
}


//-----------------------------------------------------------------------
//   get2DShapeGradsTL
//-----------------------------------------------------------------------

// VP Nguyen, 6 October 2014
// in 2D, B matrix has dimension 4x2n where n is the number of nodes
// epsilon_xy in the last row. This is needed for constitutive models
// where sigma_zz and epsilon_zz are required for 2D problems.

void              get2DShapeGradsTL

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f )

{
  JEM_PRECHECK ( b.size(0) == 4 &&
                 g.size(0) == 2 &&
                 f.size(0) == 2 &&
                 b.size(1) == 2 * g.size(1) );

  const int  nodeCount = g.size (1);

  b = 0.0;

  int    i, i1; 

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    i  = 2 * inode;

    b(0,i + 0 ) = g(0,inode) * f(0,0);  
    b(0,i + 1 ) = g(0,inode) * f(1,0);

    b(1,i + 0 ) = g(1,inode) * f(0,1);  
    b(1,i + 1 ) = g(1,inode) * f(1,1);

    b(3,i + 0 ) = g(0,inode) * f(0,1) + g(1,inode) * f(0,0);
    b(3,i + 1 ) = g(0,inode) * f(1,1) + g(1,inode) * f(1,0);
  }
}


//-----------------------------------------------------------------------
//   get3DShapeGradsTL
//-----------------------------------------------------------------------

// The strain-displacement matrix B is for a strain vector stored
// as [epsilon_xx, epsilon_yy, epsilon_zz, epsilon_xy, epsilon_yz, epsilon_zx].

void              get3DShapeGradsTL

  ( const Matrix&   b,
    const Matrix&   g,
    const Matrix&   f )

{
  JEM_PRECHECK ( b.size(0) == 6 &&
                 g.size(0) == 3 &&
                 f.size(0) == 3 &&
                 b.size(1) == 3 * g.size(1) );

  const int  nodeCount = g.size (1);

  b = 0.0;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    int  i = 3 * inode;

    b(0,i + 0) = g(0,inode) * f(0,0);
    b(0,i + 1) = g(0,inode) * f(1,0);	
    b(0,i + 2) = g(0,inode) * f(2,0);

    b(1,i + 1) = g(1,inode) * f(0,1);
    b(1,i + 1) = g(1,inode) * f(1,1);
    b(1,i + 2) = g(1,inode) * f(2,1);

    b(2,i + 2) = g(2,inode) * f(0,2);
    b(2,i + 1) = g(2,inode) * f(1,2);
    b(2,i + 2) = g(2,inode) * f(2,2);

    b(3,i + 0) = g(1,inode) * f(0,0) + g(0,inode) * f(0,1);
    b(3,i + 1) = g(1,inode) * f(1,0) + g(0,inode) * f(1,1);
    b(3,i + 2) = g(1,inode) * f(2,0) + g(0,inode) * f(2,1);

    b(4,i + 0) = g(2,inode) * f(0,1) + g(1,inode) * f(0,2);
    b(4,i + 1) = g(2,inode) * f(1,1) + g(1,inode) * f(1,2);
    b(4,i + 2) = g(2,inode) * f(2,1) + g(1,inode) * f(2,2);

    b(5,i + 0) = g(0,inode) * f(0,2) + g(2,inode) * f(0,0);
    b(5,i + 1) = g(0,inode) * f(1,2) + g(2,inode) * f(1,0);
    b(5,i + 2) = g(0,inode) * f(2,2) + g(2,inode) * f(2,0);

  }
}


//-----------------------------------------------------------------------
//   getShapeGradsTLFunc
//-----------------------------------------------------------------------


ShapeGradsTLFunc getShapeGradsTLFunc ( int rank )
{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );

  if      ( rank == 1 )
  {
    return & get1DShapeGradsTL;
  }
  else if ( rank == 2 )
  {
    return & get2DShapeGradsTL;
  }
  else
  {
    return & get3DShapeGradsTL;
  }
}


//-----------------------------------------------------------------------
// evaluate geometric stiffness matrix
//-----------------------------------------------------------------------

void evalGeometricStiffness 
  (   Matrix& tens,
    const Matrix& grads,
    const Vector& stress )
{
  int rank = grads.size(0);
  int nodeCount = grads.size(1);
  int dispCount = tens.size(0);

  MChain3     mc3;

  Matrix   H    (nodeCount, nodeCount);

  Matrix  stressTensor (rank, rank);


  if(rank == 2)
  {
    Matrix  tens1(3, 3);

    tens1 = voigtToTensorStress( stress ); //get S in tensor form

    for ( int i = 0; i < rank; i++ )
    {
      for ( int j = 0; j < rank; j++ )
      {
      stressTensor(i,j) = tens1(i,j);
      }
    }
  }
  else if (rank == 3)
  {
    stressTensor = voigtToTensorStress( stress ); //get S in tensor form
  }


  H = mc3.matmul ( grads.transpose(), stressTensor, grads ); //H = (grads)^T * S * grads

  //System::out() << "\n H \n" << H << "\n";

  for ( int inode = 0; inode < nodeCount; inode++ )
  {

    int idof = inode * rank; 

    for ( int jnode = 0; jnode < nodeCount; jnode++ )
    {

        int jdof = jnode * rank;

        tens( idof, jdof ) = tens( idof+1, jdof+1 ) = H(inode, jnode);
      }
    }
}


//-----------------------------------------------------------------------
// get Cauchy-Green tensor and its inverse (in Voigt notation) from given Lagrange strain
//-----------------------------------------------------------------------
  //3D: [Exx Eyy Ezz 2Exy 2Eyz 2Ezx] => [Cxx Cyy Czz Cxy Cyz Czx]
  //2D: [Exx Eyy Ezz 2Exy] => [Cxx Cyy Czz Cxy]
  		//C = [Cxx Cxy 0; Cxy Cyy 0; 0 0 Czz]
  // same approch for inv(C)
//-------------------------------------------------------------------

void getCauchyTensor 
	( const Vector& strain,
   		Vector&  C,
   		Vector&  invC,
   		Matrix&  invCtensor,
   		double&  I1,
   		double&  I2,
   		double&  I3 )
{
  invCtensor = 0.0; invC = 0.0;  

  I1 = I2 = I3 = 0.0;

  int strCount_ = strain.size(0);

  //-------------------------------------------------------------------
  // Green-C tensor in Voigt notation:: C_{ij} = 2 * E_{ij} + I_{ii}
  //-------------------------------------------------------------------

  for ( int i = 0; i < strCount_; ++i )
  {
  	if (i < 3)
    {
      	C[i] = 2.*strain[i] + 1.0; 

      	//if(strCount_ == 4) C[2] = 1.0;//plane strain

      	I1 += C[i]; //trace(C)
      	I2 += C[i] * C[i]; //trace(C^2)
    } 
    else
    {
    	C[i] = strain[i]; //2Exy --> Cxy
    }
  }

  I2 = 0.5 * ( I1*I1 - I2 );// 0.5*{ (tr.C)^2 - tr.(C^2)} 


  // I3 = det(C)
  //-------------------------------------

  if ( strCount_ == 6 )
  {
    I3 =  C[0] * C[1] * C[2] 
    + 2.0*C[2] * C[3] * C[4] 
    - C[0] * C[4] * C[4] 
    - C[1] * C[5] * C[5]
    - C[2] * C[3] * C[3];
  }

  else if ( strCount_ == 4 )
  {  	
    I3 = C[2] * ( C[0] * C[1] - C[3] * C[3] ); 
  }


  // inverse(C)
  //-------------------------------------

  if ( strCount_ == 6 )
  {
    invC[0] = C[1]*C[2] - C[4]*C[4]; 
    invC[1] = C[0]*C[2] - C[5]*C[5];
    invC[2] = C[0]*C[1] - C[3]*C[3];

    invC[3] = C[4]*C[5] - C[2]*C[3]; 
    invC[4] = C[3]*C[5] - C[0]*C[4];
    invC[5] = C[3]*C[4] - C[1]*C[5];
  }

  else if ( strCount_ == 4 )
  {
     invC[0] = C[1]*C[2]; 
     invC[1] = C[0]*C[2];
     invC[2] = C[0]*C[1] - C[3]*C[3]; 
     invC[3] = -C[2]*C[3]; 
  }

  invC = invC / I3; //inv(C) = cofactor/det(C);

  getInvCtensor(invC, invCtensor);
}



//-----------------------------------------------------------------------
// evaluate C_ijkl from given inv(C)
//-----------------------------------------------------------------------

void getInvCtensor 
  (  const Vector&   invC,
           Matrix&   invCtensor )
{
  int strCount_ = invC.size(0);

  invCtensor = 0.0;

  // C_{ijkl} tensor
  //-------------------------------------

  if ( strCount_ == 6)
  {
    invCtensor(0,0) = invC[0]*invC[0];  
    invCtensor(1,1) = invC[1]*invC[1]; 
    invCtensor(2,2) = invC[2]*invC[2]; 

    invCtensor(3,3) = 0.5*invC[3]*invC[3] + 0.5*invC[0]*invC[1];
    invCtensor(4,4) = 0.5*invC[4]*invC[4] + 0.5*invC[1]*invC[2]; 
    invCtensor(5,5) = 0.5*invC[5]*invC[5] + 0.5*invC[0]*invC[2];

    invCtensor(0,1) = invC[3]*invC[3]; 
    invCtensor(0,2) = invC[5]*invC[5]; 
    invCtensor(0,3) = invC[0]*invC[3]; 
    invCtensor(0,4) = invC[3]*invC[5]; 
    invCtensor(0,5) = invC[0]*invC[5];

    invCtensor(1,2) = invC[4]*invC[4];
    invCtensor(1,3) = invC[1]*invC[3]; 
    invCtensor(1,4) = invC[1]*invC[4]; 
    invCtensor(1,5) = invC[3]*invC[4];

    invCtensor(2,3) = invC[4]*invC[5]; 
    invCtensor(2,4) = invC[2]*invC[4]; 
    invCtensor(2,5) = invC[2]*invC[5];

    invCtensor(3,4) = 0.5*invC[1]*invC[5] + 0.5*invC[3]*invC[4]; 
    invCtensor(3,5) = 0.5*invC[0]*invC[4] + 0.5*invC[3]*invC[5];

    invCtensor(4,5) = 0.5*invC[2]*invC[3] + 0.5*invC[4]*invC[5];


    //generate the lower triangular matrix from symmetry
    for ( int j = 0; j < strCount_; ++j )
    {
      for ( int i = j+1; i < strCount_; ++i )
      {
        invCtensor(i,j) = invCtensor(j,i);
      }
    }
  }

  else if ( strCount_ == 4)
  {
    invCtensor(0,0) = invC[0]*invC[0];  
    invCtensor(1,1) = invC[1]*invC[1]; 
    invCtensor(2,2) = invC[2]*invC[2];

    invCtensor(3,3) = 0.5*invC[3]*invC[3] + 0.5*invC[0]*invC[1]; 

    invCtensor(0,1) = invC[3]*invC[3];
    invCtensor(0,3) = invC[0]*invC[3];
    invCtensor(1,3) = invC[1]*invC[3];

    //generate the lower triangular matrix from symmetry
    for ( int j = 0; j < strCount_; ++j )
    {
      for ( int i = j+1; i < strCount_; ++i )
      {
        invCtensor(i,j) = invCtensor(j,i);
      }
    }
  }

  // provide the minus sign

  invCtensor = - invCtensor;

}


//-----------------------------------------------------------------------
// compute Cauchy stress from given pk2 stress
// sigma = 1/J * F * S * F'
//-----------------------------------------------------------------------

void getCauchyStress 
  (       Vector&   cstress,
    const Vector&   pk2,
    const Matrix&   f )
{

  int rank_ = f.size(0);

  MChain3     mc3;

  double J;

  getDeterminant(J, f);


  if (rank_ == 2)
  {
    cstress[0] = f(0,0)*( f(0,0)*pk2[0] + f(0,1)*pk2[3] ) + f(0,1)*( f(0,1)*pk2[1] + f(0,0)*pk2[3] );
    cstress[1] = f(1,0)*( f(1,0)*pk2[0] + f(1,1)*pk2[3] ) + f(1,1)*( f(1,1)*pk2[1] + f(1,0)*pk2[3] );
    cstress[3] = f(1,0)*( f(0,0)*pk2[0] + f(0,1)*pk2[3] ) + f(1,1)*( f(0,1)*pk2[1] + f(0,0)*pk2[3] );
    cstress[2] = pk2[2];
  }
  else if (rank_ == 3)
  {
    Matrix  tens(rank_, rank_ ); Matrix  ss(rank_, rank_ );

    tens = voigtToTensorStress( pk2 );

    ss = mc3.matmul( f, tens, f.transpose() );   

    tensorToVoigtStress ( cstress, ss );
  }

  cstress = cstress/J; 
}



//-----------------------------------------------------------------------
// compute determinant of a Matrix
//-----------------------------------------------------------------------

void getDeterminant 
  (       double&  J,
    const Matrix&  f)
  {

    int rank_ = f.size(0);

    if (rank_ == 2)
    {

      J = f(0,0)*f(1,1) - f(1,0)*f(0,1);
    }
    else if (rank_ == 3)
    {

      J = f(0,0)*f(1,1)*f(2,2) + f(0,1)*f(1,2)*f(2,0) + f(0,2)*f(1,0)*f(2,1)
      - f(0,0)*f(1,2)*f(2,1) - f(0,1)*f(1,0)*f(2,2) - f(0,2)*f(1,1)*f(2,0);
    }

  }


//-----------------------------------------------------------------------
// get inverse Cauchy-Green tensor (in Voigt notation) from given eigenvalues/vector of "C"
//-----------------------------------------------------------------------

void getCauchyTensorEigen 

  ( const Vector&  evals,
    const Matrix&  evecs,
    	Vector&  invCp,
    	Vector&  invCn,
      	Matrix&  invCtensorP,
      	Matrix&  invCtensorN )
{
  int strCount_ = invCp.size(0);

  Vector      beta(3);

  invCp = invCn = 0.0;

  invCtensorP = invCtensorN = 0.0;

  //----------------------------
  // get invC+/-
  //----------------------------

 for ( int n = 0; n < 3; ++n )
  {
    beta = evecs( ALL, n );

    if ( evals[n] > 1. )
    {
    	for ( int i = 0; i < 3; ++i )
	    {
	      invCp[i] += 1./evals[n] * beta[i] * beta[i];
	    }

	    if( strCount_ == 4 )
	    {
	      invCp[3] += 1./evals[n] * beta[0] * beta[1];
	    }
	    else if( strCount_ == 6 )
	    {
	      invCp[3] += 1./evals[n] * beta[0] * beta[1];
	      invCp[4] += 1./evals[n] * beta[1] * beta[2];
	      invCp[5] += 1./evals[n] * beta[0] * beta[2];
	    }
	    else
	    {
	      System::out() << "Not supported " << endl;
	    }
  	}
  	else
    {
    	for ( int i = 0; i < 3; ++i )
	    {
	      invCn[i] += 1./evals[n] * beta[i] * beta[i];
	    }

	    if( strCount_ == 4 )
	    {
	      invCn[3] += 1./evals[n] * beta[0] * beta[1];
	    }
	    else if( strCount_ == 6 )
	    {
	      invCn[3] += 1./evals[n] * beta[0] * beta[1];
	      invCn[4] += 1./evals[n] * beta[1] * beta[2];
	      invCn[5] += 1./evals[n] * beta[0] * beta[2];
	    }
	    else
	    {
	      System::out() << "Not supported " << endl;
	    }
  	}
  }

    //----------------------------
    // get C_{ijkl}+/-
    //---------------------------- 

	getInvCtensor(invCp, invCtensorP);

	getInvCtensor(invCn, invCtensorN);    
}


//-----------------------------------------------------------------------
// get Beta tensors:beta_n_ij and beta_n_ijkl (in Voigt notation) from given eigenvector of "C"
//-----------------------------------------------------------------------

void getCauchyTensorEigen 

  ( const 	Matrix&  evecs,
      		Matrix&  Beta,
      		Cubix&  BetaTensor )
{
  int strCount_ = Beta.size(1);

  Vector      beta(3);

  //----------------------------
  // get Beta_n
  //----------------------------

 for ( int n = 0; n < 3; ++n )
  {
    beta = evecs( ALL, n );

	for ( int i = 0; i < 3; ++i )
    {
      Beta(n,i) = beta[i] * beta[i];
    }

    if( strCount_ == 4 )
    {
      Beta(n,3) = beta[0] * beta[1];
    }
    else if( strCount_ == 6 )
    {
      Beta(n,3) = beta[0] * beta[1];
      Beta(n,4) = beta[1] * beta[2];
      Beta(n,5) = beta[0] * beta[2];
    }
    else
    {
      System::out() << "Not supported " << strCount_ << endl;
    }

	//----------------------------
	// get Beta_n_ijkl
	//----------------------------

    for ( int i = 0; i < 3; ++i )
    {
    	for ( int j = i; j < 3; ++j )
    	{
    		BetaTensor(n,i,j) = beta[i] * beta[i] * beta[j] * beta[j];
    	}
    }

	if( strCount_ == 4 )
	{
  		BetaTensor(n,0,3) = beta[0] * beta[0] * beta[0] * beta[1];

  		BetaTensor(n,1,3) = beta[0] * beta[1] * beta[1] * beta[1];

  		BetaTensor(n,2,3) = beta[0] * beta[1] * beta[2] * beta[2];

  		BetaTensor(n,3,3) = beta[0] * beta[0] * beta[1] * beta[1];

  		//generate the lower triangular matrix from symmetry

	    for ( int j = 0; j < strCount_; ++j )
	    {
	      for ( int i = j+1; i < strCount_; ++i )
	      {
	        BetaTensor(n,i,j) = BetaTensor(n,j,i);
	      }
	    }
	}
	else if( strCount_ == 6 )
	{
		//b0^3*b1, b0^2*b1*b2,    b0^3*b2
  		BetaTensor(n,0,3) = beta[0] * beta[0] * beta[0] * beta[1];
  		BetaTensor(n,0,4) = beta[0] * beta[0] * beta[1] * beta[2];
  		BetaTensor(n,0,5) = beta[0] * beta[0] * beta[0] * beta[2];

  		//b0*b1^3,    b1^3*b2, b0*b1^2*b2
  		BetaTensor(n,1,3) = beta[0] * beta[1] * beta[1] * beta[1];
  		BetaTensor(n,1,4) = beta[1] * beta[1] * beta[1] * beta[2];
  		BetaTensor(n,1,5) = beta[0] * beta[1] * beta[1] * beta[2];

  		//b0*b1*b2^2,    b1*b2^3,    b0*b2^3
  		BetaTensor(n,2,3) = beta[0] * beta[1] * beta[2] * beta[2];
  		BetaTensor(n,2,4) = beta[1] * beta[2] * beta[2] * beta[2];
  		BetaTensor(n,2,5) = beta[0] * beta[2] * beta[2] * beta[2];

  		//b0^2*b1^2, b0*b1^2*b2, b0^2*b1*b2
  		BetaTensor(n,3,3) = beta[0] * beta[0] * beta[1] * beta[1];
  		BetaTensor(n,3,4) = beta[0] * beta[1] * beta[1] * beta[2];
  		BetaTensor(n,3,5) = beta[0] * beta[0] * beta[1] * beta[2];

  		//b1^2*b2^2, b0*b1*b2^2
  		BetaTensor(n,4,4) = beta[1] * beta[1] * beta[2] * beta[2];
  		BetaTensor(n,4,5) = beta[0] * beta[1] * beta[2] * beta[2];

  		//b0^2*b2^2
  		BetaTensor(n,5,5) = beta[0] * beta[0] * beta[2] * beta[2];

  		//generate the lower triangular matrix from symmetry
  		
	    for ( int j = 0; j < strCount_; ++j )
	    {
	      for ( int i = j+1; i < strCount_; ++i )
	      {
	        BetaTensor(n,i,j) = BetaTensor(n,j,i);
	      }
	    }
	}
	else
    {
      System::out() << "Not supported " << strCount_ << endl;
    }

  }


}



//-----------------------------------------------------------------------
// update stiff matrix after local iteration for plane stress state
//-----------------------------------------------------------------------

void updateStiffPlaneStress 
  (       Matrix&   stiffS )

{

  // variables for plane stress local iterations
  Matrix    stiff_mm(3, 3);
  Matrix    stiff_m3(3,1);
  Matrix    stiff_3m(1,3);
  double    stiff_33(0.0); 

  int strCount_ = 4;

// push z-components to last row and last colm
  // swap row-3 & row-4 and col.3 & col.4

  swapRows ( stiffS, 2, 3 ); //swap row-3 & row-4
  swapColms ( stiffS, 2, 3 ); //swap col.3 & col.4


  // get the blocks of matrices

  stiff_33 = stiffS(3,3); //last element after re-arrange

  for ( int i = 0; i < strCount_; ++i )
  {
    for ( int j = 0; j < strCount_; ++j )
    {
      if ( i < 3 && j < 3 ) stiff_mm(i,j) = stiffS(i,j); 

      if ( i < 3 && j >= 3 ) stiff_m3(i,0) = stiffS(i,j); 

      if ( i >= 3 && j < 3 ) stiff_3m(0,j) = stiffS(i,j); 
    }
  }


  // update stiff_mm

  stiff_mm = stiff_mm - matmul(stiff_m3, stiff_3m) * 1./stiff_33;

  // update stiffS

  for ( int i = 0; i < 3; ++i )
  {
    for ( int j = 0; j < 3; ++j )
    {
      stiffS(i,j) = stiff_mm(i,j); 
    }
  }

  // swap back to original positions

  swapColms ( stiffS, 2, 3 ); //swap col.3 & col.4
  swapRows  ( stiffS, 2, 3 ); //swap row-3 & row-4
}