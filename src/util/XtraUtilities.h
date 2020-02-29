#ifndef XTRA_UTILITIES_H
#define XTRA_UTILITIES_H

#include <jive/Array.h>
using namespace jem;

using jive::Vector;
using jive::Matrix;
using jive::IdxVector;

//-----------------------------------------------------------------------
//   public functions
//-----------------------------------------------------------------------

template <class E>

  int                     findItem

  ( IdxVector&              index,
    const Array<bool,1,E>&  expr )

{
  const int  n = expr.size ();
  int        k = count ( expr );

  index.resize ( k );

  if ( k )
  {
    idx_t*  ixp = index.addr ();

    k = 0;

    if ( expr.isContiguous() )
    {
      for ( int i = 0; i < n; i++ )
      {
        if ( expr.getFast(i) )
        {
          ixp[k++] = i;
        }
      }
    }
    else
    {
      for ( int i = 0; i < n; i++ )
      {
        if ( expr[i] )
        {
          ixp[k++] = i;
        }
      }
    }
  }

  return k;
}


  IdxVector                 SortAsc

( const  Vector&            A           );


  IdxVector                 SortDsc

( const  Vector&            A           );


//-----------------------------------------------------------------------
// Multiply Voigt stress with normal vector
//-----------------------------------------------------------------------

Vector multiplyVoigtNormal

  ( const Vector& sig,
    const Vector& n );

//-----------------------------------------------------------------------
// Multiply Voigt stiffness (DB) with normal vector
//-----------------------------------------------------------------------

Matrix multiplyVoigtNormal

  ( const Matrix& db,
    const Vector& n );


//-----------------------------------------------------------------------
//  voigt2ind
//-----------------------------------------------------------------------

Matrix        voigt2ind

  ( const Vector&       A );

//-----------------------------------------------------------------------
//  ind2voigt
//-----------------------------------------------------------------------

Vector        ind2voigt

  ( const Matrix&        A );

//-----------------------------------------------------------------------
//  Incircle
//-----------------------------------------------------------------------

bool          Incircle

 ( const Vector&    center,
   const double&    R,
   const Vector&    pt);

//-----------------------------------------------------------------------
//  swapRows
//-----------------------------------------------------------------------
void  swapRows

  (  Matrix&  f,
    const int& m,
    const int& n );

//-----------------------------------------------------------------------
//  swapColms
//-----------------------------------------------------------------------
void  swapColms

  (  Matrix&  f,
    const int& m,
    const int& n );

#endif

