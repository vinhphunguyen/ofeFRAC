
#include "XtraUtilities.h"
#include <jem/base/System.h>
#include <jem/numeric/algebra/LUSolver.h>
#include <jem/numeric/utilities.h>

using namespace jem::io;

using jem::numeric::det;
//using jem::numeric::abs;

//-----------------------------------------------------------------------
//   find()
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
//   SortAsc     Sorts a vector in asscending order and returns
//               array of indices
//-----------------------------------------------------------------------

IdxVector      SortAsc ( const  Vector&   A )

{
    int         e       = A.size();  // no. of entities
    int         count   = 0;
    double      min;

    Vector      tmp   ( e );
    IdxVector   order ( e );
    Vector      live  ( e );

    tmp         = A;
    live        = 1;
 
            
    for (int i = 0; i < e; i++)
    {
        for (int j = 0; j < e; j++)
        {
            if ( live[j] > 0 )
            {
                count += 1;

                if (count == 1)
                {
                    min = tmp[j];
                }

                if (tmp[j] <= min)
                {
                    min         = tmp[j];
                    order[i]    = j;
                }
            }
            else
            {
                continue;
            }
        }
        live[order[i]]  = 0;
        count           = 0;
    }

    return  order;
}


//-----------------------------------------------------------------------
//   SortDsc     Sorts a vector in Descending order and returns
//               array of indices
//-----------------------------------------------------------------------

IdxVector      SortDsc ( const  Vector&   A )

{
    int         e       = A.size();  // no. of entities
    int         count   = 0;
    double      max;

    Vector      tmp   ( e );
    IdxVector   order ( e );
    Vector      live  ( e );

    tmp         = A;
    live        = 1;
 
            
    for (int i = 0; i < e; i++)
    {
        for (int j = 0; j < e; j++)
        {
            if ( live[j] > 0 )
            {
                count += 1;

                if (count == 1)
                {
                    max = tmp[j];
                }

                if (tmp[j] >= max)
                {
                    max         = tmp[j];
                    order[i]    = j;
                }
            }
            else
            {
                continue;
            }
        }
        live[order[i]]  = 0;
        count           = 0;
    }

    return  order;
}

//-----------------------------------------------------------------------
//    multiplyVoigtNormal
//-----------------------------------------------------------------------

Vector multiplyVoigtNormal

  ( const Vector& sig,
    const Vector& n )

{
  int rank     =   n.size();

  Vector sigN( rank );

  if ( rank == 2 )
  {
    JEM_ASSERT( sig.size() == 3 );

    sigN[0] = sig[0]*n[0] + sig[2]*n[1];
    sigN[1] = sig[2]*n[0] + sig[1]*n[1];
  }
  else
  {
    JEM_ASSERT( sig.size() == 6 );

    sigN[0] = sig[0]*n[0] + sig[3]*n[1] + sig[5]*n[2];
    sigN[1] = sig[3]*n[0] + sig[1]*n[1] + sig[4]*n[2];
    sigN[2] = sig[5]*n[0] + sig[4]*n[1] + sig[2]*n[2];
  }

  return sigN;
}

//-----------------------------------------------------------------------
//    overloaded : multiplyVoigtNormal
//-----------------------------------------------------------------------

Matrix multiplyVoigtNormal

  ( const Matrix& db,
    const Vector& n )

{
  int dofCount = db.size(1);
  int rank     =  n.size();

  Matrix  dbN ( rank, dofCount );

  for ( int i = 0; i < dofCount; ++i )
  {
    dbN(ALL,i) = multiplyVoigtNormal( db(ALL,i), n );
  }

  return dbN;
}


//-----------------------------------------------------------------------
//  voigt2ind_
//-----------------------------------------------------------------------

Matrix        voigt2ind    ( const Vector&       A)

{
  
  Matrix    B;
  //System::out() << "A.Size()[0] is " << A.size()[0] << endl;
  System::out() << "A.size(0) is << " << A.size() << endl;

  int     rows  = A.size(0);
  
  if ( rows == 6 )
  {
    // rank_ = 3
    B.resize ( 3 , 3 );

    B(0,0) = A[0];
    B(1,1) = A[1];
    B(2,2) = A[2];
    B(0,1) = B(1,0) = A[3];
    B(0,2) = B(2,0) = A[5];
    B(1,2) = B(2,1) = A[4];

  }
  else if ( rows == 3 )
  {
    B.resize ( 2 , 2 );

    B(0,0) = A[0];
    B(1,1) = A[1];
    B(0,1) = B(1,0) = A[2];
  }
    return B;
}

//-----------------------------------------------------------------------
//  ind2voigt_
//-----------------------------------------------------------------------

Vector      ind2voigt    (const Matrix&        A)
{

  Vector    B;

  int rank = A.size(0);

  if ( rank == 3 )
  {
    B.resize ( 6 );

    B[0] = A(0,0);
    B[1] = A(1,1);
    B[2] = A(2,2);
    B[3] = A(0,1);
    B[4] = A(1,2);
    B[5] = A(0,2);

  }
  else if ( rank == 2 )
  {
    B.resize ( 3 );

    B[0] = A(0,0);
    B[1] = A(1,1);
    B[2] = A(0,1);

  }
    return B;    
}

//-----------------------------------------------------------------------
//  Incircle (NB: perform check in 2D)
//-----------------------------------------------------------------------
bool             Incircle

   ( const Vector&       center,
     const double&       R,
     const Vector&       pt )

{
 double x0, y0, x, y, Inside;
 double ax, bx, cx, dx, ay, by, cy, dy;
 bool flag;
 Matrix mat(3,3);

 x0 = center[0]; y0=center[1];
 x=pt[0];  y=pt[1];

 //Incircle test
 ax=x0-R;    ay=y0;
 bx=x0+R;    by=y0;
 cx=x0;      cy=y0+R;
 dx=x;       dy=y;

 mat(0,0) = ax-dx;
 mat(0,1) = ay-dy;
 mat(0,2) = (ax-dx) * (ax-dx) + (ay-dy) * (ay-dy);

 mat(1,0) = bx-dx;
 mat(1,1) = by-dy;
 mat(1,2) = (bx-dx) * (bx-dx) + (by-dy) * (by-dy);

 mat(2,0) = cx-dx;
 mat(2,1) = cy-dy;
 mat(2,2) = (cx-dx) * (cx-dx) + (cy-dy) * (cy-dy);

 Inside  =det(mat);
 Inside /= abs(Inside);

 flag = !(Inside < 0.);

 return flag;
}
                                                  


//-----------------------------------------------------------------------
//  swapRows
//-----------------------------------------------------------------------
void  swapRows

  (  Matrix&  f,
    const int& m,
    const int& n )
  {
    int mm = f.size(0); 
    int nn = f.size(1);

    double tt; //temp variable

    for ( int j = 0; j < nn; ++j ) //loop over colms
      {
        tt = f(m, j); // stote elements of m-th row 
        f(m, j) = f(n, j); // replace row-m with row-n
        f(n, j) = tt; // replace row-n with tempVar
      }

  }


//-----------------------------------------------------------------------
//  swapColms
//-----------------------------------------------------------------------
void  swapColms

  (  Matrix&  f,
    const int& m,
    const int& n )
  {
    int mm = f.size(0); 
    int nn = f.size(1);

    double tt; //temp variable

    for ( int i = 0; i < mm; ++i ) //loop over colms
      {
        tt = f(i, m); // stote elements of m-th row 
        f(i, m) = f(i, n); // replace row-m with row-n
        f(i, n) = tt; // replace row-n with tempVar
      }

  }