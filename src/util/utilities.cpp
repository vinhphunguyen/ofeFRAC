
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>

#include <jem/io/PrintWriter.h>
#include <jem/base/PrecheckException.h>
#include <jem/base/System.h>
#include <jem/base/Error.h>
#include <jem/base/Float.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/LUSolver.h>
#include <jem/util/Properties.h>
#include <jem/util/StringUtils.h>
#include <jem/numeric/algebra/utilities.h>

#include <jive/geom/error.h>

#include "utilities.h"

#define tolerance 0.1e-20

extern "C"
{
  #include  <math.h>
}

using jem::System;
using jem::ALL;
using jem::END;
using jem::TensorIndex;
//using jem::numeric::norm2;
using jem::io::endl;
using jem::idx_t;
using jem::Array;
//using jem::numeric::det;

//-----------------------------------------------------------------------
//   constants
//-----------------------------------------------------------------------

// number of strain components (Voigt notation)
// corresponding to the rank of the problem
//                                  1D  2D  3D
const int    STRAIN_COUNTS[4] = { 0, 1, 4, 6 };

const double PI               = 3.14159265;
const double PI3              = 0.3333333333 * PI;
const double ONE_THIRD        = 0.33333333333;
const double TWO_THIRD        = 0.66666666667;
const double EPS              = 1.e-16;


//-----------------------------------------------------------------------
//   get1DShapeGrads
//-----------------------------------------------------------------------


void              get1DShapeGrads

  ( const Matrix&   b,
    const Matrix&   g )

{
  JEM_PRECHECK ( b.size(0) == 1 &&
                 g.size(0) == 1 &&
                 b.size(1) == g.size(1) );

  b = g;
}


//-----------------------------------------------------------------------
//   get2DShapeGrads
//-----------------------------------------------------------------------

// VP Nguyen, 6 October 2014
// in 2D, B matrix has dimension 4x2n where n is the number of nodes
// epsilon_xy in the last row. This is needed for constitutive models
// where sigma_zz and epsilon_zz are required for 2D problems.

void              get2DShapeGrads

  ( const Matrix&   b,
    const Matrix&   g )

{
  JEM_PRECHECK ( b.size(0) == 4 &&
                 g.size(0) == 2 &&
                 b.size(1) == 2 * g.size(1) );

  const int  nodeCount = g.size (1);

  b = 0.0;

  int    i, i1;
  double Nix, Niy;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    i  = 2 * inode;
    i1 = i + 1;

    Nix = g(0,inode);
    Niy = g(1,inode);

    b(0,i ) = Nix;
    b(1,i1) = Niy;

    b(3,i ) = Niy;
    b(3,i1) = Nix;
  }
}


//-----------------------------------------------------------------------
//   get3DShapeGrads
//-----------------------------------------------------------------------

// The strain-displacement matrix B is for a strain vector stored
// as [epsilon_xx, epsilon_yy, epsilon_zz, epsilon_xy, epsilon_yz, epsilon_zx].

void              get3DShapeGrads

  ( const Matrix&   b,
    const Matrix&   g )

{
  JEM_PRECHECK ( b.size(0) == 6 &&
                 g.size(0) == 3 &&
                 b.size(1) == 3 * g.size(1) );

  const int  nodeCount = g.size (1);

  b = 0.0;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    int  i = 3 * inode;

    b(0,i + 0) = g(0,inode);
    b(1,i + 1) = g(1,inode);
    b(2,i + 2) = g(2,inode);

    b(3,i + 0) = g(1,inode);
    b(3,i + 1) = g(0,inode);

    b(4,i + 1) = g(2,inode);
    b(4,i + 2) = g(1,inode);

    b(5,i + 2) = g(0,inode);
    b(5,i + 0) = g(2,inode);
  }
}


//-----------------------------------------------------------------------
//   getShapeGradsFunc
//-----------------------------------------------------------------------


ShapeGradsFunc getShapeGradsFunc ( int rank )
{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );

  if      ( rank == 1 )
  {
    return & get1DShapeGrads;
  }
  else if ( rank == 2 )
  {
    return & get2DShapeGrads;
  }
  else
  {
    return & get3DShapeGrads;
  }
}

// ----------------------------------------------------------------------
//   compute hydrostatic stress
// ----------------------------------------------------------------------

double (computePStressHydroStress)

( double nu,
  const Vector& sigma )
{
  return 0.5 * ( sigma[0]  + sigma[1] );
}

double (computePStrainHydroStress)

( double nu,
  const Vector& sigma )
{
  return .333333333 * (1. + nu) * ( sigma[0] + sigma[1] );
}

double (compute3DHydroStress)

( double nu,
  const Vector& sigma )
{
  return .333333333 * ( sigma[0] + sigma[1] + sigma[2] );
}

ComputeHydroStress    getHydroStressFunc

       ( ProblemType prob )
{
  switch(prob)
  {
    case PlaneStrain  :
      return &computePStrainHydroStress;
      break; //optional
    case PlaneStress  :
      return &computePStressHydroStress;
      break; //optional
    default : //Optional
       return &compute3DHydroStress;
  }
}

// --------------------------------------------------------------------
//  get1DShapeFuncs
// --------------------------------------------------------------------

void                  get1DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n )
{
  sfuncs(0,0) = n[0];
}

// --------------------------------------------------------------------
//  get2DShapeFuncs
// --------------------------------------------------------------------

void                  get2DShapeFuncs

  ( const Matrix&       s,
    const Vector&       n )
{
  JEM_PRECHECK ( s.size(0) == 2 &&
                 s.size(1) == 2 * n.size() );

  const int  nodeCount = n.size ();

  s = 0.0;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    int  i = 2 * inode;

    s(0,i + 0) = n[inode];
    s(1,i + 1) = n[inode];
  }
}

// --------------------------------------------------------------------
//  get3DShapeFuncs
// --------------------------------------------------------------------

void                  get3DShapeFuncs

  ( const Matrix&       s,
    const Vector&       n )
{
  JEM_PRECHECK ( s.size(0) == 3 &&
                 s.size(1) == 3 * n.size() );

  const int  nodeCount = n.size ();

  s = 0.0;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    int  i = 3 * inode;

    s(0,i + 0) = n[inode];
    s(1,i + 1) = n[inode];
    s(2,i + 2) = n[inode];
  }
}

// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.

ShapeFunc        getShapeFunc

  ( int                 rank )

{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );


  if      ( rank == 1 )
  {
    return & get1DShapeFuncs;
  }
  else if ( rank == 2 )
  {
    return & get2DShapeFuncs;
  }
  else
  {
    return & get3DShapeFuncs;
  }
}

// --------------------------------------------------------------------
//  get2DFShapeFuncs
// --------------------------------------------------------------------

void                  get2DFShapeFuncs

  ( const Matrix&       s,
    const Vector&       n )
{
  JEM_PRECHECK ( s.size(0) == 2 &&
                 s.size(1) == 4 * n.size() );

  const int  nodeCount    = n.size ();
  const int  nodeCount2   = 2 * nodeCount;
  const int  nodeCount2p1 = nodeCount2 + 1;

  double ni;
  int    i;

  s = 0.0;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    i     = 2 * inode;
    ni    =   n[inode];
    s(0,i             )   =  ni;
    s(0,i + nodeCount2)   = -ni;
    s(1,i + 1)            =  ni;
    s(1,i + nodeCount2p1) = -ni;
  }
}

FShapeFunc        getFShapeFunc

  ( int                 rank )

{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );


  if ( rank == 2 )
  {
    return & get2DFShapeFuncs;
  }
}


// --------------------------------------------------------------------
//  get2DNormalFuncs
// --------------------------------------------------------------------

void                  get2DNormalFuncs

  ( const Matrix&       normal,
    const Vector&       normalV )
{
   normal(0,0) = normalV[0]; normal(0,2) = normalV[1];
   normal(1,1) = normalV[1]; normal(1,2) = normalV[0];

}

// --------------------------------------------------------------------
//  get3DNormalFuncs
// --------------------------------------------------------------------

void                  get3DNormalFuncs

  ( const Matrix&       normal,
    const Vector&       normalV )
{
   normal(0,0) = normalV[0]; normal(0,3) = normalV[1]; normal(0,5) = normalV[2];
   normal(1,1) = normalV[1]; normal(1,3) = normalV[0]; normal(1,4) = normalV[2];
   normal(2,2) = normalV[2]; normal(2,4) = normalV[1]; normal(2,5) = normalV[0];

}

// A function that returns a pointer to a function that computes the
// matrix of normal functions given the number of spatial dimensions.

NormalMatrixFunc        getNormalFunc

  ( int                 rank )

{
  if ( rank == 2 )
  {
    return & get2DNormalFuncs;
  }
  else
  {
    return & get3DNormalFuncs;
  }
}

void                  get2DNormalInMatrixFuncs

  ( const Matrix&       normal,
    const Vector&       normalV )
{
   normal(0,0) = normalV[0]; normal(0,3) = normalV[1];
   normal(1,1) = normalV[1]; normal(1,3) = normalV[0];

}

// A function that returns a pointer to a function that computes the
// matrix of normal functions given the number of spatial dimensions.

NormalInMatrixFunc    getNormalInMatrixFunc

  ( int                 rank )

{
  if ( rank == 2 )
  {
    return & get2DNormalInMatrixFuncs;
  }
  else
  {
    return & get3DNormalFuncs;
  }
}

void                  get2DTransMatrixFunc

  (       Matrix&       Q,
    const Matrix&       grads,
    const Matrix&       coord )
{
  double x0 = coord(0,1) - coord(0,0);
  double y0 = coord(1,1) - coord(1,0);

  double alpha     = ::atan2 ( y0, x0 );
  double sinAlpha  = ::sin (alpha);
  double cosAlpha  = ::cos (alpha);

  Q(0,0) = - sinAlpha; Q(0,1) = cosAlpha;
  Q(1,0) =   cosAlpha; Q(1,1) = sinAlpha;
}

void                  get3DTransMatrixFunc

  (       Matrix&       Q,
    const Matrix&       xgrads,
    const Matrix&       coord )
{
  using namespace jem;
  using jem::numeric::crossProduct;
  using jem::numeric::dotProduct;
  using namespace jive::geom;

  Tuple<double,3>     n, s1, s2;

  Vector              p(3);
  double              xg0, xg1;

  s1 = 0.0;
  s2 = 0.0;
  
  const int  nodeCount = coord.size (1);

  for ( int i = 0; i < nodeCount; i++ )
  {
    xg0   = xgrads(0,i);
    xg1   = xgrads(1,i);

    s1[0] += coord(0,i) * xg0;
    s1[1] += coord(1,i) * xg0;
    s1[2] += coord(2,i) * xg0;

    s2[0] += coord(0,i) * xg1;
    s2[1] += coord(1,i) * xg1;
    s2[2] += coord(2,i) * xg1;
  }

  n        = crossProduct ( s1, s2 );
  double a = ::sqrt ( dotProduct( n, n   ) );
  double b = ::sqrt ( dotProduct( s1, s1 ) );

  if ( jem::Float::isTiny( a ) || jem::Float::isTiny( b ) )
  {
    zeroVectorError ( "normal vector", "normal" );
  }

  // normalize n and s1

  n   = (1.0 / a) * n;
  s1  = (1.0 / b) * s1;

  // make s2 orthogonal to n and s1

  s2  = crossProduct ( n, s1 );

  // build the transformation matrix Q
  /*
   * Q=[
         e1.n e1.s1 e1.s2
         e2.n e2.s1 e2.s2
         e3.n e3.s1 e3.s2
       ]  

    v_loc = Q' * v_global
    I implemented Q' as Q in the following.
   */
  
  Q(0,0) = n[0]; Q(0,1) = s1[0]; Q(0,2) = s2[0];
  Q(1,0) = n[1]; Q(1,1) = s1[1]; Q(1,2) = s2[1];
  Q(2,0) = n[2]; Q(2,1) = s1[2]; Q(2,2) = s2[2]; 
  
  //System::out() << Q << "\n";
}


TransformationMatrixFunc getTransMatrixFunc

  ( int                 rank )
{
  if ( rank == 2 )
  {
    return & get2DTransMatrixFunc;
  }
  else
  {
    return & get3DTransMatrixFunc;
  }
}

//-----------------------------------------------------------------------
//   getDoubleShearComponentsFunc
//-----------------------------------------------------------------------

void get2DShearComponents ( Vector& v )
{
    v[3] *= 2.;
}

void get3DShearComponents ( Vector& v )
{
    v[3] *= 2.;
    v[4] *= 2.;
    v[5] *= 2.;
}

DoubleShearComponentsFunc getDoubleShearComponentsFunc

  ( int                 rank )
{
  if ( rank == 2 )
  {
    return & get2DShearComponents;
  }
  else
  {
    return & get3DShearComponents;
  }
}

//-----------------------------------------------------------------------
//   check acoustic tensor stuff...
//-----------------------------------------------------------------------

// Ortiz's algorithm
// tangent is the consistent material tangent which is a 4x4 matrix for plane
// strain problems.
// det(A)=det(n*CTO*n)= quartic equation in terms of tan(theta) by dividing
// det(A)=0 with cos^4(theta).

bool                      checkAcousticTensor2D 

  (       Vector&       normal,
    const Matrix&       tangent )
{
    bool localised (false);

    double D1111 = tangent(0,0);
    double D1212 = tangent(3,3);
    double D1112 = tangent(0,3);
    double D1211 = tangent(3,0);
    double D1222 = tangent(3,1);
    double D2212 = tangent(1,3);
    double D2211 = tangent(1,0);
    double D1122 = tangent(0,1);
    double D2222 = tangent(1,1);

    double a0 = D1111*D1212 - D1112*D1211; 
    double a1 = D1111*D1222 + D1111*D2212 - D1112*D2211 - D1122*D1211; 
    double a2 = D1111*D2222 + D1112*D1222 + D1211*D2212 - D1122*D1212 - D1122*D2211 - D1212*D2211; 
    double a3 = D1112*D2222 + D1211*D2222 - D1122*D2212 - D1222*D2211;
    double a4 = D1212*D2222 - D2212*D1222; 

    // finding minimum of a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0
    // by finding roots of a cubic equation, which is the derivative of the
    // quartic function.

    Vector ans;
    solveCubicEqua ( ans, 4.*a4, 3.*a3, 2.*a2, a1 );

     
    const int noRoots = ans.size ( );

    double xm, fm(0.);
        
    System::out() << "coefs: " << a4 << " " << a3 << " " << a2 << " "  << a1 << " " << a0 << "\n";
    System::out() << "roots: " << ans << "\n";

    for ( int i = 0; i < noRoots; i++ )
    {
        double xi = ans[i];
        double fi = a4*::pow(xi,4) + a3*::pow(xi,3) + a2*::pow(xi,2) + a1*xi + a0;
            
        System::out() << "fi: " << fi << "\n";

        // roots = [-x0,xx,x0] 
        // where -x0 and x0 both give negative fm 
        // only choose x0
        if ( ( xi > 0. ) && ( fi < fm ) )
        {
            localised = true;
            xm        = xi;
            fm        = fi;
            System::out() << "localised :)\n";
            System::out() << (180/3.14)*atan(xm) <<"\n";
        }
    }

    //xm = tan(theta);

    double theta = atan ( xm );

    normal[0] = cos ( theta );
    normal[1] = sin ( theta );

    return localised;
}

bool                      checkAcousticTensor3D 

  (       Vector&       normal,
    const Matrix&       tangent )
{
    return false;
}

CheckAcousticTensorFunc   getCheckAcousticTensorFunc

  ( int                  rank )
{
  if ( rank == 2 )
  {
    return & checkAcousticTensor2D;
  }
  else
  {
    return & checkAcousticTensor3D;
  }
}

//-----------------------------------------------------------------------
//   solveLinearEqua
//-----------------------------------------------------------------------

void                  solveLinearEqua

  ( double&       ans,
    const double  a,
    const double  b )
{
  ans = -b/a;
}

//-----------------------------------------------------------------------
//   solveQuadEqua
//-----------------------------------------------------------------------

void                  solveQuadEqua

  ( Vector&       ans,
    const double  a,
    const double  b,
    const double  c )
{
  if ( jem::numeric::abs(a) < EPS )
  {
    ans.resize( 1 );
    ans[0] = - c / b;
  }
  else
  {
    double d = b * b - 4.0 * a * c;
    
    if ( d < 0 )
    {
      using namespace jem;
      throw Error (
      JEM_FUNC,
      " imaginary principal values !!! "
      );  
    }
    else
    {
      ans.resize( 2 );
      ans[0] = 0.5 * (-b + ::sqrt( d )) / a;
      ans[1] = 0.5 * (-b - ::sqrt( d )) / a;
    }
  }
}

//-----------------------------------------------------------------------
//   solveCubicEqua
//-----------------------------------------------------------------------

void                  solveCubicEqua

  ( Vector&       ans,
    const double  a,
    const double  b,
    const double  c,
    const double  d)

{  
  // a=0: reduced to quadratic equation  

  if ( jem::numeric::abs(a) < EPS ) 
  {
    // b = 0: reduced to linear equation

    if ( jem::numeric::abs(b) < EPS )
    {
      ans.resize ( 1 );
      
      ans[0] = - d / c;

      return;
    }
    else
    {
      double D = c * c - 4.0 * b * d;
      
      if ( D < 0.0 )
      {	
	    return;
      }
      else
      {
	    D      = ::sqrt(D);

	    ans.resize ( 2 );
	    
	    ans[0] = 0.5 * (-c + D) / b;
	    ans[1] = 0.5 * (-c - D) / b;
       
	    return;
      }
    }
  }
  else
  {
    // solving cubic equation using Cardano's method.  
    // normalize to have: x^3 + aa*x^2 + bb*x + cc   
    double kk  = 1.0 / a;
    
    double aa  = b * kk;
    double bb  = c * kk;
    double cc  = d * kk;

    // change variable x=y-aa/3 to eliminate quadratic term
    // x^3 + px + q = 0

    double aa2 = aa*aa;
    double p   = ( -aa2 + 3.0 * bb ) / 9.0;
    double q   = ( 2.0 * aa * aa2 - 9.0 * aa * bb + 27.0 * cc ) / 54.0;
    double p3  = p * p * p;
    double D   = p3 + q * q;
    double aa3 = ONE_THIRD*aa;
 
    if ( jem::numeric::abs(D) < EPS )
    {
      if ( jem::numeric::abs(q) < EPS )
      {
          // one triple solution
          ans.resize(1);
          ans[0] = -aa3;
          return;
      }
      else
      {
          // one single and one double solution

          double u = pow(-q,ONE_THIRD);
          ans.resize(2);
          ans[0] = 2. * u - aa3;
          ans[1] =    - u - aa3;
          return;
      }
    }
    else
    {
        if ( D < 0. ){
            // three real solutions

            double phi = ONE_THIRD * acos ( -q / sqrt ( -p3  ) );
            double t   = 2. * sqrt ( -p );

            ans.resize(3);
            ans[0] =  t * cos ( phi       ) - aa3;
            ans[1] = -t * cos ( phi + PI3 ) - aa3;
            ans[2] = -t * cos ( phi - PI3 ) - aa3;
            return;
        }
        else{
            // one real solution

            ans.resize(1);
            double sqrtD = sqrt ( D );
            double u     = pow  ( sqrtD + jem::numeric::abs(q), ONE_THIRD );
            if ( q > 0. ){
                ans[0] = -u + p / u - aa3;
            }
            else{
                ans[0] = u - p / u - aa3;;
            }
            return;
        }
    }
  }
}

// ---------------------------------------------------------------------
//   readCracks
// ---------------------------------------------------------------------
// The file is something like this
// Cracks
// 20 (number of cracks)
// 2  (2D or 3 for 3D)
// x1 y1 [z1] x2 y2 [z2]
// x1 y1 x2 y2
// ...

void            readCracks

   (       Matrix&                   cracks,
     const String&                   fName )
{
  using  boost::lexical_cast;
  using  jem::String; 
  using  namespace jem;

  std::ifstream file ( std::string(fName.begin(), fName.end()).c_str(), 
                       std::ios::in );
  
  if ( !file ) 
  {
    std::cerr << "Unable to open crack file!!!\n\n";
    exit(1);
  }

  std::string               line;
  std::vector<std::string>  splitLine;

  getline ( file, line ); // read line "Cracks"
  getline ( file, line ); // read no of cracks
  
  const idx_t crackCount = lexical_cast<int> ( line );

  getline ( file, line ); // read number of dimensions
  const idx_t dim        = lexical_cast<int> ( line );

  cracks.resize ( crackCount, dim*2 );

  // read cracks, one by one

  for ( idx_t ie = 0; ie < crackCount; ie++ )
  {
    getline      ( file, line ); 
    boost::split ( splitLine, line, boost::is_any_of("\t ") );

    cracks(ie,0) = lexical_cast<double> ( splitLine[0] );
    cracks(ie,1) = lexical_cast<double> ( splitLine[1] );
    cracks(ie,2) = lexical_cast<double> ( splitLine[2] );
    cracks(ie,3) = lexical_cast<double> ( splitLine[3] );
  }
}

// ---------------------------------------------------------------------
//   readInterfaceMesh
// ---------------------------------------------------------------------


void readInterfaceMesh 

  ( int&              elemCount,
    int&              nodeCount,
    IntMatrix&        ielems,
    std::map<int,std::vector<int> >& matMap,
    std::vector<std::vector<int> >&  dupNodes,
    const Properties& prop,
    const Properties& conf,
    std::vector<int>&        oppositeVertices )
{
  using  boost::lexical_cast;
  using  jem::String; using namespace jem;

  String   meshFile;
  
  prop. get ( meshFile,  "meshFile" );
  conf. set ( "meshFile", meshFile  );

  std::ifstream file ( std::string(meshFile.begin(), meshFile.end()).c_str(), 
                       std::ios::in );
  
  if ( !file ) 
  {
    std::cerr << "Unable to open mesh file!!!\n\n";
    exit(1);
  }
  else
  {
    std::cout << "Reading interface elements from file " << std::string(meshFile.begin(), meshFile.end()).c_str() << "\n";
  }

  std::string               line;
  std::vector<std::string>  splitLine;
  std::vector<int>          connectivity;

  getline ( file, line ); // read line "Element"
  getline ( file, line ); // read no of elements
  elemCount = lexical_cast<int> ( line );
  getline ( file, line ); // read no of nodes/element
  nodeCount = lexical_cast<int> ( line );

  ielems.          resize ( elemCount, nodeCount );
  
  IntVector                 bulk1s(elemCount);
  IntVector                 bulk2s(elemCount);
 
  idx_t  matId;

  // read element

  for ( idx_t ie = 0; ie < elemCount; ie++ )
  {
    getline      ( file, line ); 
    boost::split ( splitLine, line, boost::is_any_of("\t ") );

    //id         = lexical_cast<int> ( splitLine[0] );
    matId      = lexical_cast<int> ( splitLine[1] );
    
    bulk1s[ie] = lexical_cast<int> ( splitLine[2] );
    bulk2s[ie] = lexical_cast<int> ( splitLine[3] );
    

    std::transform ( splitLine.begin()+4, 
	                   splitLine.end  ()-1, 
	                   std::back_inserter ( connectivity ), 
	    	             Str2IntFunctor() );

    for ( idx_t in = 0; in < connectivity.size(); in++ )
    {
      ielems(ie,in) = connectivity[in];
    }

    connectivity.clear ();

    matMap[matId].push_back ( ie );
  }

  // read node

  getline ( file, line );
  getline ( file, line );

  int noOfNodes = lexical_cast<int> ( line );

  dupNodes.resize ( noOfNodes );

  for ( idx_t in = 0; in < noOfNodes; in++ )
  {
    getline ( file, line ); 
    boost::split ( splitLine, line, boost::is_any_of("\t ") );

    std::transform ( splitLine.begin(), 
	                 splitLine.end  ()-1, 
	                 std::back_inserter ( dupNodes[in] ), 
	    	         Str2IntFunctor() );
  }

  // read opposite vertices if any (3D)

  getline ( file, line );

  if ( line == "OppositeVertices")
  {
    int oppVertex;

    for ( int in = 0; in < elemCount; in++ )
    {
      file >> oppVertex;

      oppositeVertices.push_back ( oppVertex );
    }
  }
  
  conf. set ( "bulk1s", bulk1s );
  conf. set ( "bulk2s", bulk2s );
    
  std::cout << "Reading interface elements: DONE\n ";
}
  
// ---------------------------------------------------------------------
//   readGeneralDGInput
// ---------------------------------------------------------------------


void readGeneralDGInput

  ( IntMatrix&        ielems,
    const Properties& prop,
    const Properties& conf )
{
  using  boost::lexical_cast;
  using  jem::String; using namespace jem;

  String   meshFile;
  
  prop. get ( meshFile,  "meshFile" );
  conf. set ( "meshFile", meshFile  );

  std::ifstream file ( std::string(meshFile.begin(), meshFile.end()).c_str(), 
                       std::ios::in );
  
  if ( !file ) 
  {
    std::cerr << "Unable to open mesh file!!!\n\n";
    exit(1);
  }
  else
  {
    std::cout << "Reading general DG input from file " << std::string(meshFile.begin(), meshFile.end()).c_str() << "\n";
  }

  std::string               line;
  std::vector<std::string>  splitLine;
  std::vector<int>          connectivity;

  getline ( file, line ); // read line "Element"
  getline ( file, line ); // read no of elements
  int elemCount = lexical_cast<int> ( line );
  getline ( file, line ); // read no of nodes/element
  int nodeCount = lexical_cast<int> ( line );

  ielems.          resize ( elemCount, nodeCount );
  
  IntVector                bulks(elemCount);
 
  //int id, matId;

  // read element

  for ( idx_t ie = 0; ie < elemCount; ie++ )
  {
    getline      ( file, line ); 
    boost::split ( splitLine, line, boost::is_any_of("\t ") );

    bulks[ie] = lexical_cast<int> ( splitLine[0] );

    std::transform ( splitLine.begin()+1, 
	                 splitLine.end  ()-1, 
	                 std::back_inserter ( connectivity ), 
	   	             Str2IntFunctor() );

    for ( idx_t in = 0; in < connectivity.size(); in++ )
    {
      ielems(ie,in) = connectivity[in];
    }

    connectivity.clear ();
  }

  conf. set ( "bulks",      bulks     );
  conf. set ( "inodeCount", nodeCount );

}

// ---------------------------------------------------------------------
//   readDiscreteInterfaceMesh
// ---------------------------------------------------------------------


void readDiscreteInterfaceMesh 

  ( int&                             elemCount,
    int&                             nodeCount,
    IntMatrix&                       ielems,
    std::map<int,std::vector<int> >& matMap,
    std::vector<std::vector<int> >&  dupNodes,
    const Properties&                prop,
    const Properties&                conf,
    std::vector<int>&                oppositeVertices )
{
  using  boost::lexical_cast;
  using  jem::String; 
  using namespace jem;

  String   meshFile;
  
  prop. get ( meshFile,  "meshFile" );
  conf. set ( "meshFile", meshFile  );

  std::ifstream file ( std::string(meshFile.begin(), meshFile.end()).c_str(), 
                       std::ios::in );
  
  if ( !file ) 
  {
    std::cerr << "Unable to open mesh file!!!\n\n";
    exit(1);
  }

  std::string               line;
  std::vector<std::string>  splitLine;
  std::vector<int>          connectivity;

  getline ( file, line );
  getline ( file, line );
  elemCount = lexical_cast<int> ( line );
  getline ( file, line );
  nodeCount = lexical_cast<int> ( line );

  ielems.          resize ( elemCount, nodeCount );
  oppositeVertices.resize ( elemCount );
 
  int id;
  int matId;

  // read element

  for ( int ie = 0; ie < elemCount; ie++ )
  {
    getline      ( file, line ); 
    boost::split ( splitLine, line, boost::is_any_of("\t ") );

    id    = lexical_cast<int> ( splitLine[0] );
    matId = lexical_cast<int> ( splitLine[1] );
    
    connectivity.push_back (lexical_cast<int>(splitLine[2]));
    connectivity.push_back (lexical_cast<int>(splitLine[3]));

    for ( int in = 0; in < connectivity.size(); in++ )
    {
      ielems(ie,in) = connectivity[in];
    }

    connectivity.clear ();

    matMap[matId].push_back ( ie );
  }

  // read node

  getline ( file, line );
  getline ( file, line );

  int noOfNodes = lexical_cast<int> ( line );

  dupNodes.resize ( noOfNodes );

  for ( int in = 0; in < noOfNodes; in++ )
  {
    getline ( file, line ); 
    boost::split ( splitLine, line, boost::is_any_of("\t ") );

    std::transform ( splitLine.begin(), 
	             splitLine.end  ()-1, 
	             std::back_inserter ( dupNodes[in] ), 
	   	     Str2IntFunctor() );
  }

  // read opposite vertices if any (3D)

  getline ( file, line );

  if ( line == "OppositeVertices")
  {
    int oppVertex;

    for ( int in = 0; in < elemCount; in++ )
    {
      file >> oppVertex;

      oppositeVertices.push_back ( oppVertex );
    }
  }
}

// -------------------------------------------------------
//   updateCoord2D
// -------------------------------------------------------

void            updateCoord2D 

  ( Matrix&       uCoord,
    const Matrix& coord, 
    const Vector& disp )
{
  JEM_PRECHECK ( uCoord.size(0) == 2 &&
                 coord.size(0)  == 2 );

  uCoord(0,ALL) = coord(0,ALL) + disp[slice(0,END,2)];
  uCoord(1,ALL) = coord(1,ALL) + disp[slice(1,END,2)];
  
  //const int nnode = uCoord.size(1);
  //for (int i = 0; i < nnode; ++i )
  //{
  //  uCoord(0,i) = coord(0,i) + disp[2*i];
  //  uCoord(1,i) = coord(1,i) + disp[2*i+1];
  //}
}


// -------------------------------------------------------
//   updateCoord3D
// -------------------------------------------------------

void            updateCoord3D 

  ( Matrix&       uCoord,
    const Matrix& coord, 
    const Vector& disp )
{
  JEM_PRECHECK ( uCoord.size(0) == 3 &&
                 coord.size(0)  == 3 );

  uCoord(0,ALL) = coord(0,ALL) + disp[slice(0,END,3)];
  uCoord(1,ALL) = coord(1,ALL) + disp[slice(1,END,3)];
  uCoord(2,ALL) = coord(2,ALL) + disp[slice(2,END,3)];
}


UpdateCoordFunc     getUpdateCoordFunc

  ( int                 rank )
{
  if     ( rank == 2 )
  {
    return & updateCoord2D;
  }
  else
  {
    return & updateCoord3D;
  }
}

//----------------------------------------------------------------------------
//   matrixFromTable
//----------------------------------------------------------------------------

void  matrixFromTable

  ( Matrix&              mat,
    const Table&         t,
    const StringVector&  colNames )

{
     Array<idx_t>    jcols    ( colNames.size() );
     Array<idx_t>    irows    ( t.rowCount() );
     idx_t           i; 

     const int rowCount = irows.size();
     const int colCount = jcols.size();
     
     mat.resize ( rowCount, colCount );
     mat = 0.0;

     for ( i = 0; i < rowCount; i++ )
     {
       irows[i] = i;
     }

     for ( i = 0; i < colCount; i++ )
     {
       jcols[i] = t.getColumnIndex ( colNames[i] );
     }

     t.findBlock ( mat, irows, jcols );

     mat.transpose();
}

void  matrixFromTable

  ( IntMatrix&     mat,
    const Table&   t,
    const String&  cols )

{
     using jem::util::StringUtils;

     StringVector colNames ( StringUtils::split( cols ) );
     Array<idx_t> jcols    ( colNames.size() );
     Array<idx_t> irows    ( t.rowCount() );
     Matrix       temp     ( irows.size(), jcols.size() );
     int          i;

     temp = 0.0;
     mat.resize ( irows.size(), jcols.size() );

     for ( i = 0; i < irows.size(); i++ )
     {
        irows[i] = i;
     }

     for ( i = 0; i < jcols.size(); i++ )
     {
        jcols[i] = t.getColumnIndex ( colNames[i] );
     }

     t.getBlock ( temp, irows, jcols );

     temp.transpose();

     //mat = jem::castTo<int> ( temp );
}

  void matrixFromTable

    ( Vector&        vec,
          const Table&   t,
              const String&  col )

{
     int  rowCount = t.rowCount();
     int  index;
     int  i;

     vec.resize ( rowCount );
     index = t.getColumnIndex ( col );

     vec = 0.0;

     for ( i = 0; i < rowCount; i++ )
     {
        vec[i] = t.getValue(i,index);
     }
}

void matrixFromTable

  ( IntVector&     vec,
    const Table&   t,
    const String&  col )

{
     int          rowCount = t.rowCount();
     int          index;
     int          i;

     vec.resize ( rowCount );
     index = t.getColumnIndex ( col );

     vec = 0;

     for ( i = 0; i < rowCount; i++ )
     {
        vec[i] = (int) t.getValue(i,index);
     }
}

// Vector cross (const Vector& a, const Vector& b)
// {
//    Vector v(3);

//    v[0] = a[1]*b[2] - a[2]*b[1];
//    v[1] = a[2]*b[0] - a[0]*b[2];
//    v[2] = a[0]*b[1] - a[1]*b[0];

//    return v;
// }

// -----------------------------------------------------------------
//  writeMeshToVTKFile_ 
// -----------------------------------------------------------------


void  writeMeshToVTKFile_    
       
     ( Ref<PrintWriter> vtuFile,
       const Matrix&    coords,
       const IntMatrix& connec )
{

   using jem::io::endl;

   const int nodeCount = coords.size(1);
   const int elemCount = connec.size(0);
   const int nodeElem  = connec.size(1);

   const String str1 = String::format("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"> \n"
                                      "<UnstructuredGrid> \n"
                                      "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\"> \n", 
                                      nodeCount, elemCount );

   *vtuFile << str1;

   // write node coordinates

   *vtuFile << "<Points> \n";
   *vtuFile << "<DataArray  type=\"Float64\"  NumberOfComponents=\"3\" format=\"ascii\" >\n";

   for ( int i = 0; i < nodeCount; ++i )
   {
     *vtuFile << coords(0,i) << " " << coords(1,i) << " " << coords(2,i) << '\n';
   }

   *vtuFile << "</DataArray>\n";
   *vtuFile << "</Points>\n";

   // write element connectivity

   IntVector inodes(nodeElem);

   *vtuFile << "<Cells> \n";
   *vtuFile << "<DataArray  type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

   for ( int ie = 0; ie < elemCount; ie++ )
   { 
       for ( int in = 0; in < nodeElem; in++ ){
          *vtuFile << connec(ie,in) << " "; 
       }

       *vtuFile << endl;
   }
   *vtuFile << "</DataArray>\n";

   // write cell offset

   *vtuFile << "<DataArray  type=\"Int32\"  Name=\"offsets\"  format=\"ascii\"> \n";

   int offset = 0;
   for ( int ie = 0; ie < elemCount; ie++ )
   {
      offset += 4;
      *vtuFile <<  offset << endl;
   }

   *vtuFile << "</DataArray>\n";

   // Print cell types
                                     
   *vtuFile << "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\"> \n";
                                       
   for ( int ie = 0; ie < elemCount; ie++ )
   {
     *vtuFile <<  9 << endl;
   }

  *vtuFile << "</DataArray> \n </Cells> \n";
}


void  writeMatrixToVTKFile_    
       
     ( Ref<PrintWriter> vtuFile,
       const Matrix&    matrix,
       const String&    name )
{
   const int colCount = matrix.size(1);
   //const int rowCount = matrix.size(0);

   *vtuFile << String::format(" <DataArray  type=\"Float64\"  Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\"> \n", name);

   for ( int i = 0; i < colCount; ++i )
   {
     *vtuFile << matrix(0,i) << " " << matrix(1,i) << " " << matrix(2,i) << '\n';
   }

   *vtuFile << "</DataArray> \n";
}

double mcauleyP ( double x ) 
{
  return 0.5*( x + std::abs (x) );
}

double mcauleyM ( double x ) 
{
  return 0.5*( x - std::abs (x) );
}

int sign ( double x ) 
{
  return (x < 0.) ? 0 : (x > 0.);
}

GetIandPFunc      getIandPFunc
  
  ( int                rank )
{
  if      ( rank == 2 ) 
  {
    return &getIandP2DFunc;
  }
  else if ( rank == 3 ) 
  {
    return &getIandP3DFunc;
  }
}

void                  getIandP2DFunc

  ( const Matrix&      P,
    const Vector&      I )
{
  I[0] = 1.; 
  I[1] = 1.; 
  //I[2] = 1.; 

  P(0,0) =  0.5;
  P(0,1) = -0.5;
  P(1,0) = -0.5;
  P(1,1) =  0.5;
  //P(2,0) = -1/3.;
  //P(2,1) = -1/3.;
  //P(2,2) = 2/3.;
  P(3,3) = 0.5;
}

void                  getIandP3DFunc

  ( const Matrix&      P,
    const Vector&      I )
{
  I[0] = 1.; 
  I[1] = 1.; 
  I[2] = 1.; 

  P(0,0) =  TWO_THIRD; P(0,1) = -ONE_THIRD; P(0,2) = -ONE_THIRD;
  P(1,0) = -ONE_THIRD; P(1,1) =  TWO_THIRD; P(1,2) =  -ONE_THIRD;
  P(2,0) = -ONE_THIRD; P(2,1) = -ONE_THIRD; P(2,2) =  TWO_THIRD;
  P(3,3) = 0.5;
  P(4,4) = 0.5;
  P(5,5) = 0.5;
}

void   computeGeometricFunctionMarigo

      ( double&  func,
        double&  firstDer,
        double&  secondDer,
        double   value )
{
  func      = value;
  firstDer  = 1.0;
  secondDer = 0.0;
}

void   computeGeometricFunctionAT

      ( double&  func,
        double&  firstDer,
        double&  secondDer,
        double   value )
{
  func      = value * value;
  firstDer  = 2.0 * value;
  secondDer = 2.0;
}

void   computeGeometricFunctionWu

      ( double&  func,
        double&  firstDer,
        double&  secondDer,
        double   value )
{
  func      = 2.0 * value - value*value;
  firstDer  =  2.0 - 2.0 * value;
  secondDer = -2.0;
}

GeometricFunction getGeometricFunc

  ( const String& type )
{
  if      ( type == "AT" ){
    return & computeGeometricFunctionAT;
  }
  else if ( type == "Marigo" ){
    return & computeGeometricFunctionMarigo;
  }
  else if ( type == "Wu" ){
    return & computeGeometricFunctionWu;
  }
  else
  {
      using namespace jem;
      throw Error (
      JEM_FUNC,
      " invalid geometric function !!! "
      );  
  }
}

/*
void solve3x3Eigen 

  ( const Vector& evals,
    const Matrix& evecs,
    const Vector& strain )

{
  const idx_t rank = strain.size ( );

  double exx, eyy, ezz, exy, exz(0.), ezy(0.);

  if ( rank == 6 )
  {
    exx = strain[0];
    eyy = strain[1];
    ezz = strain[2];
    exy = strain[3] / 2.;
    exz = strain[4] / 2.;
    ezy = strain[5] / 2.;
  }
  else 
  {
    exx = strain[0];
    eyy = strain[1];
    ezz = strain[2];
    exy = strain[3] / 2.;
  };

  // deviatoric matrix, only the diagonal terms

  double dev11 = TWO_THIRD * exx - ONE_THIRD * ( eyy + ezz );
  double dev22 = TWO_THIRD * eyy - ONE_THIRD * ( exx + ezz );
  double dev33 = TWO_THIRD * ezz - ONE_THIRD * ( eyy + exx );

  // (minus) second and third invariants of the deviatoric matrix

  double J2 = 0.5 * ( dev11 * dev11 + dev22 * dev22 + dev33 * dev33 + 2. * ( exy * exy +   exz * exz + ezy * ezy ) ); 
  double J3 = dev11 * dev22 * dev33 + 2. * exy * ezy * xz - exz * eyy * exz - ezy * e zy * exx - ezz * exy * exy;

  // solve for alpha

  double cos3alpha = 0.5 * J3 * pow (  3./J2 , 1.5  );
  double alpha     = ONE_THIRD * acos ( cos3alpha );

  // find the distinct eigenvalue

  double fac = 2. * sqrt ( J2/3 );
  
  if ( alpha < PI6 )
  {
    eta1 = fac * cos ( alpha );
  }
  else
  {
  }
}
  */

void computeEigenValues2D 

( const Vector& stressP,
  const Vector& stressN )
{
  double sigmaxx = stressN[0];
  double sigmayy = stressN[1];
  double sigmaxy = stressN[3];

  double delta   = pow ( sigmaxx - sigmayy, 2 ) + 4. * sigmaxy * sigmaxy;

  stressP[0]     = 0.5 * ( sigmaxx + sigmayy + sqrt ( delta ) );
  stressP[1]     = 0.5 * ( sigmaxx + sigmayy - sqrt ( delta ) );
  stressP[2]     = 0.0;
}


// --------------------------------------------------------------
//  computeMaxPrincipleStress (1D, 2D and 3D versions)
// --------------------------------------------------------------

double computeMaxPrincipleStress1D 

( const Vector& stress )
{
  return stress[0];
}

double computeMaxPrincipleStress2D 

( const Vector& stress )
{
  double sigmaxx = stress[0];
  double sigmayy = stress[1];
  double sigmaxy = stress[3];

  double delta   = sqrt ( pow ( sigmaxx - sigmayy, 2 ) + 4. * sigmaxy * sigmaxy );

  return 0.5 * ( sigmaxx + sigmayy + delta );
}

double computeMaxPrincipleStress3D 

( const Vector& stress )
{
  double sigmaxx = stress[0];
  double sigmayy = stress[1];
  double sigmazz = stress[2];
  double sigmaxy = stress[3];
  double sigmayz = stress[4];
  double sigmaxz = stress[5];

  double q  =  ( sigmaxx + sigmayy + sigmazz ) * .3333333333333;
  double p1 =  sigmaxy*sigmaxy + sigmaxz*sigmaxz + sigmayz*sigmayz;
  double p2 = pow(sigmaxx - q,2) + pow(sigmayy - q,2) + pow(sigmazz - q,2) + 2*p1;
  double p  =  sqrt(p2/6);

  // B  = (1/p)*(sigma_bar - q*Identity(3))

  Tuple<double,3,3> B;
  B(0,0) = sigmaxx-q; B(0,1) = B(1,0) = sigmaxy;
  B(1,1) = sigmayy-q; B(0,2) = B(2,0) = sigmaxz;
  B(2,2) = sigmazz-q; B(1,2) = B(2,1) = sigmayz;

  if ( abs(p) < 1e-15 ) 
  {
    B = 0.;
  } 
  else
  {
    B /= p; 
  }

  double r  = jem::numeric::det(B)/2;

  double phi;
    
  // In exact arithmetic for a symmetric matrix  -1 <= r <= 1
  // but computation error can leave it slightly outside this range.
  if      (r <= -1) 
      phi = 1.0471975511965976;
  else if (r >= 1)
      phi = 0;
  else
      phi = acos(r)/3;
    
  // the maximum eigenvalue
  
  return  q + 2*p*cos(phi);
  //ps2 = q + 2*p*cos(phi + 2*math.pi/3)
  //ps3 = I1 - ps1 - ps2

}


// --------------------------------------------------------------
//  getMaxPrincipalStressFunc
// --------------------------------------------------------------

ComputeMaxPrincipleStress        getMaxPrincipalStressFunc

( int                 rank )
{
  if           ( rank == 1 ) return &computeMaxPrincipleStress1D;
  else if      ( rank == 2 ) return &computeMaxPrincipleStress2D;
  else                       return &computeMaxPrincipleStress3D;
}

// --------------------------------------------------------------
//  computeNormOfSym2ndTensor (1D, 2D, 3D versions)
// --------------------------------------------------------------

double computeNormOfSym2ndTensor1D

  ( const Vector& stress )
{
  double v1 = stress[0];

  return sqrt ( v1 * v1 );
} 

double computeNormOfSym2ndTensor2D

  ( const Vector& stress )
{
  double v1 = stress[0];
  double v2 = stress[1];
  double v3 = stress[2];
  double v4 = stress[3];

  return sqrt ( v1 * v1 + v2 * v2 + v3 * v3 + 2.0 * v4 * v4 );
} 

double computeNormOfSym2ndTensor3D

  ( const Vector& stress )
{
  double v1 = stress[0];
  double v2 = stress[1];
  double v3 = stress[2];
  double v4 = stress[3];
  double v5 = stress[4];
  double v6 = stress[5];

  return sqrt ( v1 * v1 + v2 * v2 + v3 * v3 + 2.0 * ( v4 * v4 + v5 * v5 + v6 * v6) );
}    


// --------------------------------------------------------------
//  computeNormOfSym2ndTensor
// --------------------------------------------------------------
// 
ComputeNormOfSym2ndTensor        getNormOfSym2ndTensorFunc

( int                 rank )

{
  
  if        ( rank == 1 )
  {
    return &computeNormOfSym2ndTensor1D;
  }
  else if   ( rank == 2 ) 
  {
    return &computeNormOfSym2ndTensor2D;
  }
  else
  {
    return &computeNormOfSym2ndTensor3D;
  }
}


