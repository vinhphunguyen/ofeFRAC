
#ifndef UTILITIES_H
#define UTILITIES_H



#include <jem/base/Array.h>
#include <jem/base/Tuple.h>
#include <jem/base/String.h>
#include <jem/base/Ref.h>

#include <jive/Array.h>
#include <jive/util/Table.h>

#include <functional>
#include <string>
#include <vector>
#include <map>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

namespace jem
{
  namespace util
  {
    class Properties;
  }

  namespace io
  {
    class PrintWriter;
  }
}

using jem::Ref;
using jive::Vector;
using jive::Matrix;
using jive::IntVector;
using jive::IntMatrix;
using jem::Tuple;
using jem::String;
using jem::util::Properties;
using jem::io::PrintWriter;
using jive::util::Table;
using jive::StringVector;



enum ProblemType {
       PlaneStrain,
       PlaneStress,
       AxiSymmetric
};


//-----------------------------------------------------------------------
//   typedefs
//-----------------------------------------------------------------------

typedef Tuple<double,6>   Vector3D;
typedef Tuple<double,3>   Vector2D;
typedef Tuple<double,2>   Vector1D;

// A pointer to a function that computes the spatial derivatives of
// the interpolation matrix. This is the so-called B-matrix.
// it points to the corresponding function for 1D, 2D and 3D case.

typedef void        (*ShapeGradsFunc)

  ( const Matrix&       b,
    const Matrix&       g );

// a pointer to a function that computes the shape function matrix N
// it points to the corresponding function for 1D, 2D and 3D case.  

typedef void        (*ShapeFunc)

  ( const Matrix&       sfuncs,
    const Vector&       n );

typedef void        (*FShapeFunc)

  ( const Matrix&       sfuncs,
    const Vector&       n );

typedef void        (*NormalMatrixFunc)

  ( const Matrix&       normalM,
    const Vector&       normalV );

typedef void        (*NormalInMatrixFunc)

  ( const Matrix&       normalM,
    const Vector&       normalV );

typedef void        (*UpdateCoordFunc)

  ( Matrix&              uCoord,
    const Matrix&        coord, 
    const Vector&        disp );

typedef void       (*TransformationMatrixFunc)

  ( Matrix&             Q,
    const Matrix&       grads,
    const Matrix&       coords );


typedef void       (*DoubleShearComponentsFunc)

  ( Vector&            v );

typedef bool       (*CheckAcousticTensorFunc)

  (       Vector&        n,
    const Matrix&        tangent  );  

typedef void       (*GetIandPFunc)
  
  ( const Matrix&       P,
    const Vector&       I );

// a pointer to a function that computes the value. 1st and 2nd derivative of
// the so-called crack geometric function alpha(d)
// it points to a concrete  function. Efficient than if else in run time

typedef void       (*GeometricFunction)

      ( double&  func,
        double&  firstDer,
        double&  secondDer,
        double   value );

typedef double (*ComputeHydroStress)

      ( double nu,
        const Vector& sigma );

// pointer to a function that compute the maximum principal stress
// it points to the corresponding function for 2D and 3D case.

typedef double (*ComputeMaxPrincipleStress)

      ( const Vector& stress );      

// pointer to a function that computes the norm of a symmetric second order tensor
// stored as a vector in Voigt notation

typedef double (*ComputeNormOfSym2ndTensor)    

      ( const Vector& vec );  

//-----------------------------------------------------------------------
//   constants
//-----------------------------------------------------------------------

// An integer array that maps the number of spatial dimensions (1, 2,
// or 3) to the number of strain/stress components.

extern const int      STRAIN_COUNTS[4];


//-----------------------------------------------------------------------
//   public functions
//-----------------------------------------------------------------------

// These functions compute the B-matrix given the matrix of shape function 1st 
// gradients (in global coordinates). For 1D, 2D (plane stress/strain) and 3D problems.

void                  get1DShapeGrads

  ( const Matrix&       b,
    const Matrix&       g );

void                  get2DShapeGrads

  ( const Matrix&       b,
    const Matrix&       g );

void                  get3DShapeGrads

  ( const Matrix&       b,
    const Vector&       g );

// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.

ShapeGradsFunc        getShapeGradsFunc

  ( int                 rank );

// ----------------------------------------------------------------------
//   compute hydrostatic stress
// ----------------------------------------------------------------------

double (computePStressHydroStress)

( double nu,
  const Vector& sigma );

double (computePStrainHydroStress)

( double nu,
  const Vector& sigma );

double (compute3DHydroStress)

( double nu,
  const Vector& sigma );

ComputeHydroStress    getHydroStressFunc

 ( ProblemType prob );

// -----------------------------------------------------------------------
//   matrix of shape functions  N
// -----------------------------------------------------------------------

void                  get1DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n );

void                  get2DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n );

void                  get3DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n );

void                  get2DFShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n );

// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.

ShapeFunc             getShapeFunc

  ( int                 rank );

FShapeFunc            getFShapeFunc

  ( int                 rank );

GetIandPFunc          getIandPFunc
  
  ( int                rank );

void                  getIandP2DFunc

  ( const Matrix&      P,
    const Vector&      I );

void                  getIandP3DFunc

  ( const Matrix&      P,
    const Vector&      I );

// -----------------------------------------------------------------------
//   matrix of normal vectors
// -----------------------------------------------------------------------
  
 void                  get2DNormalFuncs

  ( const Matrix&       normalM,
    const Vector&       normalV );
  
 void                  get3DNormalFuncs

  ( const Matrix&       normalM,
    const Vector&       normalV );

// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.

NormalMatrixFunc        getNormalFunc

  ( int                 rank );

// -----------------------------------------------------------------------
//   matrix of normal vectors (new version 2x4 matrix for 2D)
// -----------------------------------------------------------------------
  
 void                  get2DNormalInMatrixFuncs

  ( const Matrix&       normalM,
    const Vector&       normalV );
  
// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.

NormalInMatrixFunc      getNormalInMatrixFunc

  ( int                 rank );

void                  get2DTransMatrixFunc

  (       Matrix&       Q,
    const Matrix&       grads,
    const Matrix&       coords );

void                  get3DTransMatrixFunc

  (       Matrix&       Q,
    const Matrix&       grads,
    const Matrix&       coords );

TransformationMatrixFunc getTransMatrixFunc

  ( int                 rank );

// These functions compute the B-matrix given an interpolation matrix.

void                  get2DShearComponents

  ( Vector&             g );

void                  get3DShearComponents

  ( Vector&             g );

// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.

DoubleShearComponentsFunc        getDoubleShearComponentsFunc

  ( int                 rank );

// check determinant of acoustic tensor for discontinuous birfurcation
// analysis

bool                  checkAcousticTensor2D 

  (       Vector&       normal,
    const Matrix&       tangent );

bool                  checkAcousticTensor3D 

  (       Vector&       normal,
    const Matrix&       tangent );

CheckAcousticTensorFunc          getCheckAcousticTensorFunc

  ( int                  rank );
  

// Solving linear equation of one variable
// a*x + b = 0

void                  solveLinearEqua

  ( double&       ans,
    const double a,
    const double b );

// Solving quadratic equation
// a*x^2 + b*x + c = 0

void                  solveQuadEqua

  ( Vector&       ans,
    const double a,
    const double b,
    const double c );

// Solving cubic equation
// a*x^3 + b*x^2 + c*x + c = 0

void                  solveCubicEqua

  ( Vector&       ans,
    const double a,
    const double b,
    const double c,
    const double d);

// Vector                cross

//   ( const Vector& a,
//     const Vector& b );

// a functor to convert a string to int - 1
// used in std algorithm

struct Str2IntFunctor : public std::unary_function<std::string,int>
{
  int operator () ( const std::string& str  ) const
  {
    return boost::lexical_cast<int> ( str ) - 1;
  }
};

// read a mesh of interface elements
// both 2D and 3D are supported

void            readInterfaceMesh 

  ( int&              elemCount,
    int&              nodeCount,
    IntMatrix&        ielems,
    std::map<int,std::vector<int> >& matMap,
    std::vector<std::vector<int> >&  dupNodes,
    const Properties& prop,
    const Properties& conf,
    std::vector<int>& oppVertices );

void            readGeneralDGInput 

  ( IntMatrix&        ielems,
    const Properties& prop,
    const Properties& conf );

// read a mesh of discrete interface elements 
// i.e. spring elements 

void            readDiscreteInterfaceMesh 

  ( int&                             elemCount,
    int&                             nodeCount,
    IntMatrix&                       ielems,
    std::map<int,std::vector<int> >& matMap,
    std::vector<std::vector<int> >&  dupNodes,
    const Properties&                prop,
    const Properties&                conf,
    std::vector<int>&                oppVertices );

// read cracks from a file 
// this is for phase field modelling, define initial cracks via H variable

void            readCracks

   (       Matrix&                   cracks,
     const String&                   file );

void            updateCoord2D 

  ( Matrix&       uCoord,
    const Matrix& coord, 
    const Vector& disp );

void            updateCoord3D 

  ( Matrix&       uCoord,
    const Matrix& coord, 
    const Vector& disp );


UpdateCoordFunc     getUpdateCoordFunc

  ( int                 rank );

void  matrixFromTable

  ( Matrix&                       mat,
    const Table&                  t,
    const StringVector&           colNames );

void  matrixFromTable

  ( IntMatrix&                    mat,
    const Table&                  t,
    const String&                 cols );

void matrixFromTable

  ( Vector&                       vec,
    const Table&                  t,
    const String&                 col );

void matrixFromTable

  ( IntVector&                    vec,
    const Table&                  t,
    const String&                 col );

void  writeMeshToVTKFile_    
       
     ( Ref<PrintWriter>  vtuFile,
       const Matrix&     coords,
       const IntMatrix&  con );

void  writeMatrixToVTKFile_    
       
     ( Ref<PrintWriter>  vtuFile,
       const Matrix&     coords,
       const String&     name );

double mcauleyP ( double x );
double mcauleyM ( double x );
int    sign     ( double x );

// common crack geometric functions alpha(d)
// used in phase-field models

void   computeGeometricFunctionMarigo // -> alpha(d) = d

      ( double&  func,
        double&  firstDer,
        double&  secondDer,
        double   value );

void   computeGeometricFunctionAT     // -> alpha(d) = d^2

      ( double&  func,
        double&  firstDer,
        double&  secondDer,
        double   value );

void   computeGeometricFunctionWu     // -> alpha(d) = 2d-d^2

      ( double&  func,
        double&  firstDer,
        double&  secondDer,
        double   value );

GeometricFunction getGeometricFunc

  ( const String& type );

void   computeEigenValues2D 

  ( const Vector& stressP,
    const Vector& stressN );

// functions to compute the maximum principal stress
// for 1D, 2D and 3D  

double computeMaxPrincipleStress1D

  ( const Vector& stress ); 

double computeMaxPrincipleStress2D

  ( const Vector& stress );  

double computeMaxPrincipleStress3D

  ( const Vector& stress );  

// and the function pointer to select the correct one based on
// the problem dimension (rank)  

ComputeMaxPrincipleStress        getMaxPrincipalStressFunc

  ( int                 rank );  

double computeNormOfSym2ndTensor1D

  ( const Vector& stress ); 

double computeNormOfSym2ndTensor2D

  ( const Vector& stress );  

double computeNormOfSym2ndTensor3D

  ( const Vector& stress );      

ComputeNormOfSym2ndTensor        getNormOfSym2ndTensorFunc   

  ( int                 rank );  

#endif

