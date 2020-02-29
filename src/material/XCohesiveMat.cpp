
#include <jem/util/Properties.h>
#include <jem/base/Error.h>

#include "XCohesiveMat.h"
#include "TuronXCohesiveMat.h"
#include "TuronXCohesiveWeibulMat.h"
#include "ElasticExponentialCohesiveMat.h"
#include "RigidExpCohesiveMat.h"

using namespace jem;

//=======================================================================
//   class XCohesiveMat
//=======================================================================

XCohesiveMat::XCohesiveMat 

  ( const int          rank,
    const Properties&  globdat )

  : CohesiveMaterial ( rank, globdat )

{}

XCohesiveMat::~XCohesiveMat()
{}

//-----------------------------------------------------------------------
//   elasticUpdate
//-----------------------------------------------------------------------

void XCohesiveMat::elasticUpdate

  ( Vector&               traction,
    Matrix&               stiff,
    const Vector&         jump )

{
  throw Error ( JEM_FUNC, 
    "pure virtuality of elasticUpdate is removed in XCohesiveMat" );
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void XCohesiveMat:: update

    ( Vector&             traction,
      Matrix&             stiff,
      const Vector&       jump0,
      const Vector&       djump,
      int                 mpoint )
{ }
  
//-----------------------------------------------------------------------
//   initShift
//-----------------------------------------------------------------------

void XCohesiveMat::initShift 

    ( const int             ipoint,
      const Vector&         traction )
{ }

//-----------------------------------------------------------------------
//   giveTensileStrength
//-----------------------------------------------------------------------

double XCohesiveMat::giveTensileStrength () const 
{
    return 0.;
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newXCohesiveMat
//-----------------------------------------------------------------------


Ref<XCohesiveMat>   newXCohesiveMat

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  Properties     matProps = props.getProps ( name );
  Properties     matConf  = conf.makeProps ( name );

  Ref<XCohesiveMat> mat;
  String                 type;
  int                    dim;

  matProps.get ( type, "type" );
  matConf .set ( "type", type );

  matProps.get ( dim, "dim"   );
  matConf .set ( "dim", dim   );
  
  if      ( type == "TuronX" )
  {
    mat = newInstance<TuronXCohesiveMat>             ( dim, globdat );
  }
  else if ( type == "TuronXWeibul" )
  {
    mat = newInstance<TuronXCohesiveWeibulMat>       ( dim, globdat );
  }
  else if ( type == "ElasticExponential" )
  {
    mat = newInstance<ElasticExponentialCohesiveMat> ( dim, globdat );
  }
  else if ( type == "RigidExponential" )
  {
    mat = newInstance<RigidExpCohesiveMat>          ( dim, globdat );
  }
  else
  {
    matProps.propertyError (
      name,
      "invalid xcohesive material type: " + type
    );
  }

  return mat;
}
