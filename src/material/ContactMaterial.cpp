/*
 */

#include <jem/base/System.h>
#include <jem/base/array/utilities.h>
#include <jem/base/array/select.h>
#include <jem/base/array/operators.h>
#include <jem/base/array/tensor.h>
#include <jem/base/array/intrinsics.h>
#include <jem/io/PrintWriter.h>
#include <jem/util/Properties.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>

#include "ContactMaterial.h"

extern "C"
{
  #include <math.h>
}

using namespace jem;
using namespace jem::io;

using jem::numeric::matmul;
using jem::numeric::norm2;
using jem::numeric::MatmulChain;
using jem::util::Properties;
using jive::Vector;

//=======================================================================
//   class ContactMaterial
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  ContactMaterial::DUMMY_PROP    = "dummy";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------

ContactMaterial::ContactMaterial

  ( const int          rank,
    const Properties&  globdat )

  : CohesiveMaterial ( rank, globdat )

{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );
  dummy_       = 1.e6;
}


ContactMaterial::~ContactMaterial ()
{} 

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void ContactMaterial::configure

  ( const Properties& props,
    const Properties& globdat )

{
  using jem::maxOf;

  props.find ( dummy_, DUMMY_PROP,0., maxOf( dummy_) );
}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void ContactMaterial::getConfig

  ( const Properties& conf , 
    const Properties& globdat ) const

{
  conf.set ( DUMMY_PROP,  dummy_ );
}

//-----------------------------------------------------------------------
//   update (regular)
//-----------------------------------------------------------------------

void  ContactMaterial::update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump,
      int                   mpoint )
{
  bool   tension = jump[0] >= 0;

  // compute secant stiffness and traction

  stiff    = 0.;
  traction = 0.;

  System::out() << jump[0] << "\n";

  if (tension == false)
  {
    traction[0] = dummy_ * jump[0];
    stiff(0,0)  = dummy_;
  }
} 


void ContactMaterial::elasticUpdate

( Vector& traction, Matrix& stiff, const Vector& jump )
{
}

// end: update (regular)

//-----------------------------------------------------------------------
//   update (moonen material)
//-----------------------------------------------------------------------

void  ContactMaterial::update

    ( Vector&               traction,
      Matrix&               stiff,
      const Matrix&         jumpStiff,
      const Vector&         tEff,
      int                   mpoint )

{ }

//-----------------------------------------------------------------------
//   update (regular)
//-----------------------------------------------------------------------

void  ContactMaterial::update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump0,
      const Vector&         djump,
      int                   mpoint )
{
    Vector jump ( rank_ ); jump = jump0 + djump;
    update ( traction, stiff, jump, mpoint );
}

// --------------------------------------------------------------------
//  commit
// --------------------------------------------------------------------

void  ContactMaterial::commit()

{
}


// --------------------------------------------------------------------
//  allocPoints
// --------------------------------------------------------------------

void  ContactMaterial::allocPoints 

    ( int     count, 
      double  dam   )

{
}

// --------------------------------------------------------------------
//  allocPoints_
// --------------------------------------------------------------------

void  ContactMaterial::allocPoints_ 

    ( const int     count, 
      const double  dam,
      const int     loading)

{
}

// --------------------------------------------------------------------
//  deallocPoints
// --------------------------------------------------------------------

  void ContactMaterial::deallocPoints ( int count )
{
}

//-----------------------------------------------------------------------
//   justTractionFree_ 
//-----------------------------------------------------------------------

  bool  ContactMaterial::justTractionFree

  ( const int   i ) const

{
  return true;
}
