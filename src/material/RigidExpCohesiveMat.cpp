/*
 *
 */

#include <jem/base/Array.h>
#include <jem/base/array/operators.h>
#include <jem/base/array/tensor.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/algebra/matmul.h>

#include "RigidExpCohesiveMat.h"

using namespace jem;

using jem::numeric::norm2;
using jem::numeric::matmul;

//=======================================================================
//   class RigidExpCohesiveMat
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* RigidExpCohesiveMat::TENSILE_STRENTH_PROP = "ft";
const char* RigidExpCohesiveMat::FRACTURE_ENERGY_PROP = "gf";
const char* RigidExpCohesiveMat::SHEAR_STIFFNESS_PROP = "dint";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

RigidExpCohesiveMat::RigidExpCohesiveMat

  ( const int          rank,
    const Properties&  globdat )

  : CohesiveMaterial ( rank, globdat ),
    XCohesiveMat     ( rank, globdat )

{}

RigidExpCohesiveMat::~RigidExpCohesiveMat ()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void  RigidExpCohesiveMat::configure

  ( const Properties&     props,
    const Properties&     globdat )

{
  props.get ( ft0_ , TENSILE_STRENTH_PROP );
  props.get ( gf_  , FRACTURE_ENERGY_PROP );
  props.get ( dint_, SHEAR_STIFFNESS_PROP );
}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------

void  RigidExpCohesiveMat::getConfig

    ( const Properties&     conf,
      const Properties&     globdat ) const

{
  conf.set ( TENSILE_STRENTH_PROP, ft0_   );
  conf.set ( FRACTURE_ENERGY_PROP, gf_    );
  conf.set ( SHEAR_STIFFNESS_PROP, dint_  );
}


//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void  RigidExpCohesiveMat::update

  ( Vector&               traction,
    Matrix&               stiff,
    const Vector&         jump,
    int                   mpoint )

{
    double hisJump = preHist_.opening[mpoint];
    double curJump = jump[0];
    double sheJump = jump[1];
    double f       = curJump - hisJump;

    double loading;
    double jumpN ( hisJump );
    double temp;

    double ft        = ft_[mpoint];
    double ftOverGf  = ft/gf_;   
  
    if      ( jem::numeric::abs ( f ) < 1.e-10 )
    {
        loading = preHist_.loading[mpoint];
    }
    else if ( f > 0. )
    {
        loading = 1;
        jumpN    = curJump;
    }
    else
    {
        loading = 0;
        jumpN    = hisJump;
    }

    stiff = 0.;

    if ( loading )
    {
        temp        = ft * exp ( -ftOverGf * curJump );
        traction[0] = temp;
        traction[1] = dint_*sheJump;

        stiff(0,0)  = temp * ( - ftOverGf ); 
        stiff(1,1)  = dint_;

        damage_[mpoint] = 1. - temp / ft;
    }
    else
    {
        double tMax = ft * exp ( -ftOverGf * hisJump );
        traction[0] = (tMax / hisJump) * curJump;
        traction[1] = dint_*sheJump;
        
        stiff(0,0)  = (tMax / hisJump);
        stiff(1,1)  = dint_;
    }

    if ( curJump < 0. )
    {
        traction[0] = 1e6*curJump;
        stiff(0,0)  =1e6;
        System::out() << "Negative jump\n";
    }

    newHist_.opening = jumpN;
    newHist_.loading = loading;

    /*
    System::out() << "opening : " << curJump << "\n";
    System::out() << "loading : " << loading << "\n";
    System::out() << "traction: " << traction << "\n";
    System::out() << "stiff   : " << stiff << "\n\n";
    */
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void  RigidExpCohesiveMat::update

  ( Vector&               traction,
    const Vector&         jump,
    int                   mpoint )

{
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void  RigidExpCohesiveMat::update

  ( Vector&               traction,
    Matrix&               stiff,
    const Vector&         jump0,
    const Vector&         djump,
    int                   mpoint )

{
  Vector tJump ( jump0 + djump );
  update ( traction, stiff, tJump, mpoint );
}


//-----------------------------------------------------------------------
//   elasticUpdate
//-----------------------------------------------------------------------

void  RigidExpCohesiveMat::elasticUpdate

  ( Vector&               traction,
    Matrix&               stiff,
    const Vector&         jump )

{
}

//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------

void  RigidExpCohesiveMat::allocPoints

  ( int      count,
    double   dam )

{
  // not yet allocated

  if ( preHist_.opening.size ( ) == 0 ) 
  {
    preHist_.opening     .resize ( count );  
    preHist_.loading     .resize ( count );
    preHist_.dissipation .resize ( count );

    newHist_.opening     .resize ( count );
    newHist_.loading     .resize ( count );
    newHist_.dissipation .resize ( count );

    preHist_.opening     = 0.;
    preHist_.loading     = 1.  ;
    preHist_.dissipation = 0. ;

    newHist_.opening     = 0.;
    newHist_.loading     = 1.  ;
    newHist_.dissipation = 0. ;

    damage_.resize ( count );
    ft_    .resize ( count );
  }

  // appending at the end of the vector

  else
  {
    for ( int i = 0; i < count; i++ )
    {
      preHist_.opening    .pushBack ( 0. );  
      preHist_.loading    .pushBack (  0.  );  
      preHist_.dissipation.pushBack (  0. );

      newHist_.opening    .pushBack ( 0. );  
      newHist_.loading    .pushBack (  0.  );  
      newHist_.dissipation.pushBack (  0. );
    }
  }
}

// --------------------------------------------------------------------
//  scaleDissipationIncrement
// --------------------------------------------------------------------

void  RigidExpCohesiveMat::scaleDissipationIncrement

  ( const double     factor,
    const int        ipoint )

{
  double dG  = newHist_.dissipation[ipoint] - preHist_.dissipation[ipoint];
         dG *= factor;

  if ( factor > 1. ) 
  {
    System::warn() << "scaleDissipationIncrement is scaling up" << endl;
  }

  newHist_.dissipation[ipoint] = preHist_.dissipation[ipoint] + dG ;
}

//-----------------------------------------------------------------------
//   evalFailure 
//-----------------------------------------------------------------------

double  RigidExpCohesiveMat::evalFailure

    ( const Vector&       sigmaN,
      int                 mpoint )  const

{
  double criterion = 1.;

  //System::out() << "B-K failure criterion: " << Fb2 << endl;

  return  criterion;
}

//-----------------------------------------------------------------------
//   deallocPoints
//-----------------------------------------------------------------------

void  RigidExpCohesiveMat::deallocPoints

  ( int      count )

{
}

//-----------------------------------------------------------------------
//   deallocPoints
//-----------------------------------------------------------------------

void RigidExpCohesiveMat:: initShift 

    ( const int             i, 
      const Vector&         traction )
{
    ft_[i] = traction[0];
}

//-----------------------------------------------------------------------
//   justTractionFree_ 
//-----------------------------------------------------------------------

bool  RigidExpCohesiveMat::justTractionFree

  ( const int   i ) const

{
  return false;
}

