
#include <jem/base/Error.h>

#include "CohesiveMaterial.h"

using jem::Error;

CohesiveMaterial::CohesiveMaterial 

  ( const int          rank,
    const Properties&  globdat )

  : rank_ ( rank ), desparateMode_ (false)

{}

CohesiveMaterial::~CohesiveMaterial()
{}

void   CohesiveMaterial::commit()
{}

void   CohesiveMaterial::cancel()
{}

void   CohesiveMaterial::configure

  ( const Properties&     props,
    const Properties&     globdat )

{}

void   CohesiveMaterial::getConfig

  ( const Properties& props , 
    const Properties& globdat ) const

{}

void CohesiveMaterial::allocPoints ( int count, double dam )
{}

void CohesiveMaterial::deallocPoints ( int count )
{}
  
void CohesiveMaterial:: update

    ( Vector&             traction,
      Matrix&             stiff,
      const Vector&       jump0,
      const Vector&       djump,
      int                 mpoint )
{ }

//--------------------------------------------------------------------
//   isFailure
//--------------------------------------------------------------------

//bool CohesiveMaterial::isFailure 
//
//    ( const Vector& stress , 
//      const Ref<Material> bulkMat )  const
//{
//  return false;
//}

//--------------------------------------------------------------------
//   evalFailure
//--------------------------------------------------------------------

//double CohesiveMaterial::evalFailure 
//
//    ( const Vector& stress , 
//      const Ref<Material> bulkMat )  const
//{
//  // if evalFailure in not implemented, call isFailure
//
//  return double( isFailure ( stress, bulkMat ) );
//}
//
//bool CohesiveMaterial::evalFailure 
//
//    ( const Vector& traction)  const
//{
//  return 1;
//}

//--------------------------------------------------------------------
//   ipointCount
//--------------------------------------------------------------------

int  CohesiveMaterial::ipointCount () const
{
  return 0;
}

//--------------------------------------------------------------------
//   despair
//--------------------------------------------------------------------

bool CohesiveMaterial::despair ()
{
  desparateMode_ = true;

  int np = ipointCount();

  hasSwitched_.resize ( np );
  useSecant_  .resize ( np );

  useSecant_ = false;

  for ( int ip = 0; ip < np; ++ip )
  {
    hasSwitched_[ip] = ( isLoading(ip) != wasLoading(ip) );
  }

  return np > 0;
}

//--------------------------------------------------------------------
//   endDespair
//--------------------------------------------------------------------

bool CohesiveMaterial::endDespair ()
{
  desparateMode_ = false;

  return true;
}

