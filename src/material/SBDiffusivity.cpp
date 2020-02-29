/*
 *  Copyright (C) 2016-2017 Monash University. All rights reserved.
 *
 *  Author: V.P. Nguyen, phu.nguyen@monash.edu
 *  Date started: 6 September 2016
 *  Status:
 *  -- 9 December 2017: 
 */

#include <jem/base/System.h>
#include <jem/base/limits.h>
#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/sparse/utilities.h>

#include "SBDiffusivity.h"

using namespace jem;
using namespace jem::io;
using jem::numeric::matmul;
using jem::numeric::dotProduct;


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


SBDiffusivity::SBDiffusivity 

    ( const Properties&     conf,
      const Properties&     props,
      const Properties&     globdat )

    : Diffusivity ( conf, props, globdat )
{
  props.get ( ft_, "ft" );
  conf .set ( "ft", ft_ );

  fac_  = c0_ / ft_ / ft_;
  fac2_ = 2. * fac_;
}


SBDiffusivity::~SBDiffusivity ()
{}


//-----------------------------------------------------------------------
//   getValues
//-----------------------------------------------------------------------

void SBDiffusivity::getValues

    ( const Properties&       params )

{
  Vector stress, cVec;

  params.get ( cVec,      "cVec"   );
  params.get ( stress,    "stress" );
  
  // 2D problems: stress = [sigmaXX sigmaYY sigmaZZ sigmaXY]

  double sigmaXX = stress[0];
  double sigmaYY = stress[1];
  double sigmaXY = stress[3];

  cVec[0] = fac_ * ( sigmaXX * sigmaXX   + sigmaXY * sigmaXY  );
  cVec[1] = fac_ * ( sigmaYY * sigmaYY   + sigmaXY * sigmaXY  );
  cVec[2] = fac_ * ( sigmaXX + sigmaYY ) * sigmaXY;

}



//-----------------------------------------------------------------------
//   getValuesAndDers
//-----------------------------------------------------------------------

void SBDiffusivity::getValuesAndDers

    ( const Properties&       params )

{
  Vector stress, cVec;
  Vector dcxx, dcxy, dcyy;

  params.get ( cVec,      "cVec"   );
  params.get ( stress,    "stress" );
  params.get ( dcxx,      "dcxx"   );
  params.get ( dcxy,      "dcxy"   );
  params.get ( dcyy,      "dcyy"   );

  double sigmaXX = stress[0];
  double sigmaYY = stress[1];
  double sigmaXY = stress[3];

  cVec[0] = fac_ * ( sigmaXX * sigmaXX   + sigmaXY * sigmaXY  );
  cVec[1] = fac_ * ( sigmaYY * sigmaYY   + sigmaXY * sigmaXY  );
  cVec[2] = fac_ * ( sigmaXX + sigmaYY ) * sigmaXY;

  // derivatives of c_{ij} w.r.t stress 
  // 2D problems: stress = [sigmaXX sigmaYY sigmaZZ sigmaXY]

  dcxx[0] = fac2_ * sigmaXX;
  dcxx[1] = 0.;
  dcxx[2] = 0.;
  dcxx[3] = fac2_ * sigmaXY;

  dcyy[0] = 0.;
  dcyy[1] = fac2_ * sigmaYY;
  dcyy[2] = 0.;
  dcyy[3] = fac2_ * sigmaXY;
  
  dcxy[0] = fac_ * sigmaXY;
  dcxy[1] = fac_ * sigmaXY;
  dcxy[2] = 0.;
  dcxy[3] = fac_ * ( sigmaXX + sigmaYY );
}

