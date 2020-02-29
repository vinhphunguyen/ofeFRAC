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

#include "ConstantDiffusivity.h"

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


ConstantDiffusivity::ConstantDiffusivity 

    ( const Properties&     conf,
      const Properties&     props,
      const Properties&     globdat )

    : Diffusivity ( conf, props, globdat )
{
}


ConstantDiffusivity::~ConstantDiffusivity ()
{}


//-----------------------------------------------------------------------
//   getValues
//-----------------------------------------------------------------------

void ConstantDiffusivity::getValues

    ( const Properties&       params )

{
  Vector  cVec;

  params.get ( cVec,      "cVec"   );

  cVec[0] = c0_;
  cVec[1] = c0_;
  cVec[2] = 0.;
}



//-----------------------------------------------------------------------
//   getValuesAndDers
//-----------------------------------------------------------------------

void ConstantDiffusivity::getValuesAndDers

    ( const Properties&       params )

{
  Vector stress, cVec;
  Vector dcxx, dcxy, dcyy;

  params.get ( cVec,      "cVec"   );
  params.get ( dcxx,      "dcxx"   );
  params.get ( dcxy,      "dcxy"   );
  params.get ( dcyy,      "dcyy"   );

  cVec[0] = c0_;
  cVec[1] = c0_;
  cVec[2] = 0.;

  dcxx = 0.;
  dcyy = 0.;
  dcxy = 0.;
}

