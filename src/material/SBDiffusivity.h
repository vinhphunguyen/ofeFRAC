/*
 * 
 *  Copyright (C) 2017 Monash Uni. All rights reserved.
 *  
 *  This class implemens the isotropic phase field material.
 *  Isotropic elastic-damage material without  strain 
 *  energy split.
 *  
 *  Author: V.P. Nguyen
 *  Date: 15 March 2017
 *
 */

#ifndef STRESS_BASED_DIFFUSITIVITY_H
#define STRESS_BASED_DIFFUSITIVITY_H

#include <jem/base/String.h>

#include "util/utilities.h"
#include "Diffusivity.h"


using namespace jem;
using jem::String;


// =======================================================
//  class SBDiffusitivity
// =======================================================


class SBDiffusivity : public Diffusivity
{
 public:

                          SBDiffusivity

    ( const Properties&     conf,
      const Properties&     props,
      const Properties&     globdat );

  virtual void            getValues

    ( const Properties&     params );
  
  virtual void            getValuesAndDers

    ( const Properties&     params );

 protected:

  virtual                ~SBDiffusivity   ();

 private:

  double                   ft_;
  double                   fac_;
  double                   fac2_;
};

#endif 
