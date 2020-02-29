/*
 * 
 *  Copyright (C) 2017 Monash University. All rights reserved.
 *  
 *  This class implements a one dimensional J2 plasticity material. To be used with 
 *  TrussModel to model reinforcement steel bars.
 *  
 *  Author: V.P. Nguyen
 *  Date: 13 Dec 2017
 *
 */

#ifndef ONE_D_J2_MATERIAL_H
#define ONE_D_J2_MATERIAL_H


#include <jive/Array.h>
#include <jem/base/Object.h>

namespace jem
{
  namespace util
  {
    class Properties;
  }
}

using jem::Ref;
using jem::String;
using jive::Vector;
using jive::Matrix;
using jem::Object;
using jem::util::Properties;

//-----------------------------------------------------------------------
//   class OneDimJ2PlasticityMat
//-----------------------------------------------------------------------

class OneDimJ2PlasticityMat : public Object
{
 public:

                         OneDimJ2PlasticityMat

    ( const Properties&     conf,
      const Properties&     props,
      const Properties&     globdat );

  /*
   * Update is the most important function of a Material class.
   * It is used to compute the stress (and tangent) for a given
   * total strain provided. This is performed at the integration point
   * "ip".
   * All the required parameters are stored in the dictonary 'params'.
   * Output also stored there as well, except tangent matrix as jive:;Properties
   * does not allow to store a matrix.
   */

  virtual void            update

    ( const Properties&      params );
  
  virtual void            commit ();

  virtual void            allocPoints

    ( jem::idx_t             count );


 protected:

  virtual                ~OneDimJ2PlasticityMat ();

 protected:

  double                    young_;
  double                    k_;               // plastic modulus
  double                    sigmaY_;          // initial yield stress

  Vector                    epsilonP_;       // plastic strain
  Vector                    epsilonP0_;      // converged plastic strain
  Vector                    alpha_;          // hardening variable alpha
  Vector                    alpha0_;         //  converged value of alpha
};


//-----------------------------------------------------------------------
//   related functions
//-----------------------------------------------------------------------


Ref<OneDimJ2PlasticityMat>          newOneDimJ2Mat

  ( const Properties&       conf,
    const Properties&       props,
    const Properties&       globdat );

#endif 
