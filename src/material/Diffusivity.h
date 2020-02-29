/*
 * 
 *  Copyright (C) 2018 Monash University. All rights reserved.
 *  
 *  This class implements a common interface for the diffusitivity matrix c
 *  in stress based gradient damage models.
 *  
 *  Author: V.P. Nguyen
 *  Date: 13 March 2018
 *
 */

#ifndef DIFFUSIVITY_H
#define DIFFUSIVITY_H


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
//   class Diffusivity
//-----------------------------------------------------------------------


class Diffusivity : public Object
{
 public:

                         Diffusivity

    ( const Properties&     conf,
      const Properties&     props,
      const Properties&     globdat );

  /*
   * calculate c matrix based on the stresses
   */

  virtual void            getValues

    ( const Properties&      params ) = 0;
  
  /*
   * calculate c matrix and its derivatives w.r.t stresses 
   * based on the stresses
   */

  virtual void            getValuesAndDers

    ( const Properties&      params ) = 0;
  
 protected:

  virtual                ~Diffusivity   ();

 protected:

  double                    c0_;
};


//-----------------------------------------------------------------------
//   related functions
//-----------------------------------------------------------------------


Ref<Diffusivity>             newDiffusivity

  ( const Properties&       conf,
    const Properties&       props,
    const Properties&       globdat );

#endif 
