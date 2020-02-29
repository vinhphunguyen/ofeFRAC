/**
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens a common interface for material models. 
 *  
 *  \author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *
 */

#ifndef MATERIAL_H
#define MATERIAL_H


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
//   class Material
//-----------------------------------------------------------------------

// The Material class represents a material model. Its main task is to
// calculate the strain-stress stiffness matrix.

class Material : public Object
{
 public:

  explicit                Material

    ( const int             rank,
      const Properties&     globdat );

  virtual void            configure

    ( const Properties&     props,
      const Properties&     globdat );

  virtual void            getConfig

    ( const Properties&     conf, 
      const Properties&     globdat ) const ; 

  /**
   * @brief      { Stress update for GP 'ip' given total strain }
   *
   * @param      stress  The stress vector
   * @param      stiff   The tangent stiffness
   * @param[in]  strain  The total strain vector
   * @param[in]  ip      The ID of Gauss point 
   */

  virtual void            update

    ( Vector&         stress,
      Matrix&         stiff,
      const Vector&   strain,
      int             ip )       = 0;

  /**
   * @brief     
   * Used for crack band approach: adjust the softening modulus according to
   * the element size "he". Note that default (empty) implementation is provided
   * in this base class so that material models do not need this will not have
   * to care about this function.
   *
   * @param      stress  The stress
   * @param      stiff   The stiff
   * @param[in]  strain  The strain
   * @param[in]  ipoint  The ipoint
   * @param[in]  he      The element characteristic size
   */
  
  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      int                   ipoint,
      double                he );

  /**
   * @brief      { stress update given strain increment }
   *
   * @param      stress   The stress
   * @param      stiff    The stiff
   * @param[in]  strain0  The strain0
   * @param[in]  dstrain  The dstrain
   * @param[in]  ipoint   The ipoint
   * @param[in]  he       The element characteristic size
   */
  
  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain0,
      const Vector&         dstrain,
      int                   ipoint,
      double                he );

  virtual void            commit ();

  virtual void            allocPoints

    ( int                   count );

  virtual double          giveHistory ( int point ) const = 0;
  inline  double          giveRho     ( ) const ;

  virtual double          checkLocalisation 

    ( Vector&               normal,
      const Vector&         stress,
      const Matrix&         tangent,
            int             ipoint ) const;

  virtual Vector          giveStress ( int ip ) const;


 protected:

  virtual                ~Material      ();

 protected:

  int                    rank_;
  double                 rho_;
};

double Material::giveRho () const { return rho_;}

//-----------------------------------------------------------------------
//   related functions
//-----------------------------------------------------------------------


Ref<Material>             newMaterial

  ( const String&           name,
    const Properties&       conf,
    const Properties&       props,
    const Properties&       globdat );

#endif 
