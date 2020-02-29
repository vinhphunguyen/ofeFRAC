/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the traction-seperation law
 *  for cohesive crack models. This is used in XFEMModel
 *  to model cohesive crack propagation.
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 21 December 2007
 *
 */


#ifndef XU_NEEDLEMAN_MATERIAL_H
#define XU_NEEDLEMAN_MATERIAL_H

#include <jem/base/String.h>
#include <jem/util/Flex.h>

#include "CohesiveMaterial.h"


using jem::String;
using jem::util::Flex;


// =======================================================
//  class XuNeedlemanMaterial
// =======================================================


class XuNeedlemanMaterial : public CohesiveMaterial
{
 public:

  static const char*      BETA_PROP;
  static const char*      TENSILE_PROP;
  static const char*      SHEAR_PROP;
  static const char*      OPEN_DISP_PROP;
  static const char*      SHEAR_DISP_PROP;
  static const char*      R_PROP;
  static const char*      Q_PROP;
  static const char*      ENERGY_PROP;
  
  /*
   *  constructor
   */

  explicit                XuNeedlemanMaterial

    ( const int             rank,
      const Properties&     globdat );

  /*
   *  configure and getConfig (from input file)
   */

  virtual void            configure

    ( const Properties&     props,
      const Properties&     globdat );

  virtual void            getConfig

    ( const Properties&     conf,
      const Properties&     globdat  )         const;

  /*
   *  compute the traction (t) and cohesive tangent stiffness matrix
   *  (stiff) at material point mpoint given the displacement jump (jump)
   *   jump[0] = crack opening displacement
   *   jump[1] = crack sliding displacement
   */

  virtual void            update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump,
      int                   mpoint );

  /*
   * update for moonen material
   */

  virtual void            update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         tEff,
      const Matrix&         jumpStiff,
      int                   mpoint );

  virtual void            elasticUpdate

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump );


  /*
   *  Called when the Newton Raphson converged
   *  Swap newHist_ to oldHist_ 
   */

  virtual void            commit           ();


  /*
   *  Allocate for history variables if it is not yet 
   *  initialized. Otherwise, extend it by appending at the
   *  end 
   */

  void                    allocPoints

    ( int                   count,
      double                dam = 0. );


  /*
   *  Return history variables at material point ipoint
   */

  inline double           giveHistory       ( int ipoint ) const;

  /*
   *  Return the tensile stress of the material, ft
   */

  inline double           giveTensileStress ( )            const;


  //virtual bool            isFailure

  //  ( const Vector&       stress,
  //    const Ref<Material> bulkMat )                        const;

  //double                  evalFailure

  //  ( const Vector&       stress,
  //    const Ref<Material> bulkMat )                        const;

  //bool                    evalFailure 

  //  ( const Vector& traction ) const;

  /*
   * Return number of integration points
   */

  virtual int             ipointCount  (  ) const
  {
    return preHist_.loading.size ( ) ;
  }

  bool                    justTractionFree

    (const int                    ip ) const;


 protected:

  virtual                ~XuNeedlemanMaterial   ();

 private:

  virtual void            getTraction_

    ( const Vector& traction,
      const Vector& separation,
      double        eqvDisp, 
      double        eqvDispMax, 
      int           loading ) const;

  virtual void            getTangentMatrix_

    ( const Matrix& tangent,
      const Vector& separation,
      double        eqvDisp, 
      double        eqvDispMax, 
      int           loading ) const;


 private:

  int                     rank_;

  double                  beta_;
  double                  r_;
  double                  shearStress_; 
  double                  tensileStress_;
  double                  critOpen_;
  double                  critSlid_;
  double                  critOpenInv_; // inverse of critOpen_;
  double                  critSlidInv_;

  double                  q_;
  double                  phin_;
  //double                  phit_;

  // history variable (equivalent crack opening), involve in time

  struct                  hist_
  {
    Flex<double>            eqvOpenMax ;
    Flex<int>               loading;
  };

  hist_                   preHist_;      // history of the previous load step
  hist_                   newHist_;      // history of the current iteration

};

// -------------------------------------------------------------------
//  giveHistory
// -------------------------------------------------------------------

inline double  XuNeedlemanMaterial::giveHistory ( int ipoint ) const
{
  return newHist_.eqvOpenMax[ipoint];
}

// -------------------------------------------------------------------
//  giveTensileStress
// -------------------------------------------------------------------

inline double  XuNeedlemanMaterial::giveTensileStress (  ) const
{
  return tensileStress_;
}

#endif 
