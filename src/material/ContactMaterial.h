/*
 *
 */


#ifndef CONTACT_MATERIAL_H 
#define CONTACT_MATERIAL_H

#include <jem/base/String.h>
#include <jem/util/Flex.h>

#include "CohesiveMaterial.h"

using jem::String;
using jem::util::Flex;


// =======================================================
//  class ContactMaterial
// =======================================================


class ContactMaterial : public virtual CohesiveMaterial
{
 public:

  static const char*      DUMMY_PROP;

  /*
   *  constructor
   */

  explicit                ContactMaterial

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
      const Properties&     globdat )         const;

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
  
  virtual void            update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump0,
      const Vector&         djump,
      int                   mpoint );

  virtual void            elasticUpdate

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump );

  /*
   * update for moonen material
   */

  virtual void            update

    ( Vector&               traction,
      Matrix&               stiff,
      const Matrix&         jumpStiff,
      const Vector&         tEff,
      int                   mpoint );

  /*
   *  Called when the Newton Raphson converged
   *  Swap newHist_ to oldHist_ 
   */

  virtual void            commit           ();
  
  virtual bool            justTractionFree

    ( const int             i ) const;


  /**
   *  Allocate for history variables if it is not yet 
   *  initialized. Otherwise, extend it by appending at the end.
   *  Initial value of damage is optional
   */

  void                    allocPoints

    ( int                   count,
      double                dam = 0. );

  virtual void            deallocPoints

    ( int                   count );

  /*
   *  Return history variables at material point ipoint
   */

  inline double           giveDissipation   ( int ipoint ) const;
  
  inline virtual double   giveTensileStrength  (  ) const;

  inline double           giveHistory       ( int ipoint ) const;

  inline int              isLoading         ( int ipoint ) const;

  virtual int             wasLoading        ( int ipoint ) const;

  /*
   * Return number of integration points
   */

  virtual int             ipointCount  (  ) const
  {
    return 1;
  }

  double                  getDummy () const {return dummy_;}

 protected:

  virtual                ~ContactMaterial   ();

  void                    allocPoints_

     ( const int             count,
       const double          dam,
       const int             loading );


 protected:

  double                  dummy_;
};

// -------------------------------------------------------------------
//  giveHistory
// -------------------------------------------------------------------

inline double  ContactMaterial::giveHistory ( int ipoint ) const
{
  return 1.;
}

// -------------------------------------------------------------------
//  giveTensileStrength
// -------------------------------------------------------------------

inline double  ContactMaterial::giveTensileStrength () const
{
  return 0.;
}

// -------------------------------------------------------------------
//  isLoading
// -------------------------------------------------------------------

inline int     ContactMaterial::isLoading ( int ipoint ) const
{
  return 1;
}

// -------------------------------------------------------------------
//  wasLoading
// -------------------------------------------------------------------

inline int     ContactMaterial::wasLoading ( int ipoint )  const
{
  return 1;
}

// -------------------------------------------------------------------
//  giveDissipation
// -------------------------------------------------------------------

inline double  ContactMaterial::giveDissipation ( int ipoint ) const
{
  return 0;
}


#endif 
