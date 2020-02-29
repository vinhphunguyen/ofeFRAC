#ifndef COHESIVE_MATERIAL_H
#define COHESIVE_MATERIAL_H

#include <jive/Array.h>
#include <jem/base/Error.h>
#include <jem/base/Object.h>
#include <jem/util/Properties.h>
#include <jive/util/utilities.h>

#include "Material.h"
//#include "DiscontinuousShape.h"

using jem::Ref;
using jem::String;
using jive::Vector;
using jive::Matrix;
using jem::Object;
using jem::util::Properties;
using jive::BoolVector;

//-----------------------------------------------------------------------
//   class CohesiveMaterial
//-----------------------------------------------------------------------

// The CohesiveMaterial class represents a cohesive material model. 
// Its main tasks are to calculate the traction and the stiffness matrix.

class CohesiveMaterial : public Object
{
 public:

  explicit                CohesiveMaterial

    ( const int             rank,
      const Properties&     globdat );

  virtual void            configure

    ( const Properties&     props,
      const Properties&     globdat );

  virtual void            getConfig

    ( const Properties&     conf,
      const Properties&     globdat )         const;

  virtual void            update

    ( Vector&             traction,
      Matrix&             stiff,
      const Vector&       jump,
      int                 mpoint )       = 0;
  
  virtual void            update

    ( Vector&             traction,
      Matrix&             stiff,
      const Vector&       jump0,
      const Vector&       djump,
      int                 mpoint );

   virtual void            elasticUpdate

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump )      = 0;
  
   virtual bool            justTractionFree

    ( const int             ipoint ) const = 0;


  virtual void            commit ();

  virtual void            cancel ();

  virtual bool            despair ();

  virtual bool            endDespair ();

  virtual void            allocPoints

    ( int                   count,
      double                dam = 0. );

  virtual void            deallocPoints

    ( int                   count );

  virtual double          giveHistory       ( int point  ) const = 0;

  virtual inline double   giveDissipation   ( int point  ) const;

  virtual inline int      isLoading         ( int point  ) const;

  virtual inline int      wasLoading        ( int point  ) const;

  //virtual bool            isFailure
  //  
  //  ( const Vector&       stress,
  //    const Ref<Material> bulkMat )      const;

  //virtual double          evalFailure
  //  
  //  ( const Vector&       stress,
  //    const Ref<Material> bulkMat )      const;

  //virtual bool            evalFailure 

  //  ( const Vector& traction )           const;

  virtual int             ipointCount () const;


 protected:

  virtual                ~CohesiveMaterial      ();

 protected:

  int                    rank_;

  bool                   desparateMode_;
  BoolVector             hasSwitched_;
  BoolVector             useSecant_;
};


//-----------------------------------------------------------------------
//   implementation
//-----------------------------------------------------------------------

inline double  CohesiveMaterial::giveDissipation ( int ipoint ) const
{
  return 0.;
}

inline int CohesiveMaterial::isLoading ( int ipoint ) const
{
  return false;
}

inline int CohesiveMaterial::wasLoading ( int ipoint ) const
{
  return false;
}

//-----------------------------------------------------------------------
//   related functions
//-----------------------------------------------------------------------


Ref<CohesiveMaterial>     newCohesiveMaterial

  ( const String&           name,
    const Properties&       conf,
    const Properties&       props,
    const Properties&       globdat );

#endif 
