#ifndef XCOHESIVE_MATERIAL_H
#define XCOHESIVE_MATERIAL_H

#include "CohesiveMaterial.h"
#include "OrthotropicMaterial.h"

using jem::Ref;
using jem::String;
using jive::Vector;
using jive::Matrix;
using jem::Object;
using jem::util::Properties;


//-----------------------------------------------------------------------
//   class XCohesiveMat
//-----------------------------------------------------------------------

// The XCohesiveMat class adds specific functionality to the 
// CohesiveMat class for use in XFEMModel
// Its main tasks are to calculate the traction and the stiffness matrix.

class XCohesiveMat : public virtual CohesiveMaterial
{
 public:

  explicit                XCohesiveMat

    ( const int             rank,
      const Properties&     globdat );

  virtual void            elasticUpdate

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump );
  
  virtual void            update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump,
      int                   mpoint ) = 0;

  
  virtual void            update

    ( Vector&               traction,
      const Vector&         jump,
      int                   mpoint ) = 0;

  // This is for cohesive models in which traction update 
  // requires incremental displacement jump.
  // VP Nguyen.

  virtual void            update

    ( Vector&             traction,
      Matrix&             stiff,
      const Vector&       jump0,
      const Vector&       djump,
      int                 mpoint );

  virtual double          evalFailure
    
    ( const Vector&         sigmaN,
      int                   mpoint ) const = 0;

  virtual void            scaleDissipationIncrement

    ( const double          factor,
      const int             ipoint ) = 0;

  // this is for initially elastic cohesive model to be used
  // in an initially rigid cohesive formulation such as XFEM, TwoScale
  // formulation. Only initially elastic cohesive model needs to implement this
  // method. VP Nguyen.
  
  virtual void            initShift 

    ( const int             ipoint,
      const Vector&         traction );


  virtual double          giveTensileStrength () const;
    
 protected:

  virtual                ~XCohesiveMat      ();

};

//-----------------------------------------------------------------------
//   related functions
//-----------------------------------------------------------------------


Ref<XCohesiveMat>    newXCohesiveMat

  ( const String&           name,
    const Properties&       conf,
    const Properties&       props,
    const Properties&       globdat );

#endif 
