#ifndef IGA_FINITE_DEFORMATION_MODEL_H
#define IGA_FINITE_DEFORMATION_MODEL_H

#include <jive/model/Model.h>
#include <jem/util/SparseArray.h>
#include "model_import.h"
#include "util/LargeDispUtilities.h"

using jem::util::SparseArray;

//-----------------------------------------------------------------------
//   class IGAElasticityModel
//-----------------------------------------------------------------------


class IGAFiniteDeformationModel : public Model
{
 public:

  static const char*        TYPE_NAME;
  static const char*        MATERIAL_PROP;
  static const char*        SHAPE_PROP;
  static const char*        DOF_NAMES[3];

                            IGAFiniteDeformationModel  ();

                            IGAFiniteDeformationModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );

  static Ref<Model>         makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  static void               declare       ();

 protected:

  virtual                  ~IGAFiniteDeformationModel  ();


 private:

  void                      calcMatrix_

    ( MBuilder&               mb,
      const Vector&           fint,
      const Vector&           state )        const;

  bool                      getTable_

    ( const Properties&       params,
      const Properties&       globdat )      const;

 private:

  Ref<XDofSpace>            dofs_;
  Ref<BezierElement>        element_;
  IdxVector                 dofTypes_;

  idx_t                     elemCount_;
  idx_t                     rank_;
  idx_t                     strCount_; // 2D:3 and 3D:6

  IdxVector                 ielems_;

  ShapeGradsFunc            getShapeGrads_;
  ShapeFunc                 getShapeFuncs_;
  BMatrixLinFunc            getBMatrixLinFunc_;
  
  Array< Ref<Material> >    materials_;
  SparseArray <int, 1>      elemMatMap_; // hold the mat index for all elements

};


#endif
