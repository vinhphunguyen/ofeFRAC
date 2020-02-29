#ifndef IGA_ELASTICITY_MODEL_H
#define IGA_ELASTICITY_MODEL_H

#include <jive/model/Model.h>
#include <jem/util/SparseArray.h>
#include "model_import.h"

namespace jive
{
   namespace util
   {
      class Constraints;
   }
}

using jive::util::Constraints;
using jem::util::SparseArray;
using jive::fem::IElement;

//-----------------------------------------------------------------------
//   class IGAElasticityModel
//-----------------------------------------------------------------------


class IGAElasticityModel : public Model
{
 public:

  static const char*        TYPE_NAME;
  static const char*        MATERIAL_PROP;
  static const char*        SHAPE_PROP;
  static const char*        DOF_NAMES[3];

                            IGAElasticityModel  ();

                            IGAElasticityModel

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

  virtual                  ~IGAElasticityModel  ();


 private:

  void                      calcMatrix_

    ( MBuilder&               mb,
      const Vector&           fint,
      const Vector&           state )        const;

  bool                      getTable_

    ( const Properties&       params,
      const Properties&       globdat )      const;

  //void                      calcFlux_

    //( XTable&                 table,
     // const Vector&           state )        const;


 private:

  Ref<XDofSpace>            dofs_;
  Ref<Constraints>          cons_;
  Ref<IElement>             element_;
  IdxVector                 dofTypes_;

  idx_t                     elemCount_;
  idx_t                     rank_;
  idx_t                     strCount_; // 2D:3 and 3D:6

  IdxVector                 ielems_;

  ShapeGradsFunc            getShapeGrads_;
  ShapeFunc                 getShapeFuncs_;
  
  Array< Ref<Material> >    materials_;
  //SparseArray <int, 1>      elemMatMap_; // hold the mat index for all elements
  IdxVector                  elemMatMap_; // hold the mat index for all elements
                                          // use of SparseArray cannot be stored in globdat

};


#endif
