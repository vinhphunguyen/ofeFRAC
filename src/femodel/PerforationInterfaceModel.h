/*
 * 
 *  Copyright (C) 2016 Monash University. All rights reserved.
 *  
 *
 *
 */


#ifndef PERFO_INTERFACE_MODEL_H
#define PERFO_INTERFACE_MODEL_H

#include <vector>

#include <jem/base/Array.h>
#include <jem/util/Flex.h>
#include <jem/util/ArrayBuffer.h>
#include <jem/util/SparseArray.h>

#include <jive/util/Assignable.h>
#include <jive/util/utilities.h>
#include <jive/util/XTable.h>
#include <jive/model/Model.h>
#include <jive/fem/NodeSet.h>
#include <jive/fem/ElementGroup.h>
#include <jive/util/Constraints.h>

namespace jem
{
  namespace util
  {
    class Properties;
  }
}

namespace jive
{
  namespace util
  {
    class XDofSpace;
  }
}

namespace jive
{
  namespace algebra
  {
    class MatrixBuilder;
  }
}

namespace jive
{
  namespace fem
  {
    class Globdat;
  }
}

namespace jive
{
  namespace geom
  {
    class InterfaceShape;
    class InternalShape;
    class BoundaryShape;
  }
}


#include "FemModels.h"
#include "material/CohesiveMaterial.h"
#include "util/utilities.h"


using jem::Array;
using jem::util::Properties;
using jem::util::ArrayBuffer;
using jem::util::Flex;
using jem::util::SparseArray;
using jive::Vector;
using jive::IdxVector;
using jive::util::XDofSpace;
using jive::util::XTable;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::algebra::MatrixBuilder;
using jive::model::Model;
using jive::geom::FShape;
using jive::geom::IShape;
using jive::fem::Globdat;
using jive::fem::NodeSet;
using jive::fem::ElementSet;
using jive::fem::ElementGroup;
using jive::IntMatrix;

//=======================================================================
//   typedefs
//=======================================================================



//=======================================================================
//   class InterfaceElementModel
//=======================================================================


class PerforationInterfaceModel : public Model
{
 public:

  typedef PerforationInterfaceModel     Self;
  typedef Model                         Super;

  // static data

  static const char*        SHAPE_PROP;
  static const char*        MATERIAL_PROP;
  static const char*        MESH_PROP;
  static const char*        THICKNESS_PROP;

  static const char*        BIOT_COEF_PROP;
  static const char*        VISCOSITY_PROP;

  // Constructor
                            PerforationInterfaceModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  // Configure itself by reading data from properties file

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );


  virtual void              getConfig

    ( const Properties&       conf )             const;

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );


 protected:

  virtual                  ~PerforationInterfaceModel       ();


 private:

  // getMatrix_: compute the tangent stiffness and internal force vector
  //  Compute K and fint for the updated displacement disp

  void                      getIntForce_

    ( const Vector&           force,
      const Vector&           disp0,
      const Vector&           disp );

  void                      getMatrix_

    ( MatrixBuilder&          mbuilder,
      const Vector&           force,
      const Vector&           disp0,
      const Vector&           disp );

  void                       getTransformationMatrix_

    ( const Matrix&           Q,
      const Matrix&           coord );
  
  void                       initPressureConstraints_ ();
  
  void                       initDispConstraints_     ();
  
  void                       removeDispConstraints_   ();

  void                       computeForceForInjection_ 

    ( const Vector&    force );
  
  void                       writeWellboreOutput_ 
    
    ( const Properties& globdat ); 

 private:
 
  Assignable<ElementGroup>  egroup_;
  Assignable<ElementSet>    elems_;
  Assignable<NodeSet>       nodes_;
   
  int                       rank_;
  int                       ielemCount_;  // num of total elements
  int                       nodeCountF_;  // no.of. all nodes per element = const
  int                       nodeCount_;   // no.of. solid nodes per element = const
  int                       nodeCount2_;  // no.of.nodes per element = const
  int                       dofCount_;   // no.of.dofs per element
  int                       ipCount_;    // num of GP per element

  Ref<FShape>               shape_;
  Ref<IShape>               sshape_;
  Ref<XDofSpace>            dofs_;
  Ref<Constraints>          cons_;

  IdxVector                 aaDofTypes_; // displacement dof types
  IdxVector                 pwDofTypes_; // pore pressure dof types
  IdxVector                 pfDofTypes_; // fracturing fluid pressure dof types
  IdxVector                 ielems_;    // connectivity of group interface elements
  
  double                    thickness_;  // thickness for plane stress problems
  
  int                       loc2Glob_;
  FShapeFunc                getFShapeFuncs_;
  
  std::vector<std::vector<int> > 
                            dupNodes_; 

  double                    permea_; // longitudial (constant) permeability
  double                    alpha_;
  double                    wf_;
  double                    dtime_;
  double                    mu_;     // 1/(12*\mu_f)
  double                    Q0_;     // injection rate [m3/s]
  IdxVector                 q0Dofs_;
  
  // for writing pressure fields to file
  // for upper nodes, lower nodes and mid nodes

  IdxVector                 inodesA_;
  IdxVector                 inodesB_;
  IdxVector                 inodesM_;
  String                    fileName_;
  int                       interval_;

  bool                      first_;
  bool                      write_;
  bool                      isConstraint_;
  bool                      writeWellbore_;

  IdxVector                 wellboreNodes_; // nodes at the wellbore
  IdxVector                 wellboreUDofA_; // nodes at the wellbore
  IdxVector                 wellborePDofA_; // nodes at the wellbore
  IdxVector                 wellboreUDofB_; // nodes at the wellbore
  IdxVector                 wellborePDofB_; // nodes at the wellbore
  IdxVector                 wellborePDofM_; // nodes at the wellbore
  Ref<PrintWriter>          wellboreOut_;
  IdxVector                 damagedElems_; // 1 => damaged, 0 => undamaged

  Matrix                    coords1D_;
};

#endif
