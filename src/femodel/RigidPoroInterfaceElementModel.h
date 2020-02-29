/*
 * 
 *  Copyright (C) 2016 Monash University. All rights reserved.
 *  
 *  This class implements a model of two dimensional poro-cohesive interface elements. 
 *  To be used with the HygroMechanicalModel for the bulk (fully saturated medium)
 *  
 *  - bottom and upper faces: displacement and pore pressure dofs.
 *  - mid-plane: fracutring fluid pressure 
 *  The so-called triple-noded interface elements.
 *
 *  Before damage, nodes are constrained to avoid large dummy stiffness.
 *
 *  Author: V.P. Nguyen, phu.nguyen@monash.edu
 *  Date: 27 October 2016
 *
 *
 */


#ifndef RG_PORO_INTERFACE_ELEMENT_MODEL_H
#define RG_PORO_INTERFACE_ELEMENT_MODEL_H

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

class XCohesiveMat;

#include "FemModels.h"
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


typedef  ElementSet     ElemSet;

//=======================================================================
//   class RigidPoroInterfaceElementModel
//=======================================================================


class RGPoroInterfaceElementModel : public Model
{
 public:

  typedef RGPoroInterfaceElementModel     Self;
  typedef Model                         Super;

  // static data

  static const char*        SHAPE_PROP;
  static const char*        MATERIAL_PROP;
  static const char*        MESH_PROP;
  static const char*        THICKNESS_PROP;
  static const char*        DISP_NODES_PROP;
  static const char*        DIRICHLET_NODES_PROP;

  static const char*        BIOT_COEF_PROP;
  static const char*        VISCOSITY_PROP;
  static const char*        DTIME_PROP;

  // Constructor
                            RGPoroInterfaceElementModel

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

  virtual                  ~RGPoroInterfaceElementModel       ();


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

  void                       initPressureConstraints_     ( );
  
  void                       zeroPressure_                ( );
  
  void                       initDisplacementConstraints_ ( );

  void                       removeDisplacementConstraints_
    
    ( int                     ie );

  void                       constraintDirichletNodes_    ( );

  void                       updateBcConstraints_         ( );

  inline void                addConstraint_ 

    ( int xdof,  int ydof,
      int xdofM, int ydofM );
  
  void                       checkCommit_

    ( const Properties& globdat, 
      const Properties& params, 
      const Vector&     disp );

  void                       initializeIPMPMap_ ( );

  void                       computeForceForInjection_ 

    ( const Vector&    force );
  
  void                       writeNodalOutput_ 
    
    ( const Properties& globdat ); 
  
  void                       writeWellboreOutput_ 
    
    ( const Properties& globdat ); 
  
  double                     evalElementFailure_

    ( int               ie ,
      int&              ig ,
      const Vector&     disp,
      const Properties& globdat,
      const Matrix&     tractions );

 private:
 
  Assignable<ElementGroup>  egroup_;
  Assignable<ElemSet>       allElems_;    // all elements
  Assignable<ElementSet>    elems_;       // interface elements
  Assignable<NodeSet>       nodes_;
   
  int                       rank_;
  int                       ielemCount_;  // num of total elements
  int                       nodeCountF_;  // no.of. all nodes per element = const
  int                       nodeCount_;   // no.of. solid nodes per element = const
  int                       nodeCount2_;  // no.of.nodes per element = const
  int                       dofCount_;   // no.of.dofs per element
  int                       ipCount_;    // num of GP per element
  
  // to compute averaged stress to check failure

  IntMatrix                 ibulkElems_;     // indices of bulk elements
  Ref<IShape>               ishape_;    // shape functions of volumetric elements
  int                       inodeCount_;   
  int                       idofCount_;   
  int                       checkInterval_;
  IdxVector                 belemMatMap_; // elem material map for bulk elems
  ShapeGradsFunc            getShapeGrads_;
  IdxVector                 damagedElems_; // 1: damaged

  Ref<FShape>               shape_;
  Ref<IShape>               sshape_;
  Ref<XDofSpace>            dofs_;
  Ref<Constraints>          cons_;

  Array<Ref<XCohesiveMat> >     
    
                            materials_; // multi materials
  
  IdxVector                 aaDofTypes_; // displacement dof types
  IdxVector                 pwDofTypes_; // pore pressure dof types
  IdxVector                 pfDofTypes_; // fracturing fluid pressure dof types
  IdxVector                 dispNodes_; // indices of nodes with prescribed disp
  IdxVector                 diriNodes_; // indices of nodes with dirichlet Bcs 
  
  IdxVector                 ielems_;    // connectivity of group interface elements
  
  SparseArray <int, 1>      elemMatMap_; // hold the mat index for all elements

  SparseArray <int, 2>      ipMpMap_;  // mapping between integration point and 

  std::vector<std::vector<int> > 
                            dupNodes_; 

  double                    thickness_;  // thickness for plane stress problems
  double                    loc2Glob_;

  bool                      initConsDone1_;
  bool                      initConsDone2_;

  FShapeFunc                getFShapeFuncs_;

  bool                      isPrsConstraint_; // true: enforce pressure continuity
  bool                      isDspConstraint_; // true: enforce displacements before damage
  bool                      isConstant_;   // constant longitudianl permeability

  double                    w0_;     // initial aperture
  double                    permea_; // longitudial (constant) permeability
  double                    kt_;     // transversal (constant) permeability, k_b=k_t
  double                    alpha_;
  double                    wf_;
  double                    dtime_;
  double                    mu_;     // 1/(12*\mu_f)
  double                    Q0_;     // injection rate [m3/s]
  int                       q0Node_; // node where injection rate is applied
  
  // for writing pressure fields to file
  // for upper nodes, lower nodes and mid nodes

  IdxVector                 inodesA_;
  IdxVector                 inodesB_;
  IdxVector                 inodesM_;
  String                    fileName_;
  int                       interval_;

  bool                      resolve_;
  bool                      tractionFree_;
  bool                      write_;
  bool                      writeWellbore_;

  IdxVector                 wellboreNodes_; // nodes at the wellbore
  IdxVector                 wellboreUDofA_; // nodes at the wellbore
  IdxVector                 wellborePDofA_; // nodes at the wellbore
  IdxVector                 wellboreUDofB_; // nodes at the wellbore
  IdxVector                 wellborePDofB_; // nodes at the wellbore
  IdxVector                 wellborePDofM_; // nodes at the wellbore
  Ref<PrintWriter>          wellboreOut_;

  Matrix                    coords1D_;
  IntMatrix                 ipStatus_;     // elastic or damaged, ipStatus(elemId,ipId)={0,1}
};

// =========================================================
//   implementation of inline functions
// =========================================================

inline void RGPoroInterfaceElementModel::addConstraint_

 ( int xdof,  int ydof,
   int xdofM, int ydofM )
{
  cons_->addConstraint ( xdof, xdofM, 1.0 );
  cons_->addConstraint ( ydof, ydofM, 1.0 );
}


#endif
