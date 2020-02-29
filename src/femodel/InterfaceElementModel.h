/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implements a model of two dimensional interface elements. 
 *  The overall procedure of using interface elements is 
 *   (1) Build a standard FE mesh
 *   (2) Using a small program to separate this mesh and insert interface
 *       elements in between interelement boundaries. The modified mesh
 *       is saved to a file named *-solid.mesh and the interface elements
 *       are written to a file named *-interface.mesh.
 *   (3) The problem is then modelled using two models:
 *
 *      (a) Continuum model using the *-solid.mesh
 *      (b) Interface model (this class) using the *-interface.mesh
 *  
 *
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 5 July 2009
 *
 *  October 2014: change getMatrix_ so that it provides [u]_0 and d[u]
 *  to the constitutive model. Before it worked on the total displacement jump
 *  [u].
 *
 */


#ifndef INTERFACE_ELEMENT_MODEL_H
#define INTERFACE_ELEMENT_MODEL_H

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
using jive::IntVector;
using jive::IdxVector;
using jive::util::XDofSpace;
using jive::util::XTable;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::algebra::MatrixBuilder;
using jive::model::Model;
using jive::geom::FShape;
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


class InterfaceElementModel : public Model
{
 public:

  typedef InterfaceElementModel     Self;
  typedef Model                     Super;

  // static data

  static const char*        SHAPE_PROP;
  static const char*        MATERIAL_PROP;
  static const char*        MESH_PROP;
  static const char*        THICKNESS_PROP;
  static const char*        DISP_NODES_PROP;
  static const char*        DIRICHLET_NODES_PROP;


  // Constructor
                            InterfaceElementModel

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

  virtual                  ~InterfaceElementModel       ();


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

  void                       initConstraints_ ();

  void                       constraintDirichletNodes_ ();

  void                       updateBcConstraints_ ();

  inline void                addConstraint_ 

    ( int xdof,  int ydof,
      int xdofM, int ydofM );

  void                       checkCommit_

    ( const Properties& params, 
      const Properties& globat );

  void                       initializeIPMPMap_ ( );

 private:
 
  Assignable<ElementGroup>  egroup_;
  Assignable<ElementSet>    elems_;
  Assignable<NodeSet>       nodes_;
   
  int                       rank_;
  int                       ielemCount_;  // num of total elements
  int                       nodeCount_;  // no.of.nodes per element = const
  int                       nodeCount2_;  // no.of.nodes per element = const
  int                       dofCount_;   // no.of.dofs per element
  int                       ipCount_;    // num of GP per element

  Ref<FShape>               shape_;
  Ref<XDofSpace>            dofs_;
  Ref<Constraints>          cons_;

  Array<Ref<CohesiveMaterial> >     
    
                            materials_; // multi materials
  
  IdxVector                 dofTypes_;
  IdxVector                 activeElems_; // 0: inactive 1: active
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

  ShapeFunc                 getShapeFuncs_;
};

// =========================================================
//   implementation of inline functions
// =========================================================

inline void InterfaceElementModel::addConstraint_

 ( int xdof,  int ydof,
   int xdofM, int ydofM )
{
  cons_->addConstraint ( xdof, xdofM, 1.0 );
  cons_->addConstraint ( ydof, ydofM, 1.0 );
}


#endif
