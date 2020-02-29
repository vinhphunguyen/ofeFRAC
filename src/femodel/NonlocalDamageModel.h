/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the nonlocal integral damage model
 *  using finite elements.
 *
 *  This models applies for element group. This design allows
 *  for example, nonlocal damage model and elastic model exist
 *  at the same time. The element group can contain subgroups
 *  for different materials. 
 *
 *  For efficiency, only mesh of one element type is allowed.
 *  
 *  Remark: In the input file, must use
 *  matrix.type = "Sparse"; for non-standard node by node assembly

 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *
 */


#ifndef NONLOCAL_DAMAGE_MODEL_H
#define NONLOCAL_DAMAGE_MODEL_H

#include <jem/numeric/func/Function.h>
#include <jive/util/Assignable.h>

#include "model_import.h"
#include "util/IntegrationPointPair.h"

using jive::util::Assignable;
using jive::IntMatrix;
using jem::util::SparseArray;

//=======================================================================
//   typedefs
//=======================================================================

// Some handy aliases to avoid some typing

typedef MatmulChain<double,1>  MChain1;
typedef MatmulChain<double,2>  MChain2;
typedef MatmulChain<double,3>  MChain3;

typedef ElementSet             ElemSet;
typedef ElementGroup           ElemGroup;



//=======================================================================
//   class NonlocalDamageModel
//=======================================================================


class NonlocalDamageModel : public Model
{
 public:

  typedef NonlocalDamageModel       Self;
  typedef Model                     Super;

  static const char*        DOF_NAMES[3];

  static const char*        SHAPE_PROP;
  static const char*        MATERIAL_PROP;
  static const char*        WEIGHT_FUNC_PROP;
  static const char*        RADIUS_PROP;
  static const char*        AVERAGING_PROP;


  // Constructor
                            NonlocalDamageModel

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

  virtual                  ~NonlocalDamageModel       ();


 private:

  // getMatrix_: compute the tangent stiffness and internal force vector
  //  Compute K and fint for the updated displacement disp

  void                      getMatrix_

    ( MBuilder&               mbuilder,
      const Vector&           force,
      const Vector&           disp );


  //  Given the displacement @disp, compute the normal strains, then
  //  the local equivalent strain and finally the nonlocal one by weighted
  //  averaging.

  void                      calcStrain_

    ( const Vector&           disp );

  // Generate the list of pairs of integration points

  void                      initIPPairs_      ();

  // Make an integration point pair (i,j)

  void                      findIPPairs_

    ( double                  radius,
      int                     ipoint,
      int                     jpoint,
      const Matrix&           coords1,
      const Matrix&           coords2,
      const Vector&           weights1,
      const Vector&           weights2 );


  /* Non symmetric averaging operator 
   * Compute the nonlocal equivalent strain from the local 
   * equivalent strain by renormalized nonsymmetric formula
   */

  void                      getNonlocalStrainNSym_

    ( const Vector& localStrain );

  /* Symmetric averaging operator 
   * Compute the nonlocal equivalent strain from the local 
   * equivalent strain by symmetric formula
   */

  void                      getNonlocalStrainSym_

    ( const Vector& localStrain );

  /* params.set ( "accept", false ) => keep same load and resolve
   * params.set ( "accept", true )  => advance to new load step (normal case)
   */

  void                      checkCommit_

    ( const Properties&       params );

  // Make stress, strain and damage tables for visualisation

  bool                      getTable_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      getStrain_

    ( XTable&                 table,
      const Vector&           weights,
      const Vector&           disp );

  void                      getStress_

    ( XTable&                 table,
      const Vector&           weights );
  
  void                      getHistory_

    ( XTable&                 table,
      const Vector&           weights );


  /**
   * @brief      This function is called by getTable_.
   *
   * @param[out] table     The table
   * @param[in]  weights   The weights
   * @param[in]  contents  Filter for which data to be written to table
   * @param[in]  state     The nodal value vector
   */

  void                      getOutputData_

    ( Ref<XTable>             table,
      const Vector&           weights,
      const String&           contents,
      const Vector&           state );      

  /* Initialize the mapping between integration points
   * and material points 
   */

  void                      initializeIPMPMap_ ();

 private:

  Assignable<ElemGroup>     egroup_;
  Assignable<ElemSet>       elems_;
  Assignable<NodeSet>       nodes_;

  int                       rank_;
  int                       numElem_;    // num of total elements
  int                       nodeCount_;  // no.of.nodes per element = const
  int                       dofCount_;   // no.of.dofs per element
  int                       strCount_;   // equals 3 (2D) and 6(3D)
  int                       ipCount_;     // num of GP per element

  Ref<IShape>               shape_;

  Ref<XDofSpace>            dofs_;
  IdxVector                 dofTypes_;
  IdxVector                 ielems_;     // indices of elements in the group

  /* strains (including nonlocal equivalent one) at integration points */

  Matrix                    strain_;

  // data structure for non-local integration

  Flex<IPointPair>          ipPairs_;
  Vector                    eqWeights_;
  Array< Flex<IPointPair> > interactMap_;

  // materials: it is possible to treat multi-material problem
  // note: material used in damage model must be a Damage material

  Array< Ref<DamageMaterial> >       materials_;

  SparseArray <int, 1>               elemMatMap_; // hold the mat index for all elements

  SparseArray <int, 2>      ipMpMap_;  // mapping between integration point and 
                                       // material point (where MaterialModel stores history)
                                       // ip counts from 0 to number of Gauss points/element
                                       // mp counts from 0 to number of Gauss points/element 
                                       // multiplied with the number of elements of this mat.

  // averaging stuffs

  Ref<Function>             weightFunc_;
  double                    radius_;
  String                    averaging_;

  // pointers to correct functions to compute B and N matrices

  ShapeGradsFunc            getShapeGrads_;
  ShapeFunc                 getShapeFuncs_;

  bool                      updated_;
  bool                      isStep_;

  /* 
   * Integer vector holds 1 or zero to indicate if element ID is 
   * active or not. Used to remove fully damaged element.
   */

  IntVector                 isActive_;
};

#endif
