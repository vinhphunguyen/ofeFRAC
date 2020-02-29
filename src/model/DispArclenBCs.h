/*
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *
 * This class implements a model to define imposed displacements
 * on certain node groups in a computation with displacement control
 * and arclength control. It is used as a child of class 
 * DispArclenModel and is able to switch between the two different
 * strategies.
 * Multiple u=0 boundaries can be defined, and ONE boundary with prescribed
 * nonzero displacements. For this loaded boundary, the boundary 
 * conditions are adapted when switching.
 *
 * Author(s):
 *
 *   Vinh Phu Nguyen, V.P.Nguyen@tudelft.nl
 *   Frans van der Meer, F.P.vanderMeer@tudelft.nl
 *
 */

#ifndef DISP_CONTROL_BCS_H
#define DISP_CONTROL_BCS_H


#include <jem/base/array/utilities.h>
#include <jem/base/array/select.h>
#include <jem/base/System.h>
#include <jem/util/Properties.h>
#include <jem/base/Error.h>

#include <jive/fem/NodeGroup.h>
#include <jive/model/Actions.h>
#include <jive/util/Assignable.h>

namespace jive
{
  namespace util
  {
    class Constraints;
    class XDofSpace;
  }
}

using namespace jem;

using jem::Error;
using jem::util::Properties;
using jem::io::endl;
using jive::Vector;
using jive::IntVector;
using jive::StringVector;
using jive::fem::NodeGroup;
using jive::fem::NodeSet;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::util::XDofSpace;

// ======================================================
//     class DispArclenBCs
// ======================================================

class DispArclenBCs : public Object
{
 public:

  typedef DispArclenBCs Self;

  static const char*         TYPE_NAME;
  static const char*         NODES_PROP;
  static const char*         DOF_PROP;
  static const char*         LOADED_PROP;

                       DispArclenBCs ();

  virtual void         configure

    ( const Properties&   myProps,
      const Properties&   globdat );

  virtual void         getConfig

    ( const Properties&   myConf,
      const Properties&   globdat )      const;

  void                 toDispControl

    ( const Properties&  params,
      const Properties&  globdat );

  void                 initConstraints

    ( const Properties&   globdat );

  void                 toArclControl

    ( const Properties&  params,
      const Properties&  globdat );

  void                 applyConstraints

    ( const Properties&  params,
      const Properties&  globdat );

  void                 getUnitLoad

    ( const Properties&  params,
      const Properties&  globdat );

  void                 storeLoadScale

    ( const Properties&   globdat,
      const Vector&       fint );

  void                 commit

    ( const bool         isDispControl,
      const Properties&  globdat );

  int                  getLoadedDofs

    ( IntVector&          idofs )        const;

  double               reduceDispIncr

    ( const double        reduction );

  inline  void         setDispValue

    ( double prescribedDisp );

  inline  void         setDispIncr

    ( double val );

  inline double        getLoadScale ()  const;

  inline double        getLoadScale0()  const;

  inline double        getDispValue ()  const;

  inline double        getDispIncr  ()  const;

  inline void          disincrement ();


 protected:

  virtual              ~DispArclenBCs ();

 private:

  Assignable<NodeSet>     nodes_;

  Ref<XDofSpace>          dofs_;
  Ref<Constraints>        cons_;

  int                     master_;

  double                  dispIncr_;
  double                  dispVal_;
  double                  dispVal0_;
  double                  loadScale_; 
  double                  loadScale0_;

  int                     ngroups_;

  /* the following members are constant input variables 
   *
   * nodeGroups   node groups for boundary conditions
   * dofTypes     dof types of boundary condition
   * loaded       index of boundary where nonzero disp is applied
   * 
   * nodeGroups.size() == dofTypes.size()
   * 0 <= loaded < nodeGroups.size()
   */

  StringVector            nodeGroups_;
  StringVector            dofTypes_;
  int                     loaded_;
};

// ======================================================
//     implementation of inline functions
// ======================================================


inline void DispArclenBCs::setDispValue

 ( double amount )
{
  dispVal_ = amount;
}

inline void DispArclenBCs::setDispIncr

 ( double val )
{
  dispIncr_ = val;
}

inline double DispArclenBCs::getLoadScale  () const 
{ 
  return loadScale_; 
}

inline double DispArclenBCs::getLoadScale0 () const 
{ 
  return loadScale0_; 
}

inline double DispArclenBCs::getDispValue  () const 
{ 
  return dispVal_; 
}

inline double DispArclenBCs::getDispIncr   () const 
{ 
  return dispIncr_;  
}

inline void DispArclenBCs::disincrement  () 
{ 
  dispVal_ = dispVal0_;
}

#endif



