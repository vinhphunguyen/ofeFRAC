/*
 *
 * Model for line constraint. Designed for use in EnergyArclenModel, 
 * such that force control can be used with a uniform displacement along
 * a boundary. 
 * And for Dirichlet boundaries in combination with LaminateMeshModule
 *
 *   - definition of master-slave relations
 *   - evaluation of external force vector
 *
 * 
 * Frans van der Meer, August 2008
 *
 */

#include <jem/base/array/utilities.h>
#include <jem/base/array/select.h>
#include <jem/base/System.h>
#include <jem/util/Properties.h>
#include <jem/base/Error.h>

#include <jive/fem/NodeGroup.h>
#include <jive/model/Actions.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/util/Assignable.h>
#include <jive/util/Constraints.h>
#include <jive/util/XDofSpace.h>

#include "UserModels.h"
#include "module/SolverNames.h"

using namespace jem;

using jem::Error;
using jem::util::Properties;
using jem::io::endl;
using jive::Vector;
using jive::IntVector;
using jive::StringVector;
using jive::fem::NodeGroup;
using jive::fem::NodeSet;
using jive::model::Model;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::util::XDofSpace;

class DirichletModel : public Model
{
 public:

  typedef DirichletModel Self;
  typedef Model          Super;

  static const char*         NODES_PROP;
  static const char*         DOF_PROP;
  static const char*         LOADED_PROP;
  static const char*         SIGN_PROP;

                       DirichletModel

    ( const String&       name,
      const Properties&   conf,
      const Properties&   props,
      const Properties&   globdat );

  virtual void         configure

    ( const Properties&   props,
      const Properties&   globdat );

  virtual void         getConfig

    ( const Properties&   conf,
      const Properties&   globdat )      const;

  virtual bool         takeAction

    ( const String&       action,
      const Properties&   params,
      const Properties&   globdat );

 protected:

  virtual              ~DirichletModel ();

 private: // functions

  void                 initConstraints_
    ( const Properties&   globdat );

 private: // variables

  StringVector            nodeGroups_;
  StringVector            dofTypes_;
  Assignable<NodeSet>     nodes_;
  Ref<XDofSpace>          dofs_;
  Ref<Constraints>        cons_;
  IntVector               slaves_;
  int                     master_;
  int                     loaded_;
  int                     ngroups_;
  int                     sign_;
  String                  loadedName_;
};

//-----------------------------------------------------------------------
//   static constants
//-----------------------------------------------------------------------

const char* DirichletModel::NODES_PROP    = "nodeGroups";
const char* DirichletModel::DOF_PROP      = "dofs";
const char* DirichletModel::LOADED_PROP   = "loaded";
const char* DirichletModel::SIGN_PROP     = "sign";

//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------


DirichletModel::DirichletModel

   ( const String&       name,
     const Properties&   conf,
     const Properties&   props,
     const Properties&   globdat ) : Super(name)
{
  Properties  myProps = props.getProps ( myName_ );
  Properties  myConf  = conf.makeProps ( myName_ );

  const String context = getContext();

  // get names of nodegroups

  myProps.get( nodeGroups_, NODES_PROP );
  myConf.set ( NODES_PROP, nodeGroups_ );

  ngroups_ = nodeGroups_.size ( );

  // get names of dof

  // can have many node groups with the same dof type

  myProps.get( dofTypes_, DOF_PROP );

  if ( dofTypes_.size() == 1 && ngroups_ > 1 )
  {
    dofTypes_.reshape( ngroups_ );

    for ( int ig = 0; ig < ngroups_; ++ig )
    {
      dofTypes_[ig] = dofTypes_[0];
    }
  }

  myConf.set ( DOF_PROP, dofTypes_ );

  // get index of boundary which is loaded

  loaded_ = -1;

  myProps.find ( loaded_, LOADED_PROP );
  myConf.set   ( LOADED_PROP, loaded_ );

  sign_ = 1.;

  myProps.find ( sign_, SIGN_PROP );
  myConf.set   ( SIGN_PROP, sign_ );
}

DirichletModel::~DirichletModel()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void DirichletModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void DirichletModel::getConfig 

  ( const Properties& conf,
    const Properties& globdat ) const

{
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool DirichletModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;

  if ( action == Actions::INIT )
  {
    initConstraints_ ( globdat );

    return true;
  }

  if ( action == Actions::GET_EXT_VECTOR )
  {
    Vector     fext;

    if ( !  params.find ( fext,  ActionParams::EXT_VECTOR   ) )
    {
      // sometimes this happens,, viz. when BasicRunData invokes 
      // takeAction( Actions::NEW_MATRIX_0, ... )

      return false;
    }

    double scale = 1.;

    globdat.find ( scale, ActionParams::SCALE_FACTOR );

    fext[master_] = sign_ * scale; 

    System::out() << "external force1: " << fext[master_] << "\n";

    return true;
  }

  if ( action == SolverNames::GET_UNIT_FEXT )
  {
    Vector     fext;

    params.get ( fext,  ActionParams::EXT_VECTOR );

    fext[master_] = sign_; 

    System::out() << "unit external force: " << fext[master_] << "\n";

    return true;
  }

  return false;
}

//-----------------------------------------------------------------------
//   initConstraints
//-----------------------------------------------------------------------

void DirichletModel::initConstraints_

  ( const Properties&  globdat )

{
  Assignable<NodeGroup> group;

  int        itype;
  int        idof;

  IntVector  inodes;

  nodes_ = NodeSet::find( globdat );
  dofs_  = XDofSpace::get ( nodes_.getData(), globdat );
  cons_  = Constraints::get ( dofs_, globdat );

  const int noGroup = nodeGroups_.size ( );

  for ( int ig = 0; ig < noGroup; ig++ )
  {
    // get information for this node group

    group  = NodeGroup::find ( nodeGroups_[ig], nodes_, globdat );

    if ( group == NIL )
    {
      System::err() << "DirichletModel : group == NIL!" << endl;
      // inodes = stdBoundaryNodes_ ( nodeGroups_[ig] );
    }

    int nn     = group.size ();

    inodes . resize ( nn );

    inodes = group.getIndices ();

    // get index of specified dof type

    itype  = dofs_->findType ( dofTypes_[ig] );

    // get indices of master and slave dofs and add constraints

    if ( ig == loaded_ )
    {
      // the first node of the group is defined as master where the force
      // is applied.

      master_ = dofs_->getDofIndex ( inodes[0], itype );

      // the remain nodes are slave nodes

      slaves_ . resize ( nn - 1 );

      for ( int in = 1; in < nn; ++in )
      {
        idof = dofs_->getDofIndex ( inodes[in], itype );

        cons_->addConstraint ( idof, master_, 1.0 );

        slaves_[in-1] = idof;

        // System::out() << "adding slave constraint on node " <<
           // inodes[in] << ", dof " << dofTypes_[ig] << endl;
      }
    }
    else
    {
      for ( int in = 0; in < nn; ++in )
      {
        idof = dofs_->getDofIndex ( inodes[in], itype );

        cons_->addConstraint ( idof, 0.0 );

        // System::out() << "adding constraint on node " <<
           // inodes[in] << ", dof " << dofTypes_[ig] << endl;
      }
    }
  }

  // compress for more efficient storage

  cons_->compress();
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newDirichletModel
//-----------------------------------------------------------------------


static Ref<Model>     newDirichletModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<DirichletModel> ( name, conf, props, globdat );
}
//-----------------------------------------------------------------------
//   declareDirichletModel
//-----------------------------------------------------------------------

void declareDirichletModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "Dirichlet", & newDirichletModel );
}
