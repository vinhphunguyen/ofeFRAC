/*
 * 
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


class DirichletInterfaceModel : public Model
{
 public:

  typedef DirichletInterfaceModel Self;
  typedef Model                   Super;

  static const char*         NODES_PROP;
  static const char*         DOF_PROP;
  static const char*         LOADED_PROP;
  static const char*         SIGN_PROP;

                       DirichletInterfaceModel

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

  virtual              ~DirichletInterfaceModel ();

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

//  int                     master_;
  int                     loaded_;
  int                     ngroups_;

  int                     sign_;

  double                  dispIncr_;
  double                  val_;
};


//-----------------------------------------------------------------------
//   static constants
//-----------------------------------------------------------------------

const char* DirichletInterfaceModel::NODES_PROP    = "nodeGroups";
const char* DirichletInterfaceModel::DOF_PROP      = "dofs";
const char* DirichletInterfaceModel::LOADED_PROP   = "loaded";
const char* DirichletInterfaceModel::SIGN_PROP     = "sign";

//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------


DirichletInterfaceModel::DirichletInterfaceModel

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
  
  myProps.find ( dispIncr_, "dispIncr");
  
  val_   = dispIncr_;
  
  nodes_ = NodeSet::find    ( globdat                   );
  dofs_  = XDofSpace::get   ( nodes_.getData(), globdat );
  cons_  = Constraints::get ( dofs_, globdat            );
}

DirichletInterfaceModel::~DirichletInterfaceModel()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void DirichletInterfaceModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void DirichletInterfaceModel::getConfig 

  ( const Properties& conf,
    const Properties& globdat ) const

{
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool DirichletInterfaceModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;

  if ( action == Actions::GET_EXT_VECTOR )
  {
    Vector     fext;

    if ( !  params.find ( fext,  ActionParams::EXT_VECTOR   ) )
    {
      // sometimes this happens,, viz. when BasicRunData invokes 
      // takeAction( Actions::NEW_MATRIX_0, ... )

      return false;
    }

    initConstraints_ ( globdat );
    
    val_ += dispIncr_;
    
    return true;
  }



  return false;
}

//-----------------------------------------------------------------------
//   initConstraints
//-----------------------------------------------------------------------

void DirichletInterfaceModel::initConstraints_

  ( const Properties&  globdat )

{
  Assignable<NodeGroup> group;

  int        itype;
  int        idof;

  Array<idx_t>  inodes;
  Array<idx_t>  jdofs;

  const int noGroup = nodeGroups_.size ( );

  for ( int ig = 0; ig < noGroup; ig++ )
  {
    // get information for this node group

    group  = NodeGroup::find ( nodeGroups_[ig], nodes_, globdat );

    if ( group == NIL )
    {
      System::err() << "DirichletInterfaceModel : group == NIL!" << endl;
    }

    int nn     = group.size ();

    inodes . resize ( nn );

    inodes = group.getIndices (); System::out() << inodes << "\n";

    // get index of specified dof type

    itype  = dofs_->findType ( dofTypes_[ig] );

    for ( int in = 0; in < nn; in++ )
    {
      idof = dofs_->getDofIndex ( inodes[in], itype );
      
      if ( cons_->isSlaveDof ( idof ) )
      {
        cons_->getMasterDofs ( jdofs, idof );
        
        System::out() << cons_->masterDofCount(idof) << "\n";
        
        if ( jdofs.size () > 1 ) System::err() << "There are more than 1 master dof!" << endl;
        
        cons_->addConstraint ( idof,     val_ );
        cons_->addConstraint ( jdofs[0], val_ );
      }
      else
      {
        cons_->addConstraint ( idof, val_ );
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
//   newDirichletInterfaceModel
//-----------------------------------------------------------------------


static Ref<Model>     newDirichletInterfaceModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<DirichletInterfaceModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareDirichletInterfaceModel
//-----------------------------------------------------------------------


void declareDirichletInterfaceModel ()

{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "DirichletInterface", & newDirichletInterfaceModel );
}

