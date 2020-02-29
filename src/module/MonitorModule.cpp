
#include <jem/base/System.h>
#include <jem/util/Properties.h>
#include <jive/util/Globdat.h>
#include <jive/util/ItemSet.h>
#include <jive/util/DofSpace.h>
#include <jive/model/StateVector.h>
#include <jive/model/Actions.h>

#include "MonitorModule.h"


//=======================================================================
//   class MonitorModule
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  MonitorModule::NODE_PROP = "node";

const char*  MonitorModule::FORCE_NAMES[3] = { "fx", "fy", "fz" };


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


MonitorModule::MonitorModule ( const String& name ) :

  Super ( name )

{
  nodeID_ = 0;
}


MonitorModule::~MonitorModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status MonitorModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  // Nothing to do.

  return OK;
}


//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------

// store the displacement at the monitored node or 
// store the internal force at this node => external force applied at this node
// in case of direct displacement control

Module::Status MonitorModule::run ( const Properties& globdat )
{
  using jem::Ref;
  using jem::System;
  using jem::io::endl;
  using jive::Vector;
  using jive::IdxVector;
  using jive::util::Globdat;
  using jive::util::ItemSet;
  using jive::util::DofSpace;
  using jive::model::StateVector;
  using jive::model::ActionParams;

  // A context string for error messages.

  const String   context = getContext ();

  // Get the default DofSpace.

  Ref<DofSpace>  dofs    = DofSpace::get  ( globdat, context );

  // Get the ItemSet with which the DofSpace is associated.

  Ref<ItemSet>   items   = dofs->getItems ();

  // Translate the node ID to its index.

  int            inode   = items->findItem ( nodeID_ );

  if ( inode < 0 )
  {
    System::warn() << myName_
		   << " : invalid node: " << nodeID_ << endl;

    return OK;
  }

  // Get all DOFs attached to the node.

  const int numOfDofTypes = dofs->typeCount();

  IdxVector  idofs  ( numOfDofTypes );
  IdxVector  jtypes ( numOfDofTypes );
  
  int  dofCount = dofs->getDofsForItem ( idofs, jtypes, inode );

 
  Vector  state;
  Vector  force;

  // Get the state vector.

  StateVector::get ( state, dofs, globdat );

  // Get the internal force vector

  //globdat.get ( force, ActionParams::INT_VECTOR); 

  // Get the Properties set in which the monitor values will be
  // stored.

  Properties  myVars = Globdat::getVariables ( myName_, globdat );

  // Store the monitor values.

  for ( int i = 0; i < dofCount; i++ )
  {
    double  value1 = state[idofs[i]];    // displacement   at i direction
 //   double  value2 = force[idofs[i]];    // internal force at i direction

    myVars.set ( dofs->getTypeName( jtypes[i] ), value1 );
 //   myVars.set ( FORCE_NAMES[i], -value2 );  
  }

  return OK;
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void MonitorModule::shutdown ( const Properties& globdat )
{
  // Nothing to do at the moment.
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void MonitorModule::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps = props.findProps ( myName_ );

  myProps.find ( nodeID_, NODE_PROP );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void MonitorModule::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );

  myConf.set ( NODE_PROP, nodeID_ );
}
