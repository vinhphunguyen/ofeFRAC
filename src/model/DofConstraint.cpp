#include <cassert>
#include <jem/base/Array.h>
#include <jem/base/Tuple.h>
#include <jem/base/System.h>
#include <jem/base/ClassTemplate.h>
#include <jem/io/FileWriter.h>
#include <jem/io/PrintWriter.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/util/Properties.h>

#include <jive/util/Constraints.h>
#include <jive/util/utilities.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/Globdat.h>
#include <jive/util/Printer.h>
#include <jive/model/Actions.h>
#include <jive/model/ModelFactory.h>
#include <jive/fem/NodeSet.h>
#include <jive/fem/NodeGroup.h>


#include "DofConstraint.h"

using jem::io::FileWriter;
using jem::io::PrintWriter;
using jive::util::Printer;


//=======================================================================
//   class DofConstraintModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  DofConstraintModel::CONS_PROP     = "constraints";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


DofConstraintModel::DofConstraintModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Model ( name )

{
  using jive::util::joinNames;
  using jive::fem::NodeGroup;

  Properties  myProps = props.getProps ( myName_ );
  Properties  myConf  = conf.makeProps ( myName_ );

  String     context  = getContext ();

  elems_   = ElementSet::get ( globdat, context );
  nodes_   = elems_.getNodes ();
  dofs_    = XDofSpace::get   ( nodes_.getData(), globdat );
  cons_    = Constraints::get ( dofs_, globdat );

  String        groupName;

  myProps.get ( groupName, "nodeGrp" );
  myProps.get ( idofs_,    "dofs" );

  myConf.set ( "nodeGrp", groupName );
  myConf.set ( "dofs"   , idofs_    );

  NodeGroup nGroup = NodeGroup::get ( groupName, nodes_, globdat, context );
  nn_ = nGroup.size ( );
  inodes_.resize ( nn_ );
  inodes_ = nGroup.getIndices ( );

  System::out() << inodes_ << "\n";
   
  doneCons_ = false;
}


DofConstraintModel::~DofConstraintModel ()
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool DofConstraintModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::util::Globdat;

  if ( action == Actions::GET_CONSTRAINTS )
  {
     if ( !doneCons_ ) 
     {
        doConstraints_ ( );
     }

     //cons_->printTo ( Printer::get() );
     //Printer::get() << "\n\n";
     //Printer::flush ();
      
    return true;
  }
    
  return false;
}

//-----------------------------------------------------------------------
//   doConstraints_
//-----------------------------------------------------------------------

void DofConstraintModel::doConstraints_ ()

{
   //cons_->printTo ( Printer::get() );
   //Printer::get() << "\n\n";
   //Printer::flush ();
   
   int       inodeM, inode2;
   int       idofM,  idof2;

   inodeM = inodes_[0];  // master node

   for (int i = 1; i < nn_; i++)
   {
      inode2 = inodes_[i];  // slave node
      for (int id = 0; id < idofs_.size(); id++)
      {
         idofM = dofs_->getDofIndex ( inodeM, idofs_[id] );
         idof2 = dofs_->getDofIndex ( inode2, idofs_[id] );

         cons_->addConstraint ( idof2, idofM, 1. ); 
      }
   }

   doneCons_ = true;

   //System::out() << "Constrainting the dofs ...done\n"; 
   cons_->printTo ( Printer::get() );
   Printer::get() << "\n\n";
   Printer::flush ();
}

//=======================================================================
//   related functions
//=======================================================================


static Ref<Model>     newDofonstraintModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<DofConstraintModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareSolidModel
//-----------------------------------------------------------------------


void declareDofConstraintModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "DofConstraints", & newDofonstraintModel );
}


