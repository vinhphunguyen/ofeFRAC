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


#include "NodeConstraint.h"

using jem::io::FileWriter;
using jem::io::PrintWriter;
using jive::util::Printer;


//=======================================================================
//   class NodeConstraintModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  NodeConstraintModel::CONS_PROP     = "constraints";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


NodeConstraintModel::NodeConstraintModel

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

  // constraints={
  //   edges = ["edge1", "edge2"];
  //   edge1 = {
  //      nGrp1 = "leftNodes1"; 
  //      nGrp2 = "leftNodes2"; 
  //      dofs  = [0,1];
  //   };
  // };

  Properties  symProps, edgeProps;
  
  doneCons_ = false;
  
  Array<idx_t>               dofs;
  StringVector            edges;
  String                  grp1, grp2;
  Assignable<NodeGroup>   edge;

  myProps.get  ( edges,  "edges" );
  myConf .set  ( "edges", edges  );

  const int edgeCount = edges.size ();
  
  nodeGrpPairs_.resize ( edgeCount );

  for ( int i = 0; i < edgeCount; i++ )
  {
     myProps.get  ( edgeProps, edges[i] );
     //System::out() << edgeProps;
  
     edgeProps.get  ( grp1,  "nGrp1" );
     edgeProps.get  ( grp2,  "nGrp2" );
     edgeProps.get  ( dofs,  "dofs"  );
     
     edge = NodeGroup::find ( grp1, nodes_, globdat );
     Array<idx_t> inodes1 = edge.getIndices ();

     nodeGrpPairs_[i].grp1.resize ( inodes1.size() );
     nodeGrpPairs_[i].grp1 = inodes1; 
     
     System::out() << dofs << "\n";

     edge = NodeGroup::find ( grp2, nodes_, globdat );
     Array<idx_t> inodes2 = edge.getIndices ();
     
     nodeGrpPairs_[i].grp2.resize ( inodes2.size() );
     nodeGrpPairs_[i].grp2 = inodes2;

     nodeGrpPairs_[i].dofs.resize ( dofs.size() );
     nodeGrpPairs_[i].dofs = dofs;
  }
}


NodeConstraintModel::~NodeConstraintModel ()
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool NodeConstraintModel::takeAction

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

void NodeConstraintModel::doConstraints_ ()

{
   cons_->printTo ( Printer::get() );
   Printer::get() << "\n\n";
   Printer::flush ();

   System::out() << "Constrainting the dofs ...\n"; 
   
   int       inode1, inode2;
   int       idof1,  idof2;

   const int edgeCount = nodeGrpPairs_.size ();

   for (int i = 0; i < edgeCount; i++)
   {
      NodeGrpPair_ npair = nodeGrpPairs_[i];
      const int nodeEdgeCount = npair.grp1.size();

      //System::out() << npair.grp1 << "\n";
      //System::out() << npair.grp2 << "\n";
         
      for (int in = 0; in < nodeEdgeCount; in++)
      {
         inode1 = npair.grp1[in];
         inode2 = npair.grp2[in];
      
         for (int id = 0; id < npair.dofs.size(); id++)
         {
            idof1 = dofs_->getDofIndex ( inode1, npair.dofs[id] );
            idof2 = dofs_->getDofIndex ( inode2, npair.dofs[id] );

                 //(!cons_->isSlaveDof ( idof2 ) ) )
            cons_->addConstraint ( idof1, idof2, 1. ); 
         }
      }
   }

   doneCons_ = true;
   System::out() << "Constrainting the dofs ...done\n"; 

   cons_->printTo ( Printer::get() );
   Printer::get() << "\n\n";
   Printer::flush ();
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newIGAThinShellModel
//-----------------------------------------------------------------------


static Ref<Model>     newNodeConstraintModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<NodeConstraintModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareSolidModel
//-----------------------------------------------------------------------


void declareNodeConstraintModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "NodeConstraints", & newNodeConstraintModel );
}


