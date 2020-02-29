//=======================================================================
//
// Model that prints averaged load-displacement data for node group or
// nodes to file
// FPM, 29-1-2008
// Modified:
//  VP Nguyen: add evalMasterDofs(fint,cons_) so that when the nodes on which
//  external forces are applied are master nodes, their internal forces are 
//  correctly defined.
//
// Data is computed for GET_MATRIX0 and GET_INT_VECTOR
//         and printed for CHECK_COMMIT
//
//=======================================================================

#include <jem/io/PrintWriter.h>
#include <jem/io/FileWriter.h>
#include <jem/base/array/select.h>
#include <jem/base/array/utilities.h>
#include <jem/base/System.h>
#include <jem/util/Properties.h>
#include <jem/mp/utilities.h>

#include <jive/util/Globdat.h>
#include <jive/util/ItemSet.h>
#include <jive/util/DofSpace.h>
#include <jive/util/Constraints.h>
#include <jive/util/Assignable.h>
#include <jive/mp/Globdat.h>
#include <jive/mp/ItemMask.h>
#include <jive/model/Actions.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/StateVector.h>
#include <jive/fem/ElementSet.h>
#include <jive/fem/NodeGroup.h>

#include "UserModels.h"

using namespace jem;
using jem::io::endl;

using jem::util::Properties;
using jem::io::PrintWriter;
using jem::io::FileWriter;
using jem::mp::allsum;
using jem::mp::Context;
using jive::Vector;
using jive::IdxVector;
using jive::IdxMatrix;
using jive::StringVector;
using jive::util::DofSpace;
using jive::util::Assignable;
using jive::util::evalMasterDofs;
using jive::util::Constraints;
using jive::model::Model;
using jive::fem::NodeSet;
using jive::fem::NodeGroup;
using jive::fem::ElementSet;


//=======================================================================
//   class LodiModel
//=======================================================================

class LodiModel : public Model
{
 public:

  static const char*        NODES;
  static const char*        GROUP;
  static const char*        FILE_NAME;


                            LodiModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       conf,
      const Properties&       globdat )      const;

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );


 protected:

  virtual                  ~LodiModel  ();


 private:

  Ref<Context>              mpx_;
  Ref<DofSpace>             dofs_;
  Ref<Constraints>          cons_;
  Ref<PrintWriter>          out_;
  Assignable<ElementSet>    elems_;
  Assignable<NodeSet>       nodes_;

  IdxVector                 inodes_;
  IdxMatrix                 idofs_;

  idx_t                     nn_;
  idx_t                     ndof_;

  String                    file_;

  bool                      toFile_;
  bool                      evalMasterDof_;
};


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  LodiModel::NODES       = "nodes";
const char*  LodiModel::GROUP       = "group";
const char*  LodiModel::FILE_NAME   = "file";


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


LodiModel::LodiModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Model ( name ), nn_(-1)

{
  using jive::mp::Globdat;

  const String   context = getContext ();

  mpx_   = Globdat::getMPContext ( globdat );
  elems_ = ElementSet::get ( globdat, context );
  nodes_ = elems_.getNodes ( );
  dofs_  = DofSpace::get   ( nodes_.getData(), globdat, context );
  cons_  = Constraints::get ( dofs_, globdat );

  ndof_  = dofs_->typeCount();

  toFile_         = false;
  evalMasterDof_  = false;
}


LodiModel::~LodiModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void LodiModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  using jive::util::Globdat;
  using jive::mp::ItemMask;

  Properties    myProps = props.findProps ( myName_ );

  const String  context = getContext ();

  IdxVector     nodeIDs;
  IdxVector     jtypes ( ndof_ );

  String        groupName;

  Properties    myVars = Globdat::getVariables ( myName_, globdat );

  bool nodes = myProps.find ( nodeIDs,   NODES );
  bool group = myProps.find ( groupName, GROUP );
  
  myProps.find ( evalMasterDof_, "evalMasterDof" );

  if ( nodes )
  {
    // Get node and dof numbers specified with 'nodes = IntArray'

    nn_ = nodeIDs.size ( );

    if ( allsum( *mpx_, nn_ ) == 0 )
    {
      myProps.propertyError ( context, "there is no node specified!!!" );
    }

    inodes_.resize ( nn_ );

    for ( int i = 0; i < nn_; i++ )
    {
      inodes_[i] = nodes_.findNode ( nodeIDs[i] );
    }
  }

  if ( group )
  {
    // Get node and dof numbers specified with 'group = String'

    NodeGroup nGroup = NodeGroup::get ( groupName, nodes_, globdat, context );

    nn_ = nGroup.size ( );

    if ( allsum( *mpx_, nn_ ) == 0 )
    {
      myProps.propertyError ( context, "there is no node in group!!!" );
    }

    inodes_.resize ( nn_ );

    inodes_ = nGroup.getIndices ( );
  }

  if ( !nodes & !group )
  {
     myProps.propertyError ( context,
			     "must specify nodes or node group!!!" );
  }

  if ( nodes || group )
  {
    Ref<ItemMask>  imask = ItemMask::get ( nodes_.getData(), globdat );

    // Remove nodes that are "owned" by another process.

    int  j = 0;

    for ( int i = 0; i < nn_; i++ )
    {
      int  inode = inodes_[i];

      if ( imask->getValue( inode ) )
      {
        inodes_[j++] = inode;
      }
    }

    if ( j < nn_ )
    {
      nn_ = j;

      inodes_.reshape ( nn_ );
    }

    for ( int i = 0; i < ndof_; i++ )
    {
      jtypes[i] = i;
    }

    idofs_.resize ( ndof_, nn_ );

    idofs_ = -1;

    dofs_->findDofIndices ( flatten( idofs_ ), inodes_, jtypes );
  }

  if (  myProps.find( file_, FILE_NAME ) )
  {
    // Only process 0 writes to the output file.

    if ( mpx_->myRank() == 0 )
    {
      toFile_ = true;
      out_    = newInstance<PrintWriter> ( newInstance<FileWriter> ( file_ ) );
    }
  }
}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void LodiModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties    myConf  = conf.makeProps ( myName_ );
  
  myConf. set ( "evalMasterDof", evalMasterDof_ );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool LodiModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::util::Globdat;
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;
  using jive::model::STATE2;

  Vector      load ( ndof_ ) ;
  Vector      disp ( ndof_ ) ;

  if ( action == Actions::GET_MATRIX0 ||
       action == Actions::GET_INT_VECTOR )
  {
    Properties  myVars = Globdat::getVariables ( myName_, globdat );

    // get global data

    Vector      fint;
    Vector      state;
    Vector      a;

    params.get       ( fint, ActionParams::INT_VECTOR );
    StateVector::get ( state, dofs_, globdat          );
    StateVector::get ( a, STATE2, dofs_, globdat      );

    Vector fint1 = fint.clone();

    
    Vector      mass;
    bool hasMass = globdat.find ( mass, "massVector" );

    // if ( !hasMass )
    // {
    //   mass.resize ( state.size() );
    //   mass = 0.;
    // }
    

    // do not forget the internal force at independent
    // dofs is increased by an amount = C^T * fint(dependent dof)

    if ( evalMasterDof_ )  evalMasterDofs   ( fint1, *cons_ );

    // compute cumulative load and average displacement

    idx_t  nnGlobal = allsum ( *mpx_, nn_ );

    load = 0.0;
    disp = 0.0;

    for ( idx_t i = 0; i < ndof_; i++ )
    {
      for ( idx_t j = 0; j < nn_; j++ )
      {
        idx_t  idof = idofs_(i,j);

        if ( idof >= 0 )
        {
          load[i] += fint1[idof];
          //System::out() << mass[idof] << " " << a[idof] << "\n";
          disp[i] += state[idof];
        }
      }

      // Sum over all processes. Note that all processes have the same
      // dofs.

      load[i] = allsum ( *mpx_, load[i] );
      disp[i] = allsum ( *mpx_, disp[i] ) / nnGlobal;
    }

    // store in globdat

    myVars.set ( "load" , load );
    myVars.set ( "disp" , disp );

    return true;
  }

  // when converged, write to file

  if ( action == "COMMIT" )
  {
    Properties  myVars = Globdat::getVariables ( myName_, globdat );

    // get internal force vector from myVars

    if ( toFile_ )
    {
      myVars.get  ( load , "load" );
      myVars.get  ( disp , "disp" );
      
      double time(0.);

      if ( globdat.find ( time, Globdat::TIME ) ){
        print ( *out_, time  );
        print ( *out_, "   " );
      }

      for ( int i = 0; i < ndof_; i++ )
      {
        *out_ << String::format( "%12.8f  ", disp[i] );
        *out_ << String::format( "%12.8f  ", load[i] );
      }

      out_->printLine();
      out_->flush();
    }
    else
    {
      // System::out() << myVars << endl;
    }

    return true;
  }

  return false;
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newLodiModel
//-----------------------------------------------------------------------


static Ref<Model>     newLodiModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<LodiModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareLodiModel
//-----------------------------------------------------------------------


void declareLodiModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "Lodi", & newLodiModel );
}
