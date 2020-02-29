/*
 *  Copyright (C) 2017 Monash University. All rights reserved.
 *
 *  This class implements a model to store the displacement at 1 node in the globdal
 *  databse, which can be retrieved at anytime by any model/module.
 *
 *  This model is working in parallel. Using allsum (...).
 *
 *  Input data block of this model:
 *
 *     monitor =
       {
         type     = "Monitor";
         node     = 4;
       };
 *
 *  Author: V.P. Nguyen, phu.nguyen@monash.edu
 *  Date: 30 November 2017
 *
 *  Updates (what and who):
 *
 */

#include <jem/base/System.h>
#include <jem/util/Properties.h>
#include <jem/mp/utilities.h>
#include <jive/util/Globdat.h>
#include <jive/util/ItemSet.h>
#include <jive/util/DofSpace.h>
#include <jive/util/Assignable.h>
#include <jive/mp/Globdat.h>
#include <jive/mp/ItemMask.h>
#include <jive/model/Actions.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/StateVector.h>
#include <jive/fem/ElementSet.h>

#include "model/UserModels.h"


using namespace jem;

using jem::util::Properties;
using jem::mp::allsum;
using jem::mp::Context;
using jive::Vector;
using jive::IdxVector;
using jive::StringVector;
using jive::util::DofSpace;
using jive::util::Assignable;
using jive::model::Model;
using jive::fem::NodeSet;
using jive::fem::ElementSet;


//=======================================================================
//   typedefs
//=======================================================================

// Some handy aliases to avoid some typing

typedef ElementSet          ElemSet;


//=======================================================================
//   class MonitorModel
//=======================================================================


class MonitorModel : public Model
{
 public:

  typedef MonitorModel      Self;
  typedef Model             Super;

  static const char*        NODE_PROP;


                            MonitorModel

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

  virtual                  ~MonitorModel  ();


 private:

  Ref<Context>              mpx_;
  Ref<DofSpace>             dofs_;
  Assignable<ElemSet>       elems_;
  Assignable<NodeSet>       nodes_;

  int                       inode_;
  IdxVector                 idofs_;
  StringVector              dofTypes_;

};


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  MonitorModel::NODE_PROP = "node";


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


MonitorModel::MonitorModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super ( name )

{
  using jive::mp::Globdat;
  
  const String  context = getContext ();
  
  mpx_   = Globdat::getMPContext ( globdat );

  elems_ = ElemSet::get  ( globdat, context );
  nodes_ = elems_.getNodes ();
  dofs_  = DofSpace::get ( nodes_.getData(), globdat, context );
  inode_ = -1;
}


MonitorModel::~MonitorModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void MonitorModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps = props.findProps ( myName_ );

  idx_t nodeID;

  if ( myProps.find( nodeID, NODE_PROP ) )
  {
    inode_ = nodes_.findNode ( nodeID );
    
    double  count = 0.0;

    if ( inode_ >= 0 )
    {
      count = 1.0;
    }

    count = allsum ( *mpx_, count );

    if ( count <= 0.0 ) 
    {
      throw jem::IllegalInputException (
        JEM_FUNC,
        "invalid node specified (not found)"
      );
    }
  }
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void MonitorModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool MonitorModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::util::Globdat;
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;


  if ( action == Actions::GET_MATRIX0 ||
       action == Actions::GET_INT_VECTOR )
  {
    Vector      state;

    StateVector::get ( state, dofs_, globdat          );
    
    double  ux    = 0.0;
    double  uy    = 0.0;

    if ( inode_ >= 0 )
    {
      idx_t  idof0 = dofs_->getDofIndex ( inode_, 0 );
      idx_t  idof1 = dofs_->getDofIndex ( inode_, 1 );

      ux    = state[idof0];  // x-displacement of crack mouth node 1
      uy    = state[idof1];  // y-displacement of crack mouth node 1
    }

    ux    = allsum ( *mpx_, ux );
    uy    = allsum ( *mpx_, uy );

    Vector disp( 2 );

    disp[0] = ux;
    disp[1] = uy;

    // store in globdat

    Properties  myVars = Globdat::getVariables ( myName_, globdat );

    myVars.set ( "disp" , disp );

    return true;
  }

  if ( action == Actions::COMMIT )
  {
    Properties  myVars = Globdat::getVariables ( myName_, globdat );

    return true;
  }

  return false;
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newMonitorModel
//-----------------------------------------------------------------------


static Ref<Model>     newMonitorModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<MonitorModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareMonitorModel
//-----------------------------------------------------------------------


void declareMonitorModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "Monitor", & newMonitorModel );
}
