#include <jem/base/Exception.h>
#include <jem/base/System.h>
#include <jem/util/Properties.h>
#include <jive/util/utilities.h>
#include <jive/util/Globdat.h>
#include <jive/model/Actions.h>
#include <jive/model/ModelFactory.h>
#include <jem/base/IllegalOperationException.h>
#include <jive/implict/SolverInfo.h>

#include "AdaptLoadScaleModel.h"


using namespace jem;
using jive::Matrix;
using jive::fem::NodeGroup;
using jive::fem::NodeSet;



//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  AdaptLoadScaleModel::TYPE_NAME       = "AdaptLoadScale";
const char*  AdaptLoadScaleModel::START_DISP      = "startDisp";
const char*  AdaptLoadScaleModel::STEP_SIZE       = "stepSize";
const char*  AdaptLoadScaleModel::MAX_STEP_SIZE   = "maxStepSize";
const char*  AdaptLoadScaleModel::MIN_STEP_SIZE   = "minStepSize";
const char*  AdaptLoadScaleModel::RED             = "reductionFactor";
const char*  AdaptLoadScaleModel::LARGER          = "enlargment_factor";
const char*  AdaptLoadScaleModel::STEPS_TO_ENLARGE = "steps_to_enlarge";
const char*  AdaptLoadScaleModel::ENLARGE_AFTER_STEPS = "enlarge_after_steps";
const char*  AdaptLoadScaleModel::MAX_DISP        = "maxDisp"; 

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


AdaptLoadScaleModel::AdaptLoadScaleModel
(  const String&         name,
   const Properties&     conf,
   const Properties&     props,
   const Properties&     globdat )

  : Model(name), start_disp(0), step_size(0.01),
    app_disp(0), app_disp_old(0), conv_count(0),
    enlarge_thres(2), larger(2), max_step_size(1), 
    enlarge_after_steps(0),min_step_size(0.0000001),
    max_disp(1.0), so_far_app_disp(0)
{
  Properties    myConf  = conf.makeProps ( myName_ );
  Properties    myProps = props.getProps ( myName_ );
  const String  context = getContext ();
 

  //get the direction the nodes are constrained in

  myProps.get ( dofTypes_, "dofs"    );
  myConf. set ( "dofs"   , dofTypes_ );

  
  myProps.get ( nodeGroups_, "nodeGroups" );
  myConf. set ( "nodeGroups", nodeGroups_ );

  //get start displacement
 
  myProps.get ( start_disp, START_DISP );
  myConf. set ( START_DISP, start_disp );

  //set starting displacement to value

  app_disp = app_disp_old = start_disp;

  //get the step size
  
  myProps.get ( step_size, STEP_SIZE );
  myConf. set ( STEP_SIZE, step_size );

  //get reduction factor

  myProps.get ( red_, RED  );
  myConf. set ( RED,  red_ );

  //get threshold when stepsize is enlarged

  myProps.get ( enlarge_thres,    STEPS_TO_ENLARGE );
  myConf. set ( STEPS_TO_ENLARGE, enlarge_thres    );

  //get enlargment factor

  myProps.get ( larger, LARGER );
  myConf. set ( LARGER, larger );

  //get max step size factor

  myProps.get ( max_step_size, MAX_STEP_SIZE );
  myConf. set ( MAX_STEP_SIZE, max_step_size );

  //get min step size factor
  
  myProps.get ( min_step_size, MIN_STEP_SIZE );
  myConf. set ( MIN_STEP_SIZE, min_step_size );

  //get min step size factor
  
  myProps.find( enlarge_after_steps, ENLARGE_AFTER_STEPS );
  myConf. set ( ENLARGE_AFTER_STEPS, enlarge_after_steps );

  //get maximum displacement to be imposed
  
  myProps.find ( max_disp, MAX_DISP );
  myConf. set  ( MAX_DISP, max_disp );
}

AdaptLoadScaleModel::~AdaptLoadScaleModel ()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void AdaptLoadScaleModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{ } 

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void AdaptLoadScaleModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{ }


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool AdaptLoadScaleModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jem::Exception;
  using jive::util::Globdat;
  using jive::implict::SolverInfo;

  if ( action == Actions::INIT )
  {
    initialize_ ( globdat );

    return true;
  }


  if ( action == Actions::GET_EXT_VECTOR ||
       action == Actions::GET_CONSTRAINTS )
  {
      int ss = conDofs1_.size ();
      for(int i = 0; i < ss; ++i)
      {
        cons_->addConstraint ( conDofs1_[i], app_disp );
      }
      
      ss = conDofs2_.size ();
      for(int i = 0; i < ss; ++i)
      {
        cons_->addConstraint ( conDofs2_[i], app_disp );
      }

      return true;
  }
  

  if(action == "REDUCE_STEP")
  {
      if( step_size >= min_step_size)
      {
          step_size *= red_;
          app_disp   = step_size;
          System::out() << "Stepsize reduced. " << app_disp << "\n\n";
      }
      else
      {
          bool ramp = true;
          params.set("RAMP", ramp);
          //throw IllegalOperationException("AdaptLoadScaleModel","Solver did not converge. No ");     
      }

      return true;
  }

  if(action == Actions::COMMIT)
  {
      //update old disp
      //update total applied displacement
      app_disp_old     = app_disp; 
      so_far_app_disp += app_disp; 

      //enlarge applied displacement
      if(conv_count == 0)
      {
         app_disp = start_disp;
      }
      else
      {
         app_disp = step_size;
      }
      
      //possibiity for step enlarge???
      Properties  info = SolverInfo::get ( globdat );

      int iterCount = 0;
      info.get(iterCount,"iterCount");

      //if(iterCount < enlarge_thres && step_size <= max_step_size 
      //                             && conv_count > enlarge_after_steps)
      //{
      //    step_size *= larger;
      //    System::out() << "Stepsize enlarged. "<< step_size << "\n\n";
      //}

      //count one up
      conv_count++;
      return true;
    }

  if(action == "CHECK_MAX_DISP")
  {
      //now check if max disp is reached
      if(max_disp < so_far_app_disp)
      {
          bool max_disp_reached = true;
          params.set("max_disp_reached", max_disp_reached);
      }
  }

  return false;
}

// -------------------------------------------------------------------
//    initialize_
// -------------------------------------------------------------------

void AdaptLoadScaleModel::initialize_ 

      ( const Properties& globdat )
{
    // get the (all) nodes
    // get the dof space
    // get the constraints associated with the DOF space.
    nodes_  = NodeSet::find    ( globdat );
    dofs_   = XDofSpace::get   ( nodes_.getData(), globdat );
    cons_   = Constraints::get ( dofs_, globdat );

    // to fix: assume there are two node groups

    Assignable<NodeGroup>  group;
    Array<idx_t>           inodes;
    idx_t                  nn, itype;

    group = NodeGroup::find ( nodeGroups_[0], nodes_, globdat );
    nn    = group.size ();
    inodes.resize ( nn );
    inodes = group.getIndices ();
    itype  = dofs_->findType ( dofTypes_[0] );

    conDofs1_.resize ( nn );
    dofs_->getDofIndices( conDofs1_, inodes, itype );
    
    group = NodeGroup::find ( nodeGroups_[1], nodes_, globdat );
    nn    = group.size ();
    inodes.resize ( nn );
    inodes = group.getIndices ();
    itype  = dofs_->findType ( dofTypes_[1] );

    conDofs2_.resize ( nn );
    dofs_->getDofIndices( conDofs2_, inodes, itype );
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newAdaptLoadScaleModel
//-----------------------------------------------------------------------

Ref<Model>            newAdaptLoadScaleModel

( const String&       name,
  const Properties&   conf,
  const Properties&   props,
  const Properties&   globdat )

{
  return newInstance<AdaptLoadScaleModel> ( name, conf, props, globdat );
}

//-----------------------------------------------------------------------
//   declare
//-----------------------------------------------------------------------


void declareAdaptLoadScaleModel ()
{
  using jive::model::ModelFactory;
  ModelFactory::declare ( "AdaptLoadScale", newAdaptLoadScaleModel );
}



