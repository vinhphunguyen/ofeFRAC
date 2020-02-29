/*
 *   Class AdaptLoadScaleModel: used in concert with module StepModule2
 *   in the sense that StepModule2 throws actions reduce step size
 *   when divergence of the NonlinModule (Newton-Raphson) was detected.
 *   Originally written by Frank Radtke and maintained by Vinh Phu Nguyen.
 *
 *   State:
 *   1. 2D only, extension to 3D is easy.
 *   2. Only constraints (for a node group) in one direction. For problems with
 *      BCs in more than one direction, use more than one model of this type,
 *      each for one direction.
 */

#ifndef ADAPT_LOADSCALE_MODEL_H
#define ADAPT_LOADSCALE_MODEL_H


#include <iostream>

#include <jem/base/String.h>

#include <jive/model/Model.h>
#include <jive/util/Assignable.h>
#include <jive/fem/NodeGroup.h>
#include <jive/fem/NodeSet.h>
#include <jive/util/Constraints.h>
#include <jem/base/Array.h>
#include <jive/util/XDofSpace.h>


using namespace jem;
using namespace jem::io;
using jem::idx_t;
using jem::Array;
using jem::util::Properties; 

using jive::model::Model;
using jive::util::Assignable;
using jive::fem::NodeSet;
using jive::fem::NodeGroup;
using jive::util::Constraints;
using jive::IntVector;
using jive::StringVector;
using jive::util::XDofSpace;


class AdaptLoadScaleModel : public Model
{
 public:

  static const char*        TYPE_NAME;
  static const char*        START_DISP;
  static const char*        STEP_SIZE;
  static const char*        RED;
  static const char*        STEPS_TO_ENLARGE;
  static const char*        LARGER;
  static const char*        MAX_STEP_SIZE;
  static const char*        MIN_STEP_SIZE;
  static const char*        ENLARGE_AFTER_STEPS;
  static const char*        MAX_DISP;


                             AdaptLoadScaleModel

    (  const String&          name,
       const Properties&      conf,
       const Properties&      props,
       const Properties&      globdat );


 protected:

  virtual                   ~AdaptLoadScaleModel  ();

 public:

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );


  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       props,
      const Properties&       globdat )        const;

 private:

   void                     initialize_

       ( const Properties&   globdat );

 private:

  //to store the nodes that are constrained
  
  Assignable<NodeSet>       nodes_;

  //first displacement applied
 
  double                    start_disp;

  //step size
  
  double                    step_size;

  //disp that is applied to constrained nodes
  
  double                    app_disp;

 //disp that is applied to constrained nodes of previous converged step
  
  double                    app_disp_old;

  //counter to count converged steps
  
  int                       conv_count;

  //dof space and constraints needed
  
  Ref<XDofSpace>            dofs_;
  Ref<Constraints>          cons_;

  //dofs to be constrained
  
  Array<idx_t>                 conDofs1_;
  Array<idx_t>                 conDofs2_;

  StringVector              nodeGroups_; //vector to store the doftype numbers
  StringVector              dofTypes_; //vector to store the doftype numbers
                                      //x - 0, y - 1, z-2
  //reduction factor
  
  double                    red_;

  //max step nr below that step enlarge
  //means:
  //if there are less iterations needed to converge than enlarge stepsize
  
  int                       enlarge_thres;
  int                       enlarge_after_steps; //start enlarging after number of steps

  double                    larger;
  double                    max_step_size;
  double                    min_step_size;
  double                    max_disp; //maximum enforced displacement to load
  double                    so_far_app_disp; //so far apllied displacement
};



#endif
