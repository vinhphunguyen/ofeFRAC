/*
 * Class DofConstraint.
 * Used to constraint nodes together i.e. enforce them to have
 * the same displacements.
 *
 * Usage:
 *
 * cons={
 *   type    = "DofConstraints";
 *   nodeGrp = "bottom";
 *   dofs    = [0,1]; // dx and dy 
 * };
 *
 * Vinh Phu Nguyen
 * Monash University 2016.
 */

#ifndef DOF_CONSTRAINT_MODEL_H
#define DOF_CONSTRAINT_MODEL_H


#include <jive/util/Assignable.h>
#include <jive/model/Model.h>
#include <jive/fem/ElementSet.h>
#include <jive/fem/NodeSet.h>

#include "femodel/model_import.h"


namespace jive
{
   namespace util
   {
      class Constraints;
   }
}

using jive::util::Assignable;
using jive::util::Constraints;


class DofConstraintModel: public Model
{
   public:
  
   static const char*        CONS_PROP;

                            DofConstraintModel  ();

                            DofConstraintModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );

  static Ref<Model>         makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  static void               declare       ();

 protected:

  virtual                  ~DofConstraintModel  ();

 private:

     void                   doConstraints_ ();

     
 private:
  
     Assignable<ElementSet>    elems_;
     Assignable<NodeSet>       nodes_;
     Ref<XDofSpace>            dofs_;
     Ref<Constraints>          cons_;

     bool                      doneCons_;

     IntVector                 inodes_;
     IntVector                 idofs_;

     int                       nn_;
};

#endif
