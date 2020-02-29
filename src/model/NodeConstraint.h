/*
 * Class NodeConstraintModel.
 * Used to constraint nodes together i.e. enforce them to have
 * the same displacements.
 * Originally written for usage in Isogeometric shell formulation where
 * constraints between control points are needed.
 *
 * Vinh Phu Nguyen
 * Cardiff 2013.
 */

#ifndef NODE_CONSTRAINT_MODEL_H
#define NODE_CONSTRAINT_MODEL_H


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


class NodeConstraintModel: public Model
{
   public:
  
   static const char*        CONS_PROP;

                            NodeConstraintModel  ();

                            NodeConstraintModel

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

  virtual                  ~NodeConstraintModel  ();

   private:

     void                   doConstraints_ ();


   private:

     struct                    NodeGrpPair_
     {
        IntVector                grp1;
        IntVector                grp2;
        IntVector                dofs;
     };
     
   private:
  
     Assignable<ElementSet>    elems_;
     Assignable<NodeSet>       nodes_;
     Ref<XDofSpace>            dofs_;
     Ref<Constraints>          cons_;

     bool                      doneCons_;
  
     Array<NodeGrpPair_>       nodeGrpPairs_;
};

#endif
