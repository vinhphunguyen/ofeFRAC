
/*
 * This class implements a model to define imposed displacements
 * on certain node groups in a computation with displacement control
 * and arclength control. It is used as a child of class 
 * DispArclenModel and is able to switch between the two different
 * strategies.
 * Multiple u=0 boundaries can be defined, and 1 with prescribed
 * nonzero displacements. For this loaded boundary, the boundary 
 * conditions are adapted when switching.
 *
 * Author(s):
 *
 *   Vinh Phu Nguyen, V.P.Nguyen@tudelft.nl
 *   Frans van der Meer, F.P.vanderMeer@tudelft.nl
 *
 */

#include <jem/util/ArrayBuffer.h>
#include <jem/numeric/utilities.h>
#include <jem/base/IllegalInputException.h>
#include <jive/util/Globdat.h>
#include <jive/util/Printer.h>
#include <jive/util/Constraints.h>
#include <jive/util/XDofSpace.h>
#include <jive/model/StateVector.h>


#include "DispArclenBCs.h"
#include "module/SolverNames.h"


using jem::util::ArrayBuffer;
using jive::model::StateVector;
using jive::model::ActionParams;

//-----------------------------------------------------------------------
//   static constants
//-----------------------------------------------------------------------

const char* DispArclenBCs::TYPE_NAME     = "DispControl";
const char* DispArclenBCs::NODES_PROP    = "nodeGroups";
const char* DispArclenBCs::DOF_PROP      = "dofs";
const char* DispArclenBCs::LOADED_PROP   = "loaded";

//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------


DispArclenBCs::DispArclenBCs ( )
{
  // initialize variables

  master_        = 0;
  loadScale_     = 0.;
  loadScale0_    = 0.;
  dispIncr_      = 0.;
  dispVal0_      = 0.;
  dispVal_       = 0.;
  loaded_        = 0;
}

DispArclenBCs::~DispArclenBCs()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void DispArclenBCs::configure

  ( const Properties&  myProps,
    const Properties&  globdat )

{
  myProps.get( nodeGroups_, NODES_PROP );

  ngroups_ = nodeGroups_.size ( );

  // get names of dof

  // can have many node groups with the same dof type

  myProps.get( dofTypes_, DOF_PROP );

  if ( dofTypes_.size() != ngroups_ )
  {
    throw IllegalInputException ( JEM_FUNC,
          "dofTypes must have the same length as nodeGroups" );
  }

  // index of boundary with nonzero displacements

  myProps.find ( loaded_, LOADED_PROP, 0, ngroups_-1 );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void DispArclenBCs::getConfig 

  ( const Properties& myConf,
    const Properties& globdat ) const

{
  myConf.set ( NODES_PROP,  nodeGroups_ );
  myConf.set ( DOF_PROP,    dofTypes_   );
  myConf.set ( LOADED_PROP, loaded_     );
}

//-----------------------------------------------------------------------
//   toDispControl
//-----------------------------------------------------------------------

void DispArclenBCs::toDispControl

  ( const Properties&  params,
    const Properties&  globdat )

{
  // erase linear constraints defined during arclength control

  IntVector  inodes;
  Assignable<NodeGroup> group;

  int ig = loaded_;
  group  = NodeGroup::find ( nodeGroups_[ig], nodes_, globdat );

  int nn = group.size ();
  inodes . resize ( nn );
  inodes = group.getIndices ();

  int itype  = dofs_->findType ( dofTypes_[ig] );

  for ( int in = 1; in < nn; in++ )
  {
    int idof = dofs_->getDofIndex ( inodes[in], itype );

    cons_->eraseConstraint ( idof );
  }

  // set correct value for prescribed displacement 
  // there are two cases: 
  // (1) to disp control after a converged load increment: converged = true
  // (2) to disp control after a diverged load increment: converged = false

  bool converged;

  params.get ( converged, SolverNames::CONVERGED );

  if ( converged )
  {
    // fix the current displacement 

    Vector  disp;

    StateVector::get ( disp,  dofs_, globdat );

    //dispVal_ = jem::numeric::abs ( disp[master_] ); // why abs here???
    dispVal_ = disp[master_] ;
  }
  else
  {
    // apply increment to old displacement

    dispVal_ = dispVal0_ + dispIncr_;
  }

  //using jive::util::Printer;
  //cons_->printTo ( Printer::get() );
  //Printer::flush ();  System::out() << "\n\n";
}

//-----------------------------------------------------------------------
//   getLoadedDofs
//-----------------------------------------------------------------------

int  DispArclenBCs::getLoadedDofs

  ( IntVector&          idofs )        const

{
  // get indices of degrees of freedom with nonzero entries in 
  // external force vector (currently always a single dof)

  ArrayBuffer<int> itmp;

  itmp.pushBack ( master_ );

  idofs.ref ( itmp.toArray() );

  return idofs.size();
}

//-----------------------------------------------------------------------
//   reduceDispIncr
//-----------------------------------------------------------------------

double DispArclenBCs::reduceDispIncr

  ( const double       reduction )

{
  // reduce the displacement increment and use the new increment
  // to update the prescribed value for this time step

  dispIncr_ *= reduction;

  dispVal_   = dispVal0_ + dispIncr_;

  // return the new value to check for bounds in DispArclenModel

  return dispIncr_;
}


//-----------------------------------------------------------------------
//   initConstraints
//-----------------------------------------------------------------------

void DispArclenBCs::initConstraints

  ( const Properties&  globdat )

{
  Assignable<NodeGroup> group;
  IntVector             inodes;

  // Get nodes, then dofs of nodes, and constraints of dofs

  nodes_ = NodeSet::find    ( globdat );
  dofs_  = XDofSpace::get   ( nodes_.getData(), globdat );
  cons_  = Constraints::get ( dofs_, globdat );

  // loop over node groups

  for ( int ig = 0; ig < ngroups_; ig++ )
  {
    group  = NodeGroup::find ( nodeGroups_[ig], nodes_, globdat );

    int nn = group.size ();
    inodes . resize ( nn );
    inodes = group.getIndices ();

    int itype  = dofs_->findType ( dofTypes_[ig] );

    if ( ig != loaded_ )  
    {
      // apply u=0 boundary conditions

      for ( int in = 0; in < nn; in++ )
      {
        int idof = dofs_->getDofIndex ( inodes[in], itype );

        cons_->addConstraint ( idof, 0. );
      }
    }
    else
    {
      // get master dof index

      master_ = dofs_->getDofIndex ( inodes[0], itype );
    }
  }

  // compress for more efficient storage

  cons_->compress();
}

//-----------------------------------------------------------------------
//   toArclControl
//-----------------------------------------------------------------------

// remove constraints already defined
// if a node group has only one node, then apply a force there
// otherwise, tie nodes together and apply a force on one 
// master node.

void DispArclenBCs::toArclControl

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::util::Printer;
  using jive::model::ActionParams;

  // System::out() << "Resetting constraints to apply force ... \n";
  
  // loop over node groups

  Array<idx_t>  inodes, idofs;
  Assignable<NodeGroup> group;

  Array<idx_t> itypes(1);

  int ig = loaded_;

  group  = NodeGroup::find ( nodeGroups_[ig], nodes_, globdat );

  int nn = group.size ();
  inodes . resize ( nn );
  idofs  . resize ( nn );
  inodes = group.getIndices ();

  int itype  = itypes[0] = dofs_->findType ( dofTypes_[ig] );

  // get dofs associated with inodes
  // then erase their constraints

  dofs_->getDofIndices    ( idofs, inodes, itypes );
  cons_->eraseConstraints ( idofs );

  master_ = dofs_->getDofIndex ( inodes[0], itype );

  for ( int in = 1; in < nn; in++ )
  {
    int idof = dofs_->getDofIndex ( inodes[in], itype );

    cons_->addConstraint ( idof, master_, 1.0 );
  }

  // when accepted and converged, this function is called after commit
  // so loadScale0_ has been updated already, 
  // when switching due to non-convergence, the old value should be used

  params.set( "OldLoadScale", loadScale0_ );

  //cons_->printTo ( Printer::get() ); System::out() << "\n\n";
  //Printer::flush ();
}

//-----------------------------------------------------------------------
//   applyConstraints
//-----------------------------------------------------------------------

void DispArclenBCs::applyConstraints

  ( const Properties&  params,
    const Properties&  globdat )

{
  // Update nonzero displacement on loaded boundary

  // System::out() <<" Imposing displacement over nodes ... \n";

  IntVector             inodes;
  Assignable<NodeGroup> group;

  // the displacement increment is adjusted outside
  // precisely in the DispArclenModel

  int        itype, idof;

  int ig = loaded_;
  group  = NodeGroup::find ( nodeGroups_[ig], nodes_, globdat );

  int nn = group.size ();
  inodes . resize ( nn );
  inodes = group.getIndices ();

  itype  = dofs_->findType ( dofTypes_[ig] );

  // get amount of displacement

  for ( int in = 0; in < nn; in++ )
  {
    idof = dofs_->getDofIndex ( inodes[in], itype );

    cons_->addConstraint ( idof, dispVal_ );
  }

  // compress for more efficient storage

  cons_->compress();

  //using jive::util::Printer;
  //cons_->printTo ( Printer::get() ); Printer::flush (); System::out() << "\n\n";
}

//-----------------------------------------------------------------------
//   storeLoadScale
//-----------------------------------------------------------------------


void DispArclenBCs::storeLoadScale

  ( const Properties&  globdat,
    const Vector&      fint )

{
  // compute the total load on the loaded boundary
  
  Assignable<NodeGroup> group;
  IntVector             inodes;

  int ig = loaded_;
  group  = NodeGroup::find ( nodeGroups_[ig], nodes_, globdat );

  int nn = group.size ();
  inodes . resize ( nn );
  inodes = group.getIndices ();

  int itype = dofs_->findType ( dofTypes_[ig] );

  loadScale_ = 0.;

  for ( int in = 0; in < nn; in++ )
  {
    int idof = dofs_->getDofIndex ( inodes[in], itype );

    loadScale_ += fint[idof];
  }
}

//-----------------------------------------------------------------------
//   getUnitLoad
//-----------------------------------------------------------------------

void DispArclenBCs::getUnitLoad

  ( const Properties&  params,
    const Properties&  globdat )

{
  // get unit force vector

  Vector       fext;

  params.get ( fext, ActionParams::EXT_VECTOR );

  fext[master_] = 1.;
}

//-----------------------------------------------------------------------
//   commit
//-----------------------------------------------------------------------

// store the converged loadscale in globdat during disp control
// and the converged displacement during arclength control

void DispArclenBCs::commit

  ( const bool         isDispControl,
    const Properties&  globdat )

{
  loadScale0_ = loadScale_;

  if ( isDispControl )
  {
    dispVal0_ = dispVal_;

    // increment the displacement value for next step

    dispVal_  = dispVal0_ + dispIncr_;
  }
  else
  {
    // get master displacement from StateVector

    Vector  disp;

    StateVector::get ( disp,  dofs_, globdat );

    dispVal0_ = jem::numeric::abs ( disp[master_] );
  } 
}
