/*
 *  Copyright (C) 2012, Ton Duc Thang University.
 *
 * This class implements a model to define imposed external forces
 * on certain nodes in a computation with a combined load control
 * and arclength control. It is used as a child of class 
 * LoadArclenModel and is able to switch between the two different strategies.
 *
 * Author(s):
 *
 *   Vinh Phu Nguyen, nvinhphu@gmail.com
 *   Saigon, Vietnam October 2012
 *
 */

#include <jem/util/ArrayBuffer.h>
#include <jem/base/Array.h>
#include <jem/numeric/utilities.h>
#include <jem/base/IllegalInputException.h>
#include <jive/util/Globdat.h>
#include <jive/util/Printer.h>
#include <jive/util/Constraints.h>
#include <jive/util/XDofSpace.h>
#include <jive/model/StateVector.h>


#include "LoadArclenBCs.h"
#include "module/SolverNames.h"


using jem::util::ArrayBuffer;
using jive::model::StateVector;
using jive::model::ActionParams;

//-----------------------------------------------------------------------
//   static constants
//-----------------------------------------------------------------------

const char* LoadArclenBCs::TYPE_NAME     = "LoadControl";
const char* LoadArclenBCs::NODES_PROP    = "nodes";
const char* LoadArclenBCs::DOF_PROP      = "dofs";
const char* LoadArclenBCs::LOAD_PROP     = "refForce";
const char* LoadArclenBCs::INCR_PROP     = "loadIncr";

//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------


LoadArclenBCs::LoadArclenBCs ( )
{
  // initialize variables

  loadScale_     = 0.;
  loadScale0_    = 0.;
  loadIncr_      = 0.;
}

LoadArclenBCs::~LoadArclenBCs()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void LoadArclenBCs::configure

  ( const Properties&  myProps,
    const Properties&  globdat )

{
  /*
   loadObj =
   {
     nodes     = [317,283];
     refForce  = [1.,-2.1436];
     dofs      = ["uy","uy"];
     loadIncr  = 1.5; 
   };
   */   

  myProps.get( nodes_,    NODES_PROP );
  myProps.get( refLoad_,  LOAD_PROP  );
  myProps.get( dofTypes_, DOF_PROP   );
  myProps.get( loadIncr_, INCR_PROP  );

  nnodes_ = nodes_.size ( );

  index_.resize      ( nnodes_ );
  initializeRefLoad_ ( globdat );

  if ( dofTypes_.size() != nnodes_ )
  {
    throw IllegalInputException ( JEM_FUNC,
          "dofTypes must have the same length as nodeGroups" );
  }
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void LoadArclenBCs::getConfig 

  ( const Properties& myConf,
    const Properties& globdat ) const

{
  myConf.set ( NODES_PROP,  nodes_    );
  myConf.set ( LOAD_PROP,   refLoad_  );
  myConf.set ( DOF_PROP,    dofTypes_ );
  myConf.set ( INCR_PROP,   loadIncr_ );
}

//-----------------------------------------------------------------------
//   toLoadControl
//-----------------------------------------------------------------------

void LoadArclenBCs::toLoadControl

  ( const Properties&  params,
    const Properties&  globdat )

{
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
  }
  else
  {
    // apply increment to old displacement

  }

  //using jive::util::Printer;
  //cons_->printTo ( Printer::get() );
  //Printer::flush ();  System::out() << "\n\n";
}

//-----------------------------------------------------------------------
//   initializeRefLoad_
//-----------------------------------------------------------------------

void LoadArclenBCs::initializeRefLoad_

  ( const Properties&  globdat )

{
  IntVector             inodes;

  // Get nodes, then dofs of nodes, and constraints of dofs

  nodeSet_ = NodeSet::find    ( globdat );
  dofs_    = XDofSpace::get   ( nodeSet_.getData(), globdat );

  // loop over node groups

  for ( int i = 0; i < nnodes_; i++ )
  {
    int nn     = nodes_[i];
    int itype  = dofs_->findType    ( dofTypes_[i] );
    int idof   = dofs_->getDofIndex ( nn, itype );
    index_[i]  = idof;
  }
}

//-----------------------------------------------------------------------
//   toArclControl
//-----------------------------------------------------------------------

void LoadArclenBCs::toArclControl

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::util::Printer;
  using jive::model::ActionParams;

  params.set( "OldLoadScale", loadScale0_ );
}

//-----------------------------------------------------------------------
//   computeExternalForce
//-----------------------------------------------------------------------

void LoadArclenBCs::computeExternalForce

  ( const Properties&  params,
    const Properties&  globdat )

{
  // get unit force vector

  Vector       fext;
  //Vector       temp;

  params.get ( fext, ActionParams::EXT_VECTOR );

  //fext[index_] = refLoad_ ;
  //temp  = loadScale_ * refLoad_;

  select ( fext, index_ ) = loadScale_ * refLoad_;
}

//-----------------------------------------------------------------------
//   getUnitLoad
//-----------------------------------------------------------------------

void LoadArclenBCs::getUnitLoad

  ( const Properties&  params,
    const Properties&  globdat )

{
  // get unit force vector

  Vector       fext;

  params.get ( fext, ActionParams::EXT_VECTOR );

  select ( fext, index_ ) = refLoad_;
}

//-----------------------------------------------------------------------
//   commit
//-----------------------------------------------------------------------


void LoadArclenBCs::commit

  ( const bool         isLoadControl,
    const Properties&  globdat )

{
  //globdat.get ( ActionParams::LOAD_SCALE, loadScale_ );


  if ( isLoadControl )
  {
    loadScale0_  = loadScale_;
    loadScale_  += loadIncr_;
  }
  else
  {
  } 
}
