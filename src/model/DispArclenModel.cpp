/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implements a model that combines
 *  displacement and arclen control in one model. This model
 *  must be used with FlexArclenModule to have a flexible path following solver.
 *
 *  Basic ideas:
 *
 *    1. Disp control: The simulation starts with displacement control
 *       where some nodes are constrained. During
 *       this stage, in FlexArclenModule, the NonlinModule is being
 *       used. This is the case until divergence occurs (snapback).
 *       Then, one node is defined as master and others as its slave.
 *       A force is then applied on this master node so that the
 *       unit external force required by ArclenModule can be
 *       constructed. 
 *
 *   2.  Arclen control (based on energy released). Now, in FlexArclenModule,
 *       the ArclenModule is active so that lambda is now an unknown.
 *       This is the case until hardening branch is detected, then
 *       switch back to (1), Load control
 *
 *   Constrained dofs are managed by a child model which is an instance
 *   of class DispControlModel.
 *
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 31 January 2009
 *
 */


#include <jem/base/array/operators.h>
#include <jem/base/System.h>
#include <jem/base/Float.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/util/Event.h>
#include <jem/util/StringUtils.h>
#include <jive/util/error.h>
#include <jive/util/Printer.h>
#include <jive/util/utilities.h>
#include <jive/util/ItemSet.h>
#include <jive/util/DofSpace.h>
#include <jive/util/Constraints.h>
#include <jive/util/VectorManager.h>
#include <jive/util/Globdat.h>
#include <jive/algebra/VectorSpace.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/model/ModelFactory.h>
#include <jive/implict/ArclenActions.h>
#include <jive/implict/SolverInfo.h>


#include "ThermNames.h"
#include "DispArclenModel.h"
#include "module/SolverNames.h"


extern "C"
{
  #include <math.h>
}

using jem::io::endl;
using jem::numeric::axpy;
using jem::util::StringUtils;
using jive::IntVector;
using jive::util::VectorManager;
using jive::model::Actions;
using jive::model::ActionParams;
using jive::model::StateVector;
using jive::implict::ArclenActions;
using jive::implict::ArclenParams;
using jive::implict::SolverInfo;


//=======================================================================
//   class DispArclenModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  DispArclenModel::TYPE_NAME       = "DispArclen";

const char*  DispArclenModel::CONSTRAINT_PROP = "constraints";
const char*  DispArclenModel::OPT_ITER_PROP   = "optIter";
const char*  DispArclenModel::SWT_ITER_PROP   = "swtIter";
const char*  DispArclenModel::SWT_ENER_PROP   = "swtEnergy";
const char*  DispArclenModel::MIN_INCR_PROP   = "minEnergyIncr";
const char*  DispArclenModel::MAX_INCR_PROP   = "maxEnergyIncr";
const char*  DispArclenModel::DISP_INCR_PROP  = "dispIncr";
const char*  DispArclenModel::MIN_DISP_PROP   = "minDispIncr";
const char*  DispArclenModel::REDUCTION_PROP  = "reduction";
const char*  DispArclenModel::PREFER_DC_PROP  = "preferDisp";

const char*  DispArclenModel::DELTA_STATE_    = "deltaState0";

const int    DispArclenModel::U_LOAD_         = 1 << 0;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


DispArclenModel::DispArclenModel

  ( const String&      name,
    const Ref<Model>&  child ) :

    Super  ( name  ),
    child_ ( newInstance<DispArclenBCs>() ),
    out_   ( System::info( name ) )

{
  updated_   = 0;
  optIter_   = 4;
  minIncr_   = 1.0e-3;
  maxIncr_   = 1.0e+1;
  dispIncr0_ = 1.0;
  arcLength_ = 0.0;
  lastArcl_  = 0.0;
  reduction_ = 0.2;
  swtEner_   = Float::MAX_VALUE;
  gC_        = Float::MAX_VALUE;
  totalDiss_ = 0.0;

  isDispControl_ = true; 
  onceDown_      = false;
  preferDisp_    = false;

  nformat_.setScientific     (   );
  nformat_.setFractionDigits ( 4 );
}


DispArclenModel::~DispArclenModel ()
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool DispArclenModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  // initialization

  if ( action == Actions::INIT )
  {
    child_->initConstraints( globdat );

    init_ ( globdat );

    return true;
  }

  // compute the external force vector
  //
  // [this is not really needed, a unit force is placed on a 
  // constrained node, no effect, except that now in the nonlinSolver
  // the residual scale factor has a minimum value 1.0, a trick to 
  // prevent tiny residual scale factor in nonlinSover when the structure
  // is completely broken]

  if ( action == Actions::GET_EXT_VECTOR  )
  {
    child_->getUnitLoad ( params, globdat );

    return false;
  }

  // apply displacement increment

  if ( action == Actions::GET_CONSTRAINTS )
  {
    if ( isDispControl_ )
    {
      child_->applyConstraints ( params, globdat );
    }

    return true;
  }

  // compute arclength function

  if ( action == ArclenActions::GET_ARC_FUNC )
  {
    evalArcFunc_ ( params, globdat );

    return true;
  }

  if ( action == ArclenActions::GET_UNIT_LOAD )
  {
    getUnitLoad_ ( params, globdat );

    return true;
  }

  if ( action == SolverNames::CHECK_COMMIT )
  {
    checkCommit_ ( params, globdat );

    checkSwitch_ ( params, globdat );

    return true;
  }

  if ( action == Actions::COMMIT )
  {
    commit_ ( params, globdat );

    child_->commit( isDispControl_, globdat );

    return true;
  }

  // this REDUCE_STEP action is thrown by FlexArclenModule

  if ( action == SolverNames::REDUCE_STEP )
  {
    reduceStep_ ( params, globdat );

    return true;
  }

  // this action is called in FlexArclenModule,
  // when switch from displacement control to arclen control

  if ( action == SolverNames::TO_ARCL )
  {
    toArclControl_ ( params, globdat );

    child_->toArclControl ( params, globdat );

    return true;
  }

  if ( action == Actions::GET_MATRIX0 ||
       action == Actions::GET_INT_VECTOR )
  {
    Vector       fint;
    params.get ( fint, ActionParams::INT_VECTOR );

    child_->storeLoadScale ( globdat, fint );

    return true;
  }

  if ( action == SolverNames::STOP_ARCL )
  {
    params.set ( SolverNames::DONE, true );

    return true;
  }

  if ( action == SolverNames::TO_DISP )
  {
    toDispControl_ ( params, globdat );

    child_->toDispControl ( params, globdat );

    return true;
  }

  if ( action == ThermNames::DISINCREMENT )
  {
    child_->disincrement ( );
  }

  if ( action == "GET_DISSIPATION" )
  {
    StringVector names ( StringUtils::split( myName_, '.' ) );
    params.set ( names[names.size()-1], totalDiss_ );
  }

  return false;
}


//-----------------------------------------------------------------------
//   findModel
//-----------------------------------------------------------------------


Model* DispArclenModel::findModel ( const String& name ) const
{
  return const_cast<Self*> ( this );
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void DispArclenModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps = props.findProps ( myName_ );

  double maxD = Float::MAX_VALUE;

  myProps.find ( optIter_, OPT_ITER_PROP, 1,        1000 );
  myProps.find ( swtIter_, SWT_ITER_PROP, 1,        1000 );
  myProps.find ( swtEner_, SWT_ENER_PROP, 0.0,      maxD );
  myProps.find ( minIncr_, MIN_INCR_PROP, 0.0,      maxD );
  myProps.find ( maxIncr_, MAX_INCR_PROP, minIncr_, maxD );

  myProps.find ( reduction_, REDUCTION_PROP,   0.0, 1.0 );
  myProps.find ( gC_,        "fractureEnergy" );

  // read the initial displacement increment

  myProps.get  ( dispIncr0_, DISP_INCR_PROP );
  child_->setDispIncr ( dispIncr0_ );

  minDispIncr_ = dispIncr0_;
  myProps.find ( minDispIncr_, MIN_DISP_PROP, dispIncr0_, maxD );

  myProps.find ( preferDisp_, PREFER_DC_PROP );

  // configuring the child object of class DispArclenBCs

  Properties childProps = myProps.findProps ( CONSTRAINT_PROP );
  child_->configure ( childProps, globdat );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void DispArclenModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );

  myConf.set ( OPT_ITER_PROP,    optIter_      );
  myConf.set ( SWT_ITER_PROP,    swtIter_      );
  myConf.set ( MIN_INCR_PROP,    minIncr_      );
  myConf.set ( MAX_INCR_PROP,    maxIncr_      );
  myConf.set ( DISP_INCR_PROP,   dispIncr0_    );
  myConf.set ( MIN_DISP_PROP,    minDispIncr_  );
  myConf.set ( PREFER_DC_PROP,   preferDisp_   );

  Properties childConf = myConf.makeProps ( CONSTRAINT_PROP );
  child_->getConfig ( childConf, globdat );
}


//-----------------------------------------------------------------------
//   setMaxIter
//-----------------------------------------------------------------------


void DispArclenModel::setMaxIter ( int count )
{
  JEM_PRECHECK ( count > 0 );

  optIter_ = count;
}


//-----------------------------------------------------------------------
//   setIncrRange
//-----------------------------------------------------------------------


void DispArclenModel::setIncrRange

  ( double  minIncr,
    double  maxIncr )

{
  JEM_PRECHECK ( minIncr >= 0.0 && minIncr < maxIncr );

  minIncr_ = minIncr;
  maxIncr_ = maxIncr;
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Model> DispArclenModel::makeNew

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  return newInstance<Self> ( name );
}


//-----------------------------------------------------------------------
//   init_
//-----------------------------------------------------------------------


void DispArclenModel::init_ ( const Properties& globdat )
{
  dofs_    = DofSpace    :: get ( globdat, getContext() );
  updated_ = 0;

  connect_ ();
}


//-----------------------------------------------------------------------
//   initLoad_
//-----------------------------------------------------------------------


void DispArclenModel::initLoad_ ( const Properties& globdat )
{
  // In the current implementation, this is absolutely ridiculous:
  // a vector with size dofCount and all 0.0's except one entry with 1.0
  // and then in getReleasedEnergy, this one entry is selected 
  // and used for multiplication
  // 
  // However, if you want to generalize to more complex load vectors,
  // this is the way to go

  Properties  params;

  load_.resize ( dofs_->dofCount() );

  load_ = 0.0;

  params.set ( ActionParams::EXT_VECTOR, load_ );

  child_->getUnitLoad ( params, globdat );

  dofs_ ->resetEvents ();

  updated_ |= U_LOAD_;
}


//-----------------------------------------------------------------------
//   evalArcFunc_
//-----------------------------------------------------------------------


void DispArclenModel::evalArcFunc_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jem::max;
  using jem::isTiny;

  if ( vspace_ == NIL )
  {
    vspace_ = VectorSpace::get ( dofs_, globdat );
  }

  if ( ! (updated_ & U_LOAD_) )
  {
    initLoad_    ( globdat );
  }

  double dEnergy = getReleasedEnergy_ ( params, globdat );
  double fvalue  = dEnergy - arcLength_ ;

  params.set ( ArclenParams::ARC_FUNC,   fvalue );
}


//-----------------------------------------------------------------------
//   getUnitLoad_
//-----------------------------------------------------------------------


void DispArclenModel::getUnitLoad_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::util::sizeError;

  Vector  f;

  if ( ! (updated_ & U_LOAD_) )
  {
    initLoad_ ( globdat );
  }

  params.get ( f, ArclenParams::UNIT_LOAD );

  if ( f.size() != load_.size() )
  {
    sizeError ( getContext(),
                "unit load vector", f.size(), load_.size() );
  }

  f = load_;
}


//-----------------------------------------------------------------------
//   commit_
//-----------------------------------------------------------------------


void DispArclenModel::commit_

  ( const Properties&  params,
    const Properties&  globdat )

{
  Properties info = SolverInfo::get ( globdat );

  int   iiter;
  info.get ( iiter, SolverInfo::ITER_COUNT );

  // store total dissipated energy (for postprocessing)

  totalDiss_ += lastArcl_;

  if ( isDispControl_ )
  {
    // adjust increment according to number of iterations
    // do not do for the initial elastic branch since no help

    double oldIncr = child_->getDispIncr();
    double newIncr = 0.;

    if ( onceDown_ )
    {
      double n       = ( iiter - optIter_ ) / 4.0;
             newIncr = oldIncr * ::pow ( 0.5, n );
    }
    else
    {
      newIncr = oldIncr;
    }

    // check if it falls within allowable interval (defined by user)

    //newIncr = max ( newIncr, minDispIncr_ );
    //newIncr = min ( newIncr, maxDispIncr_ );

    child_->setDispIncr ( newIncr );

    out_ << "DispControl, old disp increment " << oldIncr << "\n"
         << "             new disp increment " << newIncr << "\n\n";

    out_.flush();
  }
  else
  {
    // adjust increment according to number of iterations

    params.get ( iiter, ArclenParams::IITER );

    double n    = ( iiter - optIter_ ) / 4.0;
    arcLength_ *= ::pow ( 0.5, n );

    // check if it falls within allowable interval (defined by user)

    if        ( arcLength_ < minIncr_ )
    {
      arcLength_ = minIncr_;
    }
    else if   ( arcLength_ > maxIncr_ )
    {
      arcLength_ = maxIncr_;
    }

/*
    // check for switch if preferDisp_ option is set

    if ( preferDisp_ && arcLength_ < swtEner_ )
    {
      // if positive load increment: switch back to load control
      // add abs so that it works for negative force

      double lambda0 = abs ( child_->getLoadScale0() );
      double dlambda = abs ( child_->getLoadScale () ) - lambda0;

      onceDown_ |= dlambda < 0.;

      if ( ( dlambda > 0. ) & ( onceDown_ ) )
      {
        out_ << "Loads increases, go back to displacement control...\n";

        params.set ( SolverNames::DO_SWITCH, "please" );

        onceDown_ = false;
      }
    }
*/
    out_ << "EnergyArclen, old dEnergy " << lastArcl_ << "\n"
         << "              new dEnergy " << arcLength_ << ", energy0: " << energy0_ << "\n\n";

    out_.flush();
  }
}

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void DispArclenModel::checkCommit_

  ( const Properties&  params,
    const Properties&  globdat )

{
  if ( ! isDispControl_ )
  {
    // don't accept solution when the load changed sign

    if ( child_->getLoadScale() * child_->getLoadScale0() < 0. )
    {
      System::warn() << "Load changed sign." << endl;

      reduceStep_ ( params, globdat );

      params.set ( SolverNames::ACCEPT, false );
    }
  }
}

//-----------------------------------------------------------------------
//   checkSwitch_
//-----------------------------------------------------------------------

void DispArclenModel::checkSwitch_

  ( const Properties&  params,
    const Properties&  globdat )

{
  int    iiter;

  if ( isDispControl_ )
  {
    lastArcl_ = getReleasedEnergy_( params, globdat );

    params.get ( iiter, "iterCount" );

    out_ << "DispControl: dEnergy " << lastArcl_ << 
            ", loadScale " << child_->getLoadScale() <<  
            ", energy0:"   << energy0_ << endl;

    // store value anyway

    if ( ( lastArcl_ > swtEner_ ) || 
         ( iiter >= swtIter_ && lastArcl_ > minIncr_ ) ||
         ( energy0_ > gC_ )
       ) 
    {
      params.set ( SolverNames::DO_SWITCH, "please" );
    }

    // check if it's going down 

    double lambda0 = ::fabs ( child_->getLoadScale0() );
    double dlambda = ::fabs ( child_->getLoadScale () ) - lambda0;

    onceDown_ |= dlambda < 0.;

    if (onceDown_) out_ << "Going down !!! \n";
  }
  else
  {
    lastArcl_ = arcLength_;

    if ( preferDisp_ && arcLength_ < swtEner_ )
    {
      // if positive load increment: switch back to load control
      // add abs so that it works for negative force

      double lambda0 = ::fabs ( child_->getLoadScale0() );
      double dlambda = ::fabs ( child_->getLoadScale () ) - lambda0;

      onceDown_ |= dlambda < 0.;

      if ( ( dlambda > 0. ) & ( onceDown_ ) )
      {
        out_ << "Loads increases, go back to displacement control...\n";

        params.set ( SolverNames::DO_SWITCH, "please" );

        onceDown_ = false;
      }
    }
  }
}

//-----------------------------------------------------------------------
//   reduceStep_
//-----------------------------------------------------------------------

// For now, only adjust the load increment for Arclen control
// by either (1) reduce the amount of released energy
// or (2) switch back to load control since the load is going up again

void DispArclenModel::reduceStep_

  ( const Properties&  params,
    const Properties&  globdat )

{
  if ( isDispControl_ )
  {
    double newIncr = child_->reduceDispIncr ( reduction_ );

    if ( newIncr > minDispIncr_ )
    {
      out_ << "Displacement increment reduced to "
           << newIncr << "\n";

      params.set ( SolverNames::DONE, false );
    }
  }
  else
  {
    if ( arcLength_ > minIncr_ )
    {
      arcLength_ *= reduction_; // reduce arc-length

      out_ << "Path following parameter reduced to "
           << arcLength_ << "\n";

      params.set ( SolverNames::DONE, true );
    }
  }
}

//-----------------------------------------------------------------------
//   toArclControl_
//-----------------------------------------------------------------------

void DispArclenModel::toArclControl_

  ( const Properties&  params,
    const Properties&  globdat )

{
  isDispControl_ = false;
  arcLength_     = lastArcl_;

  if        ( arcLength_ < minIncr_ )
  {
    arcLength_ = minIncr_;
  }
  else if   ( arcLength_ > maxIncr_ )
  {
    arcLength_ = maxIncr_;
  }

  out_ << "Switching to arclength control with dEnergy: " 
       << arcLength_ << endl;
}

//-----------------------------------------------------------------------
//   toDispControl_
//-----------------------------------------------------------------------

void DispArclenModel::toDispControl_

  ( const Properties&  params,
    const Properties&  globdat )

{
  isDispControl_ = true;

  bool converged;

  params.get ( converged, SolverNames::CONVERGED );

  if ( converged )
  {
    out_ << "Switching temporarily to disp control "
         << "fixing the displacement at " << child_->getDispValue() << "\n";
  }
  else
  {
    out_ << "Switching back to disp control "
         << "with disp incr: " << dispIncr0_ << "\n";
  }
}


//-----------------------------------------------------------------------
//   connect_
//-----------------------------------------------------------------------


void DispArclenModel::connect_ ()
{
  using jem::util::connect;

  connect ( dofs_->newSizeEvent,   this, & Self::dofsChanged_ );
  connect ( dofs_->newOrderEvent,  this, & Self::dofsChanged_ );

  dofs_->resetEvents ();
}


//-----------------------------------------------------------------------
//   dofsChanged_
//-----------------------------------------------------------------------


void DispArclenModel::dofsChanged_ ()
{
  updated_ = 0;
}


//-----------------------------------------------------------------------
//   consChanged_
//-----------------------------------------------------------------------


void DispArclenModel::consChanged_ ()
{
}

//-----------------------------------------------------------------------
//   getReleasedEnergy_
//-----------------------------------------------------------------------

// compute the energy released during one load step
// used to check for switch to use ArclenModule 

double DispArclenModel::getReleasedEnergy_

 ( const Properties&  params,
   const Properties&  globdat )

{
  if ( ! (updated_ & U_LOAD_) )
  {
    initLoad_ ( globdat );
  }

  // get master dof index with nonzero displacement
  // [currently DispControlModel deals only with one, 
  //  but this function is more general]

  IntVector idofs;

  int nLoaded = child_->getLoadedDofs ( idofs );

  double jac11, lambda, dlambda, lambda0, dEnergy;
  Vector jac10, u, u0, du, tmp, fth0;

  StateVector::get    ( u,  dofs_, globdat );
  StateVector::getOld ( u0, dofs_, globdat );

  Vector loadB ( select ( load_, idofs ) );
  Vector uB    ( select (     u, idofs ) );
  Vector u0B   ( select (    u0, idofs ) );

  du. resize ( nLoaded );
  tmp.resize ( nLoaded );

  if ( params.find ( lambda,  ArclenParams::LOAD_SCALE ) )
  {
    // when called in evalArcFunc 

    params.get ( lambda0, ArclenParams::OLD_LOAD_SCALE );
    params.get ( jac10,   ArclenParams::JACOBIAN10 );

    dlambda = lambda - lambda0;
    jac11   = -0.5 * dot ( u0B, loadB );
    jac10   = 0.;

    select( jac10, idofs ) = 0.5 * lambda0 * loadB;

    params.set ( ArclenParams::JACOBIAN11, jac11 );
  }
  else
  {
    // when called in checkSwitch

    lambda0 = child_->getLoadScale0();
    dlambda = child_->getLoadScale () - lambda0;
  }

  axpy ( du, uB, -1.0, u0B );

  energy0_  = 0.5 * dot ( u0B, loadB );
  energy0_ *= lambda0;

  tmp       = lambda0 * du - dlambda * u0B;
  dEnergy   = 0.5 * dot ( tmp, loadB );

  // double power = lambda0 * dot ( du, loadB );

  if ( globdat.find ( fth0, "var.thermForce" ) )
  {
    // correct energy increment for thermal strain contribution
    // to elastic energy

    int ndofs = dofs_->dofCount();
    int size0 = fth0.size();

    if ( ndofs > size0 )
    {
      // this is not done anymore?

      fth0.reshape ( ndofs );

      for ( int i = size0; i < ndofs; ++i )
      {
        fth0[i] = 0.;
      }
    }

    VectorManager::getVector ( du, DELTA_STATE_, dofs_, globdat );

    axpy ( du, u, -1.0, u0 );

    JEM_ASSERT ( du.size() == fth0.size() );

    double duFth = dot ( du, fth0 );

    dEnergy += 0.5 * duFth;

    if ( jac10.size() > 0 )
    {
      jac10 += 0.5 * fth0;
    }

    // NB: elastic energy (energy0_) is not adapted,
    // because this would require another vector assembly
  }

  return dEnergy;
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareDispArclenModel
//-----------------------------------------------------------------------


void declareDispArclenModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( DispArclenModel::TYPE_NAME,
                          & DispArclenModel::makeNew );
}
