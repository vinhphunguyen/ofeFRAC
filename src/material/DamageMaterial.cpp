/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the isotropic elastic-damage material
 *  This represents the material at a point in space.
 *  It is implemented in such a way that can be used for any
 *  discretisation methods, say finite elements, EFG and so on.
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *
 */

#include <jem/base/System.h>
#include <jem/base/limits.h>
#include <jem/base/Float.h>
#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/base/Array.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/utilities.h>
#include <jem/io/NumberFormat.h>
#include <algorithm>

#include "util/utilities.h"
#include "DamageMaterial.h"


using namespace jem;

using namespace jem::io;


using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jem::util::Properties;
using jive::Vector;


const double EPS = 1e-10;




//=======================================================================
//   class DamageMaterial
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  DamageMaterial::SOFTENING_PROP       = "softening";
const char*  DamageMaterial::EQUISTRAIN_PROP      = "equistrain";
const char*  DamageMaterial::KAPPAI_PROP          = "kappaI";
const char*  DamageMaterial::KAPPAC_PROP          = "kappaC";
const char*  DamageMaterial::TENSILE_PROP         = "ft";
const char*  DamageMaterial::FRACTURE_ENERGY_PROP = "gf";
const char*  DamageMaterial::ALPHA_PROP           = "alpha";
const char*  DamageMaterial::BETA_PROP            = "beta";
const char*  DamageMaterial::ETA_PROP             = "eta";
const char*  DamageMaterial::B_PROP               = "b";
const char*  DamageMaterial::LENGTH_PROP          = "lengthscale";
const char*  DamageMaterial::REMOVE_DAMAGE_PROP   = "damThreshold";

const char*  DamageMaterial::CRACK_WIDTH_PROP     = "crackWidth";

const char*  DamageMaterial::MAZARS_EQUI_STRAIN   = "Mazars";
const char*  DamageMaterial::MISES_EQUI_STRAIN    = "vonMises";
const char*  DamageMaterial::RANKINE_EQUI_STRAIN  = "Rankine";
const char*  DamageMaterial::ENERGY_STRAIN        = "Energy";

const char*  DamageMaterial::PERFECT_SOFTENING    = "perfect";
const char*  DamageMaterial::LINEAR_SOFTENING     = "linear";
const char*  DamageMaterial::EXPONENT1_SOFTENING  = "exponential1"; // ki,alpha, beta
const char*  DamageMaterial::EXPONENT2_SOFTENING  = "exponential2"; // ki, kc
const char*  DamageMaterial::EXPONENT3_SOFTENING  = "exponential2Reg"; // Jirasek regularised of 2
const char*  DamageMaterial::POWER_SOFTENING      = "power";
const char*  DamageMaterial::HYPERBOLIC_SOFTENING = "hyperbolic";
const char*  DamageMaterial::EXPOENERGY_SOFTENING = "expoEnergy";

const double DamageMaterial::CRITICAL_DAMAGE      = 0.999;

vector<String> DamageMaterial::equiStrainDefs   (4);
vector<String> DamageMaterial::softeningLawDefs (8);


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


DamageMaterial::DamageMaterial

  ( int rank, const Properties& globdat )
    : Material ( rank, globdat )
{
  kappaI_     = 0.0;
  kappaC_     = 0.0;
  alpha_      = 0.0;
  beta_       = 0.0;
  eta_        = 0.0;
  b_          = 0.0;
  c_          = 0.0;
  threshold_  = 0.1;

  softening_  = Linear;
  equiStrain_ = VonMises ;

  elasticMat_ = newInstance<HookeMaterial> ( rank, globdat );

  elasticMod_ .resize ( STRAIN_COUNTS[rank], STRAIN_COUNTS[rank] );

  equiStrainDefs[0] = MAZARS_EQUI_STRAIN  ;
  equiStrainDefs[1] = MISES_EQUI_STRAIN   ;
  equiStrainDefs[2] = RANKINE_EQUI_STRAIN ;
  equiStrainDefs[3] = ENERGY_STRAIN ;

  softeningLawDefs[0] = PERFECT_SOFTENING ;
  softeningLawDefs[1] = LINEAR_SOFTENING ;
  softeningLawDefs[2] = EXPONENT1_SOFTENING ;
  softeningLawDefs[3] = EXPONENT2_SOFTENING ;
  softeningLawDefs[4] = POWER_SOFTENING ;
  softeningLawDefs[5] = HYPERBOLIC_SOFTENING ;
  softeningLawDefs[6] = EXPONENT3_SOFTENING ;
  softeningLawDefs[7] = EXPOENERGY_SOFTENING ; // German, IJNME 07
}


DamageMaterial::~DamageMaterial ()
{} 

//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------

void DamageMaterial::allocPoints ( int count )
{
  preHist_.eqveps. resize ( count );
  preHist_.loading.resize ( count );

  newHist_.eqveps. resize ( count );
  newHist_.loading.resize ( count );

  preHist_.eqveps  = 0.0;
  preHist_.loading = 0;

  newHist_.eqveps  = 0.0;
  newHist_.loading = 0;
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void DamageMaterial::configure ( const Properties& props,
                                 const Properties& globdat )
{
  using jem::maxOf;
  using std::find;

  elasticMat_->configure ( props, globdat );

  bool ki = props.find ( kappaI_, KAPPAI_PROP, 0.0, maxOf( kappaI_) );

  if ( !ki )
  {
    props.get ( ft_, TENSILE_PROP, 0.0, maxOf ( ft_ ) );

    kappaI_ = ft_ / elasticMat_->giveYoung ();
  }
  else
  {
    ft_ = kappaI_ * elasticMat_->giveYoung ();
  }

  props.find ( gf_,     FRACTURE_ENERGY_PROP, 0.0, maxOf ( gf_) );
  props.find ( lambda_, CRACK_WIDTH_PROP,     0.0, 1000.0       );
  props.find ( c_,      LENGTH_PROP,          0.0, maxOf ( c_)  );

  // read type of softening law used

  String      softening;

  props.get ( softening, SOFTENING_PROP  );

  if ( find ( softeningLawDefs.begin (),
              softeningLawDefs.end   (),
              softening ) == softeningLawDefs.end () )
  {
    throw Error (
      JEM_FUNC,
      String("unexpected definition of softening law!!!\n") +
      String("Supported softening laws include: \n") +
      PERFECT_SOFTENING   + String(", ") + LINEAR_SOFTENING + String(", ") +
      EXPONENT1_SOFTENING + String(", ") + EXPONENT2_SOFTENING + String(", ") +
      POWER_SOFTENING     + String(", ") + HYPERBOLIC_SOFTENING
    );
  }

  // define correctly the parameters for each softening law

  if      ( softening == EXPONENT1_SOFTENING )
  {
    props.get ( alpha_,  ALPHA_PROP,  0.0, maxOf ( alpha_  ) );
    props.get ( beta_,   BETA_PROP,   0.0, maxOf ( beta_   ) );

    softening_ = Exponential3Params;
  }
  else if ( softening == EXPONENT2_SOFTENING )
  {
    props.get ( kappaC_, KAPPAC_PROP, 0.0, maxOf ( kappaC_ ) );

    softening_ = Exponential2Params;
  }
  else if ( softening == EXPONENT3_SOFTENING )
  {
    softening_ = Exponential2ParamsReg;
  }
  else if ( softening == POWER_SOFTENING )
  {
    props.get ( kappaC_, KAPPAC_PROP, 0.0, maxOf ( kappaC_ ) );
    props.get ( alpha_,  ALPHA_PROP,  0.0, maxOf ( alpha_  ) );
    props.get ( beta_,   BETA_PROP,   0.0, maxOf ( beta_   ) );

    softening_ = Power;
  }
  else if ( softening == HYPERBOLIC_SOFTENING )
  {
    props.get ( b_, B_PROP, 0.0, maxOf ( b_ ) );

    softening_ = Hyperbolic;
  }
  else if ( softening == LINEAR_SOFTENING )
  {
    softening_ = Linear;

    props.get ( kappaC_, KAPPAC_PROP, 0.0, maxOf ( kappaC_ ) );
  }
  else if ( softening == EXPOENERGY_SOFTENING )
  {
    props.get ( kappaI_,  KAPPAI_PROP,  0.0, maxOf ( kappaI_  ) );
    props.get ( kappaC_,  KAPPAC_PROP,  0.0, maxOf ( kappaC_  ) );
    props.get ( s_,       "s",          0.0, maxOf ( s_       ) );

    softening_ = ExpoEnergy;
  }

  // read type of equivalent strain definition

  String equiStrain;

  props.get ( equiStrain, EQUISTRAIN_PROP );

  if ( find ( equiStrainDefs.begin (), 
              equiStrainDefs.end   (), 
              equiStrain ) == equiStrainDefs.end () )
  {
    throw Error (
      JEM_FUNC,
      String ("unexpected definition of equivalent strain!!!\n") +
      String ("Available ones are: ") + MAZARS_EQUI_STRAIN + String(", ") +
      MISES_EQUI_STRAIN + String(", ") + RANKINE_EQUI_STRAIN
    );
  }

  if      ( equiStrain == MISES_EQUI_STRAIN )
  {
    props.get ( eta_, ETA_PROP, 0.0, maxOf ( eta_ ) );

    equiStrain_ = VonMises;
  }
  else if ( equiStrain == MAZARS_EQUI_STRAIN )
  {
    equiStrain_ = Mazars;
  }
  else if ( equiStrain == RANKINE_EQUI_STRAIN )
  {
    equiStrain_ = Rankine;
  }
  else if ( equiStrain == ENERGY_STRAIN )
  {
    equiStrain_ = StrainEnergy;
  }

  vonMisesEqv_.init ( eta_, elasticMat_->givePoisson ( ) );

  elasticMod_ = elasticMat_->getStiffMat();
    
  props.find ( threshold_, REMOVE_DAMAGE_PROP, 0.0, 1.0 );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void DamageMaterial::getConfig ( const Properties& conf,
                                 const Properties& globdat ) const
{
  elasticMat_->getConfig ( conf , globdat);

  conf.set ( SOFTENING_PROP,  softening_  );
  conf.set ( EQUISTRAIN_PROP, equiStrain_ );

  conf.set ( KAPPAI_PROP, kappaI_ );
  conf.set ( KAPPAC_PROP, kappaC_ );
  conf.set ( ALPHA_PROP,  alpha_  );
  conf.set ( BETA_PROP,   beta_   );
  conf.set ( ETA_PROP,    eta_    );
  conf.set ( REMOVE_DAMAGE_PROP, threshold_ );

  if ( c_ != 0 ) conf.set ( LENGTH_PROP, c_ );
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

// Remark: applicable to both local and nonlocal damage model (integral or gradient
// enhanced)

// Nonlinear solver procedure of jem/jive:
//  - beginning of every load step:
//      + compute K and fint with ZERO strain increment 
//  - iteration i:
//      + solve for displacement correction
//      + update the displacement
//      + compute K and fint for this new updated displacement
//      + check for convergence: if no, do again 
//
// Attention for the recomputation of K and fint at the beginning of load step
// 

void  DamageMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      int                   ipoint )
{
  // Get the nonlocal equivalent strain:
  //
  // The nonlocal equi. strain is computed at the global level and stored
  // in the final component of the strain vector, hence
  // equiStrain = strain[n];

  const int n          = stress.size ();
  const int m          = strain.size ();

  double    equiStrain;

  if ( n == m ) // local damage models
  {
    equiStrain = getEquiStrain ( strain );
  }
  else
  {
    equiStrain = strain[n];
  }

  // Get the history variable of the previous converged load step

  double    kappa0     = preHist_.eqveps[ipoint];

  // Compute the loading function f

  double    f          = equiStrain - kappa0;

  double kappa;
  int    load;

  // f == 0 (exactly equals to zero), at the beginning of every 
  // loading step, this condition holds for loading integration point

  if ( f == 0.0 )
  {
    kappa = kappa0;
    load  = preHist_.loading[ipoint];

    // System::out() << " kappa0             " << kappa0      <<"\n";
    // System::out() << " equivalent strain  " << equiStrain  <<"\n";
    // System::out() << " loading function    " << f          <<"\n";
  }  
  else 
  {
    // avoid spurious loading/unloading in 1d softening bar

    //  if ( abs (f) < EPS )
    //  {
    //    f = 0.0;
    //  } 
  
    // f > 0: loading, update history variable

    if ( f > 0.0 )
    {  
      kappa = equiStrain;  
      load  = ( kappa < kappaI_ ) ? 0 : 1;
      // System::out() << kappa0 << " "  << kappa << " "  << ipoint  <<"\n\n";
    }
 
    // f < 0: unloading, keep the same history variable

    else 
    {
      kappa = kappa0;
      load  = 0;
    }
  }

  // System::out() << " kappa0             " << kappa0      <<"\n";
  // System::out() << " equivalent strain  " << equiStrain  <<"\n";
  // System::out() << " kappa              " << kappa       <<"\n";

  // compute the damage variable

  double omega  = damageEvolution_ ( kappa );

  // secant stiffness matrix only
  // the second term will be computed by non-local integral  model

  stiff         = ( 1.0 - omega ) * elasticMod_;

  // Compute stress vector

  MatmulChain<double,1>     mc1;

  stress = ( 1.0 - omega ) * mc1.matmul ( elasticMod_, strain[slice(BEGIN,n)] );

  // update history variables

  newHist_.eqveps[ipoint]  = kappa;
  newHist_.loading[ipoint] = load;
      
  //System::out() << omega << "\n";

  if ( jem::Float::isNaN( sum(  stiff  ) ) ||
       jem::Float::isNaN( sum(  stress ) ) )
  {
      System::out() << omega << " " << strain << "\n";
     System::out() << "Invalid matrix element" <<
        "\n stress: " << stress <<
         "\n stiff   :\n" << stiff << endl;
  }

}

// ----------------------------------------------------------------
//  update (overloaded version)
// ----------------------------------------------------------------

void  DamageMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      int                   ipoint,
      double                he )
{
  MatmulChain<double,1>     mc1;
  MatmulChain<double,2>     mc2;

  double    equiStrain = getEquiStrain ( strain ); // Get the local equivalent strain:
  double    kappa0     = preHist_.eqveps[ipoint];  // Get the history variable of the previous converged load step
  double    f          = equiStrain - kappa0;      // Compute the loading function f

  double kappa;
  int    load;

  // f == 0 (exactly equals to zero), at the beginning of every 
  // loading step, this condition holds for loading integration point

  if ( f == 0.0 )
  {
    kappa = kappa0;
    load  = preHist_.loading[ipoint];
  }  
  else 
  {
    if ( f > 0.0 )
    {
      kappa = equiStrain;  
      load  = ( kappa < kappaI_ ) ? 0 : 1;
    }
 
    // f < 0: unloading, keep the same history variable

    else 
    {
      kappa = kappa0;
      load  = 0;
    }
  }

  // compute the damage variable

  double omega  = damageEvolution_ ( kappa );

  // secant stiffness matrix only
  // the second term will be computed by non-local integral  model

  stiff         = ( 1.0 - omega ) * elasticMod_;
        
  if ( load )
  {
     double  dOmega      = getdOmegadKappa ( kappa   );
     Vector  dEpsBardEps = getdEpsBardEps  ( strain  );	
  
     //System::out() << " strain              " << strain       <<"\n";
     //System::out() << " depv              " << dEpsBardEps       <<"\n";

     stiff -= dOmega * mc2.matmul(elasticMod_,matmul(strain,dEpsBardEps));	
  }

  // Compute stress vector

  stress = ( 1.0 - omega ) * mc1.matmul ( elasticMod_, strain );

  // update history variables

  newHist_.eqveps[ipoint]  = kappa;
  newHist_.loading[ipoint] = load;
   
  //System::out() << " stress              " << stress       <<"\n";
  
  /*
  if ( jem::Float::isNaN( sum(  stiff  ) ) ||
       jem::Float::isNaN( sum(  stress ) ) )
  {
      System::out() << omega << " " << strain << "\n";
     System::out() << "Invalid matrix element" <<
        "\n stress: " << stress <<
         "\n stiff   :\n" << stiff << endl;
  }*/

  if ( ipoint == 0 ) {
      System::out() << strain <<  "\n";
      System::out() << load <<  "\n";
      System::out() << omega <<  "\n";
  }

}

// ----------------------------------------------------------------
//  update (overloaded version)
// ----------------------------------------------------------------

void  DamageMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain0,
      const Vector&         dstrain,
      int                   ipoint,
      double                 he )
{
    Vector tstrain (STRAIN_COUNTS[rank_]);
    tstrain = strain0 + dstrain;
    update ( stress, stiff, tstrain, ipoint, he );
}

// --------------------------------------------------------------------
//  commit
// --------------------------------------------------------------------

void  DamageMaterial::commit()
{
  newHist_.eqveps. swap ( preHist_.eqveps  );
  newHist_.loading.swap ( preHist_.loading );
}


// --------------------------------------------------------------------
//  checkLocalisation
// --------------------------------------------------------------------

double DamageMaterial :: checkLocalisation
 
  (        Vector&    normal,
    const  Vector&    stress,
    const  Matrix&    tangent,
           int        ip )   const
{
    double damage        = giveOmega ( ip );
    bool   localised     = ( damage > 0 ) ? true : false;
        
    //System::out() << damage << "\n";

    if ( localised )
    {
        Vector sigmaI ( 3 );
        Matrix dir    ( 3, 3 );
        elasticMat_ -> getPrincipalStressAndDirections ( stress, sigmaI, dir );

        int    i   = 0;
        double max = sigmaI[0];

        if ( max < sigmaI[1] ){
            max = sigmaI[1]; i = 1;
        }
        
        if ( max < sigmaI[2] ){
            max = sigmaI[2]; i = 2;
        }

        normal = dir(slice(BEGIN,rank_),i);
        //System::out() << sigmaI << "\n";
        //System::out() << damage << "\n";
        //normal[0]=1.;normal[1]=0.;normal[2]=0.;
    }

    return damage;
}

// -------------------------------------------------------------------
//  giveElasticMatrix
// -------------------------------------------------------------------

Matrix  DamageMaterial::giveElasticMatrix ( ) const
{
  return elasticMod_;
}

//-----------------------------------------------------------------------
//   damageEvolution_
//-----------------------------------------------------------------------

// Remark: The following Fortran-like implementation is ugly
// but efficient. There is an alternative using policy-based
// class by writing orthogonal classes for every softening
// law. It is, however, slower than the current implementation.

double DamageMaterial::damageEvolution_

  ( double k ) const

{
  double omega = 0.0;

  if      ( softening_ == Perfect )
  {
    omega = perfectSoftening( kappaI_, k ); 
  }
  else if ( softening_ == Linear )
  {
    omega = linearSoftening( kappaI_, kappaC_ , k );
  } 
  else if ( softening_ == Exponential3Params )
  {
    omega = exponentialSoftening1 ( kappaI_, alpha_, beta_, k );
  }
  else if ( softening_ == Exponential2Params )
  {
    omega = exponentialSoftening2 ( kappaI_, kappaC_ , k );
  }
  else if ( softening_ == Power )
  {
    omega = powerSoftening ( kappaI_, kappaC_, alpha_, beta_, k );
  }
  else if ( softening_ == Hyperbolic )
  {
    omega = hyperbolicSoftening ( kappaI_, b_, k );
  }
  else if ( softening_ == ExpoEnergy )
  {
    omega = expoEnergySoftening ( kappaI_, kappaC_, s_, k );
  }

  return omega;
}

//-----------------------------------------------------------------------
//   damageEvolution_
//-----------------------------------------------------------------------

double DamageMaterial::damageEvolution_

  ( double k, double he ) const

{
  double omega = 0.0;

  if      ( softening_ == Exponential3Params )
  {
    omega = exponentialSoftening1 ( kappaI_, alpha_, gf_, ft_ , k, he );
  }
  else if ( softening_ == Exponential2Params )
  {
    omega = exponentialSoftening2 
             ( kappaI_, kappaC_ , k, lambda_, he );
  }
  else if ( softening_ == Exponential2ParamsReg )
  {
    omega = exponentialSoftening3 ( kappaI_, k, gf_, ft_, he );
  }

  return omega;
}

//-----------------------------------------------------------------------
//   getdOmegadKappa
//-----------------------------------------------------------------------

double DamageMaterial:: getdOmegadKappa

    ( double k )  const

{
  double dOmega = 0.0;

  if           ( softening_ == Perfect  )
  {
    dOmega = getDerivPerfectSoftening( kappaI_, k );
  } 
  else if      ( softening_ ==  Linear )
  {
    dOmega = getDerivLinearSoftening( kappaI_, kappaC_ , k );
  } 
  else if      ( softening_ == Exponential3Params )
  {
    dOmega = getDerivExponentialSoftening1 ( kappaI_, alpha_, beta_, k );
  }
  else if      ( softening_ == Exponential2Params )
  {
    dOmega = getDerivExponentialSoftening2 ( kappaI_, kappaC_, k );
  }
  else if      ( softening_ == Power )
  {
    dOmega = getDerivPowerSoftening ( kappaI_, kappaC_, alpha_, beta_, k );
  }
  else if      ( softening_ == Hyperbolic )
  {
    dOmega = getDerivHyperbolicSoftening ( kappaI_, b_,k );
  }
  else if      ( softening_ == ExpoEnergy )
  {
    dOmega = getDerivExpoEnergySoftening ( kappaI_, kappaC_, s_, k );
  }

  return dOmega;
}

//-----------------------------------------------------------------------
//   getdOmegadKappa (for regularised local damage)
//-----------------------------------------------------------------------

double DamageMaterial:: getdOmegadKappa

    ( double k, double he )  const

{
  double dOmega = 0.0;

  if       ( softening_ == Exponential3Params )
  {
    dOmega = getDerivExponentialSoftening1 ( kappaI_, alpha_, gf_, ft_, k, he );
  }
  else if  ( softening_ == Exponential2Params )
  {
    dOmega = getDerivExponentialSoftening2 
                ( kappaI_, kappaC_, k, lambda_, he );
  }
  else if  ( softening_ == Exponential2ParamsReg )
  {
    dOmega = getDerivExponentialSoftening3
                ( kappaI_, k, gf_, ft_, he );
  }

  //System::out() << dOmega <<"\n";

  return dOmega;
}

//-----------------------------------------------------------------------
//   getEquiStrain
//-----------------------------------------------------------------------

double DamageMaterial::getEquiStrain

  ( const Vector& strain ) 

{
  double eqvStr = 0.0;

  if      ( equiStrain_ == VonMises  )
  {
    double I1, J2;

    elasticMat_->getI1andJ2  ( strain, I1, J2 );  

    eqvStr     = vonMisesEqv_ ( I1, J2 );
  } 
  else if ( equiStrain_ == Mazars  )
  {
    // code for 1d bar in tension !!!
    // then the only positive principal strain is the axial strain 
    // to run 1D, you have to uncomment the following

/*
    if ( strain.size() == 1 )
    {
      return  ::fabs ( strain[0] );
    }
*/

    Vector princ = elasticMat_->getPrincipalStrains ( strain );

    eqvStr  = 0.0;

    for ( int i = 0; i < 3; i++ )
    {
      double x =  princ[i] ;

      if ( x > 0 )
      {
	    eqvStr += x * x ;
      }
    }

    eqvStr = ::sqrt ( eqvStr );
  }

  // equivalent strain according to Rankine def
  // based on positive principal stress

  else if ( equiStrain_ == Rankine  )
  {
    MatmulChain<double,1>     mc1;

    double young  = elasticMat_->giveYoung ();

    Vector stress = mc1.matmul ( elasticMod_, strain );
    Vector newStr = elasticMat_->giveStressVector   ( stress );
    Vector princ  = elasticMat_->getPrincipalValues ( newStr );

    eqvStr        = equiStrainRankine ( young, princ );
  }

  else if ( equiStrain_ == StrainEnergy  )
  {
    MatmulChain<double,1>     mc1;

    eqvStr        = 0.5 * dot ( strain, mc1.matmul ( elasticMod_, strain ) );
  }
 
  return eqvStr;
}

// --------------------------------------------------------------------
//  getdEpsBardEps (derivatives of equivalent strain w.r.t strain vector)
// --------------------------------------------------------------------

Vector DamageMaterial::getdEpsBardEps

  ( const Vector& strain ) 

{
  const int strCount = strain.size ( );

  Vector  ret ( strCount );

  if      ( equiStrain_ == VonMises  )
  {
    double I1, J2;

    Vector dI1( strCount );
    Vector dJ2( strCount );

    elasticMat_->getI1J2andGrads ( strain, I1, J2, dI1, dJ2 );

    ret  = vonMisesEqv_ ( dI1, dJ2, I1, J2 );
  }

  else if ( equiStrain_ == Mazars  )
  {
    // code for 1d bar in tension !!!
    // then the only positive principal strain is the axial strain 

    if ( strCount == 1 ) 
    {
      ret = evalHeaviside ( strain[0] ); 
    }
    else
    {
      ret = getdEquiStraindEpsilon1_ ( strain );
    }
  } 

  else if ( equiStrain_ == Rankine  )
  {
    MatmulChain<double,1>     mc1;

    double young  = elasticMat_->giveYoung   ();

    Vector stress = mc1.matmul ( elasticMod_, strain );

    ret           = getdEquiStraindEpsilon3_ ( stress, young, elasticMod_ );
  }

  else if ( equiStrain_ == StrainEnergy  )
  {
    MatmulChain<double,1>     mc1;

    ret          = mc1.matmul ( elasticMod_, strain );;
  }

  return ret;
}

// ---------------------------------------------------------------
//  isFullyDamaged ( int mpoint )
// ---------------------------------------------------------------

bool DamageMaterial::isFullyDamaged ( int mpoint ) const
{
  double kappa = newHist_.eqveps[ mpoint ];
  double omega = damageEvolution_ ( kappa );

  return ( 1.0 - omega < threshold_ ) ? true : false ;
}

// ---------------------------------------------------------------
//  compute the derivatives of equivalent strain
//  w.r.t the strain vector
//  Mazars formulas
// ---------------------------------------------------------------

Vector  DamageMaterial:: getdEquiStraindEpsilon1_

  ( const Vector&   strain )

{
  const int s = strain.size();

  Vector   ret ( s );

  ret   = 0.0 ;

  // 3D problem

  if ( s == 6 )
  {
    Vector  inv    (2)  ;
    Vector  princ  (3);

    Vector    dI1  (6);
    Vector    dI2  (6);
    Vector    dI3  (6);

    Vector    dEpsIdEps (6);

    elasticMat_->getIandDIDStrain ( strain, inv, princ, dI1, dI2, dI3 );

    double epsBar = 0.0;

    for ( int i = 0; i < 3; i++ )
    {
      double x =  princ[i] ;

      if ( x > 0 )
      {
	    epsBar += x * x ;
      }
    }

    epsBar = ::sqrt ( epsBar );

    // check for case of numericallly zero strain vector

    if ( epsBar < EPS )
    {
      return ret; 
    }

    double I1 = inv[0];
    double I2 = inv[1];

    double vi, temp;
 
    for ( int i = 0 ; i < 3; i++ )
    {
      vi  = princ[i];

      if ( vi > 0 )
      {
	    temp      = 1.0 / ( 3.0 * vi * vi - 2.0 * I1 * vi + I2);
	    dEpsIdEps = temp * ( vi * vi * dI1 - vi * dI2 + dI3 );
	   ret       += vi * dEpsIdEps ;
      }
    }

    ret          *= (1.0 / epsBar) ;
  }
  else
  {
    // 2D: strain={eps_xx, eps_yy, eps_zz, eps_xy}  

    double exx =       strain[0];
    double eyy =       strain[1];
    double exy = 0.5 * strain[3];

    double exxMeyy = exx - eyy;
    double exxPeyy = exx + eyy;

    double poi  = elasticMat_->givePoisson();
    double poiD = poi / ( poi - 1.0 );

    // principal strains are roots of a quadratic equation with det = d

    double d   = ::sqrt ( exxMeyy * exxMeyy + 4.0 * exy * exy );

    double prinstr0 = 0.5 * ( exxPeyy + d );
    double prinstr1 = exxPeyy - prinstr0;
    double prinstr2 = elasticMat_->giveState() == PlaneStress ? poiD * exxPeyy : 0.0;

    prinstr0 = prinstr0 < 0. ? 0. : prinstr0;
    prinstr1 = prinstr1 < 0. ? 0. : prinstr1;
    prinstr2 = prinstr2 < 0. ? 0. : prinstr2;

    double den = ::sqrt ( prinstr0 * prinstr0 + prinstr1 * prinstr1 + prinstr2 * prinstr2 ); 

    // check for case of numericallly zero strain vector

    if ( den < 1e-08 )
    {
/*      using jem::io::NumberFormat;

      NumberFormat nformat;

      nformat.setScientific     (   );
      nformat.setFractionDigits ( 8 );

      System::out() << "mazars: so small eqv strain " 
                    << nformat.print ( den ) << "\n"; 

      System::out() << "mazars: determinant d so small " 
                    << nformat.print ( d ) << "\n"
                    << prinstr << "\n"; 
*/
      return ret; 
    }

    // derivatives of equivalent strain w.r.t principal strains

    double denInv = 1. / den;

    double dedprin0 = prinstr0 * denInv;
    double dedprin1 = prinstr1 * denInv;
    double dedprin2 = prinstr2 * denInv;

    // derivatives of principal strains w.r.t strain tensor

    // avoid division by zero

    if ( d == 0.0 )
    {
      d = 1.0;
    }

    double fac     =  0.5 / d; 

    double de1dexx = 0.5 + fac * exxMeyy;
    double de1deyy = 1.0 - de1dexx;
    double de1dexy = 2.0 * fac * exy ; // attention here, 2 = 4 (from derivation) * 0.5

    double de2dexx =  de1deyy;
    double de2deyy =  de1dexx;
    double de2dexy = -de1dexy;

    double de3dexx = poiD;
    double de3deyy = de3dexx;

    // finally, derivatives of equivalent strain w.r.t strain tensor
    // by chain rule

    ret[0] = dedprin0 * de1dexx +  dedprin1 * de2dexx +  dedprin2 * de3dexx;
    ret[1] = dedprin0 * de1deyy +  dedprin1 * de2deyy +  dedprin2 * de3deyy;
    ret[3] = dedprin0 * de1dexy +  dedprin1 * de2dexy;
  }

  return ret;
}


// ---------------------------------------------------------------
//  compute the derivatives of equivalent strain
//  w.r.t the strain vector
//  Rankine formulas
// ---------------------------------------------------------------

Vector DamageMaterial::getdEquiStraindEpsilon3_

  ( const Vector&   stress,
    const double    young,
    const Matrix&   De )

{
  MatmulChain<double,1>     mc1;

  const int s       =  stress.size();

  Vector    ret( s );

  ret = 0.0;

  // 3D problem

  if ( s == 6 )
  {
    Vector    dSigIdEps ( s );

    Vector    inv;
    Vector    princ;

    double    vi,temp;
    double    eqvStr;

    elasticMat_->getPrincipalValues ( inv, princ, stress );

    // compute the equivalent strain

    eqvStr  = equiStrainRankine ( young, princ ) ;

    if ( eqvStr < EPS )
    {
      return ret;
    }

    // compute derivatives of invariants w.r.t stress vector

    Vector dI1dSigma = elasticMat_->getDI1DStrain ( stress );
    Vector dI2dSigma = elasticMat_->getDI2DStrain ( stress );
    Vector dI3dSigma = elasticMat_->getDI3DStrain ( stress );

    for ( int i = 0; i < 3; i++)
    {
      vi = princ[i];

      if ( vi > 0 )
      {
        temp      = 1.0 / ( 3.0 * vi * vi -
                    2.0 * inv[0] * vi + inv[1] );

        dSigIdEps = temp * (  vi * vi * dI1dSigma -
                    vi * dI2dSigma + dI3dSigma );

        ret       += vi * dSigIdEps ;
      }
    }

    ret  = mc1.matmul ( De, ret );
    ret *= 1.0 / ( young * young * eqvStr);
  }

  else // plane problem
  {
    Vector     prinstr(3);      // principal stresses
    Vector     dedprin(3);      // derivatives of eqv strain w.r.t principal stresses

    double sxx = stress[0];
    double syy = stress[1];
    double sxy = stress[2];

    double poi  = elasticMat_->givePoisson();
    double eInv = 1.0 / young;

    // principal stresses are root of a quadratic equation with det = d

    double d   = ::sqrt ( ( sxx - syy ) * ( sxx - syy ) + 4.0 * sxy * sxy );

    prinstr[0] = 0.5 * ( sxx + syy + d );
    prinstr[1] = 0.5 * ( sxx + syy - d );

    //ProblemType state = elasticMat_->giveState();

    if ( elasticMat_->giveState() == PlaneStrain )
    {
      prinstr[2] = poi * ( sxx + syy );
    }
    else
    {
      prinstr[2] = 0.0;
    }

    // Rankine equivalent strain, den = E * eqv strain

    double den = 0.0;

    for ( int i = 0; i < 3; i++ )
    {
      if ( prinstr[i] < 0.0 ) 
      {
        prinstr[i] = 0.0;
      }

      den += prinstr[i] * prinstr[i];
    }

    den = ::sqrt ( den );

    // without the following, get invalid elements (NaN)

    if ( den < 1e-012 )
    {
/*      using jem::io::NumberFormat;

      NumberFormat nformat;

      nformat.setScientific     (   );
      nformat.setFractionDigits ( 8 );

      System::out() << "rankine: so small eqv strain " 
                    << nformat.print ( den ) << " " 
                    << nformat.print ( d ) << "\n"; 
*/
      return ret;
    }

    // derivatives of equivalent strain w.r.t principal stresses

    for ( int i = 0; i < 3; i++ )
    {
      dedprin[i] = prinstr[i] / den;
    }

    // derivatives of principal stresses w.r.t stress tensor

    if ( d == 0.0 )
    {
      d = 1.0;
    }

    double fac     =  0.5 / d; 

    double ds1dsxx = 0.5 + fac * ( sxx - syy );
    double ds1dsyy = 0.5 + fac * ( syy - sxx );
    double ds1dsxy = 4.0 * fac * sxy ; // 4 or 2???

    double ds2dsxx =  ds1dsyy;
    double ds2dsyy =  ds1dsxx;
    double ds2dsxy = -ds1dsxy;

    double ds3dsxx = poi;
    double ds3dsyy = poi;

    // derivatives of stress w.r.t strain

    //double muy = 0.5 * young / ( 1.0 + poi );
    //double lam = young * poi / ( 1.0 + poi ) / ( 1.0 - 2.0 * poi );


    double dsxxdexx = De(0,0);
    double dsxxdeyy = De(0,1);

    double dsyydexx = dsxxdeyy;
    double dsyydeyy = De(1,1);

    double dsxydexy = De(2,2);

/*
    if ( state == "PLANE_STRAIN" )
    {
      dsxxdexx = lam + 2.0 * muy;
      dsxxdeyy = lam;

      dsyydexx = lam;
      dsyydeyy = dsxxdexx;

      dsxydexy = 0.5 * muy;
    }
    else
    {
      double a = 2.0 * muy / ( lam + 2.0 * muy ); 

      dsxxdexx = 2.0 * ( lam + muy ) * a;
      dsxxdeyy = lam *a;

      dsyydexx = dsxxdeyy;
      dsyydeyy = dsxxdexx;

      dsxydexy = 0.5 * muy *a;
    }
*/

    double ds1dexx = ds1dsxx * dsxxdexx + ds1dsyy * dsyydexx;
    double ds1deyy = ds1dsxx * dsxxdeyy + ds1dsyy * dsyydeyy;
    double ds1dexy = ds1dsxy * dsxydexy;

    double ds2dexx = ds2dsxx * dsxxdexx + ds2dsyy * dsyydexx;
    double ds2deyy = ds2dsxx * dsxxdeyy + ds2dsyy * dsyydeyy;
    double ds2dexy = ds2dsxy * dsxydexy;

    double ds3dexx = ds3dsxx * dsxxdexx + ds3dsyy * dsyydexx;
    double ds3deyy = ds3dsxx * dsxxdeyy + ds3dsyy * dsyydeyy;

    // finally, derivatives of equivalent strain w.r.t strain tensor
    // by chain rule

    ret[0] = dedprin[0] * ds1dexx +  dedprin[1] * ds2dexx +  dedprin[2] * ds3dexx;
    ret[1] = dedprin[0] * ds1deyy +  dedprin[1] * ds2deyy +  dedprin[2] * ds3deyy;
    ret[2] = dedprin[0] * ds1dexy +  dedprin[1] * ds2dexy;

    ret *= eInv;
  }

  return ret;
}


// ======================================================================
//   Implementation of related functions
// ======================================================================


//-----------------------------------------------------------------------
//   perfectSoftening
//-----------------------------------------------------------------------

double                      perfectSoftening

  ( double Ki, double K )
{
  if ( K < Ki )
  {
    return 0.0;
  }

  double omega = 1.0 - Ki / K;

  return omega;
}

//-----------------------------------------------------------------------
//   getDerivPerfectSoftening
//-----------------------------------------------------------------------

double                      getDerivPerfectSoftening

    ( double Ki, double K )

{
  if ( K < Ki )
  {
    return 0.0;
  }

  return  Ki / ( K * K );
}


//-----------------------------------------------------------------------
//   linearSoftening
//-----------------------------------------------------------------------

double                      linearSoftening

  ( double Ki,
    double Kc,
    double K )
{
  if ( K < Ki )
  {
    return 0.0;
  }

  double omega = K > Kc ? 1.0 : Kc * (K - Ki) / (K * (Kc - Ki));
 
  if ( omega == 1.0 )
  {
    omega = DamageMaterial::CRITICAL_DAMAGE;
  }

  return omega;
}

//-----------------------------------------------------------------------
//   getDerivLinearSoftening
//-----------------------------------------------------------------------

double                      getDerivLinearSoftening

    ( double Ki,
      double Kc,
      double K )

{
  return  (K < Kc) ?  ( Ki * Kc ) / ( (Kc - Ki) * K * K ) : 0.0;
}


//-----------------------------------------------------------------------
//   exponentialSoftening1
//-----------------------------------------------------------------------

double                      exponentialSoftening1

  ( double Ki,
    double alpha,
    double beta,
    double K )
{
  if ( K < Ki ) return 0.0;

  double omega = 1.0 - (Ki / K) * ( 1.0 - alpha + alpha * exp (-beta * (K - Ki) ) );

  //return isTiny ( 1.0 - omega ) ? DamageMaterial::CRITICAL_DAMAGE : omega;

  return omega;
}

// overloaded version

double                      exponentialSoftening1

  ( double Ki,
    double alpha,
    double gf, double ft,
    double K,
    double he )
{
  if ( K < Ki ) return 0.0;

  double beta  = alpha / ( -0.5 * Ki + gf / ft /he );

  double omega = 1.0 - (Ki / K) * ( 1.0 - alpha + alpha * exp (-beta * (K - Ki) ) );

  return omega;
}

//-----------------------------------------------------------------------
//   getDerivExponentialSoftening1
//-----------------------------------------------------------------------

double                      getDerivExponentialSoftening1

  ( double Ki,
    double alpha,
    double beta,
    double K )
{
  const double h = K - Ki;	
  const double a = Ki / (K * K);
  const double b = a * K;
  const double c = alpha * exp(-beta * h );

  return a * (1.0 - alpha + c) + beta * b * c;
}

// overloaded version

double                      getDerivExponentialSoftening1

  ( double Ki,
    double alpha,
    double gf, double ft,
    double K,
    double he )
{
  double beta    = alpha / ( -0.5 * Ki + gf / ft /he );

  const double h = K - Ki;	
  const double a = Ki / (K * K);
  const double b = a * K;
  const double c = alpha * exp(-beta * h );

  return a * (1.0 - alpha + c) + beta * b * c;
}

//-----------------------------------------------------------------------
//    exponentialSoftening2
//-----------------------------------------------------------------------

double                       exponentialSoftening2

  ( double Ki,
    double Kc,
    double K )
{
  if ( K < Ki )
  {
    return 0.0;
  }

  double omega = 1.0 - (Ki / K) * exp ( (Ki - K) / (Kc - Ki) );

  return omega;
}

// overloaded version

double                       exponentialSoftening2

  ( double Ki,
    double Kc,
    double K, double lambda,
    double he )
{
  if ( K < Ki )
  {
    return 0.0;
  }

  double KC    = lambda / he * ( Kc - 0.5 * Ki ) + 0.5 * Ki;

  double omega = 1.0 - (Ki / K) * exp ( (Ki - K) / (KC - Ki) );

  return omega;  
}

//-----------------------------------------------------------------------
//    getDerivExponentialSoftening2
//-----------------------------------------------------------------------

double                       getDerivExponentialSoftening2

    ( double Ki,
      double Kc,
      double K )

{ 
  double a = Ki - K;
  double b = Kc - Ki;
  double c = exp ( a / b );
  double d = Ki / (K * K);

  return d * c + d * K * c / b;
}

// overloaded version

double                       getDerivExponentialSoftening2

    ( double Ki,
      double Kc,
      double K, double lambda,
      double he )

{
  double KC= lambda / he * ( Kc - 0.5 * Ki ) + 0.5 * Ki;

  double a = Ki - K;
  double b = KC - Ki;
  double c = exp ( a / b );
  double d = Ki / (K * K);

  return d * c + d * K * c / b;
}

// Jirasek regularised exponential softening law

double                  exponentialSoftening3

  ( double Ki, double K,
    double gf, double ft, 
    double he )
{
  if ( K < Ki )
  {
    return 0.0;
  }

  double Kc    = 0.5 * Ki + gf / ft / he;
  double omega = 1.0 - (Ki / K) * exp ( (Ki - K) / (Kc - Ki) );

  return  omega;
}

double                  getDerivExponentialSoftening3

   ( double Ki, double K,
     double gf, double ft,
     double he )
{
  double Kc = 0.5 * Ki + gf / ft / he;

  double a  = Ki - K;
  double b  = Kc - Ki;
  double c  = exp ( a / b );
  double d  = Ki / (K * K);

  return d * c + d * K * c / b;
}

//-----------------------------------------------------------------------
//    hyperbolicSoftening
//-----------------------------------------------------------------------

double                       hyperbolicSoftening

  ( double Ki,
    double b,
    double K )
{
  if ( K < Ki )
  {
    return 0.0;
  }

  double omega = 1.0 - 1.0 / ( 1.0 + b * ( K - Ki ) );

  return isTiny ( 1.0 - omega ) ? DamageMaterial::CRITICAL_DAMAGE : omega;  
}

//-----------------------------------------------------------------------
//    getDerivExponentialSoftening2
//-----------------------------------------------------------------------

double                       getDerivHyperbolicSoftening

    ( double Ki,
      double b,
      double K )

{

  return b / ( 1.0 + b * ( K - Ki ) ) / ( 1.0 + b * ( K - Ki ) ) ;
 
}


//-----------------------------------------------------------------------
//  powerSoftening
//-----------------------------------------------------------------------

// particularly implemented for the composite compact tension test

double                      powerSoftening

  ( double Ki,
    double Kc,
    double alpha,
    double beta,
    double K )
{
  double omega;

  if ( K < Ki )
  {
    omega = 0.0;
  }
  else
  {
    double Kcc    =  Kc;
    double omegac = 1.0 - ::pow ( Ki / Kcc, beta ) * 
                          ::pow ( (Kc - Kcc) / (Kc - Ki) , alpha );

    omega  = (K > Kcc) ? omegac : 1.0 -::pow ( Ki / K, beta ) * 
                                       ::pow ( (Kc - K) / (Kc - Ki) , alpha );
  }

  return omega; 

  /*

  // the damage evolution in Geers - Enhanced solution control

  double omega;

  if ( K < Ki )
  {
    omega = 0.0;
  }
  else
  {
    const double gamma = 0.01;
    const double eta   = 1.0 - gamma;
    const int    n     = -5;

    double  muy  = ( gamma - n * eta ) / ( n * eta * ::pow ( Ki, n ) );

    omega        = 1.0 - gamma * Ki / K - eta * ::pow ( K / Ki, n) *
                    ::exp ( muy * ( ::pow ( K, n ) - ::pow( Ki, n ) ) );  
  }
  
  return omega; 
  */
  
}

//-----------------------------------------------------------------------
//  getDerivPowerSoftening
//-----------------------------------------------------------------------

double                      getDerivPowerSoftening

  ( double Ki,
    double Kc,
    double alpha,
    double beta,
    double K )
{
  double ret;
  double Kcc       =  Kc;

  if ( K > Kcc )
  {
    ret =  0.0;
  }
  else
  {
    const double a = Ki / K ;
    const double b = 1.0 / ( Kc - Ki );
    const double c = ( Kc - K ) * b;
    const double d = ::pow ( a, beta  - 1.0 );
    const double e = ::pow ( c, alpha - 1.0 );

    ret  = d * e * ( (a / K) * beta * c  +  a * alpha * b );
  }

  return ret ; 

/*
  const double gamma = 0.01;
  const double eta   = 1.0 - gamma;
  const int    n     = -5;

  double  muy  = ( gamma - n * eta ) / ( n * eta * ::pow ( Ki, n ) );

  double fac1  =  exp ( muy * ( ::pow ( K, n )  - ::pow( Ki, n ) ) );

  double ret   =  gamma * Ki / K / K - eta * n * ::pow ( K / Ki, n-1) / Ki * fac1
                 - eta * pow ( K / Ki, n) * fac1 * muy * n * ::pow ( K, n-1);
 
  return ret;
*/
}

//-----------------------------------------------------------------------
//   expoEnergySoftening
//-----------------------------------------------------------------------


double                     expoEnergySoftening

 ( double kappaI, double kappaC,
   double s, double kappa )
{
  if ( kappa < kappaI ) return 0.;

  double tem = pow ( ( kappa - kappaI ) / kappaC, s );

  return ( 1. - exp ( - tem ) );
}

//-----------------------------------------------------------------------
//   getDerivExpoEnergySoftening
//-----------------------------------------------------------------------


double                      getDerivExpoEnergySoftening

  ( double kappaI, double kappaC,
    double s, double kappa )
{
  if ( kappa < kappaI ) return 0.;

  double tem1 = pow ( ( kappa - kappaI ) / kappaC, s   );
  double tem2 = pow ( ( kappa - kappaI ) / kappaC, s-1 );

  return  exp ( - tem1 ) * tem2 * s / kappaC;
}

//-----------------------------------------------------------------------
//   equiStrainRankine
//-----------------------------------------------------------------------

double                      equiStrainRankine

  ( const double  young,
    const Vector& princ )

{
  double ret   = 0.0;

  const  int s = princ.size ();

  for ( int i = 0; i < s; i++ )
  {
    double x =  princ[i] ;

    if ( x > 0 )
    {
      ret += x * x ;
    }
  }

  return ::sqrt ( ret ) / young;
}



//-----------------------------------------------------------------------
//   evalMcAuley
//-----------------------------------------------------------------------

double                evalMcAuley

  ( const double x )
{
  return (x > 0 ? x : 0.0);
}

//-----------------------------------------------------------------------
//    evalHeaviside
//-----------------------------------------------------------------------

double               evalHeaviside

  ( const double x )
{
 return (x < 0 ? -1.0 : 1.0);
} 
