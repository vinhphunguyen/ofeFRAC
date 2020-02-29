
#include "UserModels.h"


//-----------------------------------------------------------------------
//   declareModels
//-----------------------------------------------------------------------


void declareUserModels ()
{
  declareArclenModel                 ();
  declareEnergyArclenModel           ();
  declareLoadArclenModel             ();
  declareDispArclenModel             ();

  declareLineLoadModel               ();

  declareMonitorModel                ();
  declareLodiModel                   ();
  declareNodeConstraintModel         ();
  declareDofConstraintModel          ();

  declareAdaptLoadScaleModel         ();
}
