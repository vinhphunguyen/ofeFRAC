
#include "FemModels.h"


//-----------------------------------------------------------------------
//   declareModels
//-----------------------------------------------------------------------


void declareFemModels ()
{
  declareElasticityModel             ();
  declareElasticityModelNURBS        ();
  declareIGAElasticityModel          ();

  declareIGAFiniteDeformationModel   ();


  declareLocalDamageModel            ();
  declarePureLocalDamageModel        ();
  declareNonLocalDamageModel         ();
  declareGradientDamageModel         ();
  declareStressBasedGEDModel         ();

  declareInterfaceElementModel       ();
  declareDiscreteInterfaceElementModel       ();
  declareInterfaceElement3DModel     ();

  declarePoroInterfaceElementModel   ();
  declareRGPoroInterfaceElementModel ();
  declarePerforationInterfaceModel   ();


  declareTwoNodeTrussModel           ();
}
