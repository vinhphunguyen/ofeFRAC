#ifndef FEM_MODELS_H
#define FEM_MODELS_H


//-----------------------------------------------------------------------
//   public functions
//-----------------------------------------------------------------------


void  declareFemModels                   ();

void  declareElasticityModel             ();
void  declareElasticityModelNURBS        ();

// NURBS/T-splines IGA models

void  declareIGAElasticityModel          ();
void  declareIGAFiniteDeformationModel   ();
// non-local integral damage model
// and gradient-enhanced damage models

void  declareNonLocalDamageModel         ();
void  declareLocalDamageModel            (); // regularised local 
void  declarePureLocalDamageModel        ();
void  declareGradientDamageModel         ();
void  declareStressBasedGEDModel         ();

// zero-thickness cohesive interface element models

void  declareInterfaceElementModel       ();
void  declareDiscreteInterfaceElementModel ();
void  declareInterfaceElement3DModel     ();

// poro-mechanical models with CZM
// for fracking simulation

void  declarePoroInterfaceElementModel   ();
void  declareRGPoroInterfaceElementModel ();
void  declarePerforationInterfaceModel   ();

void  declareTwoNodeTrussModel               ();

#endif


