#ifndef USER_MODELS_H
#define USER_MODELS_H


//-----------------------------------------------------------------------
//   public functions
//-----------------------------------------------------------------------


void  declareUserModels                  ();

void  declareMonitorModel                ();
void  declareLodiModel                   ();
void  declareCoordsStateModel            ();
void  declareNodeConstraintModel         ();
void  declareDofConstraintModel          ();

void  declareArclenModel                 ();
void  declareLoadArclenModel             ();
void  declareDispArclenModel             ();
void  declareEnergyArclenModel           ();

void  declareDirichletModel              ();
void  declareDirichletInterfaceModel     ();

void  declareLineLoadModel               (); 

void  declareDirichletBCModel            ();
void  declareAdaptLoadScaleModel         ();


#endif


