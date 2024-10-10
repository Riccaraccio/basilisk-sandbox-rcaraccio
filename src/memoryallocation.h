/*
** TO BE REVISED **
Memory allocation for the OpenSMOKE_Interface
*/

#include "OpenSMOKE_Interface.h"

#pragma autolink -L$OPENSMOKE_INTERFACE/build -lopensmoke

char* kinfolder;

unsigned int NGS, NSS;

scalar* YGList = NULL;
scalar* YSList = NULL;
scalar* Dmix2List = NULL;
scalar T[];
scalar epsilon[];

double* gas_start; 
double* sol_start;
double* gas_MWs;
double* sol_MWs;
double Pref;
double cp1, cp2;
double mass0 = 0;
scalar f0[], fS[], fG[];
face vector fsS[], fsG[];

#ifdef SOLVE_TEMPERATURE
extern double lambda1, lambda2, cp1, cp2;

scalar TInt[];
face vector lambda1f[], lambda2f[];
#endif

event defaults (i = 0) {

  //Read the kinetic scheme
  char kinfolder_root[120];
  sprintf (kinfolder_root, "%s/kinetics/%s/kinetics",
      getenv ("OPENSMOKE_INTERFACE"), kinfolder);

  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder_root);
  OpenSMOKE_ReadSolidKinetics (kinfolder_root);

  NGS = OpenSMOKE_NumberOfSpecies();
  NSS = OpenSMOKE_NumberOfSolidSpecies(); 

  //Allocate gas species fields
  for (int jj = 0; jj<NGS; jj++){
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    YGList = list_append (YGList, a);
  } 
  
  //Allocate solid species fields
  for (int jj = 0; jj<NSS; jj++){
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s",OpenSMOKE_NamesOfSolidSpecies(jj));
    a.name = strdup (name);
    YSList = list_append (YSList, a);
  } 

  //Allocate diff coeff fields, not used for now
  for (int jj = 0; jj<NGS; jj++) {
    scalar s = new scalar;
    free (s.name);
    char name[20];
    sprintf (name, "D_%s",OpenSMOKE_NamesOfSpecies(jj));
    s.name = strdup (name);
    Dmix2List = list_append (Dmix2List, s);
  }  

  //initialize vector with initial values
  gas_start = (double *)malloc(NGS * sizeof(double));
  sol_start = (double *)malloc(NSS * sizeof(double));
  gas_MWs = (double *)malloc(NGS * sizeof(double));
  sol_MWs = (double *)malloc(NSS * sizeof(double));

  for (int jj=0; jj<NGS; jj++){
    gas_start[jj] = 0.;
    gas_MWs[jj] = OpenSMOKE_MW(jj);
  }

  for (int jj=0; jj<NSS; jj++){
    sol_start[jj] = 0.;
    sol_MWs[jj] = OpenSMOKE_MW_Solid(jj);
  }

// For adaptive meshes, not used for now
#if TREE
  for (scalar s in YGList) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
    s.dirty = true; // boundary conditions need to be updated
  }
#endif

#if TREE
  for (scalar s in YSList) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
    s.dirty = true; // boundary conditions need to be updated
  }
#endif

}

event init(i = 0){
  //initialize gas fractions fields
  foreach(){
    for (int jj = 0; jj<NGS; jj++){
      scalar YG = YGList[jj];
      YG[] = gas_start[jj];
    }
  }
  //initialize solid fractions fields
  foreach(){
    for (int jj = 0; jj<NSS; jj++){
      scalar YS = YSList[jj];
      YS[] = sol_start[jj]*f[];
    }
  }
  // calculate initial mass

  foreach(){
    if(f[]>F_ERR){
      for (int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        mass0 += YS[]*rho1*epsilon[]*f[];
      }
    }
  }
}

event cleanup (t = end)
{
  OpenSMOKE_Clean ();

  delete (YSList), free (YSList), YSList = NULL;
  delete (YGList), free (YGList), YGList = NULL;
  delete (Dmix2List), free (Dmix2List), Dmix2List = NULL;
  free(gas_start), gas_start = NULL;
  free(sol_start), sol_start = NULL;
  free(gas_MWs), gas_MWs = NULL;
  free(sol_MWs), sol_MWs = NULL;
  #ifdef SOLVE_TEMPERATURE
    delete ({T});
  #endif
}


// EVAPORATIO  TESTS
// double mEvapVal = 0.1;
scalar mEvap[];
scalar * mEvapList = {mEvap};
double omega = 0;
event phasechange (i++)
{
  // foreach() {
  //   mEvap[] = 0.;
  //   if (f[] > F_ERR && f[] < 1.-F_ERR)
  //     mEvap[] = mEvapVal;
  //   //drhodt[] = -0.1;
  // }

}