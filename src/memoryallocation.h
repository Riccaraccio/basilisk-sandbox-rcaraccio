#include "opensmoke.h"

unsigned int NGS, NSS;
scalar omega[];

scalar* YGList_G    = NULL;
scalar* YGList_S    = NULL;
scalar* YGList_Int  = NULL;
scalar* YSList      = NULL;
scalar* Dmix2List_G = NULL;
scalar* Dmix2List_S = NULL;

double* gas_start;
double* sol_start;
double* gas_MWs;
double* sol_MWs;
double Pref = 101325.;

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
  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_S",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    YGList_S = list_append (YGList_S, a);
  }
  reset (YGList_S, 0.);

  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_G",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    YGList_G = list_append (YGList_G, a);
  }
  reset (YGList_G, 0.);

  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_Int",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    YGList_Int = list_append (YGList_Int, a);
  }
  reset (YGList_Int, 0.);

  //Allocate solid species fields
  for (int jj = 0; jj<NSS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s",OpenSMOKE_NamesOfSolidSpecies(jj));
    a.name = strdup (name);
    YSList = list_append (YSList, a);
  }
  reset (YSList, 0.);

  //Allocate diff coeff fields
  for (int jj = 0; jj<NGS; jj++) {
    scalar s = new scalar;
    free (s.name);
    char name[20];
    sprintf (name, "D_%s_S",OpenSMOKE_NamesOfSpecies(jj));
    s.name = strdup (name);
    Dmix2List_S = list_append (Dmix2List_S, s);
  }
  reset (Dmix2List_S, 0.);

  for (int jj = 0; jj<NGS; jj++) {
    scalar s = new scalar;
    free (s.name);
    char name[20];
    sprintf (name, "D_%s_G",OpenSMOKE_NamesOfSpecies(jj));
    s.name = strdup (name);
    Dmix2List_G = list_append (Dmix2List_G, s);
  }
  reset (Dmix2List_G, 0.);

  //initialize vector with initial values
  gas_start = (double *)malloc(NGS * sizeof(double));
  sol_start = (double *)malloc(NSS * sizeof(double));
  gas_MWs = (double *)malloc(NGS * sizeof(double));
  sol_MWs = (double *)malloc(NSS * sizeof(double));

  for (int jj=0; jj<NGS; jj++) {
    gas_start[jj] = 0.;
    gas_MWs[jj] = OpenSMOKE_MW(jj);
  }

  for (int jj=0; jj<NSS; jj++) {
    sol_start[jj] = 0.;
    sol_MWs[jj] = OpenSMOKE_MW_Solid(jj);
  }

  for (scalar s in YGList_S)
    s.inverse = true;

  for (scalar s in YGList_G)
    s.inverse = false;

  for (scalar s in YSList)
    s.inverse = false;

  f.tracers = list_concat (f.tracers, YGList_S);
  f.tracers = list_concat (f.tracers, YGList_G);
  f.tracers = list_concat (f.tracers, YSList);

#if TREE
  for (scalar s in YGList_S) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
    s.dirty = true;
  }
#endif

#if TREE
  for (scalar s in YGList_G) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
    s.dirty = true;
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
    s.dirty = true;
  }
#endif

}

event init (i = 0){
  //initialize gas fractions fields
  foreach() {
    for (int jj = 0; jj<NGS; jj++) {
      scalar YG = YGList_S[jj];
      YG[] = gas_start[jj]*(1-f[]);
    }

    for (int jj = 0; jj<NGS; jj++) {
      scalar YG = YGList_G[jj];
      YG[] = gas_start[jj]*f[];
    }
  }

  foreach() {
    for (int jj 0; jj<NGS; jj++) {
      scalar YGInt = YGList_Int[jj];
      YGInt[] = 0.;
    }
  }

  //initialize solid fractions fields
  foreach() {
    for (int jj = 0; jj<NSS; jj++) {
      scalar YS = YSList[jj];
      YS[] = sol_start[jj]*f[];
    }
  }
}

event cleanup (t = end) {
  delete (YSList), free (YSList), YSList = NULL;
  delete (YGList_S), free (YGList_S), YGList_S = NULL;
  delete (Dmix2List_S), free (Dmix2List_S), Dmix2List_S = NULL;
  delete (YGList_G), free (YGList_G), YGList_G = NULL;
  delete (Dmix2List_G), free (Dmix2List_G), Dmix2List_G = NULL;
  delete (YGList_Int), free (YGList_Int), YGList_Int = NULL;
  free(gas_start), gas_start = NULL;
  free(sol_start), sol_start = NULL;
  free(gas_MWs), gas_MWs = NULL;
  free(sol_MWs), sol_MWs = NULL;
}
