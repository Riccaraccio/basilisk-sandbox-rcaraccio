#include "opensmoke.h"

unsigned int NGS, NSS;
scalar omega[];

scalar* YGList_G    = NULL;
scalar* YGList_S    = NULL;
scalar* YGList_Int  = NULL;
scalar* sSexpList   = NULL;
scalar* sGexpList   = NULL;
scalar* YSList      = NULL;
scalar* DmixGList_G = NULL;
scalar* DmixGList_S = NULL;

double* gas_start;
double* sol_start;
double* gas_MWs;
double* sol_MWs;
double Pref = 101325.;

#ifdef SOLVE_TEMPERATURE
scalar T[];
double lambdaS = 1.; double lambdaG = 1.;
double cpS = 1.; double cpG = 1.;

double TS0 = 300.; double TG0 = 300.;

scalar TInt[];
scalar TS, TG;
scalar sST[], sGT[];
face vector lambda1f[], lambda2f[];
vector lambda1v[], lambda2v[];
#endif

bool success;
scalar MWmixG_G[], MWmixG_S[];
scalar fG[], fS[];
face vector fsS[], fsG[];

#ifdef VARPROP //TODO: check if necessary
scalar dummy[];
scalar rhoSv0[], rhoGv0_G[], rhoGv0_S[];
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
  
  YGList_G    = NULL;
  YGList_S    = NULL;
  YGList_Int  = NULL;
  sSexpList   = NULL;
  sGexpList   = NULL;
  YSList      = NULL;
  DmixGList_G = NULL;
  DmixGList_S = NULL;

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
    DmixGList_S = list_append (DmixGList_S, s);
  }
  reset (DmixGList_S, 0.);

  for (int jj = 0; jj<NGS; jj++) {
    scalar s = new scalar;
    free (s.name);
    char name[20];
    sprintf (name, "D_%s_G",OpenSMOKE_NamesOfSpecies(jj));
    s.name = strdup (name);
    DmixGList_G = list_append (DmixGList_G, s);
  }
  reset (DmixGList_G, 0.);

// fields for the source therms
for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    scalar b = new scalar;
    free (a.name);
    free (b.name);
    char aname[20];
    char bname[20];
    sprintf (aname, "sSexp_%s", OpenSMOKE_NamesOfSpecies(jj));
    sprintf (bname, "sGexp_%s", OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (aname);
    b.name = strdup (bname);
    a.nodump = true;
    b.nodump = true;
    sSexpList = list_append (sSexpList, a);
    sGexpList = list_append (sGexpList, b);
  }
  reset (sSexpList, 0.);
  reset (sGexpList, 0.);

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
    s.inverse = false;

  for (scalar s in YGList_G)
    s.inverse = true;

  for (scalar s in YSList)
    s.inverse = false;

  f.tracers = list_concat (f.tracers, YGList_S);
  f.tracers = list_concat (f.tracers, YGList_G);
  f.tracers = list_concat (f.tracers, YSList);

  fS.nodump = true;
  fG.nodump =true;

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

#ifdef SOLVE_TEMPERATURE
  TS = new scalar;
  TG = new scalar;

  TS.inverse = false;
  TG.inverse = true;

  sST.nodump = true;
  sGT.nodump = true;

  f.tracers = list_append (f.tracers, TS);
  f.tracers = list_append (f.tracers, TG);

# if TREE
  TS.refine = refine_linear;
  TS.restriction  = restriction_volume_average;
  TS.dirty = true;

  TG.refine = refine_linear;
  TG.restriction  = restriction_volume_average;
  TG.dirty = true;
# endif
#endif

#ifdef VARPROP
  if (is_constant(rhoGv_G)) {
    scalar * l = list_copy(all);
    rhoGv_G = new scalar;
    free(all);
    all = list_concat({rhoGv_G}, l);
    free(l);
  }
  if (is_constant(rhoGv_S)) {
    scalar * l = list_copy(all);
    rhoGv_S = new scalar;
    free(all);
    all = list_concat({rhoGv_S}, l);
    free(l);
  }
  if (is_constant(rhoSv)) {
    scalar * l = list_copy(all);
    rhoSv = new scalar;
    free(all);
    all = list_concat({rhoSv}, l);
    free(l);
  }
  if (is_constant(muGv_G)) {
    scalar * l = list_copy(all);
    muGv_G = new scalar;
    free(all);
    all = list_concat({muGv_G}, l);
    free(l);
  }
  if (is_constant(muGv_S)) {
    scalar * l = list_copy(all);
    muGv_S = new scalar;
    free(all);
    all = list_concat({muGv_S}, l);
    free(l);
  }
  if (is_constant(lambdaGv_G)) {
    scalar * l = list_copy(all);
    lambdaGv_G = new scalar;
    free(all);
    all = list_concat({lambdaGv_G}, l);
    free(l);
  }
  if (is_constant(lambdaGv_S)) {
    scalar * l = list_copy(all);
    lambdaGv_S = new scalar;
    free(all);
    all = list_concat({lambdaGv_S}, l);
    free(l);
  }
  if (is_constant(lambdaSv)) {
    scalar * l = list_copy(all);
    lambdaSv = new scalar;
    free(all);
    all = list_concat({lambdaSv}, l);
    free(l);
  }
  if (is_constant(cpGv_G)) {
    scalar * l = list_copy(all);
    cpGv_G = new scalar;
    free(all);
    all = list_concat({cpGv_G}, l);
    free(l);
  }
  if (is_constant(cpGv_S)) {
    scalar * l = list_copy(all);
    cpGv_S = new scalar;
    free(all);
    all = list_concat({cpGv_S}, l);
    free(l);
  }
  if (is_constant(cpSv)) {
    scalar * l = list_copy(all);
    cpSv = new scalar;
    free(all);
    all = list_concat({cpSv}, l);
    free(l);
  }
  #endif
}

event init (i = 0) {
  //initialize gas fractions fields
  // check if the sum of the fractions is 1
  double sum = 0.;
  for (int jj=0; jj<NGS; jj++)
    sum += gas_start[jj];
  if (fabs(sum-1.) > 1e-10) {
    fprintf(stderr, "Sum of gas fractions is not 1. Exiting...\n");
    exit(1);
  }

  foreach() {
    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList_G[jj];
      YG[] = gas_start[jj]*(1-f[]);
    }

    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList_S[jj];
      YG[] = gas_start[jj]*f[];
    }
  }

  foreach() {
    for (int jj=0; jj<NGS; jj++) {
      scalar YGInt = YGList_Int[jj];
      YGInt[] = 0.;
    }
  }

  //initialize solid fractions fields
  // check if the sum of the fractions is 1
  sum = 0.;
  for (int jj=0; jj<NSS; jj++)
    sum += sol_start[jj];
  if (fabs(sum-1.) > 1e-6) {
    fprintf(stderr, "Sum of solid fractions is not 1. Exiting...\n");
    exit(1);
  }

  foreach() {
    for (int jj=0; jj<NSS; jj++) {
      scalar YS = YSList[jj];
      YS[] = sol_start[jj]*f[];
    }
  }

  foreach() {
    for (int jj=0; jj<NGS; jj++) {
      scalar sSexp = sSexpList[jj];
      scalar sGexp = sGexpList[jj];
      sSexp[] = 0.;
      sGexp[] = 0.;
    }
  }

#ifdef SOLVE_TEMPERATURE
  foreach() {
    TS[] = TS0*f[];
    TG[] = TG0*(1. - f[]);
    T[]  = TS[] + TG[];
    TInt[] = f[]<1-F_ERR && f[] > F_ERR ? TS0 : 0.;
  }
#endif
}

event cleanup (t = end) {
  //delete(f.tracers),  free(f.tracers),   f.tracers = NULL;
  delete (YSList),      free (YSList),      YSList = NULL;
  delete (YGList_S),    free (YGList_S),    YGList_S = NULL;
  delete (YGList_G),    free (YGList_G),    YGList_G = NULL;
  delete (YGList_Int),  free (YGList_Int),  YGList_Int = NULL;
  delete (sSexpList),   free (sSexpList),   sSexpList = NULL;
  delete (sGexpList),   free (sGexpList),   sGexpList = NULL;
#ifdef VARPROP
  delete (DmixGList_S), free (DmixGList_S), DmixGList_S = NULL;
  delete (DmixGList_G), free (DmixGList_G), DmixGList_G = NULL;
  // delete (CpGList_S),   free (CpGList_S),   CpGList_S = NULL;
  // delete (CpGList_G),   free (CpGList_G),   CpGList_G = NULL;
  // delete (CpSList),     free (CpSList),     CpSList = NULL;
#endif
  free(gas_start),  gas_start = NULL;
  free(sol_start),  sol_start = NULL;
  free(gas_MWs),    gas_MWs = NULL;
  free(sol_MWs),    sol_MWs = NULL;
#ifdef SOLVE_TEMPERATURE
  delete ({TS,TG});
#endif

#ifdef TEMPERATURE_PROFILE
  TemperatureProfile_Free();
#endif
}