#include "opensmoke.h"

unsigned int NGS, NSS;
scalar omega[];

scalar* YGList_G    = NULL;
scalar* YGList_S    = NULL;
scalar* YGList_Int  = NULL;
// scalar* XGList_G    = NULL;
// scalar* XGList_S    = NULL;
// scalar* XGList_Int  = NULL;
scalar* sSexpList   = NULL;
scalar* sGexpList   = NULL;
scalar* YSList      = NULL;
scalar* Dmix2List_G = NULL;
scalar* Dmix2List_S = NULL;

double* gas_start;
double* sol_start;
double* gas_MWs;
double* sol_MWs;
double Pref = 101325.;

scalar T[];
double lambdaS = 1.; double lambdaG = 1.;
double cpS = 1.; double cpG = 1.;
bool success;

#ifdef SOLVE_TEMPERATURE
double TS0 = 300.; double TG0 = 300.;

scalar TInt[];
scalar TS, TG;
scalar sST[], sGT[];
face vector lambda1f[], lambda2f[];
#endif

scalar fG[], fS[];
face vector fsS[], fsG[];
scalar fTmp[], fSpc[]; //dummy tracers

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
  // XGList_G    = NULL;
  // XGList_S    = NULL;
  // XGList_Int  = NULL;
  sSexpList     = NULL;
  sGexpList     = NULL;
  YSList      = NULL;
  Dmix2List_G = NULL;
  Dmix2List_S = NULL;

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

  // for (int jj = 0; jj<NGS; jj++) {
  //   scalar a = new scalar;
  //   free (a.name);
  //   char name[20];
  //   sprintf (name, "X_%s_S",OpenSMOKE_NamesOfSpecies(jj));
  //   a.name = strdup (name);
  //   XGList_S = list_append (XGList_S, a);
  // }
  // reset (XGList_S, 0.);

  // for (int jj = 0; jj<NGS; jj++) {
  //   scalar a = new scalar;
  //   free (a.name);
  //   char name[20];
  //   sprintf (name, "X_%s_Int",OpenSMOKE_NamesOfSpecies(jj));
  //   a.name = strdup (name);
  //   XGList_Int = list_append (XGList_Int, a);
  // }
  // reset (XGList_Int, 0.);

  // for (int jj = 0; jj<NGS; jj++) {
  //   scalar a = new scalar;
  //   free (a.name);
  //   char name[20];
  //   sprintf (name, "X_%s_G",OpenSMOKE_NamesOfSpecies(jj));
  //   a.name = strdup (name);
  //   XGList_G = list_append (XGList_G, a);
  // }
  // reset (XGList_G, 0.);

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

  // for (scalar s in XGList_S)
  //   s.inverse = false;

  // for (scalar s in XGList_G)
  //   s.inverse = true;

  fTmp.tracers = NULL;
  fSpc.tracers = NULL;

  // fSpc.tracers = list_concat (fSpc.tracers, YGList_S);
  // fSpc.tracers = list_concat (fSpc.tracers, YGList_G);
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

  // fTmp.tracers = list_append (fTmp.tracers, TS);
  // fTmp.tracers = list_append (fTmp.tracers, TG);
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
}

event init (i = 0){
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
  if (fabs(sum-1.) > 1e-10) {
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

  foreach(){
    fTmp[] = f[];
    fSpc[] = f[];
  }

  scalar* interfaces2 = {fTmp, fSpc};
#if TREE
  for (scalar c in interfaces2) {
    c.refine = c.prolongation = fraction_refine;
    c.dirty = true;
    scalar* tracers = c.tracers;
    for (scalar t in tracers) {
      t.restriction = restriction_volume_average;
      t.refine = t.prolongation = vof_concentration_refine;
      t.dirty = true;
      t.c = c;
    }
  }
#endif
  for (scalar c in interfaces2) {
    scalar* tracers = c.tracers;
    for (scalar t in tracers)
      t.depends = list_add (t.depends, c);
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
  delete (YSList), free (YSList), YSList = NULL;
  delete (YGList_S), free (YGList_S), YGList_S = NULL;
  delete (Dmix2List_S), free (Dmix2List_S), Dmix2List_S = NULL;
  delete (YGList_G), free (YGList_G), YGList_G = NULL;
  delete (Dmix2List_G), free (Dmix2List_G), Dmix2List_G = NULL;
  delete (YGList_Int), free (YGList_Int), YGList_Int = NULL;
  // delete (XGList_G), free (XGList_G), XGList_G = NULL;
  // delete (XGList_S), free (XGList_S), XGList_S = NULL;
  // delete (XGList_Int), free (XGList_Int), XGList_Int = NULL;
  delete (sSexpList), free (sSexpList), sSexpList = NULL;
  delete (sGexpList), free (sGexpList), sGexpList = NULL;
  free(gas_start), gas_start = NULL;
  free(sol_start), sol_start = NULL;
  free(gas_MWs), gas_MWs = NULL;
  free(sol_MWs), sol_MWs = NULL;
#ifdef SOLVE_TEMPERATURE
  delete ({TS,TG});
#endif
  delete (fTmp.tracers), free(fTmp.tracers), fTmp.tracers = NULL;
  delete (fSpc.tracers), free(fSpc.tracers), fSpc.tracers = NULL;

#ifdef TEMPERATURE_PROFILE
  TemperatureProfile_Free();
#endif
}
