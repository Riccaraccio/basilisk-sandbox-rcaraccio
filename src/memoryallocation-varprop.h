#include "opensmoke.h"

unsigned int NGS, NSS;
scalar omega[];

scalar* YGList    = NULL;
scalar* sSexpList   = NULL;
scalar* sGexpList   = NULL;
scalar* YSList      = NULL;
scalar* DmixGList = NULL;
#ifdef MOLAR_DIFFUSION
scalar* XGList    = NULL;
#endif
#ifdef MASS_DIFFUSION_ENTHALPY
scalar* cpGList = NULL;
#endif

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

lambda2v.n[bottom] = neumann(0.); 
lambda2v.t[bottom] = neumann(0.);
lambda2v.n[left] = neumann(0.);
lambda2v.t[left] = neumann(0.);
lambda2v.n[right] = neumann(0.);
lambda2v.t[right] = neumann(0.);
lambda2v.n[top] = neumann(0.);
lambda2v.t[top] = neumann(0.);

lambda1v.n[bottom] = neumann(0.);
lambda1v.t[bottom] = neumann(0.);
lambda1v.n[left] = neumann(0.);
lambda1v.t[left] = neumann(0.);
lambda1v.n[right] = neumann(0.);
lambda1v.t[right] = neumann(0.);
lambda1v.n[top] = neumann(0.);
lambda1v.t[top] = neumann(0.);
#endif

bool success;
scalar MWmixG[];
scalar fG[], fS[];
face vector fsS[], fsG[];

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
  
  YGList      = NULL;
  sSexpList   = NULL;
  sGexpList   = NULL;
  YSList      = NULL;
  DmixGList = NULL;
  #ifdef MOLAR_DIFFUSION
  XGList    = NULL;
  #endif
  #ifdef MASS_DIFFUSION_ENTHALPY
  cpGList = NULL;
  #endif

  //Allocate gas species fields
  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    YGList = list_append (YGList, a);
  } 
  reset (YGList, 0.);

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
    sprintf (name, "D_%s",OpenSMOKE_NamesOfSpecies(jj));
    s.name = strdup (name);
    DmixGList = list_append (DmixGList, s);
  }
  reset (DmixGList, 0.);

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

  #ifdef MOLAR_DIFFUSION
  //Allocate gas molar fraction fields
  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "XG_%s",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    XGList = list_append (XGList, a);
  }
  reset (XGList, 0.);
  #endif
  #ifdef MASS_DIFFUSION_ENTHALPY
  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "cpG_%s",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    cpGList = list_append (cpGList, a);
  }
  reset (cpGList, 0.);
  #endif

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

  for (scalar s in YSList)
    s.inverse = false;

  scalar* temp = list_concat(f.tracers, YSList);
  if (f.tracers)
    free(f.tracers);
  f.tracers = temp;

  fS.nodump = true;
  fG.nodump =true;

#ifdef SOLVE_TEMPERATURE
  TS = new scalar;
  TG = new scalar;

  TS.inverse = false;
  TG.inverse = true;

  sST.nodump = true;
  sGT.nodump = true;

  f.tracers = list_append (f.tracers, TS);
  f.tracers = list_append (f.tracers, TG);
#endif

#ifdef VARPROP
  if (is_constant(rhoGv)) {
    scalar * l = list_copy(all);
    rhoGv = new scalar;
    free(all);
    all = list_concat({rhoGv}, l);
    free(l);
  }
  if (is_constant(rhoSv)) {
    scalar * l = list_copy(all);
    rhoSv = new scalar;
    free(all);
    all = list_concat({rhoSv}, l);
    free(l);
  }
  if (is_constant(muGv)) {
    scalar * l = list_copy(all);
    muGv = new scalar;
    free(all);
    all = list_concat({muGv}, l);
    free(l);
  }
  if (is_constant(lambdaGv)) {
    scalar * l = list_copy(all);
    lambdaGv = new scalar;
    free(all);
    all = list_concat({lambdaGv}, l);
    free(l);
  }
  if (is_constant(lambdaSv)) {
    scalar * l = list_copy(all);
    lambdaSv = new scalar;
    free(all);
    all = list_concat({lambdaSv}, l);
    free(l);
  }
  if (is_constant(cpGv)) {
    scalar * l = list_copy(all);
    cpGv = new scalar;
    free(all);
    all = list_concat({cpGv}, l);
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
      scalar YG = YGList[jj];
      YG[] = gas_start[jj];
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
  delete (YSList),      free (YSList),      YSList = NULL;
  delete (YGList),    free (YGList),    YGList = NULL;
  delete (sSexpList),   free (sSexpList),   sSexpList = NULL;
  delete (sGexpList),   free (sGexpList),   sGexpList = NULL;
  delete (DmixGList), free (DmixGList), DmixGList = NULL;
  free(gas_start),  gas_start = NULL;
  free(sol_start),  sol_start = NULL;
  free(gas_MWs),    gas_MWs = NULL;
  free(sol_MWs),    sol_MWs = NULL;
  delete(f.tracers),  free(f.tracers),   f.tracers = NULL;

#ifdef MOLAR_DIFFUSION
  delete (XGList), free (XGList), XGList = NULL;
#endif

#ifdef MASS_DIFFUSION_ENTHALPY
  delete (cpGList), free (cpGList), cpGList = NULL;
#endif

#ifdef SOLVE_TEMPERATURE
  delete ({TS,TG});
#endif

#ifdef TEMPERATURE_PROFILE
  TemperatureProfile_Free();
#endif
}