#define F_ERR 1e-10

double TG0 = 700.;
int NGS, NSS, maxlevel = 5;
double *sol_MWs, *gas_MWs;

#include "run.h"
#include "opensmoke.h"
#include "reactors.h"

int main() {
  init_grid (1 << maxlevel);
  kinfolder = "biomass/dummy-solid-gas";
  run();
}

event defaults (i = 0) {
  char kinfolder_root[128];
  sprintf (kinfolder_root, "%s/kinetics/%s/kinetics",
      getenv ("OPENSMOKE_INTERFACE"), kinfolder);

  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder_root);
  OpenSMOKE_ReadSolidKinetics (kinfolder_root);
  NGS = OpenSMOKE_NumberOfSpecies(); // Number of gas species
  NSS = OpenSMOKE_NumberOfSolidSpecies(); // Number of solid species

  sol_MWs = (double *)malloc(NSS * sizeof(double));
  gas_MWs = (double *)malloc(NGS * sizeof(double));
  
  for (int jj = 0; jj < NGS; ++jj)
    gas_MWs[jj] = OpenSMOKE_MW(jj);

  for (int jj = 0; jj < NSS; ++jj)
    sol_MWs[jj] = OpenSMOKE_MW_Solid(jj);
}

event init(i = 0) {
  OpenSMOKE_InitODESolver();

  fprintf(stderr, "Number of gas species: %d\n", NGS);
  fprintf(stderr, "Number of solid species: %d\n", NSS);
}

event cleanup(t = end) {
  OpenSMOKE_CleanODESolver();
}