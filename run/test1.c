#include "opensmoke.h"

int main () {
  char* kinfolder = "biomass/dummy-solid";
  char kinfolder_root[120];
  sprintf (kinfolder_root, "%s/kinetics/%s/kinetics",
      getenv ("OPENSMOKE_INTERFACE"), kinfolder);

  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder_root);
  fprintf(stderr, "Reading solid kinetics...\n");
  OpenSMOKE_ReadSolidKinetics (kinfolder_root);
  fprintf(stderr, "Loaded successfully\n");

  int NGS = OpenSMOKE_NumberOfSpecies();
  int NSS = OpenSMOKE_NumberOfSolidSpecies();

  return 0;
}
