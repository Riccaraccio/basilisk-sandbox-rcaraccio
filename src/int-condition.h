#include "intgrad.h"
#ifndef USE_GSL
# define USE_GSL
#endif
#include "fsolve-gsl.h"

//Extern variables
extern scalar fS, fG;
extern face vector fsS, fsG;
extern scalar* Dmix2List_G;
extern scalar* Dmix2List_S;
extern scalar* YGList_S;
extern scalar* YGList_G;
extern scalar* YGList_Int;

#ifndef SOLVE_TEMPERATURE
typedef struct {
  coord c;
} UserDataNls;
#endif

void EqSpecies (const double* xdata, double* fdata, void* params) {
  UserDataNls* data = (UserDataNls*)params;

  double YGInti[NGS];
  // double XGInti[NGS];
  double jG_S[NGS];
  double jG_G[NGS];
  bool success = false;

  foreach_point(data->c.x, data->c.y, data->c.z) {
    // OpenSMOKE_GasProp_SetTemperature (TInt[]);
    // OpenSMOKE_GasProp_SetPressure (p[]+Pref);

    for (int jj=0; jj<NGS; jj++)
      YGInti[jj] = xdata[jj];

    // OpenSMOKE_MoleFractions_From_MassFractions(XGInti, gas_MWs, YGInti);

    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList_G[jj];
      scalar Dmix2 = Dmix2List_G[jj];
      double gtrgrad = ebmgrad (point, YG, fS, fG, fsS, fsG, true, YGInti[jj], &success);
      jG_G[jj] = rhoG*Dmix2[]*gtrgrad; // to be changed to molar diffusion
      //jG_G = rhoG*Dmix2[]*gas_MWs[jj]/MWmixGInt*gtrgrad;
    }

    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList_S[jj];
      scalar Dmix2 = Dmix2List_S[jj];
      double gtrgrad = ebmgrad (point, YG, fS, fG, fsS, fsG, false, YGInti[jj], &success);
      jG_S[jj] = rhoG*Dmix2[]*gtrgrad;
    }

    for (int jj=0; jj<NGS; jj++)
      fdata[jj] = jG_G[jj] + jG_S[jj];
  }
}

int EqSpeciesGsl (const gsl_vector* x, void* params, gsl_vector*f) {
  double* xdata = x->data;
  double* fdata = f->data;

  EqSpecies (xdata, fdata, params);
  return GSL_SUCCESS;
}

void intConcentration () {
  foreach() {
    if (f[]>F_ERR && f[]<1.-F_ERR) {

      Array* arrUnk = array_new();

      for (int jj=0; jj<NGS; jj++) {
        scalar YGInt = YGList_Int[jj];
        double vali = YGInt[];
        array_append (arrUnk, &vali, sizeof(double));
      }

      UserDataNls data;

      coord o = {x,y,z};
      foreach_dimension()
        data.c.x = o.x;

      fsolve (EqSpeciesGsl, arrUnk, &data);

      {
        double* unk = (double*)arrUnk->p;
        for (int jj=0; jj<NGS; jj++) {
          scalar YGInt = YGList_Int[jj];
          YGInt[] = unk[jj];
        }
      }

      array_free (arrUnk);
    }
  }
}
