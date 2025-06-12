#include "intgrad.h"
#ifndef USE_GSL
# define USE_GSL
#endif
#include "fsolve-gsl.h"

//Extern variables
extern scalar fS, fG;
extern face vector fsS, fsG;
extern scalar* DmixGList_G;
extern scalar* DmixGList_S;
extern scalar* YGList_S;
extern scalar* YGList_G;
extern scalar* YGList_Int;
extern double* gas_MWs;

#ifdef MOLAR_DIFFUSION
extern scalar* XGList_S;
extern scalar* XGList_G;
#endif


#ifndef SOLVE_TEMPERATURE
typedef struct {
  coord c;
} UserDataNls;
#endif

int EqSpecies(const gsl_vector * xdata, void * params, gsl_vector * fdata) {
  UserDataNls *data = (UserDataNls *)params;

  double YGInti[NGS];
  for (int jj = 0; jj < NGS; jj++)
    YGInti[jj] = gsl_vector_get(xdata, jj);

  double jG_S[NGS];
  double jG_G[NGS];
  bool success = false;

  Point point = locate(data->c.x, data->c.y, data->c.z);
  // foreach_point(data->c.x, data->c.y, data->c.z, serial) {

#ifdef MOLAR_DIFFUSION
  //convert mass to mole fractions
  double XGInti[NGS], MWmixInt;
  OpenSMOKE_MoleFractions_From_MassFractions(XGInti, &MWmixInt, YGInti);
#endif

  for (int jj = 0; jj < NGS; jj++) {
    scalar DmixG = DmixGList_G[jj];
    double rhoGvh_G;
#ifdef VARPROP
    rhoGvh_G = rhoGv_G[];
#else
    rhoGvh_G = rhoG;
#endif

#ifdef MOLAR_DIFFUSION
    scalar XG = XGList_G[jj];
    double gtrgrad = ebmgrad(point, XG, fS, fG, fsS, fsG, true, XGInti[jj], &success);
    jG_G[jj] = rhoGvh_G * DmixG[] * gas_MWs[jj] / MWmixInt * gtrgrad;
#else
    scalar YG = YGList_G[jj];
    double gtrgrad = ebmgrad(point, YG, fS, fG, fsS, fsG, true, YGInti[jj], &success);
    jG_G[jj] = rhoGvh_G * DmixG[] * gtrgrad;
#endif
  }

  for (int jj = 0; jj < NGS; jj++) {
    scalar DmixG = DmixGList_S[jj];
    double rhoGvh_S;
#ifdef VARPROP
    rhoGvh_S = rhoGv_S[];
#else
    rhoGvh_S = rhoG;
#endif

#ifdef MOLAR_DIFFUSION
    scalar XG = XGList_S[jj];
    double gtrgrad = ebmgrad(point, XG, fS, fG, fsS, fsG, false, XGInti[jj], &success);
    jG_S[jj] = rhoGvh_S * DmixG[] * gas_MWs[jj] / MWmixInt * gtrgrad;
#else
    scalar YG = YGList_S[jj];
    double gtrgrad = ebmgrad(point, YG, fS, fG, fsS, fsG, false, YGInti[jj], &success);
    jG_S[jj] = rhoGvh_S * DmixG[] * gtrgrad;
#endif
  }
  // }

  double jGtot_S = 0., jGtot_G = 0.;
#ifdef FICK_CORRECTED
  for (int jj = 0; jj < NGS; jj++) {
    jGtot_S += jG_S[jj];
    jGtot_G += jG_G[jj];
  }
#endif

  for (int jj = 0; jj < NGS; jj++)
    gsl_vector_set(fdata, jj, (jG_S[jj] - jGtot_S*YGInti[jj]) + 
                              (jG_G[jj] - jGtot_G*YGInti[jj]));

  return GSL_SUCCESS;
}

void intConcentration () {
  gsl_vector* unk = gsl_vector_alloc(NGS);
  foreach() {
    if (f[]>F_ERR && f[]<1.-F_ERR) {

      for (int jj=0; jj<NGS; jj++) {
        scalar YGInt = YGList_Int[jj];
        gsl_vector_set(unk, jj, YGInt[]);
      }

      UserDataNls data;

      coord o = {x,y,z};
      foreach_dimension()
        data.c.x = o.x;

      fsolve (EqSpecies, unk, &data, "EqSpecies");

      for (int jj = 0; jj < NGS; jj++) {
        scalar YGInt = YGList_Int[jj];
        YGInt[] = gsl_vector_get(unk, jj);
      }
    }
  }
  gsl_vector_free(unk);
}
