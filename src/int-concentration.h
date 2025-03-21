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

#ifndef SOLVE_TEMPERATURE
typedef struct {
  coord c;
} UserDataNls;
#endif

void EqSpecies(const double *xdata, double *fdata, void *params) {
  UserDataNls *data = (UserDataNls *)params;

  double YGInti[NGS];
  // double XGInti[NGS];
  double jG_S[NGS];
  double jG_G[NGS];
  bool success = false;

  Point point = locate(data->c.x, data->c.y, data->c.z);
  for (int jj = 0; jj < NGS; jj++)
    YGInti[jj] = xdata[jj];

  for (int jj = 0; jj < NGS; jj++) {
    scalar YG = YGList_G[jj];
    scalar DmixG = DmixGList_G[jj];
    double gtrgrad = ebmgrad(point, YG, fS, fG, fsS, fsG, true, YGInti[jj], &success);

    double rhoGvh_G;
#ifdef VARPROP
    rhoGvh_G = rhoGv_G[];
#else
    rhoGvh_G = rhoG;
#endif

    jG_G[jj] = rhoGvh_G*DmixG[]*gtrgrad;
  }

  for (int jj = 0; jj < NGS; jj++) {
    scalar YG = YGList_S[jj];
    scalar DmixG = DmixGList_S[jj];
    double gtrgrad = ebmgrad(point, YG, fS, fG, fsS, fsG, false, YGInti[jj], &success);

    double rhoGvh_S;
#ifdef VARPROP
    rhoGvh_S = rhoGv_S[];
#else
    rhoGvh_S = rhoG;
#endif

    jG_S[jj] = rhoGvh_S * DmixG[] * gtrgrad;
  }

  for (int jj = 0; jj < NGS; jj++)
    fdata[jj] = jG_G[jj] + jG_S[jj];
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
          YGInt[] = clamp(YGInt[], 0., 1.);
        }
      }

      array_free (arrUnk);
    }
  }
}
