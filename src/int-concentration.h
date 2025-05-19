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

int EqSpecies(const gsl_vector * xdata, void * params, gsl_vector * fdata) {
  double YGInti[NGS];
  for (int jj = 0; jj < NGS; jj++)
    YGInti[jj] = gsl_vector_get(xdata, jj);
  
  UserDataNls *data = (UserDataNls *)params;

  // double XGInti[NGS];
  double jG_S[NGS];
  double jG_G[NGS];
  bool success = false;

  //Point point = locate(data->c.x, data->c.y, data->c.z);
  foreach_point(data->c.x, data->c.y, data->c.z, serial) {

    for (int jj = 0; jj < NGS; jj++) {
      scalar YG = YGList_G[jj];
      scalar DmixG = DmixGList_G[jj];
      fprintf (stderr, "YGInti[jj] = %g\n", YGInti[jj]);
      double gtrgrad = ebmgrad(point, YG, fS, fG, fsS, fsG, true, YGInti[jj], &success);

      double rhoGvh_G;
#ifdef VARPROP
      rhoGvh_G = rhoGv_G[];
#else
      rhoGvh_G = rhoG;
#endif

      jG_G[jj] = rhoGvh_G * DmixG[] * gtrgrad;
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
  }

  for (int jj = 0; jj < NGS; jj++)
    gsl_vector_set(fdata, jj, jG_S[jj] + jG_G[jj]);

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
    
      fprintf (stderr, "before fsolve, f=%g\n",f[]);
      fsolve (EqSpecies, unk, &data, "EqSpecies");
      fprintf (stderr, "after fsolve\n");

      for (int jj = 0; jj < NGS; jj++) {
        scalar YGInt = YGList_Int[jj];
        YGInt[] = gsl_vector_get(unk, jj);
      }
    }
  }
  gsl_vector_free(unk);
}
