#include "intgrad.h"
#ifndef USE_GSL
# define USE_GSL
#endif
#include "fsolve-gsl.h"

//Extern variables
extern scalar fS, fG;
extern face vector fsS, fsG;
extern scalar porosity;
extern vector lambda1v, lambda2v;
extern double TG0;
extern scalar TInt, TS, TG, T;

typedef struct {
  coord c;
} UserDataNls;

#ifndef RADIATION_INTERFACE
# define RADIATION_INTERFACE 0.
#endif

double divq_rad_int (double TInti, double Tbulk = 300., double alphacorr = 1.) {
  return alphacorr*5.670373e-8*(pow(Tbulk, 4.) - pow(TInti, 4.));
}

void EqTemperature (const double* xdata, double* fdata, void* params) {
  UserDataNls* data = (UserDataNls*) params;

  foreach_point (data->c.x, data->c.y, data->c.z, serial) {
    double TInti = xdata[0];
    bool success = false;

    double gradTGn = ebmgrad (point, TG, fS, fG, fsS, fsG, true, TInti, &success);
    double gradTSn = ebmgrad (point, TS, fS, fG, fsS, fsG, false, TInti, &success);

    coord n = facet_normal (point, fS, fsS);
    double lambda1vh = lambda1v.x[]*n.x + lambda1v.y[]*n.y;
  
    //lambdaG is not dipendent on the direction, this is just to treat it the same way as lambdaS
    double lambda2vh = lambda2v.x[]*n.x + lambda2v.y[]*n.y;

    //Interface energy balance
    fdata[0] = - divq_rad_int (TInti, TG0, RADIATION_INTERFACE)
                     + lambda1vh*gradTSn
                     + lambda2vh*gradTGn;
  }
}

int EqTemperatureGsl (const gsl_vector* x, void* params, gsl_vector* f) {
  double* xdata = x->data;
  double* fdata = f->data;

  EqTemperature (xdata, fdata, params);
  return GSL_SUCCESS;
}

void ijc_CoupledTemperature() {
  foreach() {
    if (f[]>F_ERR && f[] < 1.-F_ERR) {
      Array * arrUnk = array_new();
      {
        double vali = TInt[];
        array_append (arrUnk, &vali, sizeof(double));
      }

      UserDataNls data;
      coord o = {x,y,z};
      foreach_dimension()
        data.c.x = o.x;

      fsolve (EqTemperatureGsl, arrUnk, &data);

      {
        double* unk  = (double*)arrUnk->p;
        TInt[] = unk[0];
      }
      array_free (arrUnk);
    }
  }
}
