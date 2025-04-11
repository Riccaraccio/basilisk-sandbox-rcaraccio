#define VARPROP

#ifndef F_ERR
# define F_ERR 1e-10
#endif


// IDK WHY THIS DOESN'T WORK
// (const) scalar rhoGv_G = zeroc, rhoGv_S = zeroc, rhoSv = zeroc;
// (const) scalar muGv_G = zeroc, muGv_S = zeroc;
// (const) scalar lambdaGv_G= zeroc, lambdaGv_S = zeroc, lambdaSv = zeroc;
// (const) scalar cpGv_G = zeroc, cpGv_S = zeroc, cpSv = zeroc;

//USE THIS INSTEAD FOR NOW, also seems to be much faster
scalar rhoGv_G[], rhoGv_S[], rhoSv[];
scalar muGv_G[], muGv_S[];
scalar lambdaGv_G[], lambdaGv_S[], lambdaSv[];
scalar cpGv_G[], cpGv_S[], cpSv[];

typedef struct {
  double T, P;
  double * x;
} ThermoState;

typedef struct {
  // Mixture properties
  double (* rhov) (void *);
  double (* muv) (void *);
  double (* lambdav) (void *);
  double (* cpv) (void *);
  // Species properties
  void (* diff) (void *, double *);
  double (* cps)  (void *, int);
} ThermoProps;

#define aavg(f,v1,v2) (clamp(f,0.,1.)*(v1 - v2) + v2)
#define havg(f,v1,v2) (1./(clamp(f,0,1)*(1./(v1) - 1./(v2)) + 1./(v2)))

extern scalar f;
extern face vector alphav;
extern scalar rhov;
#ifdef FILTERED
extern scalar sf;
#else
# define sf f
#endif

event properties (i++) {

  foreach_face() {
    double ff = (sf[] + sf[-1])/2.;

    double rhoGvh_S = 0.5*(rhoGv_S[] + rhoGv_S[-1]);
    double rhoGvh_G = 0.5*(rhoGv_G[] + rhoGv_G[-1]);
    
    alphav.x[] = fm.x[]/aavg (ff, rhoGvh_S, rhoGvh_G);

    {
      face vector muv = mu;
      double mu1vh = 0.5*(muGv_S[] + muGv_S[-1]);
      double mu2vh = 0.5*(muGv_G[] + muGv_G[-1]);

      muv.x[] = fm.x[]*aavg (ff, mu1vh, mu2vh);
    }
  }

  foreach() {
    double rho1vh = rhoGv_S[];
    double rho2vh = rhoGv_G[];

    rhov[] = cm[]*aavg (sf[], rho1vh, rho2vh);
  }
}

/**
## Useful functions

We define functions that are useful for variable properties
simulations.
*/

/**
### *check_termostate()*: check that the thermodynamic state is
reasonable. */

int check_thermostate (ThermoState * ts, int NS) {
  double sum = 0.;
  for (int jj=0; jj<NS; jj++)
    sum += ts->x[jj];

  int T_ok = (ts->T > 180. && ts->T < 4000.) ? true : false;
  int P_ok = (ts->P > 1e3 && ts->P < 1e7) ? true : false;
  int X_ok = (sum > 1.-1.e-3 && sum < 1.+1.e-3) ? true : false;

  return T_ok*P_ok*X_ok;
}

/**
### *print_thermostate()*: print the thermodynamic state of the mixture.
*/

void print_thermostate (ThermoState * ts, int NS, FILE * fp = stdout) {
  fprintf (fp, "Temperature = %g - Pressure = %g\n", ts->T, ts->P);
  for (int jj=0; jj<NS; jj++)
    fprintf (fp, "  Composition[%d] = %g\n", jj, ts->x[jj]);
  fprintf (fp, "\n");
}

/**
### *gasprop_thermal_expansion()*: Thermal expansion coefficient of an ideal gas
*/

double gasprop_thermal_expansion (ThermoState * ts) {
  return ts->T > 0. ? 1./ts->T : 0.;
}