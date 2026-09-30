/**
# Variable properties
This file defines the variable properties for the gas and solid phases, as well as the functions to compute them. 
*/


#define VARPROP

#ifndef F_ERR
# define F_ERR 1e-10
#endif

/**
Variable properties fields. 
Given that we consider posous media, we define separate properties for the gas phase inside the pores and for the solid matrix.
The suffixe "_S" refers to internal gas  and "_G" to the surrounding gas.
 */
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
  double (* rhov)     (void *);
  double (* muv)      (void *);
  double (* lambdav)  (void *);
  double (* cpv)      (void *);
  // Species properties
  void   (* diff)     (void *, double *);
  double (* cps)      (void *, int);
  void   (* cpvs)     (void *, double *);
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

/**
We overwrite the properties for the Navier-Stokes solver with the variable properties.
*/

/**
Caution: `rhov[]` holds `cm[]*rhomix`, so it already carries the metric. Do
not give it metric-aware tree operators (`refine_linear`,
`restriction_volume_average`). They apply the metric a second time, and near
the axis the coarse value becomes non-positive. `viscosity.h` then stops
with SIGFPE at `dt/rho[]` on a coarse level, which a leaf loop never visits.
`two-phase-generic.h` gives the same field the default operators. Keep that
choice.

This event runs in the `adapt` event of `centered.h`, through
`event ("properties")`, after `adapt_wavelet` changes the grid. The chain
calls `update_properties()` first (`qcc -events` shows the order), so the
densities here come from the prolongated state of the new grid.

`update_properties()` sets both densities to 0 and fills each one only where
its gate passes, for example `TG[] > 0`. If a cell gets `rhomix` of 0, the
run stops at `1./rhomix` below, because Basilisk arms the floating point
traps. This event is then not the cause. The measured case was a runaway of
the gas temperature in a sliver cell: see the interface heat source in
`multicomponent-varprop.h`. Fix the cause. Do not add a guard here that
hides it. */

event properties (i++) {

  scalar alphacenter[], mucenter[];
  foreach() {
    double rhomix = rhoGv_G[]*(1.-f[]) + rhoGv_S[]*f[];
    alphacenter[] = 1./rhomix;
    mucenter[] = (muGv_G[]*(1.-f[]) + muGv_S[]*f[]);
    rhov[] = cm[]*rhomix;
  }

  foreach_face() {
    alphav.x[] = fm.x[]*face_value(alphacenter, 0);
    {
      face vector muv = mu;
      muv.x[] = fm.x[]*face_value(mucenter, 0);
    }
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
