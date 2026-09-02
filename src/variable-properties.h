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
not give it metric-aware tree operators. Commit `f44028d` added

    rhov.refine = refine_linear;
    set_restriction (rhov, restriction_volume_average);

in a `defaults` event, and from that commit no case of this sandbox started
from scratch: every one stopped at step 1 with SIGFPE at `viscosity.h:167`,
which computes `dt/rho[]`. Both operators apply the metric a second time to a
field that already has it. `restriction_volume_average` sums `cm[]*s[]` over
the children and divides by the `cm[]` of the parent. Near the axis `cm` is
`y`, and `y` changes by a factor of 3 between sibling cells, so the coarse
value becomes non-positive. The block below is commented out for this reason.
Do not restore it. `two-phase-generic.h` gives the same field the default
operators, in the same axisymmetric cases. Keep that choice.

The symptom hides on the leaves. This event rewrites every leaf from
`rhoGv_G`, `rhoGv_S` and `f` on each timestep, so a `foreach()` loop finds
`rhov` positive everywhere. The measurement on the failing run gave a minimum
of 9.8e-5 and no cell at or below zero, out of 65536. But `viscosity.h`
solves on a multigrid with `foreach_level_or_leaf`, so it reads the coarse
levels too, and a leaf loop never visits those. Examine the coarse levels
before you conclude that the density is good.

This event carried three guards until 2026-09-03: a fallback that read the
density back from `rhov[]/cm[]`, a test before `1./rhomix`, and a test before
the write of `rhov[]`. They are gone, and they are not needed. `rhomix` stays
positive in every configuration of this sandbox. A guard on `rhov[]` also
does not work: it stops a good value from becoming 0, and it cannot repair a
cell whose stored value is already bad. That was measured on the from-scratch
start, and the crash did not change.

If a later change does make `rhomix` reach 0, the run stops at `1./rhomix`
in the loop below, because Basilisk arms the floating point traps. Fix the
cause in `update_properties()`, which sets both densities to 0 and refills
them only where its gates pass. Do not add a guard here that hides it.

Caution: the signal code of such a SIGFPE reads 7 (invalid), not 3 (divide by
zero). `fsolve-gsl.h` masks the traps around the GSL solve, the invalid
operations there set a sticky flag, and the kernel reads that flag first. Do
not identify the operation from the signal code after the first interface
solve of a run.
*/
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
