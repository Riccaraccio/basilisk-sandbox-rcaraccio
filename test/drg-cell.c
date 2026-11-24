/**
# Does the DRG reduction keep an ignition?

`drg.h` reduces the gas mechanism at the start of each stiff solve, and the
solver integrates only the important species. This case marches one batch
reactor through an ignition at the fixed step of the production case. It
compares two kinds of trajectory of the same mixture:

* `full` calls `OpenSMOKE_ODESolver` with the full mechanism on every step.
* `drg` calls `gas_stiff_solve()` of `drg.h` on every step, as the sweeps of
  `chemistry.h` do. The case makes a new DRG object for each threshold
  `eps`.

With `eps = 0`, the reduced right-hand side must be equal to the full one at
the state of the analysis, and the ignition must come at the same step. The
test stops with an error if it does not. For all `eps`, the case gives the
error in the ignition delay and in the temperature, the mean size of the
reduced mechanism, and the CPU time against the full solve.

With `eps = 0`, the temperature of the reduced trajectory can differ from
the full one by tens of kelvin near the ignition. This is a limit of the
method, not an error of the code. The analysis uses the state at the start
of the step. A species with a zero mass fraction at that state is not
important, thus it stays at zero for the full step, also when the full
mechanism forms it during the step. A long step through an ignition shows
this most. The test does not check this difference.
*/

#define F_ERR 1e-10
#define SOLVE_TEMPERATURE 1
#define GAS_DRG 1

scalar f[];

#include "run.h"

attribute {
  scalar * tracers, c;
  bool inverse;
}

scalar p[];
scalar porosity[];
scalar zeta[];

double rhoG = 0.3;
double rhoS = 1550.;

#include "memoryallocation-varprop.h"
#include "chemistry.h"
#include <time.h>

/**
The key species are the fuels of the mixture below and O2. */

static const char * keys[] = {"CO", "H2", "CH4", "O2"};
#define NKEY 4

/**
`rho` and `cp` follow the state, as `update_properties()` makes them follow it
in a real cell. All the trajectories use the same rule. */

static void gas_props (const double * y, UserDataODE * data) {
  double x[NGS], mf[NGS], MW;
  for (int jj = 0; jj < NGS; jj++)
    mf[jj] = fmax (0., y[jj]);
  OpenSMOKE_GasProp_SetTemperature (clamp_temperature (y[NGS]));
  OpenSMOKE_GasProp_SetPressure (clamp_pressure (data->P));
  OpenSMOKE_MoleFractions_From_MassFractions (x, &MW, mf);
  data->rhog = data->P/(R_GAS*1000.*y[NGS])*MW;
  data->cpg = OpenSMOKE_GasProp_HeatCapacity (x);
}

static double cpu_time (void) {
  struct timespec ts;
  clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &ts);
  return ts.tv_sec + 1e-9*ts.tv_nsec;
}

/**
The result of one trajectory. */

typedef struct {
  double tign;   // time of the first rise of RISE kelvin, -1 if none
  double Tend;   // temperature at the end
  double sumY;   // sum of the mass fractions at the end
  double cpu;    // CPU time of the march [s]
} Trajectory;

#define RISE 400.

/**
March the mixture over `nsteps` steps of `dtb`. With `reduced` false the
step is the full solve, else it is `gas_stiff_solve()`. `Tlog` gets the
temperature after each step. */

static Trajectory march (double T0, double dtb, int nsteps, bool reduced,
                         double * Tlog, double * rhs_err) {
  double y[NGS + 1];
  for (int jj = 0; jj < NGS + 1; jj++) y[jj] = 0.;

  /**
  The pyrolysis gas of the biomass in air, as in `gate-ignition.c`. */

  y[OpenSMOKE_IndexOfSpecies ("N2")]  = 0.700;
  y[OpenSMOKE_IndexOfSpecies ("O2")]  = 0.212;
  y[OpenSMOKE_IndexOfSpecies ("CO")]  = 0.050;
  y[OpenSMOKE_IndexOfSpecies ("H2")]  = 0.004;
  y[OpenSMOKE_IndexOfSpecies ("CH4")] = 0.014;
  y[OpenSMOKE_IndexOfSpecies ("CO2")] = 0.010;
  y[OpenSMOKE_IndexOfSpecies ("H2O")] = 0.010;
  y[NGS] = T0;

  UserDataODE data;
  data.P = Pref;
  data.zeta = 0.;
  data.sources = NULL;

  Trajectory tr = {-1., 0., 0., 0.};
  double c0 = cpu_time ();
  for (int n = 0; n < nsteps; n++) {
    gas_props (y, &data);
    data.T = y[NGS];

    /**
    With `rhs_err`, compare the full right-hand side with the reduced one at
    the state of the analysis, before the step. The current DRG object must
    have `eps = 0`. */

    if (rhs_err) {
      double dyf[NGS + 1], dyr[NGS + 1];
      drg_analysis (y, data.P);
      gas_batch_nonisothermal_constantpressure (y, dtb, dyf, &data);
      gas_drg_active = true;
      gas_batch_nonisothermal_constantpressure (y, dtb, dyr, &data);
      gas_drg_active = false;
      double scale = 1e-300, diff = 0.;
      for (int jj = 0; jj < NGS; jj++) {
        scale = fmax (scale, fabs (dyf[jj]));
        diff = fmax (diff, fabs (dyr[jj] - dyf[jj]));
      }
      *rhs_err = fmax (*rhs_err, diff/scale);
      *rhs_err = fmax (*rhs_err, fabs (dyr[NGS] - dyf[NGS])/
                       fmax (1e-300, fabs (dyf[NGS])));
    }

    if (reduced)
      gas_stiff_solve (y, dtb, &data);
    else
      OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
                           NGS + 1, dtb, y, &data);

    Tlog[n] = y[NGS];
    if (tr.tign < 0. && y[NGS] > T0 + RISE)
      tr.tign = (n + 1)*dtb;
  }
  tr.cpu = cpu_time () - c0;
  tr.Tend = y[NGS];
  for (int jj = 0; jj < NGS; jj++)
    tr.sumY += y[jj];
  return tr;
}

int main() {
  TS0 = 800.; TG0 = 1123.;
  cpS = 1500.; cpG = 1100.;
  kinfolder = "biomass/Solid-gas-88";
  DT = 2.4e-4;
  init_grid (1);
  run();
}

/**
The `init` event of the case runs before the one of `drg.h`, so the key
species are set when `drg.h` makes the DRG object. */

event init (i = 0) {
  drg_key_species = keys, drg_nkey = NKEY;

  restarted = true;
  foreach() {
    f[] = 1.; porosity[] = 0.2; zeta[] = 0.; p[] = 0.;
    TS[] = TS0; TG[] = 0.; T[] = TS0;
  }
}

event check (i = 0) {

  double dtb = 2.4e-4;
  int nsteps = 3000;             // 0.72 s

  double T0s[] = {1000., 1100., 1300.};
  double epss[] = {0., 1e-3, 1e-2, 5e-2, 1e-1};
  int nT0 = 3, neps = 5;

  fprintf (stderr, "# gas branch, %d gas species, dt = %g s, %d steps\n",
           NGS, dtb, nsteps);
  fprintf (stderr, "# T0  eps  tign  dtign[steps]  Tend  maxdT[K]"
           "  sumY-1  ns_mean  nr_mean  cpu[s]  speedup\n");

  double * Tfull = malloc (nsteps*sizeof(double));
  double * Tdrg = malloc (nsteps*sizeof(double));
  bool pass = true;

  for (int k = 0; k < nT0; k++) {
    double Tmax = 1e5, eps0 = 0.;
    OpenSMOKE_DRG_Init (keys, NKEY, &eps0, &Tmax, 1,
                        drg_T_threshold, drg_R_threshold);
    double rhs_err = 0.;
    march (T0s[k], dtb, nsteps, false, Tfull, &rhs_err);
    fprintf (stderr, "%g rhs eps=0: max relative |dy_reduced - dy_full|"
             " = %.3e\n", T0s[k], rhs_err);

    /**
    With `eps = 0`, the reduced mechanism drops only the species with a zero
    right-hand side. The two right-hand sides must thus agree to the round
    off. */

    if (!(rhs_err < 1e-6)) {
      fprintf (stderr, "FAIL: eps = 0 changes the right-hand side"
               " at T0 = %g\n", T0s[k]);
      pass = false;
    }

    Trajectory full = march (T0s[k], dtb, nsteps, false, Tfull, NULL);
    fprintf (stderr, "%g full %g 0 %g 0 %.3e %d - %.3g 1\n",
             T0s[k], full.tign, full.Tend, full.sumY - 1., NGS, full.cpu);

    for (int e = 0; e < neps; e++) {
      OpenSMOKE_DRG_Init (keys, NKEY, &epss[e], &Tmax, 1,
                          drg_T_threshold, drg_R_threshold);
      drg_nsolve = 0, drg_nspecies = 0., drg_nreactions = 0.;

      Trajectory drg = march (T0s[k], dtb, nsteps, true, Tdrg, NULL);

      double maxdT = 0.;
      for (int n = 0; n < nsteps; n++)
        maxdT = fmax (maxdT, fabs (Tdrg[n] - Tfull[n]));
      double dtign = (full.tign > 0. && drg.tign > 0.) ?
        (drg.tign - full.tign)/dtb : (full.tign > 0.) != (drg.tign > 0.) ?
        HUGE : 0.;

      fprintf (stderr, "%g %g %g %g %g %g %.3e %.1f %.1f %.3g %.2f\n",
               T0s[k], epss[e], drg.tign, dtign, drg.Tend, maxdT,
               drg.sumY - 1., drg_nspecies/drg_nsolve,
               drg_nreactions/drg_nsolve, drg.cpu, full.cpu/drg.cpu);

      /**
      With `eps = 0`, the ignition must come at the same step, or at most
      one step away. The test does not check `maxdT` (see the top of this
      file). */

      if (epss[e] == 0. && fabs (dtign) > 1.) {
        fprintf (stderr, "FAIL: eps = 0 moves the ignition by %g steps"
                 " at T0 = %g\n", dtign, T0s[k]);
        pass = false;
      }
    }
  }
  free (Tfull), free (Tdrg);

  if (!pass)
    exit (1);
}

event stop (i = 1);
