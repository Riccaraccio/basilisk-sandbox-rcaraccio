/**
# Does the frozen-cell gate change an ignition?

The gate of `chemistry.h` decides from the rate at the start of the step
whether the stiff solve can change the state. An ignition is the case that a
rate test can miss: during the induction period the species move slowly, the
radical pool grows by a factor each step, and the mixture then runs away inside
one step. This case measures whether that happens.

It marches a batch reactor at the fixed step of the production case and
compares two trajectories of the same mixture.

* `full` calls `OpenSMOKE_ODESolver` on every step.
* `gated` applies the gate of `chemistry.h` first, exactly as written there,
  and calls the solver only when the gate stays open.

The two must give the same ignition delay. The case sweeps the initial
temperature, so it covers the induction periods from a fraction of a
millisecond to longer than the run.
*/

#define F_ERR 1e-10
#define SOLVE_TEMPERATURE 1

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

/**
The gate of the gas branch, copied from `chemistry.h`. It returns `true` when
it skipped the solve. */

static bool gas_gate (double * y, double dtb, UserDataODE * data) {
  double dy[NGS + 1];
  gas_batch_nonisothermal_constantpressure (y, dtb, dy, data);

  double dYmax = 0.;
  for (int jj = 0; jj < NGS; jj++)
    dYmax = fmax (dYmax, fabs (dy[jj]));

  if (dYmax*dtb < FROZEN_CELL_YTOL && fabs (dy[NGS])*dtb < FROZEN_CELL_TTOL) {
    for (int jj = 0; jj < NGS + 1; jj++)
      y[jj] += dtb*dy[jj];
    return true;
  }
  return false;
}

/**
`rho` and `cp` follow the state, as `update_properties()` makes them follow it
in a real cell. Both trajectories use the same rule, so the comparison is
fair. */

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

/**
March the mixture and return the time at which the temperature first rises by
`RISE` kelvin. `nskip` counts the steps that the gate skipped. */

#define RISE 400.

static double ignition_delay (double T0, double dtb, int nsteps, bool gated,
                              long * nskip, double * Tend, double * Tlog,
                              int nlog) {
  double y[NGS + 1];
  for (int jj = 0; jj < NGS + 1; jj++) y[jj] = 0.;

  /**
  The pyrolysis gas of the biomass burning in air: CO, H2 and the light
  hydrocarbons that the solid releases, at an equivalence ratio near one. */

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

  *nskip = 0;
  double tign = -1.;
  for (int n = 0; n < nsteps; n++) {
    gas_props (y, &data);
    data.T = y[NGS];

    bool skipped = false;
    if (gated)
      skipped = gas_gate (y, dtb, &data);
    if (!skipped)
      OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
                           NGS + 1, dtb, y, &data);
    if (skipped)
      (*nskip)++;

    if (n < nlog && Tlog)
      Tlog[n] = y[NGS];
    if (tign < 0. && y[NGS] > T0 + RISE)
      tign = (n + 1)*dtb;
  }
  *Tend = y[NGS];
  return tign;
}

int main() {
  TS0 = 800.; TG0 = 1123.;
  cpS = 1500.; cpG = 1100.;
  kinfolder = "biomass/Solid-gas-88";
  DT = 2.4e-4;
  init_grid (1);
  run();
}

event init (i = 0) {
  restarted = true;
  foreach() {
    f[] = 1.; porosity[] = 0.2; zeta[] = 0.; p[] = 0.;
    TS[] = TS0; TG[] = 0.; T[] = TS0;
  }
}

event check (i = 0) {

  double dtb = 2.4e-4;
  int nsteps = 3000;             // 0.72 s, longer than any delay below
  int nlog = 3000;

  fprintf (stderr, "# gas branch, dt = %g s, %d steps\n", dtb, nsteps);
  fprintf (stderr, "# T0   tign_full  tign_gated  dt_err[steps]"
           "  Tend_full  Tend_gated  skipped\n");

  double * Ta = malloc (nlog*sizeof(double));
  double * Tb = malloc (nlog*sizeof(double));

  double T0s[] = {800., 900., 950., 1000., 1050., 1100., 1200., 1300.};
  for (int k = 0; k < 8; k++) {
    long na, nb;
    double Ea, Eb;
    double ta = ignition_delay (T0s[k], dtb, nsteps, false, &na, &Ea, Ta, nlog);
    double tb = ignition_delay (T0s[k], dtb, nsteps, true,  &nb, &Eb, Tb, nlog);

    double dTmax = 0.;
    for (int n = 0; n < nlog; n++)
      dTmax = fmax (dTmax, fabs (Ta[n] - Tb[n]));

    fprintf (stderr, "%g %g %g %g %g %g %ld  maxdT %g K\n",
             T0s[k], ta, tb,
             (ta > 0. && tb > 0.) ? (tb - ta)/dtb : 0.,
             Ea, Eb, nb, dTmax);
  }
  free (Ta), free (Tb);
}


/**
## Where the gate of the gas branch closes

The sweep above uses a mixture that burns. This one removes the fuel step by
step and finds the level at which the gate starts to skip the cell. It then
integrates the same mixture for the whole sweep with the solver, so that the
temperature rise says whether that mixture was inert in truth.
*/

event boundary_sweep (i = 0) {

  double dtb = 2.4e-4;
  int nsteps = 3000;

  fprintf (stderr, "\n# gas branch: how much fuel closes the gate?\n");
  fprintf (stderr, "# Yfuel  gate_closes  skipped/%d  dT_over_0.72s[K]"
           "  maxdT_gated_vs_full[K]\n", nsteps);

  double yf[] = {1e-2, 1e-4, 1e-6, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 0.};
  for (int k = 0; k < 9; k++) {

    double y[NGS + 1], ya[NGS + 1];
    UserDataODE data;
    data.P = Pref; data.zeta = 0.; data.sources = NULL;

    for (int jj = 0; jj < NGS + 1; jj++) y[jj] = 0.;
    y[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;
    y[OpenSMOKE_IndexOfSpecies ("CO")] = yf[k];
    y[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765 - yf[k];
    y[NGS] = 1123.;
    for (int jj = 0; jj < NGS + 1; jj++) ya[jj] = y[jj];

    /**
    Does the gate close on the first step? */

    double ycheck[NGS + 1];
    for (int jj = 0; jj < NGS + 1; jj++) ycheck[jj] = y[jj];
    gas_props (ycheck, &data);
    data.T = ycheck[NGS];
    bool closes = gas_gate (ycheck, dtb, &data);

    /**
    March both trajectories over the whole sweep. */

    long nskip = 0;
    double maxdT = 0.;
    for (int n = 0; n < nsteps; n++) {
      gas_props (ya, &data); data.T = ya[NGS];
      OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
                           NGS + 1, dtb, ya, &data);

      gas_props (y, &data); data.T = y[NGS];
      if (gas_gate (y, dtb, &data))
        nskip++;
      else
        OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
                             NGS + 1, dtb, y, &data);

      maxdT = fmax (maxdT, fabs (ya[NGS] - y[NGS]));
    }

    fprintf (stderr, "%g %d %ld %g %g\n",
             yf[k], closes, nskip, ya[NGS] - 1123., maxdT);
  }
}

event stop (i = 1);
