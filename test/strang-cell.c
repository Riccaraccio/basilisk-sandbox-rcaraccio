/**
# Lie and Strang split of the gas chemistry in one stirred cell

This test supports `GAS_CHEMISTRY_STRANG` in `chemistry.h` (item TL-1 of
`~/discretization-report/time-level-review.md`). It answers two questions:

1. What order in `dt` does each split give, on a problem where the reaction
   and the transport balance each other, as in a flame?
2. What does the second half cost? The stiff solver adapts its internal
   step, so two solves over `dt/2` can cost less than two solves over `dt`.

## The model

One cell of gas, a perfectly stirred reactor. Two operators act on the state
`y = (Y_1 .. Y_NGS, T)`:

* `R(h)`: the gas batch reactor of `reactors.h` over `h`, with the stiff
  solver of OpenSMOKE++. This is the call of `chemistry.h`.
* `M(h)`: the transport, a relaxation to the feed state `y_in` with the
  residence time `tau`. Its exact solution is

      y <- y_in + (y - y_in)*exp(-h/tau)

  It plays the part of the advection and the diffusion, which bring fresh
  fuel and oxygen into the reaction zone.

The feed is pyrolysis gas and air at 1123 K. The cell ignites and reaches a
steady burning state, where the reaction and the feed balance.

The two splits:

    Lie      y <- M(dt) R(dt) y          (the code at 0)
    Strang   y <- R(dt/2) M(dt) R(dt/2) y (the code at 1)

Both read `T` at the end of the step, as `OutputData` does.

The reference is the Strang split at `DT_REF/2`. The test prints the
steady `T` of each split and its error against the reference.

## What the test shows

Column `ratio` is the error at `dt` over the error at `dt/2`. First order
gives 2 and second order gives 4.

At a constant `dt`, the Strang sequence is the Lie sequence with one extra
`R(dt/2)` at the end: `R_h M R_h R_h M R_h = R_h M R M R_h`. So the Strang
state at the end of a step is `R(dt/2)` applied to the Lie state after the
transport. The column `R_h(Lie)` checks this identity. Thus at constant
`dt` the split changes the time at which the state is read, and the
distribution of the expansion over the projections, but not the sequence
itself.

Result on 2026-09-22 (`TAU` 5e-3 s, feed 3 % TAR and 2 % H2O in air, steady
burning state at 1492.17 K, the two references agree to 1.3e-3 K):

    dt        err_Lie  ratio   err_Strang  ratio   R_h(Lie) - T_Strang
    8e-4      -42.79           +1.763              -0.002
    4e-4      -21.22   2.02    +0.504      3.50    -0.000
    2e-4      -10.51   2.02    +0.187      2.70    +0.002
    1e-4       -5.20   2.02    +0.105      1.78    +0.009
    5e-5       -2.55   2.04    +0.097      1.08     0.000
    2.5e-5     -1.22   2.08    +0.082      1.18     0.000

Lie is first order. Strang is second order down to about 0.1 K, where the
tolerance of the stiff solver sets a floor.

## The cost

The test times each call of `R`. `cost ratio` is the time of the Strang step
(two solves over `dt/2`) over the time of the Lie step (one solve over
`dt`), on the same state. The rows `burn` use the steady burning state, the
rows `cold` use air at 1123 K with no fuel, which is the state of most gas
cells of the real case.

Result on 2026-09-22, on a machine with other jobs (so the times are noisy):
`burn` 1.78, 0.40, 1.33, 0.80, 1.06, 1.13 and `cold` 3.61, 1.16, 2.65,
1.58, 2.20, 2.10 for `dt` 8e-4 to 2.5e-5. A cold cell costs about 2 times,
because each call of the Gear solver has a fixed start cost. A burning cell
costs about 1.1 times.

Run with `make strang-cell.tst` from `test/`, or with `qcc` directly. */

#define F_ERR 1e-10
#define SOLVE_TEMPERATURE 1

#include <time.h>

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
double rhoS = 1500.;

/**
See `test/gas-source-cell.c`: these fields replace `variable-properties.h`,
and `VARPROP` makes the reactor recompute `rho` and `cp` at each evaluation,
as in the real case. */

#define VARPROP
scalar rhoGv_G[], rhoGv_S[], rhoSv[];
scalar muGv_G[], muGv_S[];
scalar lambdaGv_G[], lambdaGv_S[], lambdaSv[];
scalar cpGv_G[], cpGv_S[], cpSv[];

#include "memoryallocation-varprop.h"

event timestep (i++) {
  dtnext (DT);
}

#include "reactors.h"

#ifndef TAU
# define TAU 5e-3
#endif

#ifndef TEND_CELL
# define TEND_CELL 0.2
#endif

#ifndef DT_REF
# define DT_REF 2e-6
#endif

#ifndef NREP
# define NREP 1000
#endif

static double yin[16];
static double cpu_R = 0.;
static long ncall_R = 0;

static double now (void) {
  struct timespec ts;
  clock_gettime (CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + 1e-9*ts.tv_nsec;
}

static void opR (double * y, double h) {
  if (!(h > 0.))
    return;
  UserDataODE d;
  d.P = Pref; d.T = y[NGS]; d.sources = NULL;
  d.rhog = 0.3; d.cpg = 1200.;
  double t0 = now();
  OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
                       NGS + 1, h, y, &d);
  cpu_R += now() - t0;
  ncall_R++;
  for (int jj = 0; jj < NGS; jj++)
    y[jj] = y[jj] > 0. ? y[jj] : 0.;
}

static void opM (double * y, double h) {
  double e = exp (-h/TAU);
  for (int jj = 0; jj < NGS + 1; jj++)
    y[jj] = yin[jj] + (y[jj] - yin[jj])*e;
}

/**
March one split over `TEND_CELL` from the state `y`, and return the
temperature at the end of the last step. */

static double march (int strang, double dt, double * y) {
  int n = (int) (TEND_CELL/dt + 0.5);
  for (int k = 0; k < n; k++) {
    if (strang) {
      opR (y, 0.5*dt);
      opM (y, dt);
      opR (y, 0.5*dt);
    }
    else {
      opR (y, dt);
      opM (y, dt);
    }
  }
  return y[NGS];
}

int main() {
  TG0 = 1123.; TS0 = 300.;
  cpS = 1500.; cpG = 1100.;
  kinfolder = "biomass/dummy-solid-gas";
  DT = 1e-4;
  init_grid (1);
  run();
}

event init (i = 0) {
  restarted = true;
  OpenSMOKE_InitODESolver ();
  assert (NGS + 1 <= 16);

  for (int jj = 0; jj < NGS + 1; jj++)
    yin[jj] = 0.;
  int iTAR = OpenSMOKE_IndexOfSpecies ("TAR");
  int iH2O = OpenSMOKE_IndexOfSpecies ("H2O");
  int iO2  = OpenSMOKE_IndexOfSpecies ("O2");
  int iN2  = OpenSMOKE_IndexOfSpecies ("N2");
  yin[iTAR] = 0.03;
  yin[iH2O] = 0.02;
  yin[iO2]  = 0.235*0.95;
  yin[iN2]  = 1. - yin[iTAR] - yin[iH2O] - yin[iO2];
  yin[NGS]  = 1123.;

  /**
  A burning start: march the feed with a fine Strang split to the steady
  state, so each run of the ladder starts on the same burning branch. The
  reference then marches the same time again, at `DT_REF` and at
  `DT_REF/2`. The difference of the two references gives the accuracy of
  the reference, and the drift from the start gives the distance to the
  steady state. */

  double yb[NGS + 1];
  for (int jj = 0; jj < NGS + 1; jj++)
    yb[jj] = yin[jj];
  yb[NGS] = 2000.;
  double Tstart = march (1, DT_REF, yb);

  double yref[NGS + 1];
  for (int jj = 0; jj < NGS + 1; jj++)
    yref[jj] = yb[jj];
  double Tref = march (1, DT_REF, yref);
  for (int jj = 0; jj < NGS + 1; jj++)
    yref[jj] = yb[jj];
  double Tref2 = march (1, 0.5*DT_REF, yref);
  fprintf (stderr, "# tau = %g s, tend = %g s, Tin = %g K\n"
           "# reference Strang: T(dt_ref) = %.6f K, T(dt_ref/2) = %.6f K,"
           " drift from the start %.3g K, dt_ref = %g s\n",
           (double) TAU, (double) TEND_CELL, yin[NGS], Tref, Tref2,
           Tref - Tstart, (double) DT_REF);
  Tref = Tref2;
  fprintf (stderr, "# %9s %12s %10s %6s %12s %10s %6s %12s\n",
           "dt", "T_Lie", "err_Lie", "ratio", "T_Strang", "err_Str",
           "ratio", "R_h(Lie)");

  double dts[] = {8e-4, 4e-4, 2e-4, 1e-4, 5e-5, 2.5e-5};
  double eL0 = 0., eS0 = 0.;
  for (int a = 0; a < 6; a++) {
    double dt = dts[a];
    double yl[NGS + 1], ys[NGS + 1];
    for (int jj = 0; jj < NGS + 1; jj++)
      yl[jj] = ys[jj] = yb[jj];
    double TL = march (0, dt, yl);
    double TS_ = march (1, dt, ys);
    opR (yl, 0.5*dt); // R(dt/2) of the Lie state after the transport
    double eL = TL - Tref, eS = TS_ - Tref;
    fprintf (stderr, "  %9.2e %12.4f %10.4f %6.2f %12.4f %10.4f %6.2f %12.4f\n",
             dt, TL, eL, a ? eL0/eL : 0., TS_, eS, a ? eS0/eS : 0., yl[NGS]);
    eL0 = eL, eS0 = eS;
  }

  /**
  The cost of one step on a fixed state: one solve over `dt` against two
  solves over `dt/2`. Each measurement repeats the step on a copy of the
  state and takes the mean. */

  fprintf (stderr, "\n# cost of one step: %8s %12s %12s %8s\n",
           "dt", "t_R(dt) [s]", "t_2R(dt/2)", "ratio");
  double yc[NGS + 1];
  for (int jj = 0; jj < NGS + 1; jj++)
    yc[jj] = yin[jj];
  yc[iTAR] = 0.; yc[iH2O] = 0.; yc[iO2] = 0.235; yc[iN2] = 0.765;
  for (int s = 0; s < 2; s++) {
    double * ystate = s == 0 ? yb : yc;
    for (int a = 0; a < 6; a++) {
      double dt = dts[a];
      int nrep = NREP;
      double y[NGS + 1];

      cpu_R = 0.;
      for (int r = 0; r < nrep; r++) {
        for (int jj = 0; jj < NGS + 1; jj++) y[jj] = ystate[jj];
        opR (y, dt);
      }
      double t1 = cpu_R/nrep;

      cpu_R = 0.;
      for (int r = 0; r < nrep; r++) {
        for (int jj = 0; jj < NGS + 1; jj++) y[jj] = ystate[jj];
        opR (y, 0.5*dt);
        opR (y, 0.5*dt);
      }
      double t2 = cpu_R/nrep;
      fprintf (stderr, "  %-6s %8.2e %12.3e %12.3e %8.2f\n",
               s == 0 ? "burn" : "cold", dt, t1, t2, t2/t1);
    }
  }

  exit (0);
}
