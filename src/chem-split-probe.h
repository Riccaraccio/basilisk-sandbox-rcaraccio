/**
# Probe TL-1: the split of the gas chemistry from the transport

Item TL-1 of `~/discretization-report/time-level-review.md`. This header
changes nothing in the run. It measures one quantity only, and it writes one
file.

## The question

Each step runs the chemistry first, over the full `dt`, from the state of the
start of the step. The advection and the implicit diffusion run after it. The
splitting is a Lie splitting: no half step and no symmetric form exist. The
review predicts three results of this form:

1. The gap between the maximum temperature after the chemistry and the
   maximum temperature at the end of the step grows with `dt`.
2. The width of the reaction zone grows with `dt`, as `dt` or as the square
   root of `dt`.
3. The middle value of the two maxima depends much less on `dt` than the
   maximum at the end of the step does. The maximum at the end of the step
   gives -38 K per unit of `ln dt`.

Result 3 is the decisive one. If the middle value is nearly flat against
`ln dt`, the splitting carries the dependence of `Tmax` on the timestep, and
a Strang splitting is the fix.

## What the probe measures

`TG` holds the gas temperature in tracer form in this part of the step, so
the intrinsic value is `TG/(1-f)`. The probe reads it in the cells where
`1-f > F_ERR`.

* `Tmax_pre` the largest intrinsic gas temperature before the chemistry.
* `Tmax_chem` the same maximum after the chemistry.
* `Tmax_end` the largest value of `T` at the end of the same step, after the
  diffusion and before `adapt`. It is the quantity of column 3 of
  `OutputData`, `statsf(T).max`. The flame sits in pure gas cells, where `T`
  is the intrinsic gas temperature, so the three maxima are values of one
  quantity.
* `Tmid` the middle value `0.5*(Tmax_chem + Tmax_end)`.
* `dT_max` and `dT_mean` the largest and the mean increment of the intrinsic
  gas temperature over the chemistry. The mean covers the cells that the
  reactor changed, and `nreact` gives their number.
* `nzone` the number of cells whose heat release is above one half of the
  largest heat release. This is the width of the reaction zone.
* `qmax` the largest heat release per unit volume of the cell, and `Qchem`
  its integral over the domain.

The heat release of one cell is

    q = rhoGv_G*cpGv_G*(TG_chem - TG_pre)/dt

per unit volume of the cell, in W/m^3, with `TG` in tracer form. This is the
same quantity that `chemistry.h` puts into `DTDtG`, and `update_divergence()`
gives it the same weight `1-f`. The probe builds it from the temperature
increment, so it does not depend on `GAS_SOURCE_EXACT` or on
`GAS_SOURCE_AVERAGED`. A negative increment gives a negative `q`, and only a
positive `q` counts in `nzone` and in `nreact`.

Caution: `TURN_OFF_HEAT_OF_REACTION` zeroes the temperature source of the
reactor. Every column of this file then holds the transport alone.

## How to read `chemsplit.dat`

    #t(1) i(2) dt(3) Tmax_pre(4) Tmax_chem(5) Tmax_end(6) Tmid(7)
     dT_max(8) dT_mean(9) nreact(10) nzone(11) qmax(12) Qchem(13)

The file holds one line for each output time. All the columns of a line
come from one step: the step that starts at the time of column 1. Columns 4,
5 and 8 to 13 come from its chemistry, and column 6 from its end, before
`adapt`.

Caution: column 3 of `OutputData` on the line of the same time is NOT the
same number. `OutputData` reads `T` at the start of that step, thus at the
end of the step before, after `adapt`.

Plot `Tmax_end`, `Tmid` and `nzone` against the mean `dt` of each run of a
timestep ladder. Read the slope against `ln dt`. A slope of `Tmid` well
below the slope of `Tmax_end` confirms TL-1.

Caution: `dt` moves from step to step, because `dtnext()` snaps the step to
the next event time. Do not read one line as a measurement of `dt`. Compare
runs, and use the mean `dt` of a window of each run.


## Flags

* `CHEM_SPLIT_PROBE` 1 turns the probe on. The default is 0.
* `CHEM_SPLIT_PROBE_FILE` the name of the file.

## Cost

The probe works only on the step that starts at an output time, one step in
about 20 at `DT = 5e-4`. On that step it makes four sweeps of the grid. On
the other steps it tests one flag. The probe adds one field to each cell. */

#ifndef CHEM_SPLIT_PROBE
# define CHEM_SPLIT_PROBE 0
#endif

#if CHEM_SPLIT_PROBE

#include <time.h>

#ifndef CHEM_SPLIT_PROBE_FILE
# define CHEM_SPLIT_PROBE_FILE "chemsplit.dat"
#endif


/**
The intrinsic gas temperature of the start of the step, in tracer form. */

scalar csp_TGpre[];

static bool csp_armed = false, csp_snap = false;
static double csp_cpu = 0.;
static int csp_ncall = 0;
double csp_Tmax_pre = 0., csp_Tmax_chem = 0.;
double csp_dTmax = 0., csp_dTmean = 0.;
double csp_nreact = 0., csp_nzone = 0.;
double csp_qmax = 0., csp_Qchem = 0.;

event defaults (i = 0) {
  csp_TGpre.nodump = true;
}

/**
The arm. An event with a time condition and a new name runs at the start of
the step, before `stability`, so it arms the step that starts at the output
time. The interval is the one of the output events of `run/test.c`. The
event adds no new event time, so it does not change `dt`. */

event chemsplit_arm (t += 0.01) {
  csp_armed = true;
  csp_snap = false;
}

/**
This instance of `reset_sources` runs before the one of `chemistry.h`,
because same-name events run in reverse declaration order. So it runs before
the chemistry. */

event reset_sources (i++) {
  if (!csp_armed)
    return 0;
  clock_t c0 = clock();
  foreach()
    csp_TGpre[] = TG[];
  csp_snap = true;
  csp_cpu += (double)(clock() - c0)/CLOCKS_PER_SEC;
}

/**
This instance of `phasechange` runs after the `chemistry` event and before
the `phasechange` event of `shrinking.h`. `f` does not change over the
chemistry, so the same `1-f` converts both states. */

event phasechange (i++) {
  if (!csp_armed || !csp_snap || !(dt > 0.))
    return 0;
  clock_t c0 = clock();

  double Tpre = -HUGE, Tchem = -HUGE, dTmax = -HUGE;
  double dTsum = 0., nreact = 0., qmax = 0., Qchem = 0.;

  foreach (reduction(max:Tpre) reduction(max:Tchem) reduction(max:dTmax)
           reduction(+:dTsum) reduction(+:nreact)
           reduction(max:qmax) reduction(+:Qchem)) {
    double fG = 1. - f[];
    if (fG > F_ERR) {
      double Tp = csp_TGpre[]/fG, Tc = TG[]/fG, dT = Tc - Tp;
#ifdef VARPROP
      double rc = rhoGv_G[]*cpGv_G[];
#else
      double rc = rhoG*cpG;
#endif
      double q = rc*(TG[] - csp_TGpre[])/dt;

      Tpre = max (Tpre, Tp);
      Tchem = max (Tchem, Tc);
      dTmax = max (dTmax, dT);
      qmax = max (qmax, q);
      Qchem += q*dv();
      if (q > 0.) {
        dTsum += dT;
        nreact += 1.;
      }
    }
  }

  /**
  The count of the reaction zone needs the largest heat release, so it needs
  a second sweep. */

  double nzone = 0.;
  if (qmax > 0.)
    foreach (reduction(+:nzone)) {
      double fG = 1. - f[];
      if (fG > F_ERR) {
#ifdef VARPROP
        double rc = rhoGv_G[]*cpGv_G[];
#else
        double rc = rhoG*cpG;
#endif
        double q = rc*(TG[] - csp_TGpre[])/dt;
        if (q > 0.5*qmax)
          nzone += 1.;
      }
    }

  csp_Tmax_pre = (Tpre > -HUGE) ? Tpre : 0.;
  csp_Tmax_chem = (Tchem > -HUGE) ? Tchem : 0.;
  csp_dTmax = (dTmax > -HUGE) ? dTmax : 0.;
  csp_dTmean = (nreact > 0.) ? dTsum/nreact : 0.;
  csp_nreact = nreact;
  csp_nzone = nzone;
  csp_qmax = qmax;
  csp_Qchem = Qchem;
  csp_cpu += (double)(clock() - c0)/CLOCKS_PER_SEC;
}

/**
The writer. `end_timestep` runs after the diffusion and the projection, and
before `adapt`. So `T` holds the end state of the armed step. */

event end_timestep (i++) {
  if (!csp_armed)
    return 0;
  csp_armed = false;
  if (!csp_snap)
    return 0;
  clock_t c0 = clock();

  double Tend = statsf (T).max;

  if (pid() == 0) {
    static FILE * fp = NULL;
    if (!fp) {
      fp = fopen (CHEM_SPLIT_PROBE_FILE, restarted ? "a" : "w");
      if (fp == NULL) {
        fprintf (stderr, "Error opening %s\n", CHEM_SPLIT_PROBE_FILE);
        exit (1);
      }
      if (!restarted)
        fprintf (fp, "#t(1) i(2) dt(3) Tmax_pre(4) Tmax_chem(5) Tmax_end(6)"
                     " Tmid(7) dT_max(8) dT_mean(9) nreact(10) nzone(11)"
                     " qmax(12) Qchem(13)\n");
    }
    fprintf (fp, "%g %d %g %g %g %g %g %g %g %g %g %g %g\n",
             t, i, dt, csp_Tmax_pre, csp_Tmax_chem, Tend,
             0.5*(csp_Tmax_chem + Tend), csp_dTmax, csp_dTmean,
             csp_nreact, csp_nzone, csp_qmax, csp_Qchem);
    fflush (fp);
  }
  csp_cpu += (double)(clock() - c0)/CLOCKS_PER_SEC;
  csp_ncall++;
}

/**
The cost of the probe, in the log at the end of the run. */

event chemsplit_cost (t = end) {
  if (pid() == 0)
    fprintf (stderr, "# chem-split-probe: %g s CPU in %d steps\n",
             csp_cpu, csp_ncall);
}

#endif // CHEM_SPLIT_PROBE
