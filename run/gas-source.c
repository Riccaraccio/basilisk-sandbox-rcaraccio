/**
# Where does the flame-sheet jitter of fatehi-combustion come from?

The temperature oscillation of `fatehi-combustion.c` is not an oscillation of
the flame temperature. The peak temperature is steady to 3 K. The flame *sheet*
moves by about one cell, and a probe on a fixed line converts that displacement
into a 33 K swing, because the flame gradient is 42 K per cell.

The signature that separates the two is

~~~
ratio = hf_std(T on a fixed line) / hf_std(Tmax)
~~~

* below 1: the flame flickers in amplitude and stays put
* above 1: the flame sheet moves

The cluster campaign measured, at matched sampling and matched filter:

| run              | zeta policy | Tmax_hf | Tline_hf | ratio |
|------------------|-------------|---------|----------|-------|
| `expansion`      | SWELLING    | 4.21    | 1.38     | 0.33  |
| `combustion`     | SWELLING    | 0.43    | 0.17     | 0.39  |
| `gravity`        | SWELLING    | 1.69    | 0.51     | 0.30  |
| `var-prop-const` | SWELLING    | 1.64    | 0.49     | 0.30  |
| `com-shrink`     | SHRINK      | 2.17    | 2.16     | 0.99  |
| **real fatehi**  | REACTION    | 3.05    | 7.23     | 2.37  |

Every case that holds the interface still sits near 0.3. Only the real
configuration displaces the sheet. So a reduced case with `ZETA_SWELLING`
cannot be used to study this: it deletes the mechanism. This case therefore
keeps the full physics of `fatehi-combustion.c` and puts one suspect behind
each macro.

**Acceptance test: the reference case must reproduce a ratio near 2.4.** If it
does not, it is not reproducing the problem, and nothing measured on it means
anything.

## The cases

The source is one file, `gas-source.c`. Each case is a symbolic link to it that
the `Makefile` creates on demand and gives its own `CFLAGS`. That is the usual
Basilisk pattern (see `collapse-inviscid.c` in the Basilisk test suite), except
that here the links are generated rather than committed: `Makefile.defs`
declares `%.c` precious, so they survive, and `.gitignore` keeps them out of the
repository. `make gas-source-unlink` removes them.

Each case runs in its own directory, so the output names do not collide.

| case                | knob                        | hypothesis under test          |
|---------------------|-----------------------------|--------------------------------|
| `gas-source-ignite` | `TEND=12`, `IGNITION=1`     | builds the shared start state  |
| `gas-source-inst`   | none: the reference         | current behaviour              |
| `gas-source-avg`    | `GAS_SOURCE_AVERAGED=1`     | form of the gas expansion source |
| `gas-source-zmix`   | `FLAME_PRINT_TIME=0.005`    | the stale adapt criterion      |
| `gas-source-zeta`   | `ZETA_POLICY=ZETA_CONST`    | the global maximum in `set_zeta` |
| `gas-source-dt`     | `CFLNUM=0.125`              | the operator split             |
| `gas-source-exact`  | `GAS_SOURCE_EXACT=1`        | exact expansion + source filter |
| `gas-source-exact0` | `GAS_SOURCE_EXACT=1`, `GAS_SOURCE_FILTER_PASSES=0` | exact expansion, no filter |
| `gas-source-dtfix`  | `DTMAX=1.75e-4`, `CFLNUM=2` | the step-size loop             |
| `gas-source-dtvar`  | `CFLNUM=0.347`              | the control for `dtfix`        |

### The hypotheses

1. `GAS_SOURCE_AVERAGED` selects the form of the gas-phase reaction source in
   `chemistry.h`. With `0` the source is the instantaneous rate at the
   converged end-of-step state, which is exponential in `T`. It reaches the
   velocity field through `drhodt` and the pressure Poisson equation, so it
   moves the flame, which changes the rate again. With `1` the source is the
   step-averaged rate, `(state_end - state_start)/dt`.

2. `flame.h` writes `zmix` only inside `event flame (t += FLAME_PRINT_TIME)`,
   which defaults to 0.1 s. The `adapt` event below reads `zmix - zsto` on
   **every** timestep. For 0.1 s the mesh follows a frozen flame position, then
   it snaps. `run/glarborg-particle.c` already sets 0.01 for this reason.

3. `ZETA_REACTION` divides `omega[]` by a global maximum (`shrinking.h`), so
   noise in the peak reaction rate rescales the interface velocity everywhere
   at once. That fits the measured correlation of 0.87 to 0.93 between stations
   9 mm apart. `ZETA_CONST` removes the global coupling.

4. The case leaves `DT = 1`, so the advective CFL is what limits `dt`, and
   `dtnext` then re-quantises it so an integer number of steps lands on each
   event time. A first-order split at the flame depends on `dt`. Quartering
   `CFL` quarters `dt`; a cap through `DTMAX` would do nothing unless it fell
   below the CFL-limited step, which is why this case uses `CFLNUM`.

5. `gas-source-sdiag` measured a period-2 oscillation in `Qrho` with
   `corr(dQrho, d(dt)) = +0.84`. With `GAS_SOURCE_AVERAGED=1` the source is
   `rho*(y_end - y_start)/dt`. When the step consumes the fuel of the cell,
   the numerator stops to grow and the source becomes proportional to `1/dt`.
   The `dt` controller then reads a velocity that this source produced, so the
   loop closes: source, divergence, velocity, CFL, `dt`, source. A loop with a
   lag of one step oscillates with a period of two steps.

   `gas-source-dtfix` opens the loop. It puts `DTMAX` below the CFL-limited
   step, so every step has the same size. `gas-source-dtvar` is the control:
   it runs at the same mean step, but the CFL still sets it, so the step still
   varies. A period-2 mode that disappears in `dtfix` and stays in `dtvar`
   proves the loop. A mode that disappears in both only proves that a smaller
   step helps, which hypothesis 4 already covers.

   Caution: verify that column 2 of `probe.dat` is constant in `dtfix`. If
   `DTMAX` does not bind on every step, the test says nothing. Lower
   `GSAB_DTFIX` in the Makefile and run it again.

   Both cases use `GAS_SOURCE_AVERAGED=1`, because the measurement that
   started this comes from that build. Do not compare them to
   `gas-source-inst`.

## How to run

Ignition takes about 10 s of simulated time. Run it once, then branch every
case from the same state, so each case differs by one variable and by nothing
else.

~~~bash
cd run
make gas-source-ignite.tst                 # to t = 12 s, writes ignition-snapshot

for c in inst avg zmix zeta dt; do
    mkdir -p gas-source-$c
    cp gas-source-ignite/ignition-snapshot gas-source-$c/
done

make gas-source-inst.tst gas-source-avg.tst gas-source-zmix.tst \
     gas-source-zeta.tst gas-source-dt.tst

python3 ../misc/gas-source.py \
    gas-source-inst/probe.dat gas-source-avg/probe.dat  \
    gas-source-zmix/probe.dat gas-source-zeta/probe.dat \
    gas-source-dt/probe.dat --t0 13
~~~

A shell glob works too, but do not write one inside a C comment: the `*` and
the `/` of a path such as `gas-source-<star>/probe.dat` close the comment.

On restart the case prefers its own `last-snapshot`, and falls back to the
shared `ignition-snapshot`. A branch never overwrites `ignition-snapshot`, so
re-running a case always restarts it from the same state.

## Compile-time parameters

`GAS_SOURCE_AVERAGED`, `ZETA_POLICY`, `ADAPT_ZDIFF`, `DTMAX`, `CFLNUM`,
`MAXLEVEL`, `TEND`, `TCUT`, `IGNITION`. `FLAME_PRINT_TIME` belongs to
`flame.h` and is overridable in the same way.

## What the probe records

Every diagnostic runs on `i++`. Sampling at `t += 0.1`, as the real case does,
puts the Nyquist limit at 5 Hz, and an 8 mm buoyant flame flickers near 17 Hz,
so everything above 5 Hz folds into the low-frequency band. That aliasing is
why the real run appears to wander at about 1 Hz.

`r_flame` comes from reductions only, so no lattice of sample points jitters
against the grid and no `nodata` can appear. `Tp1` and `Tp2` sit on the
x = H0/2 + 2 mm line, at y = 4 mm and y = 9.6 mm, which is where the real run
shows the largest jitter. `Qrho` measures the gas expansion source directly:
watch its sign, because a reacting gas that expands must give
`Qsrc + Qrho < 0`.
*/

#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define SOLID_SOURCE_DIAG 1
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1

#ifndef MAXLEVEL
# define MAXLEVEL 10
#endif
#ifndef DTMAX
# define DTMAX 1.
#endif
#ifndef TEND
# define TEND 20.
#endif
#ifndef ZETA_POLICY
# define ZETA_POLICY ZETA_REACTION
#endif
#ifndef ADAPT_ZDIFF
# define ADAPT_ZDIFF 1
#endif
#ifndef IGNITION
# define IGNITION 0
#endif

/**
`DTMAX` caps the timestep, `CFLNUM` scales it. The real case leaves `DT = 1`,
so the advective CFL is what actually limits `dt`. A cap above the CFL-limited
step therefore changes nothing, which is why the timestep hypothesis is tested
with `CFLNUM` and not with `DTMAX`.

Caution: the `init` event applies `CFLNUM`, not `main()`. See the comment
there. To confirm that the value took effect, read the `CFL=` field in the
probe.dat header. Do not read the `CFLNUM=` field in the log. The log prints
the macro, which is correct even when the solver ignores it. */

#ifndef CFLNUM
# define CFLNUM 0.5
#endif

/**
`qcc` turns every event into `int action (const int i, const double t, Event *)`,
so the name `t` inside an event body is the parameter, not the global. The
parameter holds the time before the event runs. `restore()` writes the global,
so an event that calls `restore()` cannot see the new time through `t`. This
function is at file scope, so it reads the global.

Caution: do not replace `global_time()` with `t` inside an event. The value
becomes the time before the restore, which is 0 on the first step. */

static double global_time (void) { return t; }

/**
`TCUT` must sit between the ambient 1123 K and the peak flame temperature, on
the steep side of the profile, so that the measured radius follows the flame
sheet.

Caution: check that `Vhot` (column 6) is not zero over the analysis window. A
zero means the flame never crossed `TCUT`, and then `r_flame` carries no
information. The default suits the level 10 flame, which reaches 2330 K. */

#ifndef TCUT
# define TCUT 1800.
#endif

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "gravity.h"
#include "superquadric.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
#include "flame.h"

const double Uin = 0.13; //inlet velocity
u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

u.n[right]   = neumann (0.);
u.t[right]   = neumann (0.);
p[right]     = dirichlet (0.);
psi[right]   = neumann (0.);

int maxlevel = MAXLEVEL, minlevel = 2;
double D0 = 8e-3, H0 = 8e-3;
double solid_mass0 = 0.;

int main() {

  /**
  Caution: under MPI every rank shares this stderr. Guard every message with
  `pid() == 0`, or the log carries one copy per rank. */

  if (pid() == 0)
    fprintf (stderr, "# GAS_SOURCE_AVERAGED=%d GAS_SOURCE_RHO_MEAN=%d"
                     " GAS_SOURCE_EXACT=%d FILTER_PASSES=%d"
                     " ZETA_POLICY=%d ADAPT_ZDIFF=%d"
                     " FLAME_PRINT_TIME=%g DTMAX=%g CFLNUM=%g MAXLEVEL=%d"
                     " TEND=%g TCUT=%g IGNITION=%d nranks=%d\n",
             (int) gas_source_averaged, (int) gas_source_rho_mean,
             (int) GAS_SOURCE_EXACT,
#if GAS_SOURCE_EXACT
             gas_source_filter_passes,
#else
             0,
#endif
             (int) (ZETA_POLICY), ADAPT_ZDIFF,
             (double) FLAME_PRINT_TIME, (double) DTMAX, (double) CFLNUM,
             MAXLEVEL, (double) TEND, (double) TCUT, IGNITION, npe());

  lambdaSmodel = L_TENWOLDE;
  TS0 = 300.; TG0 = 1123.;
  rhoS = 1550;
  eps0 = 0.2; // low, compressed pellet

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_POLICY;

  DT = DTMAX;

  G.x = -9.81;

  kinfolder = "biomass/dummy-solid-gas";
  shift_prod = true;

  L0 = 20*D0;
  origin (-L0/2, 0);
  init_grid (1 << maxlevel);

  emissivity = emissivity_diblasi;

  run();
}

event init (i = 0) {
  /**
  Caution: do not set `CFL` in `main()`. The `defaults` event of
  `navier-stokes/centered.h` sets `CFL = 0.8`, and that event runs after
  `main()`. The `stability` event of `vof.h` then clamps `CFL` down to 0.5.
  An `init` event runs after every `defaults` event, so the value survives.
  `vof.h` only lowers `CFL`, so a value below 0.5 stays. */

  CFL = CFLNUM;

  scalar f0[];
  fraction (f0, superquadric (x, y, 20, 0.5*H0, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 0.935;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")]   = 0.061;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]     = 0.004;

  foreach()
    porosity[] = eps0*f0[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += f0[]*(1. - eps0)*rhoS*dv();

  TG[left]   = dirichlet (TG0);
  TG[top]    = dirichlet (TG0);
  TG[right]  = neumann (0.);
  TG[bottom] = neumann (0.);

  for (int jj = 0; jj < NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[left] = dirichlet (0.765);
      YG[top]  = dirichlet (0.765);
    } else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      YG[left] = dirichlet (0.235);
      YG[top]  = dirichlet (0.235);
    } else {
      YG[left] = dirichlet (0.);
      YG[top]  = dirichlet (0.);
    }
  }

  /**
  Prefer this case's own checkpoint, so an interrupted run resumes. Otherwise
  branch from the shared state that the ignition case wrote. A branch never
  writes `ignition-snapshot`, so a re-run always restarts from the same
  state. */

  if (restore (file = "last-snapshot", list = all)) {
    if (pid() == 0)
      fprintf (stderr, "# resumed from last-snapshot at t = %g\n",
               global_time());
    restarted = true;
  }
  else if (restore (file = "ignition-snapshot", list = all)) {
    if (pid() == 0)
      fprintf (stderr, "# branched from ignition-snapshot at t = %g\n",
               global_time());
    restarted = true;

    /**
    `restore()` reinstates `t` from the dump header. A branch that starts at
    t = 0 therefore holds no ignition state: it would run the whole ignition
    again, spend the allocation, and leave a probe.dat that answers a
    different question. Stop instead.

    A snapshot written at t = 0 usually means the ignition case was built
    without `-DIGNITION=1 -DTEND=12`, so it never reached the dump. */

    if (global_time() <= 0.) {
      if (pid() == 0)
        fprintf (stderr,
                 "ERROR: ignition-snapshot carries t = %g, so it is not an"
                 " ignition state.\n"
                 "       Rebuild the ignition case with -DIGNITION=1"
                 " -DTEND=12 and run it again.\n", global_time());
      exit (1);
    }
  }
  else {
    if (pid() == 0)
      fprintf (stderr, "# no snapshot found, starting from scratch\n");
    foreach() {
      f[] = f0[];
      porosity[] = eps0*f[];
    }
  }
}

/**
## Fixed probes

Caution: `interpolate` returns `nodata` (1e30) when no rank locates the point.
That happened twice in the real run and produced a 190 K spike, because
`print_profile` wrote the value raw and `T_H2O_weigthed_average` had no guard
either. Guard it. */

static double probe (scalar s, double xp, double yp) {
  double v = interpolate (s, xp, yp);
  return (v > 1e20) ? nodata : v;
}

/**
## Per-step diagnostics

Everything below is a reduction or a collective `interpolate`, so every rank
holds the same values, and only rank 0 writes. */

event probe_flame (i++, last) {

  /**
  Mixture temperature. `TS` and `TG` are in tracer form at the end of the step,
  so their sum is the volume-weighted mixture temperature, the same definition
  the solver gives to `T`. Build it here so the value is unambiguous at this
  point of the timestep. */

  scalar Tm[];
  foreach()
    Tm[] = TS[] + TG[];

  /**
  Divergence budget. The projection enforces

    div(uf) = -(gas_source + drhodt)

  and `drhodt` exists only under VARPROP. A probe that watches `gas_source`
  alone reports a large residual on this case for no reason. Watch both. */

  double Qsrc = 0., Qrho = 0., Qdiv = 0., resmax = 0.;
  foreach (reduction(+:Qsrc) reduction(+:Qrho) reduction(+:Qdiv)
           reduction(max:resmax)) {
    double d = 0.;
    foreach_dimension()
      d += uf.x[1] - uf.x[];
    d /= Delta;                    // = cm*div(u), same weighting as gas_source

    Qsrc += gas_source[]*sq(Delta);
    Qrho += drhodt[]*sq(Delta);
    Qdiv += d*sq(Delta);
    resmax = max (resmax, fabs (d + gas_source[] + drhodt[]));
  }

  /**
  Flame position and size, from reductions only.

    Vhot     volume inside the T = TCUT surface
    r_flame  radius of the flame sheet, weighted by (T - TCUT)
    x_flame  axial position of the same weighted centroid

  A jitter of one cell in the flame position shows up directly in `r_flame`. */

  double Vhot = 0., wsum = 0., rw = 0., xw = 0.;
  foreach (reduction(+:Vhot) reduction(+:wsum) reduction(+:rw) reduction(+:xw)) {
    if (Tm[] > TCUT) {
      double w = Tm[] - TCUT;
      Vhot += dv();
      wsum += w*dv();
      rw   += w*sqrt (sq(x) + sq(y))*dv();
      xw   += w*x*dv();
    }
  }
  double r_flame = (wsum > 0.) ? rw/wsum : 0.;
  double x_flame = (wsum > 0.) ? xw/wsum : 0.;

  stats sT = statsf (Tm);
  stats sr = statsf (drhodt);
  stats so = statsf (omega);

  /**
  The 2 mm line of the real case. `Tp2` sits where `T_profile_2mm.dat` shows
  the largest jitter, so `hf_std(Tp2)/hf_std(Tmax)` is the ratio that the
  reference case must reproduce. */

  double xp = 0.5*H0 + 2e-3;
  double Tp1 = probe (Tm, xp, 4.0e-3);
  double Tp2 = probe (Tm, xp, 9.6e-3);

  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[] - porosity[])*rhoS*dv();

  /**
  Health of the solid source. `Qsrc` comes from `omega[]`, which comes from
  the solid branch of `chemistry.h`, which `TS` drives. A smooth `TS` with a
  flickering `Qsrc` means that something switches cells in and out of the sum.
  These seven numbers rank the candidates.

    n_solid     cells with f > F_ERR
    n_skip      of those, the cells a guard skipped (solid_diag == 2)
    n_bad       of those, the cells whose solve diverged (solid_diag == 3)
    n_off       of those, the cells with omega == 0, which also counts the
                cells that the clamps of shrinking.h removed
    n_nan       cells whose gas_source is not finite; one poisons the
                Poisson solve
    rhoGvS_min  minimum of rhoGv_S over the solid cells. update_properties()
                resets it to zero and refills it only under its own
                conditions, and shrinking.h:182 divides by it with no guard
    dTS_max     largest change of TS over one step, per unit of f

  Read `n_skip`, `n_bad` and `n_off` first. If they oscillate, the flicker is
  a switch and the sampling of `omega` is secondary. Read `dTS_max` second: it
  says whether the end-state sampling of `omega` is worth a fix at all. */

  double n_solid = 0., n_skip = 0., n_bad = 0., n_off = 0., n_nan = 0.;
  double rhoGvS_min = HUGE, dTS_max = 0.;
  foreach (reduction(+:n_solid) reduction(+:n_skip) reduction(+:n_bad)
           reduction(+:n_off) reduction(+:n_nan)
           reduction(min:rhoGvS_min) reduction(max:dTS_max)) {
    if (f[] > F_ERR) {
      n_solid += 1.;
#if SOLID_SOURCE_DIAG
      if (solid_diag[] == 2.) n_skip += 1.;
      if (solid_diag[] == 3.) n_bad += 1.;
      dTS_max = max (dTS_max, dTS_step[]);
#endif
      if (omega[] == 0.) n_off += 1.;
      rhoGvS_min = min (rhoGvS_min, rhoGv_S[]);
    }
    if (!isfinite (gas_source[]))
      n_nan += 1.;
  }
  if (n_solid == 0.)
    rhoGvS_min = 0.;

  if (pid() == 0) {
    static FILE * fp = NULL;
    if (!fp) {
      fp = fopen ("probe.dat", restarted ? "a" : "w");
      if (fp == NULL) {
        fprintf (stderr, "Error opening probe.dat\n");
        exit (1);
      }
      fprintf (fp, "#1:t 2:dt 3:Tmax 4:r_flame 5:x_flame 6:Vhot 7:Qsrc 8:Qrho"
                   " 9:Qdiv 10:resmax 11:drhodt_min 12:drhodt_max 13:omega_max"
                   " 14:Tp1 15:Tp2 16:mgp_i 17:mgp_resa 18:ncells 19:mass"
                   " 20:n_solid 21:n_skip 22:n_bad 23:n_off 24:n_nan"
                   " 25:rhoGvS_min 26:dTS_max\n");
      fprintf (fp, "# averaged=%d rhomean=%d zeta=%d adapt_zdiff=%d"
                   " flame_dt=%g DT=%g"
                   " CFL=%g maxlevel=%d Tcut=%g nranks=%d\n",
               (int) gas_source_averaged, (int) gas_source_rho_mean,
               (int) (ZETA_POLICY), ADAPT_ZDIFF,
               (double) FLAME_PRINT_TIME, DT, CFL, maxlevel, (double) TCUT,
               npe());
    }
    fprintf (fp, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %g %ld %g"
                 " %g %g %g %g %g %g %g\n",
             t, dt, sT.max, r_flame, x_flame, Vhot, Qsrc, Qrho, Qdiv, resmax,
             sr.min, sr.max, so.max, Tp1, Tp2, mgp.i, mgp.resa,
             grid->tn, solid_mass/solid_mass0,
             n_solid, n_skip, n_bad, n_off, n_nan, rhoGvS_min, dTS_max);
    fflush (fp);
  }
}

#if TREE
event adapt (i++) {
  scalar oxidiser = YGList_G[OpenSMOKE_IndexOfSpecies ("O2")];

#if ADAPT_ZDIFF
  /**
  Caution: `zmix` is written only every `FLAME_PRINT_TIME`. Between those
  writes this criterion follows a frozen flame position. That is hypothesis 2. */

  scalar zdiff[];
  foreach()
    zdiff[] = zmix[] - zsto[];

  adapt_wavelet_leave_interface ({T, oxidiser, zdiff}, {f},
      (double[]){5e0, 1e-2, 1e-2}, maxlevel, minlevel, 2);
#else
  adapt_wavelet_leave_interface ({T, oxidiser}, {f},
      (double[]){5e0, 1e-2}, maxlevel, minlevel, 2);
#endif

  // Unrefine for outflow condition
  unrefine (x > L0*0.4);
}
#endif

/**
Checkpoint on `i++`, not on `t += ...`, so it does not re-quantise `dt`. */

event snapshot (i = 2000; i += 2000) {
  dump ("last-snapshot");
}

event stop (t = TEND) {
#if IGNITION
  dump ("ignition-snapshot");
  if (pid() == 0)
    fprintf (stderr, "# wrote ignition-snapshot at t = %g\n", t);
#endif
}
