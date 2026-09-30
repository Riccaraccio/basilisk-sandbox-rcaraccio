/**
# Combustion of a small biomass particle in hot air: the Glarborg case

A biomass particle burns in air at 1473 K, which enters at 1.5 m/s.
`CASE_NUMBER` selects the shape. All five shapes have a volume near that of
a sphere of 3 mm:

| `CASE_NUMBER` | shape | D0 [mm] | H0 [mm] |
|---|---|---|---|
| 0 | sphere | 3.00 | - |
| 1 | cylinder | 2.08 | 4.16 |
| 2 | cylinder | 1.65 | 6.60 |
| 3 | cylinder | 1.44 | 8.65 |
| 4 | cylinder | 1.31 | 10.48 |

## The configuration

The numerical settings are the ones of `run/fatehi-combustion.c`. Read its
header for the reasons and the known limits. They were tested on the Fatehi
case, so here they are a transfer, not a measurement on this case.

| setting | value |
|---|---|
| domain | `L0 = 20 max (D0, H0)`, axisymmetric |
| grid | `MAXLEVEL` 10, minimum level 3, adapt on `{T, O2}` with `{5, 1e-2}` |
| timestep | cap 5e-4 s, `CFL` 0.5 (the CFL binds far below the cap) |
| projection | `TOLERANCE` 1e-5, `NITERMIN` 2 |
| shrinkage | `ZETA_REACTION` |
| inlet | air at 1473 K, 1.5 m/s |
| outlet | Neumann `u`, `TG`, `YG`; Dirichlet `p`, `pf` |

The domain scales with the particle, so the cell size changes with the case.
At level 10 it is 59 um for case 0 and 205 um for case 4. The Fatehi case
runs at 78 um. Caution: case 4 has only 6 cells across the diameter at level
10. Build it with `MAXLEVEL=11` or more for a result that you quote.

Every collective call (`statsf`) runs on every rank, and only rank 0 writes a
message or a file. */

#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1
#define FLAME_PRINT_TIME 0.01

/**
The knobs below are overridable with `-D`. `SNAPSHOT_EVERY` writes
`snapshot-<t>` every that many seconds beside `last-snapshot`. Set it to 0 to
turn the snapshots off. */

#ifndef CASE_NUMBER
# define CASE_NUMBER 0
#endif

#ifndef MAXLEVEL
# define MAXLEVEL 10
#endif

#ifndef TEND
# define TEND 15.
#endif

#ifndef SNAPSHOT_EVERY
# define SNAPSHOT_EVERY 5
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
#include "view.h"
#include "flame.h"

const double Uin = 1.5; //inlet velocity
/**
Caution: give `pf` the same conditions as `p`. Basilisk does not copy them,
so without these lines every boundary of `pf` is Neumann, and the projection
of `uf` has no solution. */

u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
pf[left]     = neumann (0.);
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);
psi[right]    = neumann (0.);

const double tend = TEND; //simulation time
int maxlevel = MAXLEVEL, minlevel = 3;
double solid_mass0 = 0.;
double D0, H0;

double D0_arr[] = {3e-3, 2.08e-3, 1.65e-3, 1.44e-3, 1.31e-3};
double H0_arr[] = {0., 4.16e-3, 6.60e-3, 8.65e-3, 10.48e-3};

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

int main() {

  lambdaSmodel = L_LU;
  TS0 = 300.; TG0 = 1473.;
  rhoS = 1000;
  eps0 = 0.4;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_REACTION;

  DT = 5e-4; // the CFL binds far below this cap

  G.x = -9.81;

  kinfolder = "biomass/Solid-gas-88";
  shift_prod = true;

  if (CASE_NUMBER < 0 || CASE_NUMBER >= 5) {
    if (pid() == 0)
      fprintf (stderr, "Invalid CASE_NUMBER %d. Must be between 0 and 4.\n",
               CASE_NUMBER);
    return 1;
  } else {
    D0 = D0_arr[CASE_NUMBER];
    H0 = H0_arr[CASE_NUMBER];
  }

  L0 = 20*max (D0, H0);
  origin (-L0/2, 0);
  emissivity = emissivity_lu;

  /**
  One cell holds 1104 fields. Start coarse, so that the first adapt and the
  chemistry event of `i = 0` do not run on a uniform grid at `maxlevel`.
  `event init` refines near the particle. */

  init_grid (1 << min (maxlevel, 8));

  /**
  The projection. `project_sf()` passes `TOLERANCE/sq(dt)` to `poisson()`, so
  the default 1e-3 stops the solve after one cycle at every step. Keep both
  lines together.

  Caution: no `defaults` event resets `TOLERANCE` or `NITERMIN`, so `main()`
  is the right place for them. `CFL` is the opposite case; see `event init`. */

  TOLERANCE = 1e-5;
  NITERMIN = 2;

  /**
  Read this line before you quote a run. It prints what the build enabled,
  not what the directory name promises. */

  if (pid() == 0)
    fprintf (stderr, "# glarborg: case=%d D0=%g H0=%g maxlevel=%d DT=%g"
                     " Uin=%g tend=%g tol=%g nitermin=%d filter=%d"
                     " strang=%d frozen=%d snapevery=%d nranks=%d\n",
             CASE_NUMBER, D0, H0, MAXLEVEL, DT, Uin, (double) TEND,
             TOLERANCE, NITERMIN, gas_source_filter_passes,
             GAS_CHEMISTRY_STRANG, FROZEN_CELL_GATE, SNAPSHOT_EVERY, npe());

  run();
}

event init (i = 0) {

  /**
  Caution: `navier-stokes/centered.h` assigns `CFL = 0.8` in its `defaults`
  event, which runs after `main()`. So the value belongs here. */

  CFL = 0.5;

  /**
  Refine near the particle BEFORE `fraction()`. `fraction()` computes the
  volume fraction on the grid that exists at this point.

  Caution: do not move this `refine()` to `main()`. `run()` calls
  `init_grid (N)` again (`$BASILISK/run.h:17`), and `init_grid` of the tree
  frees the grid. A `refine()` in `main()` does nothing: the particle then
  starts on a uniform grid at level 8, the first adapt builds the finest
  cells from coarse PLIC lines, and the solid loses 0.3 % at the first step.

  The disc holds the particle and no more: the corner of a square pellet is at
  0.71 of its size. One cell holds 1104 fields, so a disc of 4 sizes at level
  11 would hold 2.3 GB, and the chemistry event of `i = 0` would run on all of
  it. The adapt of the first steps refines the gas near the particle. */

  refine (circle (x, y, 0.75*max (D0, H0)) > 0. && level < maxlevel);

  scalar f0[];
  if (CASE_NUMBER == 0) {
    fraction (f0, circle (x, y, 0.5*D0));
  } else {
    fraction (f0, superquadric (x, y, 20, 0.5*H0, 0.5*D0));
  }

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.4752;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")]  = 0.2039;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.1547;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0202;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0000;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0332;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0168;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0030;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.0930;

  /**
  Caution: the porosity follows `f0`, not `f`. The solver has not assigned
  `f` yet at this point. */

  foreach()
    porosity[] = eps0*f0[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += f0[]*(1. - eps0)*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  TG[left] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[left] = dirichlet (0.765);
      YG[top] = dirichlet (0.765);
    } else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      YG[left] = dirichlet (0.235);
      YG[top] = dirichlet (0.235);
    }
    else {
      YG[left] = dirichlet (0.);
      YG[top] = dirichlet (0.);
    }
  }

  if (restore (file = "last-snapshot", list = all)) {
    if (pid() == 0)
      fprintf (stderr, "Restart file found!\n");
    restarted = true;
  } else {
    if (pid() == 0)
      fprintf (stderr, "No restart file found, starting from scratch!\n");

    foreach() {
      f[] = f0[];
      porosity[] = eps0*f[];
    }
  }
}

/**
The scalar series.

Caution: keep the sampling at 0.01 s. A slower sampling folds the flicker of
the flame into the low band.

Caution: `statsf` is collective. It runs on every rank, before the guard on
`pid()`. Only rank 0 opens and writes the file. */

event output (t += 0.01) {

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[] - porosity[])*rhoS*dv();

  // Calaculate overall char mass
  double char_mass = 0., wood_mass = 0.;
  foreach (reduction(+:char_mass) reduction(+:wood_mass)) {
    if (f[] > F_ERR) {
      double local_char_fraction = calculate_char_fraction(point, YSList, f);
      double local_moist_fraction = calculate_moisture_fraction(point, YSList, f);
      char_mass += local_char_fraction*(f[] - porosity[])*rhoS*dv();
      wood_mass += (1. - local_char_fraction - local_moist_fraction)*(f[] - porosity[])*rhoS*dv();
    }
  }

  /**
  The early warning of a runaway of the gas temperature in a thin cell, as in
  `fatehi-combustion.c`. Nothing bounds `TG` in a sliver cell. `TGmin_gas` is the smallest `TG` over the cells with
  `f < F_ERR`, and `nTGneg` the number of cells with `TG < 0`. Both do not
  depend on the position within the step. If `nTGneg` is not 0, or
  `TGmin_gas` falls below about 250 K, stop the run. */

  double TGmin_gas = HUGE, nTGneg = 0.;
  foreach (reduction(min:TGmin_gas) reduction(+:nTGneg)) {
    if (f[] < F_ERR)
      TGmin_gas = min (TGmin_gas, TG[]);
    if (TG[] < 0.)
      nTGneg += 1.;
  }

  stats sT = statsf (T);

  if (pid() != 0)
    return 0;

  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = NULL;
  if (fp == NULL) {
    fp = fopen (name, restarted ? "a" : "w");
    if (fp == NULL) {
      fprintf (stderr, "Error opening %s\n", name);
      exit (1);
    }
  }

  if (i == 0)
    fprintf (fp, "# t(1), Ms/Ms0(2), Tmax(3), Char/Ms0(4), Wood/Ms0(5),"
                 " TGmin_gas(6), nTGneg(7), dt(8)\n");

  fprintf (fp, "%g %g %g %g %g %g %g %g\n", t, solid_mass/solid_mass0, sT.max,
           char_mass/solid_mass0, wood_mass/solid_mass0, TGmin_gas, nTGneg,
           dt);

  fflush (fp);
}

#if TREE
event adapt (i++) {
  scalar oxidiser = YGList_G[OpenSMOKE_IndexOfSpecies ("O2")];

  adapt_wavelet_leave_interface ({T, oxidiser}, {f},
    (double[]){5e0, 1.e-2}, maxlevel, minlevel, 2);

  // Unrefine for outflow condition
  unrefine (x > L0*0.4);
}
#endif

event movie (t += 0.1) {
  clear();
  view (theta=0, phi=0, psi=-pi/2., width = 1080, height = 1080);
  squares ("T", min = 300, max = 2000, spread = -1, linear = true);
  isoline ("T", val = statsf(T).max);
  isoline ("zmix - zsto", lw = 1.5, lc = {1., 1., 1.});
  draw_vof ("f", lw = 1.5);
  mirror ({0, 1}) {
    squares ("O2_G + O2_S", min = 0., max = 0.235, spread = -1, linear = true);
    isoline ("zmix - zsto", lw = 1.5, lc = {1., 0., 0.});
    draw_vof ("f", lw = 1.5);
  }
  save ("movie.mp4");
}

event dump (t = 1; t += 1) {
  dump ("last-snapshot");
}

/**
The numbered snapshots, see `SNAPSHOT_EVERY`. A restart reads
`last-snapshot`, so copy the chosen `snapshot-<t>` to `last-snapshot` before
the restart. */

#if SNAPSHOT_EVERY > 0
event snapshot (t = SNAPSHOT_EVERY; t += SNAPSHOT_EVERY) {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  dump (name);
}
#endif

/**
Caution: `return 1` is what ends the run, not the time of the event. If a
later edit adds an event with a condition such as `t <= X`, a bare
`event stop (t = tend)` does not end the run. */

event stop (t = tend) {
  return 1;
}

/**
~~~gnuplot
~~~
**/
