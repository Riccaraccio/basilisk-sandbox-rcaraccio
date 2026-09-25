/**
# Combustion of a biomass particle in hot air: the Lu case

A spherical biomass particle of 9.5 mm burns in air at 1050 K, which enters at
0.6 m/s. The case compares against the measurements in `data/lu/`:

| file | quantity |
|---|---|
| `mass` | the remaining solid mass |
| `T-particle-core` | the temperature at the centre of the particle |
| `T-particle-surf` | the temperature at the surface of the particle |

## The configuration

This case takes the configuration that the oscillation campaign of 2026-08-31
to 2026-09-10 settled for `run/fatehi-combustion.c`. Read the header of that
file for the run that supports each value. Do not change a value here without
the run that supports it. The campaign ran on the Fatehi case only, so a value
below is a transfer, not a measurement on this case.

Applied:

- `TOLERANCE = 1e-5` and `NITERMIN = 2`, because Basilisk scales the tolerance
  of the projection as `TOLERANCE/dt^2`.
- `init_grid (1 << min (maxlevel, 8))` in `main()`, and a `refine()` near the
  particle in `event init`, before `fraction()`. One cell holds 1104 fields,
  so a uniform grid at `maxlevel` allocates gigabytes before the first adapt.
  A `refine()` in `main()` does nothing: read `event init`.
- `event output (t += 0.01)`, so that the flicker of the flame does not fold
  into the band below 1 Hz.
- `FROZEN_CELL_GATE = 1` and `CORRECTIVE_CFL = 0.8`.
- `zeta_policy = ZETA_REACTION`. The shrinkage follows the local rate of
  reaction. The old case used `ZETA_CONST`.
- `INT_TEMP_VOFBC` with `INT_TEMP_PICARD`, as in the Fatehi case.
- `CFL` is set in `event init`, not in `main()`.
- A guard on `pid() == 0` for every message and every file. Every collective
  call (`statsf`, `interpolate`, `avg_interface`) runs on every rank.

`GAS_PHASE_REACTIONS` is NOT set here. No file in `src/` tests it, so it is a
dead flag. Use `TURN_OFF_GAS_REACTIONS` to switch off the gas kinetics. */

#define INT_TEMP_VOFBC 1
#define INT_TEMP_PICARD 1

#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1
#define RADIATION_TEMP 1273

/**
The knobs below are overridable from the `Makefile`, so that an A/B keeps one
source.

Caution: `CORRECTIVE_CFL` is a floating constant. The preprocessor rejects a
floating constant in an `#if`, so test it at run time, never with `#if`. */

#ifndef FROZEN_CELL_GATE
# define FROZEN_CELL_GATE 1
#endif

#ifndef CORRECTIVE_CFL
# define CORRECTIVE_CFL 0.8
#endif

#ifndef MAXLEVEL
# define MAXLEVEL 11
#endif

#ifndef TEND
# define TEND 100.
#endif


#ifndef ZETA_POLICY
# define ZETA_POLICY ZETA_REACTION
#endif


/**
`DT_VALUE` caps the step. The CFL usually binds below it, so this value only
stops a runaway. Do not push the step below 2e-4: the cost of the chemistry
per second of physical time grows as 1/dt. */

#ifndef DT_VALUE
# define DT_VALUE 5e-4
#endif

/**
`CFL_VALUE` is applied in `event init`, never in `main()`. The `defaults`
event of `navier-stokes/centered.h` sets `CFL = 0.8` and runs after `main()`.
An `init` event runs after every `defaults` event, so a value set there
survives. */

#ifndef CFL_VALUE
# define CFL_VALUE 0.5
#endif

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "gravity.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
#include "view.h"
#include "flame.h"

#if VTK_OUTPUT
#include "vtk.h"
#endif

const double Uin = 0.6; //inlet velocity
/**
Caution: give `pf` the same conditions as `p`. Basilisk does not copy them,
so without these lines every boundary of `pf` is Neumann, and the projection
of `uf` has no solution. See the comment in `fatehi-combustion.c`. */

u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
pf[left]     = neumann (0.);
psi[left]    = dirichlet (0.);

u.n[top]     = neumann (0.);
psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);
psi[right]    = neumann (0.);

const double tend = TEND; //simulation time
int maxlevel = MAXLEVEL; int minlevel = 3;

double D0 = 9.5e-3;
double solid_mass0 = 0.;

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

int main() {

  /**
  Caution: under MPI every rank shares this stderr. Guard the message with
  `pid() == 0`, or the log carries one copy per rank.

  Read this line before you quote a run. It prints what the build enabled,
  not what the directory name promises. */

  if (pid() == 0)
    fprintf (stderr, "# lu: maxlevel=%d DT=%g CFL=%g Uin=%g tend=%g"
                     " zeta=%d frozen=%d corrCFL=%g"
                     " averaged=%d exact=%d nranks=%d\n",
             MAXLEVEL, (double) DT_VALUE, (double) CFL_VALUE, Uin,
             (double) TEND, (int) ZETA_POLICY, FROZEN_CELL_GATE, (double) CORRECTIVE_CFL,
             (int) gas_source_averaged, (int) GAS_SOURCE_EXACT, npe());

  lambdaSmodel = L_LU;
  TS0 = 300.; TG0 = 1050.;
  rhoS = 1000.; // biomass density
  eps0 = 0.4;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_POLICY;

  DT = DT_VALUE;

  G.x = -9.81;

  kinfolder = "biomass/Solid-gas-88";
  shift_prod = true;

  L0 = 20*D0;
  origin (-L0/2, 0);

  /**
  One cell holds 1104 fields. Start coarse, so that the first adapt and the
  chemistry event of `i = 0` do not run on a uniform grid at `maxlevel`.
  `event init` refines near the particle. */

  init_grid (1 << min (maxlevel, 8));

  emissivity = emissivity_lu;

  /**
  The projection. `project_sf()` passes `TOLERANCE/sq(dt)` to `poisson()`, so
  the default 1e-3 stops the solve after one cycle at every step. Keep both
  lines together.

  Caution: no `defaults` event resets `TOLERANCE` or `NITERMIN`, so `main()`
  is the right place for them. `CFL` is the opposite case; see `event init`. */

  TOLERANCE = 1e-5;
  NITERMIN = 2;

  run();
}

/**
Open one output file. Only rank 0 writes, so only rank 0 opens; a handle on
another rank would truncate the same path. */

static FILE * open_output (const char * name)
{
  if (pid() != 0)
    return NULL;

  FILE * fp = fopen (name, restarted ? "a" : "w");
  if (fp == NULL) {
    fprintf (stderr, "Error opening %s\n", name);
    exit (1);
  }
  return fp;
}

event init (i = 0) {

  /**
  Caution: `navier-stokes/centered.h` assigns `CFL = 0.8` in its `defaults`
  event, which runs after `main()`. So the value belongs here. */

  CFL = CFL_VALUE;

  /**
  Refine near the particle BEFORE `fraction()`. `fraction()` computes the
  volume fraction on the grid that exists at this point.

  Caution: do not move this `refine()` to `main()`. `run()` calls
  `init_grid (N)` again (`$BASILISK/run.h:17`), and `init_grid` of the tree
  frees the grid. A `refine()` in `main()` does nothing, so the particle
  started on a uniform grid at level 8. The first adapt then built the finest
  cells from coarse PLIC lines: the solid lost 0.3 % at the first step, and the
  corner of the pellet smeared. `run/shrink-corner.c` measures this.

  The disc holds the particle and no more: the corner of a square pellet is at
  0.71 of its size. One cell holds 1104 fields, so a disc of 4 sizes at level
  11 would hold 2.3 GB, and the chemistry event of `i = 0` would run on all of
  it. The adapt of the first steps refines the gas near the particle.

  Caution: runs before this change took `solid_mass0` from the level 8 `f0`.
  Their normalized mass reads about 0.3 % lower. Compare new runs with new
  runs. */

  refine (circle (x, y, 0.75*D0) > 0. && level < maxlevel);

  scalar f0[];
  fraction (f0, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]   = 0.3104;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")]   = 0.1175;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]   = 0.1245;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]   = 0.0004;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]   = 0.0001;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]   = 0.0407;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]    = 0.0004;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]    = 0.0060;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")]  = 0.1700;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BMOIST")] = 0.2300;

  /**
  Caution: the porosity follows `f0`, not `f`. The solver has not assigned
  `f` yet at this point. */

  foreach()
    porosity[] = eps0*f0[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += f0[]*(1. - eps0)*rhoS*dv();

  TG[left] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);
  TG[right] = neumann (0.);
  TG[bottom] = neumann (0.);

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

Caution: keep the sampling at 0.01 s. At 0.1 s the Nyquist limit is 5 Hz, and
the flicker of the flame folds into the band below it.

Caution: `statsf`, `interpolate` and `avg_interface` are collective. They run
on every rank, before the guard on `pid()`. */

event output (t += 0.01) {
  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[] - porosity[])*rhoS*dv();

  double T_center = interpolate (T, 0., 0.);

  // calculate surface temperature as average over interface cells
  double T_surface = avg_interface (TInt, f, F_ERR);

  // calculate solid temperature as average over interface cells
  double TS_avg = 0.;
  int count = 0;
  foreach (reduction(+:TS_avg) reduction(+:count)) {
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      TS_avg += TS[]/f[];
      count++;
    }
  }

  if (count > 0)
    TS_avg /= count;

  stats sT = statsf (T);

  if (pid() != 0)
    return 0;

  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = open_output (name);

  fprintf (fp, "%g %g %g %g %g %g\n", t, solid_mass/solid_mass0, sT.max,
           T_center, T_surface, TS_avg);
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

event movie (t += 1) {
  scalar XH2O_G = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar XH2O_S = XGList_S[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar XH2O[];
  foreach()
    XH2O[] = XH2O_S[]*f[] + XH2O_G[]*(1. - f[]);

  clear();
  view (theta=0, phi=0, psi=-pi/2., width = 2400, height = 2400);
  draw_vof ("f", lw = 1.5);
  squares ("T", min = 300, max = 2000, spread = -1, linear = true);
  isoline ("zmix - zsto", lw = 1.5, lc = {1., 1., 1.});
  mirror ({0, 1}) {
    squares ("XH2O", min = 0, max = 1., spread = -1, linear = true);
    cells();
    draw_vof ("f", lw = 1.5);
  }
  save ("movie.mp4");
}

event dump (t = 1; t += 1) {
  dump ("last-snapshot");
}

#if VTK_OUTPUT
event vtk (t += 10; t <= 80) {

  mixture_fraction (zmix);
  scalar zdiff[];
  foreach()
    zdiff[] = zmix[] - zsto[];

  char name[120];
  sprintf (name, "fatehi-%d.vtk", (int) (t));
  FILE* fvtk = fopen (name, "w");

  // H2O
  scalar XH2O_G = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar XH2O_S = XGList_S[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar XH2O[];
  foreach()
    XH2O[] = XH2O_S[]*f[] + XH2O_G[]*(1. - f[]);

  // CO2
  scalar XCO2_G = XGList_G[OpenSMOKE_IndexOfSpecies ("CO2")];
  scalar XCO2_S = XGList_S[OpenSMOKE_IndexOfSpecies ("CO2")];
  scalar XCO2[];
  foreach()
    XCO2[] = XCO2_S[]*f[] + XCO2_G[]*(1. - f[]);

  // OH
  scalar XOH_G = XGList_G[OpenSMOKE_IndexOfSpecies ("OH")];
  scalar XOH_S = XGList_S[OpenSMOKE_IndexOfSpecies ("OH")];
  scalar XOH[];
  foreach()
    XOH[] = XOH_S[]*f[] + XOH_G[]*(1. - f[]);

  output_vtk ({f, T, XH2O, XCO2, XOH, u.x, u.y, zdiff}, n=(1<<maxlevel), fp=fvtk, linear=true);
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
~~~gnuplot Mass
reset
set terminal svg size 450, 450
set output "lu-mass.svg"
set xlabel "Time [s]"
set ylabel "Normalized Mass [-]"
set yrange [0.:1.1]
set xrange [0:100]

plot  "../../data/lu/mass" u 1:(1 - ($2 + $3)/2):(abs(($2 + $3)/2 -$2)) w yerrorbars pt 4 ps 0.8 lc "black" notitle, \
      "cluster/10/lu-combustion/OutputData-10" u 1:2 w l lw 2 lc "black" t "Lu lambda", \
      "cluster/galgano/lu-combustion/OutputData-9" u 1:2 w l lw 2 lc "red" t "Galgarno lambda",\
      "cluster/ignition/lu-combustion/OutputData-9" u 1:2 w l lw 2 lc "blue" t "Increased rho"

~~~

~~~gnuplot Temperature evolution
reset
set terminal svg size 450, 450
set output "T-lu.svg"
set xlabel "Time [s]"
set ylabel "Temperature [K]"
set yrange [200:1900]
set xrange [0:100]
set key bottom right

set dashtype 11 (2,2,2,2)

plot  "../../data/lu/T-particle-core" u 1:(($2 + $3)/2):(abs(($2 + $3)/2 -$2)) w yerrorbars pt 4 ps 0.8 lc "black" notitle, \
      "../../data/lu/T-particle-surf" u 1:(($2 + $3)/2):(abs(($2 + $3)/2 -$2)) w yerrorbars pt 4 ps 0.8 lc "red" notitle, \
      "cluster/10/lu-combustion/OutputData-10" u 1:4 w l lw 2 dt 1  lc "black" title "Lu lambda", \
      "cluster/10/lu-combustion/OutputData-10" u 1:6 w l lw 2 dt 11 lc "black" notitle, \
      "cluster/galgano/lu-combustion/OutputData-9" u 1:4 w l lw 2 dt 1  lc "red" title "Galgano lambda", \
      "cluster/galgano/lu-combustion/OutputData-9" u 1:6 w l lw 2 dt 11 lc "red" notitle, \
      "cluster/ignition/lu-combustion/OutputData-9" u 1:4 w l lw 2 dt 1 lc "blue" title "Increased rho", \
      "cluster/ignition/lu-combustion/OutputData-9" u 1:6 w l lw 2 dt 11 lc "blue" notitle
~~~

*/
