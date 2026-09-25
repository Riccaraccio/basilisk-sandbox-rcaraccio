/**
# Combustion of a biomass pellet in hot air: the Fatehi case

A cylindrical biomass pellet of 8 mm burns in air at 1123 K, which enters at
0.13 m/s. The case compares against the measurements in `data/fatehi/`:

| file | quantity | time range |
|---|---|---|
| `mass` | the remaining solid mass | 0 to 268 s |
| `T-2mm`, `T-11mm` | the temperature of the path | 0 to 66 s |
| `yH2O-2mm`, `yH2O-11mm` | the H2O mole fraction of the path | 0 to 68 s |

The experiment reads the temperature and the H2O along a horizontal beam, so
the case must give a path average, not a point value. `path_average_H2O()`
does that. `misc/WMS/WMS.py` reads the profile files below and gives the
sensor-equivalent value that the plots use.

## The configuration

The oscillation campaign of 2026-08-31 to 2026-09-10 measured every knob that
this case can turn. The block below carries what that campaign settled. Read
`run/test.c` for the solid ladder and `run/gas-source.c` for the gas ladder.
Do not change a value here without the run that supports it.

Applied, each one a removal of an artifact:

- `TOLERANCE = 1e-5` and `NITERMIN = 2`. Basilisk scales the tolerance of the
  projection as `TOLERANCE/dt^2`, so the default 1e-3 stops the solve after
  one multigrid cycle at every step. A refinement near the flame then left a
  residual 100 times larger and a spike of 30 K that lasted 0.3 s.
- `init_grid (1 << min (maxlevel, 8))` in `main()`, and a `refine()` of a disc
  of 0.75 D0 in `event init`, before `fraction()`. One cell of this case holds
  1104 fields, which is 8.8 kB, so `init_grid (1 << 10)` allocates 9.3 GB
  before the first adapt, and the chemistry event of `i = 0` runs on all of
  it. A `refine()` in `main()` does nothing: read `event init`.
- `event output (t += 0.01)`. At 0.1 s the Nyquist limit is 5 Hz, and the
  flame flickers near 17 Hz, so the jitter folds into the band below 1 Hz and
  the run appears to wander. The profile events stay at 0.1 s, because they
  write 1025 samples per line and per file.
- `FROZEN_CELL_GATE`. `biomass/Solid-gas-88` costs about 100 evaluations of
  the right-hand side in a cell that does not react. The gate spends one
  evaluation to find that out. At the default tolerance of 1e-15 it gives the
  same answer in every column of every row, and it reaches 1.28 times more
  simulated time for the same wall clock.
- `CORRECTIVE_CFL = 0.8`. The corrective velocity of `MOLAR_DIFFUSION` and
  `FICK_CORRECTED` reaches 0.333 m/s against an inflow of 0.13 m/s, and no
  event limited the step by it, so the scheme ran above a Courant number of 1
  at the front. The limit at 0.5 costs 2.1 times the step; at 0.8, which
  matches the flow CFL, it costs 1.3 times.
- `zeta_policy = ZETA_REACTION`. The shrinkage follows the local rate of
  reaction. `ZETA_CONST` splits the released volume in half everywhere; the
  A/B from one ignition snapshot gives it a worse probe temperature (50.8 K
  against 35.4 K peak to peak) and a worse peak reaction rate (34.7 % against
  18.5 %), and it does not remove the slow band.
- A guard on `nodata` in `print_profile()`. `interpolate()` gives 1e30 when no
  rank holds the point. Two rows of the archived run carried that value, and
  `WMS.py` read it as a temperature.

OPEN, and switched off here for one reason only. Read this before you decide:

- `INT_TEMP_VOFBC` is the intended replacement for `INT_TEMP_ROBIN`, and
  `multicomponent-varprop.h:78` makes the two an `#error` against each other.
  It puts the interface temperature inside the operator of the Poisson solve,
  so `plic_flux()` rebuilds the interface gradient on every relaxation sweep
  instead of freezing it at step n. It removes both sliver crashes. It is NOT
  refuted, and it is not the default of `src/`: `INT_TEMP_VOFBC` has no
  `#ifndef` block anywhere, so an undefined name reads as 0.

  Two of the three ladders run to the end with it, and neither one warns:

  | run | flags beyond the reduced case | reached | `pf` warnings |
  |---|---|---|---|
  | `test-vofbcfl` | none | t = 40 | 0 |
  | `test-vofbcm` | `MOISTURE` | t = 39.74 | 0 |
  | `test-fullvofbc` | the full set | **t = 10.34** | **84** |

  Only the full set fails, and this case carries the full set. It fails as a
  collapse of the timestep, which is the signature that `INT_TEMP_ROBIN`
  gave: over 0.05 s the step falls 3.4e-4 -> 1.7e-5 and the residual of the
  projection climbs 2.3 -> 113. The `pf` warnings start at line 5 of the log,
  which is the first step, so the projection is unhappy from the beginning.
  The build is current: `stat -c %y test-fullvofbc/warn` gives 2026-09-09
  16:38, after the fix of `src/plicbc.h` at 11:13.

  The difference between `test-vofbcm`, which passes, and `test-fullvofbc`,
  which fails, is five ingredients: `Da`, gravity, the shape, the emissivity
  of Di Blasi, and the three transport flags. The transport flags are the
  suspect, because `INT_TEMP_VOFBC` and they act on the same diffusion solves
  of `multicomponent-varprop.h`, and no other VOFBC run carries them.

  So: switch it on with `-DINT_TEMP_VOFBC=1 -DINT_TEMP_PICARD=1` when the
  interaction is understood, and repeat `test-fullvofbc` first. Do not run a
  comparison against the experiment on a build that stops at t = 10.34 s.

- `INT_TEMP_PICARD` belongs with `INT_TEMP_VOFBC`, and the rung `TLAD_VOFBC`
  carries both. On its own it is open, not refuted. The measurement that
  reads badly for it, a slow band of 23.5 % against 10.4 %, comes from a
  build that also carried the conductance at `INT_TEMP_ROBIN_SMAX = 1`, and
  that setting alone moves the surface temperature by 2 K where the loop
  moves it by 0.07 K. Every Picard run of the archive inherits that artefact
  and has to be repeated.

- `INT_TEMP_ROBIN`. The conductance biases the surface temperature by 20 to
  25 K at `SMAX = 1`, and `SMAX = 20` is a mitigation, not a cure. Prefer
  `INT_TEMP_VOFBC`, which supersedes it.

Refuted, and switched OFF. Do not switch one on without a new run:

- `GAS_SOURCE_EXACT`. The exact logarithmic form of the expansion is a no-op
  on this two-dimensional flame: it matches the averaged form to three digits.
- `PHASE_AWARE_PROPERTIES`. The zero mixture density comes from a runaway of
  the gas temperature in a thin cell, not from the prolongation after adapt.
- `PIN_SOLID_INTERIOR`. Leave it at 0 for a production run. Set it to 1 only
  for a campaign that compares two cases, so that the interior of the pellet
  coarsens the same way in both.
- `maxlevel = 11`. The slow band falls 1.4 to 1.9 times from level 10 to
  level 11, but level 11 and level 12 both meet the runaway of the gas
  temperature in a thin cell.

Taken from the defaults of `src/`, which the campaign moved:

- `GAS_SOURCE_AVERAGED = 1`. The expansion source of the gas is the step
  average, not the rate at the end state. The end-state form has the wrong
  sign in a cell that burns out inside one step.
- `CORRECTIVE_LIMITER = 1` and `MDE_INTERFACE = 1`, the other two fixes of the
  transport flags.
- `INT_TEMP_TOL = 1`. It scales the tolerance of each temperature solve to
  `theta/dt`, which halves the multigrid cycles of the gas and leaves every
  column equal to printed precision.
- `GAS_STATE_FALLBACK = 1`, which repairs a cell whose gas state is bad.

`GAS_PHASE_REACTIONS` is NOT set here. No file in `src/` tests it, so it is a
dead flag: the mechanism of the pore gas runs whatever its value is. */

#define INT_TEMP_VOFBC 1
#define INT_TEMP_PICARD 1

#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1

/**
The knobs below are overridable from the `Makefile`, so that an A/B keeps one
source. Every one of them holds the value that the campaign settled.

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
# define TEND 150.
#endif

/**
`PROFILE_TEND` stops the six profile files. The measurements of the beam end
at 68 s, so 100 s covers every comparison and keeps the files small. */

#ifndef PROFILE_TEND
# define PROFILE_TEND 100.
#endif

/**
`DT_VALUE` caps the step. The CFL usually binds below it, so this value only
stops a runaway. Do not push the step below 2e-4: the amplitude of the slow
band saturates there, and the cost of the chemistry per second of physical
time grows as 1/dt. */

#ifndef DT_VALUE
# define DT_VALUE 5e-4
#endif

/**
`CFL_VALUE` is applied in `event init`, never in `main()`. The `defaults`
event of `navier-stokes/centered.h` sets `CFL = 0.8` and runs after `main()`,
and the `stability` event of `vof.h` then clamps it to 0.5. An `init` event
runs after every `defaults` event, so a value set there survives. */

#ifndef CFL_VALUE
# define CFL_VALUE 0.5
#endif

#define PRINT_FIELDS 0
#define VTK_OUTPUT 1

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

#if VTK_OUTPUT
#include "vtk.h"
#endif

const double Uin = 0.13; //inlet velocity
/**
Caution: give `pf` the same conditions as `p`. Basilisk does not copy the
conditions of `p` to `pf`. Without the lines for `pf`, every boundary of `pf`
is Neumann. The projection of `uf` then has no solution when the expansion
term `drhodt` is not zero, because the outflow flux is fixed. `pf` drifts
without limit, the divergence error grows, and at ignition `dt` collapses. The
production runs of 2026-09-16/17 crashed in this way (`pf` near -6e5, `sum`
2e10 in the log of the solver). */

u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
pf[left]     = neumann (0.);
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

/**
## The conditions at the outlet

`OUTLET_BC` selects how the outlet (right) treats gas that flows back into the
domain. A plain Neumann outlet lets the backflow grow without limit: in the
runs of 2026-09-16/17 the axis cell of the outlet flowed back at -0.35 m/s at
t = 10 and at -58 m/s at t = 10.227, and `dt` collapsed.

- `OUTLET_BC = 0`: the old outlet. Neumann for `u`, `TG` and `YG`.
- `OUTLET_BC = 1`: strict block (default). Where the gas flows in (`u.x < 0`
  in the cell next to the outlet), `u.n` and `u.t` are 0. The outlet is a
  wall for that flow.
- `OUTLET_BC = 2`: controlled inflow. The gas can flow in. `u.n` stays
  Neumann, `u.t` is 0, and the gas that enters is ambient gas: `TG0` and air.
  Where the gas flows out, all fields stay Neumann.

The two runs from t = 0 of 2026-09-17/18 decided the default. With 2, the
inflow at the outlet started at t = 6 and grew by about 2 times every 0.5 s,
to -0.76 m/s at t = 8 and -159 m/s at t = 8.17, where the run crashed before
ignition. With 1, `uf` on the outlet stayed at or above 0, the centred `u.x`
next to the outlet stayed above -0.035 m/s, and the run went through ignition
(Tmax 1718 K at t = 11.16) at a normal `dt`.

Caution: `-DOUTLET_BC` without a value defines the flag as 1. Always give the
value.

With 1 and 2, the gas that enters is ambient gas in both cases (see
`event init`). Only the velocity differs.

Caution: the physics of the heating phase needs an inflow at the outlet. While
the moisture evaporates, a cold plume sinks from the particle toward the inlet,
and gas from above must replace it. `OUTLET_BC = 1` makes the outlet a wall
for that gas, so it changes the flow of the heating phase.

Caution: these conditions act on the centred `u`. The projection corrects the
boundary face of `uf` with the gradient of `pf`, so `uf` can still flow in.
Read the outlet columns of `dtlimits.dat` to see it. */

#ifndef OUTLET_BC
# define OUTLET_BC 1
#endif

#if OUTLET_BC == 1
u.n[right]    = u.x[] < 0. ? dirichlet (0.) : neumann (0.);
u.t[right]    = u.x[] < 0. ? dirichlet (0.) : neumann (0.);
#elif OUTLET_BC == 2
u.n[right]    = neumann (0.);
u.t[right]    = u.x[] < 0. ? dirichlet (0.) : neumann (0.);
#else
u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
#endif
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);
psi[right]    = neumann (0.);

const double tend = TEND; //simulation time
int maxlevel = MAXLEVEL; int minlevel = 2;

double D0 = 8e-3, H0 = 8e-3;
double solid_mass0 = 0.;

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

int main() {

  /**
  Caution: under MPI every rank shares this stderr. Guard the message with
  `pid() == 0`, or the log carries one copy per rank.

  Read this line before you quote a run. It prints what the build enabled,
  not what the directory name promises. */

  if (pid() == 0)
    fprintf (stderr, "# fatehi: maxlevel=%d DT=%g CFL=%g Uin=%g tend=%g"
                     " zeta=REACTION frozen=%d corrCFL=%g"
                     " averaged=%d exact=%d outlet=%d nranks=%d\n",
             MAXLEVEL, (double) DT_VALUE, (double) CFL_VALUE, Uin,
             (double) TEND, FROZEN_CELL_GATE, (double) CORRECTIVE_CFL,
             (int) gas_source_averaged, (int) GAS_SOURCE_EXACT, OUTLET_BC,
             npe());

  lambdaSmodel = L_TENWOLDE;
  TS0 = 300.; TG0 = 1123.;
  rhoS = 1550;
  eps0 = 0.2; // low, compressed pellet

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_REACTION;

  DT = DT_VALUE;

  G.x = -9.81;

  kinfolder = "biomass/Solid-gas-88";
  shift_prod = true;

  L0 = 20*D0;
  origin (-L0/2, 0);

  /**
  One cell holds 1104 fields, which is 8.8 kB. A uniform grid at `maxlevel`
  would allocate 9.3 GB before the first adapt, and the chemistry event of
  `i = 0` would run on all of it. Start coarse. `event init` refines near the
  pellet. */

  init_grid (1 << min (maxlevel, 8));

  emissivity = emissivity_diblasi;

  /**
  The projection. `project_sf()` passes `TOLERANCE/sq(dt)` to `poisson()`, so
  the default 1e-3 with `dt = 3e-4` gives about 9e3 and the solve stops after
  one cycle at every step. Keep both lines together.

  Caution: no `defaults` event resets `TOLERANCE` or `NITERMIN`, so `main()`
  is the right place for them. `CFL` is the opposite case; see `event init`.

  Caution: do not read `mgp.resa` against `TOLERANCE`. The quantity that
  Basilisk controls is `resa*dt^2`. */

  TOLERANCE = 1e-5;
  NITERMIN = 2;

  run();
}

// Output files for profiles
FILE* fTprofile_2mm, * fTprofile_11mm;
FILE* fxH2Oprofile_2mm, * fxH2Oprofile_11mm;
FILE* fxOHprofile_2mm, * fxOHprofile_11mm;

/**
Open one profile file. Only rank 0 writes them, so only rank 0 opens them; a
handle on another rank would truncate the same path. */

static FILE * open_profile (const char * name)
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
  event, which runs after `main()`, and `vof.h` then clamps it to 0.5. So the
  value belongs here. To confirm that it took effect, read the live variable,
  not the macro. */

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

  refine (circle (x, y, 0.75*max (D0, H0)) > 0. && level < maxlevel);

  scalar f0[];
  fraction (f0, superquadric (x, y, 20, 0.5*H0, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.4344;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("GMSW")]  = 0.2108;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.1347;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0786;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0178;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0167;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0419;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0041;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.0610;

  /**
  Caution: the porosity follows `f0`, not `f`. The solver has not assigned
  `f` yet at this point. */

  foreach()
    porosity[] = eps0*f0[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += f0[]*(1. - eps0)*rhoS*dv(); //Note: (1-e) = (1-ef) != (1-e)f

  TG[left] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);
#if OUTLET_BC
  TG[right] = u.x[] < 0. ? dirichlet (TG0) : neumann (0.);
#else
  TG[right] = neumann (0.);
#endif
  TG[bottom] = neumann (0.);

  /**
  With `OUTLET_BC` set, the gas that enters through the outlet is air. The
  values must be constants in each branch: a boundary condition cannot read
  a local variable of this loop. */

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[left] = dirichlet (0.765);
      YG[top] = dirichlet (0.765);
#if OUTLET_BC
      YG[right] = u.x[] < 0. ? dirichlet (0.765) : neumann (0.);
#endif
    } else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      YG[left] = dirichlet (0.235);
      YG[top] = dirichlet (0.235);
#if OUTLET_BC
      YG[right] = u.x[] < 0. ? dirichlet (0.235) : neumann (0.);
#endif
    }
    else {
      YG[left] = dirichlet (0.);
      YG[top] = dirichlet (0.);
#if OUTLET_BC
      YG[right] = u.x[] < 0. ? dirichlet (0.) : neumann (0.);
#endif
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

  fTprofile_2mm     = open_profile ("T_profile_2mm.dat");
  fTprofile_11mm    = open_profile ("T_profile_11mm.dat");
  fxH2Oprofile_2mm  = open_profile ("xH2O_profile_2mm.dat");
  fxH2Oprofile_11mm = open_profile ("xH2O_profile_11mm.dat");
  fxOHprofile_2mm   = open_profile ("xOH_profile_2mm.dat");
  fxOHprofile_11mm  = open_profile ("xOH_profile_11mm.dat");
}

/**
## Line diagnostics along the line x = x_interp, from y = 0 to y = length

Caution: this function samples a fixed lattice of spacing `length/n_samples`
that does not follow the grid. Keep that spacing at or below the cell size, or
the quadrature throws away what the grid resolves. The default gives
`length/n_samples = Delta/2`, which is two samples per finest cell. That is
enough, because `interpolate_linear` holds no information below `Delta`.

Caution: `interpolate()` is collective. Every rank must call it the same
number of times, so the loop that samples runs on every rank and only rank 0
writes.

Caution: `interpolate()` gives `nodata`, which is 1e30, when no rank holds the
point. `WMS.py` masks no such value, so it would read 1e30 as a temperature.
This function drops the sample instead. Two rows of the archived run carried
one and gave a spike of 190 K.

The position of the sample comes from an integer index, not from an
accumulated sum. A sum of `n_samples` additions drifts, and `WMS.py` sorts the
samples of one instant by position. */

void print_profile (scalar s, double x_interp, FILE * fp, double time,
                    int n_samples = 1 << maxlevel, const double length = L0/2)
{
  double * val = malloc ((n_samples + 1)*sizeof (double));

  for (int k = 0; k <= n_samples; k++)
    val[k] = interpolate (s, x_interp, k*length/n_samples);

  if (pid() == 0)
    for (int k = 0; k <= n_samples; k++)
      if (val[k] < 1e20)                       // drop nodata
        fprintf (fp, "%g %g %g\n", time, k*length/n_samples, val[k]);

  free (val);
}

/**
## The path averages that the experiment measures

The beam of the experiment crosses the plume, so both quantities are averages
along the line x = x_interp:

  `xH2O`  the plain mean of the H2O mole fraction
  `Tavg`  the temperature that the same H2O column weights, which is
          `sum(x_H2O) / sum(x_H2O/T)`

One sweep of `foreach_region` gives both. Do not split them again: each sweep
costs one collective reduction.

`foreach_region` interpolates inside the cell that holds each sample point, so
`interpolate_linear` cannot give `nodata` here. That is why this function
needs no guard and `print_profile` does. */

void path_average_H2O (double x_interp, double * xH2O, double * Tavg,
                       int n_samples = 1 << (maxlevel - 1),
                       const double length = L0/4.)
{
  scalar XH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];

  double sum = 0., count = 0., den = 0.;
  coord pos, box[2] = {{x_interp, 0.}, {x_interp, length}}, nn = {1, n_samples};
  foreach_region (pos, box, nn, reduction(+:sum) reduction(+:count)
                  reduction(+:den)) {
    double x_local = interpolate_linear (point, XH2O, pos.x, pos.y, pos.z);
    sum   += x_local;
    count += 1.;
    den   += x_local/interpolate_linear (point, T, pos.x, pos.y, pos.z);
  }

  *xH2O = count > 0. ? sum/count : 0.;
  *Tavg = den   > 0. ? sum/den   : TG0;   // avoid a division by 0
}

/**
The profiles that `misc/WMS/WMS.py` reads. They stay at 0.1 s: each call
writes 1025 samples to each of six files, so 0.01 s would give 6 million lines
for each second of physical time. The scalar series of `event output` carries
the sampling rate instead. */

event print_profile (t += 0.1; t <= PROFILE_TEND) {
  scalar XH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar XOH = XGList_G[OpenSMOKE_IndexOfSpecies ("OH")];

  // Temperature profiles
  print_profile (T, H0/2 + 2e-3, fTprofile_2mm, t);
  print_profile (T, H0/2 + 11e-3, fTprofile_11mm, t);

  // Water vapor mole fraction profile
  print_profile (XH2O, H0/2 + 2e-3, fxH2Oprofile_2mm, t);
  print_profile (XH2O, H0/2 + 11e-3, fxH2Oprofile_11mm, t);

  print_profile (XOH, H0/2 + 2e-3, fxOHprofile_2mm, t);
  print_profile (XOH, H0/2 + 11e-3, fxOHprofile_11mm, t);

  if (pid() == 0) {
    fflush (fTprofile_2mm);
    fflush (fTprofile_11mm);
    fflush (fxH2Oprofile_2mm);
    fflush (fxH2Oprofile_11mm);
    fflush (fxOHprofile_2mm);
    fflush (fxOHprofile_11mm);
  }
}

/**
The scalar series.

Caution: keep the sampling at 0.01 s. At 0.1 s the Nyquist limit is 5 Hz, the
flame flickers near 17 Hz, and everything above 5 Hz folds into the band below
it. That aliasing is what made the archived run appear to wander at about
1 Hz. */

event output (t += 0.01) {
  double xH2O[5], Tavg[5];
  double sample_points[5] = {H0/2 + 2e-3, H0/2 + 4e-3, H0/2 + 8e-3,
                             H0/2 + 11e-3, H0/2 + 15e-3};

  for (int ii = 0; ii < 5; ii++) {
    path_average_H2O (sample_points[ii], &xH2O[ii], &Tavg[ii]);

    /**
    An empty path gives a weighted mean of nothing. Report the ambient
    temperature, as the experiment does before the plume arrives. */

    if (xH2O[ii] < 1e-4)
      Tavg[ii] = TG0;
  }

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();

  /**
  `statsf` is collective, so it runs on every rank, outside the guard. */

  stats sT = statsf (T);

  if (pid() != 0)
    return 0;

  static FILE *fpxH2O = open_profile ("xH2OProfile.dat");

  static FILE * fpT = open_profile ("TemperatureProfile.dat");
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = open_profile (name);

  fprintf (fpxH2O, "%g %g %g %g %g %g\n", t, xH2O[0], xH2O[1], xH2O[2],
           xH2O[3], xH2O[4]);
  fflush (fpxH2O);

  fprintf (fpT, "%g %g %g %g %g %g\n", t, Tavg[0], Tavg[1], Tavg[2], Tavg[3],
           Tavg[4]);
  fflush (fpT);

  fprintf (fp, "%g %g %g\n", t, solid_mass/solid_mass0, sT.max);
  fflush (fp);
}

/**
## The probe of the timestep limits

`DT_PROBE` writes `dtlimits.dat` at every step. Use it to find which velocity
makes the timestep collapse at ignition. Three limits set `dt`:

- `dt_uf`, the CFL limit of the flow velocity `uf` (`centered-phasechange.h`),
- `dt_ubf`, the CFL limit of the shrinkage velocity `ubf` (`shrinking.h`),
- `dt_corr`, the limit of the corrective velocity (`multicomponent-varprop.h`).

`bind` gives the smallest of them: 0 is `DT`, 1 is `uf`, 2 is `ubf` and 3 is
the corrective velocity. `dt` can be smaller than all four, because
`timestep()` lets the step grow only slowly.

For `uf` and `ubf`, the file also gives the face with the largest
`|u|/(fm*Delta)`. It gives the position, the level, and `f`, `T`, `TG`,
`omega` and `gas_source` of the two cells on each side of that face. `L` is
the cell on the left or bottom side, `R` is the cell on the right or top side.
`TG` is the raw field, not the gas temperature. `omega` is the value of the
previous step, because the `chemistry` event comes after `stability`.

Caution: do not call `timestep()` here. It keeps a static `previous` value,
and a second call changes the step of the solver. This probe repeats its
loop instead.

This `stability` event is declared after the headers, so it runs before
their `stability` events. The fields are thus the ones that set `dt`. The
`vof` event writes the line, because `dt` is final at that point.

### The neighbourhood of the fastest `uf` face

Columns 37 to 56 describe the neighbourhood of the fastest `uf` face. Use them
to tell a checkerboard mode from a smooth jet.

- `dir` is 0 for a face normal to x (axial) and 1 for a face normal to y.
- `u0` is the signed velocity on the face. `um1` and `up1` are the velocities
  on the previous and the next face along the normal. `ub` and `ut` are the
  velocities on the two faces beside it, below and above for an x face. All
  five are `uf/fm`, thus physical velocities in m/s. In a checkerboard mode
  the signs alternate. In a jet they do not.
- `rho` is `rhov/cm`, the density that the flow solver uses.
- `drhodt` and `divsrc` are the expansion term and the full source of the
  projection (`div_source`), per unit volume. `divuf` is the divergence of
  `uf` per unit volume. After the projection `divuf + divsrc` is close to 0.
  Caution: this probe reads them in `stability`, thus after the projection of
  the previous step, and `drhodt` and `divsrc` belong to that step.
- `p` and `pf` are the two pressures of the solver.

All cell values come as a pair: `L` first, `R` second.

### Snapshots

`DTP_DUMP_T1` and `DTP_DUMP_T2` give two times for a `dump()`. The files are
`dtprobe-t<time>`. A third `dump()`, `dtprobe-dt`, happens once when `dt`
falls below `DTP_DUMP_DT`. It catches the late stage of the collapse, because
the time of the crash changes from one restart to the next. These files also
contain `p` and `pf`. They do not replace `last-snapshot`. To look at one,
restore it in a short program and write a VTK file. */

#ifndef DT_PROBE
# define DT_PROBE 1
#endif

#if DT_PROBE
#define DTP_NCELL 5  // f, T, TG, omega, gas_source
#define DTP_NFACE (4 + 2*DTP_NCELL) // rate, x, y, level, L cells, R cells

#define DTP_NNEAR 20 // dir, 5 velocities, 7 cell pairs
#define DTP_NOUT 5    // outlet: min u.x, min uf/fm, inflow, outflow, inlet

static double dtp_out[DTP_NOUT];

/**
The outlet columns. `umin` is the smallest `u.x` in the cells next to the
outlet. `ufmin` is the smallest `uf.x/fm.x` on the faces of the outlet. `Qin`
and `Qout` are the sums of `uf.x*Delta` over the faces of the outlet where
the gas flows in and out: `uf` carries the metric, so this is the volume flux
per radian. `Qinlet` is the same sum on the inlet, for reference. A value of
`Qin` that grows toward `Qinlet` or beyond is the backflow. */

static void dtp_outlet (double * out)
{
  double umin = HUGE, ufmin = HUGE;
  double flux_in = 0., flux_out = 0., flux_inlet = 0.;
  double xr = X0 + L0, xl = X0;
  foreach_face (x, reduction(min:umin) reduction(min:ufmin)
                reduction(+:flux_in) reduction(+:flux_out)
                reduction(+:flux_inlet)) {
    if (x > xr - 1e-6*L0) {
      umin = min (umin, u.x[-1]);
      if (fm.x[] > 0.)
        ufmin = min (ufmin, uf.x[]/fm.x[]);
      if (uf.x[] < 0.)
        flux_in += uf.x[]*Delta;
      else
        flux_out += uf.x[]*Delta;
    }
    else if (x < xl + 1e-6*L0)
      flux_inlet += uf.x[]*Delta;
  }
  double v[DTP_NOUT] = {umin, ufmin, flux_in, flux_out, flux_inlet};
  for (int k = 0; k < DTP_NOUT; k++)
    out[k] = v[k];
}

static double dtp_uf[DTP_NFACE], dtp_ubf[DTP_NFACE];
static double dtp_near[DTP_NNEAR];
static double dtp_corr = HUGE, dtp_DT = HUGE;

/**
Find the face with the largest `|u|/(fm*Delta)`, and fill `out` with the
data of that face. Every rank calls it, because the loops are collective. */

static void dtp_scan (face vector u, double * out)
{
  double rate = 0.;
  foreach_face (reduction(max:rate))
    if (u.x[] != 0. && fm.x[] > 0.)
      rate = max (rate, fabs (u.x[])/(fm.x[]*Delta));

  /**
  Every rank sets a value only at the face that has the maximum rate. The
  others keep `-HUGE`, so the `max` reduction returns the value of that
  face. */

  double xf = -HUGE, yf = -HUGE, lev = -HUGE;
  double fL = -HUGE, fR = -HUGE, TL = -HUGE, TR = -HUGE;
  double TGL = -HUGE, TGR = -HUGE, oL = -HUGE, oR = -HUGE;
  double sL = -HUGE, sR = -HUGE;
  if (rate > 0.)
    foreach_face (reduction(max:xf) reduction(max:yf) reduction(max:lev)
                  reduction(max:fL) reduction(max:fR)
                  reduction(max:TL) reduction(max:TR)
                  reduction(max:TGL) reduction(max:TGR)
                  reduction(max:oL) reduction(max:oR)
                  reduction(max:sL) reduction(max:sR))
      if (u.x[] != 0. && fm.x[] > 0. &&
          fabs (u.x[])/(fm.x[]*Delta) >= rate) {
        xf = x; yf = y; lev = level;
        fL = f[-1];           fR = f[];
        TL = T[-1];           TR = T[];
        TGL = TG[-1];         TGR = TG[];
        oL = omega[-1];       oR = omega[];
        sL = gas_source[-1];  sR = gas_source[];
      }

  double v[DTP_NFACE] = {rate, xf, yf, lev,
                         fL, TL, TGL, oL, sL,
                         fR, TR, TGR, oR, sR};
  for (int k = 0; k < DTP_NFACE; k++)
    out[k] = v[k];
}

/**
Fill `out` with the neighbourhood of the face of `u` that has the rate `rate`.
The array reduction uses the same `-HUGE` method as `dtp_scan()`. */

static inline double dtp_vel (double uff, double fmf)
{
  return fmf > 0. ? uff/fmf : 0.;
}

static void dtp_neighbourhood (face vector u, double rate, double * out)
{
  double nb[DTP_NNEAR];
  for (int k = 0; k < DTP_NNEAR; k++)
    nb[k] = -HUGE;
  int uxi = u.x.i; // not rotated: the loop compares it with the rotated u.x

  if (rate > 0.)
    foreach_face (reduction(max:nb[:DTP_NNEAR]))
      if (u.x[] != 0. && fm.x[] > 0. &&
          fabs (u.x[])/(fm.x[]*Delta) >= rate) {
        nb[0] = (u.x.i == uxi) ? 0. : 1.;
        nb[1] = dtp_vel (u.x[-1], fm.x[-1]);
        nb[2] = dtp_vel (u.x[], fm.x[]);
        nb[3] = dtp_vel (u.x[1], fm.x[1]);
        nb[4] = dtp_vel (u.x[0,-1], fm.x[0,-1]);
        nb[5] = dtp_vel (u.x[0,1], fm.x[0,1]);
        nb[6] = cm[-1] > 0. ? rhov[-1]/cm[-1] : 0.;
        nb[7] = cm[] > 0. ? rhov[]/cm[] : 0.;
        nb[8] = cm[-1] > 0. ? drhodt[-1]/cm[-1] : 0.;
        nb[9] = cm[] > 0. ? drhodt[]/cm[] : 0.;
        nb[10] = cm[-1] > 0. ? div_source[-1]/cm[-1] : 0.;
        nb[11] = cm[] > 0. ? div_source[]/cm[] : 0.;
        nb[12] = cm[-1] > 0. ?
          (u.x[] - u.x[-1] + u.y[-1,1] - u.y[-1])/(Delta*cm[-1]) : 0.;
        nb[13] = cm[] > 0. ?
          (u.x[1] - u.x[] + u.y[0,1] - u.y[])/(Delta*cm[]) : 0.;
        nb[14] = p[-1];   nb[15] = p[];
        nb[16] = pf[-1];  nb[17] = pf[];
        nb[18] = porosity[-1]; nb[19] = porosity[];
      }

  for (int k = 0; k < DTP_NNEAR; k++)
    out[k] = nb[k];
}

event stability (i++) {
  dtp_DT = dtmax;
  dtp_scan (uf, dtp_uf);
  dtp_neighbourhood (uf, dtp_uf[0], dtp_near);
  dtp_outlet (dtp_out);
  dtp_scan (ubf, dtp_ubf);
#ifdef FICK_CORRECTED
  dtp_corr = (CORRECTIVE_CFL > 0. && corrective_uodx > 0.) ?
    CORRECTIVE_CFL/corrective_uodx : HUGE;
#endif
}

static void dtp_print_face (FILE * fp, const double * d)
{
  for (int k = 0; k < DTP_NFACE; k++)
    fprintf (fp, " %g", d[k]);
}

event vof (i++) {
  if (pid() != 0)
    return 0;

  static FILE * fp = NULL;
  if (!fp) {
    fp = open_profile ("dtlimits.dat");

    /**
    A restart appends to the file. The file can be new all the same, for
    example after a restart in a new folder. So write the header when the
    file is empty, not when `restarted` is 0. */

    fseek (fp, 0, SEEK_END);
    if (ftell (fp) == 0)
      fprintf (fp, "#t(1) i(2) dt(3) DT(4) dt_uf(5) dt_ubf(6) dt_corr(7)"
               " bind(8)"
               " uf: rate(9) x(10) y(11) level(12)"
               " fL(13) TL(14) TGL(15) omegaL(16) gsL(17)"
               " fR(18) TR(19) TGR(20) omegaR(21) gsR(22)"
               " ubf: rate(23) x(24) y(25) level(26)"
               " fL(27) TL(28) TGL(29) omegaL(30) gsL(31)"
               " fR(32) TR(33) TGR(34) omegaR(35) gsR(36)"
               " ufnear: dir(37) um1(38) u0(39) up1(40) ub(41) ut(42)"
               " rhoL(43) rhoR(44) drhodtL(45) drhodtR(46)"
               " divsrcL(47) divsrcR(48) divufL(49) divufR(50)"
               " pL(51) pR(52) pfL(53) pfR(54) porL(55) porR(56)"
               " outlet: umin(57) ufmin(58) Qin(59) Qout(60) Qinlet(61)\n");
  }

  double lim[4] = {
    dtp_DT,
    dtp_uf[0] > 0. ? CFL/dtp_uf[0] : HUGE,
    dtp_ubf[0] > 0. ? CFL/dtp_ubf[0] : HUGE,
    dtp_corr
  };
  int bind = 0;
  for (int k = 1; k < 4; k++)
    if (lim[k] < lim[bind])
      bind = k;

  fprintf (fp, "%g %d %g %g %g %g %g %d", t, i, dt,
           lim[0], lim[1], lim[2], lim[3], bind);
  dtp_print_face (fp, dtp_uf);
  dtp_print_face (fp, dtp_ubf);
  for (int k = 0; k < DTP_NNEAR; k++)
    fprintf (fp, " %g", dtp_near[k]);
  for (int k = 0; k < DTP_NOUT; k++)
    fprintf (fp, " %g", dtp_out[k]);
  fputc ('\n', fp);
  fflush (fp);
  return 0;
}

#ifndef DTP_DUMP_T1
# define DTP_DUMP_T1 10.235
#endif
#ifndef DTP_DUMP_T2
# define DTP_DUMP_T2 10.244
#endif
#ifndef DTP_DUMP_DT
# define DTP_DUMP_DT 1e-5
#endif

/**
Caution: an event at a fixed time makes `dtnext()` shorten the step before
it. After a restart, the run thus follows a slightly different path than a
run without these events. */

/**
`centered.h` sets `nodump` on `p` and `pf`, so a normal `dump()` does not write
them. `dtp_dump()` writes them in the probe files only. It sets the flag again
after the dump, so `last-snapshot` does not change. */

static void dtp_dump (const char * name)
{
  p.nodump = pf.nodump = false;
  dump (name);
  p.nodump = pf.nodump = true;
}

event dtprobe_dump (t = {DTP_DUMP_T1, DTP_DUMP_T2}) {
  char name[80];
  sprintf (name, "dtprobe-t%g", t);
  dtp_dump (name);
}

event dtprobe_dump_dt (i++) {
  static bool done = false;
  if (!done && i > 0 && dt < DTP_DUMP_DT) {
    dtp_dump ("dtprobe-dt");
    done = true;
  }
}
#endif // DT_PROBE

#if TREE
event adapt (i++) {
  scalar oxidiser = YGList_G[OpenSMOKE_IndexOfSpecies ("O2")];

  scalar zdiff[];
  foreach()
    zdiff[] = zmix[] - zsto[];

  adapt_wavelet_leave_interface ({T, oxidiser, zdiff}, {f},
    (double[]){5e0, 1e-2, 1e-2}, maxlevel, minlevel, 2);

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
  squares ("T", min = 300, max = 2200, spread = -1, linear = true);
  isoline ("zmix - zsto", lw = 1.5, lc = {1., 1., 1.});
  draw_vof ("f", lw = 1.5);
  mirror ({0, 1}) {
    squares ("XH2O", min = 0, max = 1., spread = -1, linear = true);
    isoline ("zmix - zsto", lw = 1.5, lc = {1., 1., 1.});
    draw_vof ("f", lw = 1.5);
  }
  save ("movie.mp4");
}

#if VTK_OUTPUT
event vtk (t += 5; t <= 80) {

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

#if PRINT_FIELDS
event save_fields (t += 10) {
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

  // CO
  scalar XCO_G = XGList_G[OpenSMOKE_IndexOfSpecies ("CO")];
  scalar XCO_S = XGList_S[OpenSMOKE_IndexOfSpecies ("CO")];
  scalar XCO[];
  foreach()
    XCO[] = XCO_S[]*f[] + XCO_G[]*(1. - f[]);

  // LVG
  scalar LVG_G = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  scalar LVG_S = YGList_S[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  scalar LVG[];
  foreach()
    LVG[] = LVG_S[]*f[] + LVG_G[]*(1. - f[]);

  char name[80];
  sprintf (name, "fields-%g", t);
  FILE * fs = open_profile (name);
  output_field ({T, f, u.x, u.y, XH2O, XCO2, XOH, XCO, LVG, omega}, fs);
  fclose (fs);
}
#endif

event dump (t = 1; t += 1) {
  dump ("last-snapshot");
}

/**
Caution: `return 1` is what ends the run, not the time of the event.

`events()` in `$BASILISK/grid/events.h` keeps the loop alive while any event
still has a condition and is still alive. `event print_profile` carries
`t <= PROFILE_TEND`, so it holds the loop open, and a bare `event stop
(t = tend)` cannot end a run with `tend` below `PROFILE_TEND`. A run at
`TEND = 0.25` went on to t = 0.29 before this line existed. An action that
returns a non-zero value gives `event_stop`, which ends the run whatever the
other events ask for. */

event stop (t = tend) {
  if (pid() == 0) {
    fclose (fTprofile_2mm);
    fclose (fTprofile_11mm);
    fclose (fxH2Oprofile_2mm);
    fclose (fxH2Oprofile_11mm);
    fclose (fxOHprofile_2mm);
    fclose (fxOHprofile_11mm);
  }
}

/** 
~~~gnuplot Mass
set terminal svg size 450, 450
set output "mass-fatehi.svg"
set xlabel "Time [s]"
set ylabel "Normalized solid mass [-]"
set xrange [0:100]
set yrange [0:1.1]

set size square
error_margin = 0.025

plot "../../data/fatehi/mass" u ($1-4):2:(error_margin) w yerrorbars pt 64 ps 1 lc "black" notitle, \
     "cluster/temp/OutputData-10" u 1:2 w l lw 2 lc "black" notitle, \
     #"../../data/fatehi/mass" u 1:2 w p pt 64 ps 1 lc "black" notitle", \
~~~

~~~gnuplot H2O mole fraction
reset
set terminal svg size 450, 450
set output "xH2O-fatehi.svg"
set xlabel "Time [s]"
set ylabel "H2O Mole Fraction [-]"
set yrange [0.:0.5]
set xrange [0:70]

#folder = "cluster/new/lu/fatehi-combustion/"
folder = "cluster/temp/"

error_margin = 0.03
shift = 6

plot  "../../data/fatehi/yH2O-11mm" every 2 u 1:2:(error_margin) w yerrorbars pt 64 ps 0.8 lc "black" notitle  ,\
      "../../data/fatehi/yH2O-2mm" every 2 u 1:2:(error_margin) w yerrorbars pt 64 ps 0.8 lc "red" notitle  ,\
      sprintf("%s%s", folder, "results_2mm/effective_values.dat")  every 5 u ($1-shift):($5/100)   w l lw 2 lc "red" title "2 mm", \
      sprintf("%s%s", folder, "results_11mm/effective_values.dat") every 5 u ($1-shift):($5/100)  w l lw 2 lc "black" title "11 mm" 
      #"../../data/fatehi/yH2O-2mm" every 2 u 1:2 w lp pt 64 ps 0.8 lc "red" notitle  ,\
      #"../../data/fatehi/yH2O-11mm" every 2 u 1:2 w lp pt 64 ps 0.8 lc "black" notitle  ,\
~~~

~~~gnuplot Temperature
reset
set terminal svg size 450, 450
set output "temperature-fatehi.svg"
set xlabel "Time [s]"
set ylabel "Temperature [K]"
set yrange [800:2000]
set xrange [0:60]

folder = "cluster/temp/"

error_margin = 50
shift = 6

plot  "../../data/fatehi/T-11mm" every 2 u 1:2:(error_margin) w yerrorbars pt 4 ps 0.8 lw 1 lc "black" notitle, \
      "../../data/fatehi/T-2mm" every 2 u 1:2:(error_margin) w yerrorbars pt 4 ps 0.8 lw 1 lc "red" notitle, \
      sprintf("%s%s", folder, "results_2mm/effective_values.dat")  every 5 u ($1-shift):2  w l lw 2 lc "red" title "2 mm", \
      sprintf("%s%s", folder, "results_11mm/effective_values.dat") every 5 u ($1-shift):2  w l lw 2 lc "black" title "11 mm", \
      #sprintf("%s%s", folder, "TemperatureProfile.dat") u ($1-shift):2 w l lw 2 lc "light-coral" title "2 mm", \
      #sprintf("%s%s", folder, "TemperatureProfile.dat") u ($1-shift):5 w l lw 2 lc "gray" title "11 mm", \
      #"../../data/fatehi/T-2mm" every 2 u 1:2 w lp pt 4 ps 0.8 lw 1 lc "red" notitle, \
      #"../../data/fatehi/T-11mm" every 2 u 1:2 w lp pt 4 ps 0.8 lw 1 lc "black" notitle, \
~~~

~~~gnuplot Temperature evolution
reset
set terminal svg size 550, 550 
set output "temperature-evolution-fatehi-2mm.svg"
load "/home/rcaraccio/gnuplot-palettes/inferno.pal"

set xlabel "y position"
set ylabel "Temperature"
set key right top
set ytics 900,200,2300
set yrange [800:2300]
set xrange [-0.06:0.06]
set title "Temperature 2mm"

end_time = 100 
step_time = 20

set cbrange [0:end_time]
unset colorbox

folder = "cluster/temp/"

plot for [t=0:end_time:step_time] sprintf("%s%s", folder, "T_profile_2mm.dat") u ($1==t ? $2 : 1/0):3:(t) w l lw 3 lc palette title sprintf("%d s", t), \
     for [t=0:end_time:step_time] sprintf("%s%s", folder, "T_profile_2mm.dat") u ($1==t ? -$2 : 1/0):3:(t) w l lw 3 lc palette notitle
~~~

~~~gnuplot Temperature evolution
reset
set terminal svg size 550, 550
set output "temperature-evolution-fatehi-11mm.svg"
load "/home/rcaraccio/gnuplot-palettes/inferno.pal"

set xlabel "Axial distance [m]"
set ylabel "Temperature [K]"
set key right top
set ytics 900,200,2300
set yrange [950:2500]
set xrange [-0.06:0.06]
set title "Temperature 11mm"

folder = "cluster/latest/"

end_time = 100
step_time = 20

set cbrange [0:end_time]
unset colorbox

plot for [t=0:end_time:step_time] sprintf("%s%s", folder, "T_profile_11mm.dat") u ($1==t ? $2 : 1/0):3:(t) w l lw 3 lc palette title sprintf("%d s", t), \
     for [t=0:end_time:step_time] sprintf("%s%s", folder, "T_profile_11mm.dat") u ($1==t ? -$2 : 1/0):3:(t) w l lw 3 lc palette notitle
~~~

~~~gnuplot xH2O evolution
reset
set terminal svg size 550, 550
set output "xH2O-evolution-fatehi-2mm.svg"
load "/home/rcaraccio/gnuplot-palettes/inferno.pal"
set xlabel "y position"
set ylabel "H2O Mole Fraction"
set key right top
set yrange [0:0.5]
set xrange [-0.06:0.06]
set title "xH2O 2mm"

end_time = 100
step_time = 20

set cbrange [0:end_time]
unset colorbox 

plot for [t=0:end_time:step_time] "cluster/temp/xH2O_profile_2mm.dat" u ($1==t ? $2 : 1/0):3:(t) w l lw 2 lc palette title sprintf("%d s", t), \
     for [t=0:end_time:step_time] "cluster/temp/xH2O_profile_2mm.dat" u ($1==t ? -$2 : 1/0):3:(t) w l lw 2 lc palette notitle

~~~

~~~gnuplot xH2O evolution
reset
set terminal svg size 550, 550
set output "xH2O-evolution-fatehi-11mm.svg"
load "/home/rcaraccio/gnuplot-palettes/inferno.pal"
set xlabel "y position"
set ylabel "H2O Mole Fraction"
set key right top
set yrange [0:0.5]
set xrange [-0.06:0.06]
set title "xH2O 11mm"

end_time = 100
step_time = 20

set cbrange [0:end_time]
unset colorbox 

folder = "cluster/temp/"

plot for [t=0:end_time:step_time] sprintf("%s%s", folder, "xH2O_profile_11mm.dat") u ($1==t ? $2 : 1/0):3:(t)   w l lw 2 lc palette title sprintf("%d s", t), \
     for [t=0:end_time:step_time] sprintf("%s%s", folder, "xH2O_profile_11mm.dat") u ($1==t ? -$2 : 1/0):3:(t)  w l lw 2 lc palette notitle

~~~

~~~gnuplot Test plot
reset
set terminal svg size 450, 450
set output "test-plot-fatehi.svg"
set xlabel "Time [s]"
set ylabel "H2O Mole Fraction[-]"
set yrange [0.:0.5]
set xrange [0:60]

endtime = 70
step = 1
n = endtime/step + 1
error_margin = 0.03
shift = 6

folder_path = "cluster/temp/"

# Define the upper space limit
space_max = 0.015

# ---- 2mm case ----
system(sprintf("awk '$2 <= %f { sum[$1]+=$3; count[$1]++; if(!($1 in max) || $3>max[$1]) max[$1]=$3 } END { for(t in sum) print t, sum[t]/count[t], max[t] }' %sxH2O_profile_2mm.dat | sort -n > %savg_h2o_2mm.dat", space_max, folder_path, folder_path))
system(sprintf("awk '$2 <= %f { sum[$1]+=$3; count[$1]++; if(!($1 in max) || $3>max[$1]) max[$1]=$3 } END { for(t in sum) print t, sum[t]/count[t], max[t] }' %sxH2O_profile_11mm.dat | sort -n > %savg_h2o_11mm.dat", space_max, folder_path, folder_path))


plot  "../../data/fatehi/yH2O-11mm" every 2 u 1:2:(error_margin) w yerrorbars pt 64 ps 0.8 lc "black" notitle  ,\
      "../../data/fatehi/yH2O-2mm" every 2 u 1:2:(error_margin) w yerrorbars pt 64 ps 0.8 lc "red" notitle  ,\
      sprintf("%s%s", folder, "results_2mm/effective_values.dat")  every 5 u ($1-shift):($7/100)   w l lw 2 lc "red" title "2 mm", \
      sprintf("%s%s", folder, "results_11mm/effective_values.dat") every 5 u ($1-shift):($7/100)  w l lw 2 lc "black" title "11 mm" 
      #sprintf("%savg_h2o_11mm.dat", folder_path) u ($1-shift):3 w l lw 2 lc "black" t "max 11mm", \
      #sprintf("%savg_h2o_2mm.dat", folder_path) u ($1-shift):3 w l lw 2 lc "red" t "max 2mm", \
      #"../../data/fatehi/yH2O-2mm" every 2 u 1:2 w lp pt 64 ps 0.8 lc "red" notitle  ,\
      #"../../data/fatehi/yH2O-11mm" every 2 u 1:2 w lp pt 64 ps 0.8 lc "black" notitle  ,\

~~~

~~~gnuplot H2O mole fraction
reset
set terminal svg size 550, 350
set output "xH2O-cd.svg"
set xlabel "Time [s]"
set ylabel "Column Desnity [% m]"
set xrange [0:60]
set yrange [0.:]

folder = "cluster/temp/"
length = 0.015

error_margin = 0.03*length*100

shift = 5

plot  "../../data/fatehi/yH2O-11mm" every 2 u 1:($2*length*100):(error_margin) w yerrorbars pt 64 ps 0.8 lc "black" notitle  ,\
      "../../data/fatehi/yH2O-2mm" every 2 u 1:($2*length*100):(error_margin) w yerrorbars pt 64 ps 0.8 lc "red" notitle  ,\
      "../../data/fatehi/yH2O-2mm" every 2 u 1:($2*length*100) w lp pt 64 ps 0.8 lc "red" notitle  ,\
      "../../data/fatehi/yH2O-11mm" every 2 u 1:($2*length*100) w lp pt 64 ps 0.8 lc "black" notitle  ,\
      sprintf("%s%s", folder, "results_2mm/cd_xH2O_2mm.dat") u ($1-shift):2 w l lw 2 lc "red" title "2 mm", \
      sprintf("%s%s", folder, "results_11mm/cd_xH2O_11mm.dat") u ($1-shift):2 w l lw 2 lc "black" title "11 mm" 
~~~

**/
