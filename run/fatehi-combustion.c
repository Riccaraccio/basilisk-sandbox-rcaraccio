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

| setting | value |
|---|---|
| domain | `L0 = 20 D0`, axisymmetric, pellet of 8 mm (superquadric of exponent 20) |
| grid | `MAXLEVEL` 11, minimum level 2, adapt on `{T, O2}` with `{5, 1e-2}` |
| timestep | cap 5e-4 s, `CFL` 0.5 |
| projection | `TOLERANCE` 1e-5, `NITERMIN` 2 |
| gas chemistry | Strang split and the frozen-cell gate (`chemistry.h`) |
| expansion source | exact form, filtered in 4 passes of width `Delta_min` |
| shrinkage | `ZETA_REACTION`: it follows the local rate of reaction |
| inlet, top | air at 1123 K, 0.13 m/s |
| outlet | Neumann `u`, `TG`, `YG`; Dirichlet `p`, `pf` |

Notes on the settings:

- `TOLERANCE` and `NITERMIN`: Basilisk scales the tolerance of the projection
  as `TOLERANCE/dt^2`, so the default 1e-3 stops the solve after one
  multigrid cycle at every step.
- One cell holds 1104 fields with `biomass/Solid-gas-88`, which is 8.8 kB.
  `main()` starts on level 8 and `event init` refines a disc near the
  pellet. A uniform grid at level 10 allocates 9.3 GB before the first adapt.
- `event output` samples at 0.01 s. At 0.1 s the Nyquist limit is 5 Hz, and
  the flicker of the flame near 17 Hz folds into the band below 1 Hz.

Known limits:

- Nothing bounds `TG` in a gas sliver. The output writes `TGmin_gas` and
  `nTGneg` as the early warning. Read the caution at the interface heat
  source in `multicomponent-varprop.h`.
- The expansion source carries a spatial filter
  (`navier-stokes/centered-phasechange.h`). State it, with its width
  `sigma = Delta_min`, in a publication.
- `shift_field()` makes the solid mass drift by about -0.4 per cent of
  `solid_mass0` at level 11. It is not fixed.

`GAS_PHASE_REACTIONS` does nothing. Use `TURN_OFF_GAS_REACTIONS` to switch off
the gas kinetics. */

#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1

/**
The knobs below are overridable with `-D`.

`FATEHI_DUMMY = 1` selects `biomass/dummy-solid-gas` and a solid of BIOMASS,
MOIST and ASH. It is for a quick check of a build. That mechanism has no OH,
so the OH outputs are off with it. */

#ifndef FATEHI_DUMMY
# define FATEHI_DUMMY 0
#endif

#if FATEHI_DUMMY
# define FATEHI_KINFOLDER "biomass/dummy-solid-gas"
#else
# define FATEHI_KINFOLDER "biomass/Solid-gas-88"
#endif

/**
`SNAPSHOT_EVERY` writes a numbered snapshot `snapshot-<t>` every that many
seconds of physical time, beside `last-snapshot`. It is an integer, so that
the preprocessor can test it. Set it to 0 to turn the snapshots off.

Caution: a snapshot at level 11 with `biomass/Solid-gas-88` is large. Check
the free space of the disk before a run of 150 s. */

#ifndef SNAPSHOT_EVERY
# define SNAPSHOT_EVERY 5
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
#include "vtk.h"

const double Uin = 0.13; //inlet velocity
/**
Caution: give `pf` the same conditions as `p`. Basilisk does not copy the
conditions of `p` to `pf`. Without the lines for `pf`, every boundary of `pf`
is Neumann. The projection of `uf` then has no solution when the expansion
term `drhodt` is not zero, because the outflow flux is fixed. `pf` drifts
without limit, the divergence error grows, and at ignition `dt` collapses. */

u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
pf[left]     = neumann (0.);
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

/**
## The outlet

Neumann for `u`, `TG` and `YG`, Dirichlet for `p` and `pf`.

Caution: a Neumann outlet lets gas flow back into the domain. The physics of
the heating phase needs it: while the moisture evaporates, a cold plume sinks
from the pellet toward the inlet, and gas must replace it. But runs with
other settings saw the backflow on the axis grow by a constant factor until
`dt` collapsed. If the inflow at the outlet grows in this way, stop the run.
The tag `oscillation-campaign-2026-09` keeps two other outlets (`OUTLET_BC`
1 and 2). */

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);
psi[right]    = neumann (0.);

const double tend = TEND; //simulation time
int maxlevel = MAXLEVEL; int minlevel = 2;

double D0 = 8e-3, H0 = 8e-3;
double solid_mass0 = 0.;

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

int main() {
  lambdaSmodel = L_TENWOLDE;
  TS0 = 300.; TG0 = 1123.;
  rhoS = 1550;
  eps0 = 0.2; // low, compressed pellet

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_REACTION;

  DT = 5e-4; // the CFL usually binds below this cap

  G.x = -9.81;

  kinfolder = FATEHI_KINFOLDER;
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

  /**
  Read this line before you quote a run. It prints what the build enabled,
  not what the directory name promises. Under MPI every rank shares stderr,
  so only rank 0 writes. */

  if (pid() == 0)
    fprintf (stderr, "# fatehi: kin=%s maxlevel=%d DT=%g Uin=%g tend=%g"
                     " tol=%g nitermin=%d filter=%d snapevery=%d nranks=%d\n",
             FATEHI_KINFOLDER, MAXLEVEL, DT, Uin, (double) TEND, TOLERANCE,
             NITERMIN, gas_source_filter_passes, SNAPSHOT_EVERY, npe());

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
  value belongs here. */

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
  fraction (f0, superquadric (x, y, 20, 0.5*H0, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;

#if FATEHI_DUMMY
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 0.935;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")]   = 0.061;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]     = 0.004;
#else
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.4344;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("GMSW")]  = 0.2108;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.1347;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0786;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0178;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0167;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0419;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0041;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.0610;
#endif

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
This function drops the sample instead. Two rows of an earlier run carried
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

  // Temperature profiles
  print_profile (T, H0/2 + 2e-3, fTprofile_2mm, t);
  print_profile (T, H0/2 + 11e-3, fTprofile_11mm, t);

  // Water vapor mole fraction profile
  print_profile (XH2O, H0/2 + 2e-3, fxH2Oprofile_2mm, t);
  print_profile (XH2O, H0/2 + 11e-3, fxH2Oprofile_11mm, t);

#if !FATEHI_DUMMY
  scalar XOH = XGList_G[OpenSMOKE_IndexOfSpecies ("OH")];
  print_profile (XOH, H0/2 + 2e-3, fxOHprofile_2mm, t);
  print_profile (XOH, H0/2 + 11e-3, fxOHprofile_11mm, t);
#endif

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
it. That aliasing made an earlier run appear to wander at about
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
  The early warning of a runaway of the gas temperature in a thin cell.
  Nothing bounds `TG` in a sliver cell. The value of `TG` is `TG*(1-f)` or
  `TG` itself, which depends on the position within the step. So the event
  writes two quantities that do not depend on it:

  - `TGmin_gas`, the smallest `TG` over the cells with `f < F_ERR`, where the
    two forms are equal,
  - `nTGneg`, the number of cells with `TG < 0`, because both forms have the
    same sign.

  A runaway starts as a negative `TG` in a cell with a gas fraction of about
  1e-4, and it reaches the pure gas cells after a few steps. If `nTGneg` is
  not 0, or `TGmin_gas` falls below about 250 K, stop the run. */

  double TGmin_gas = HUGE, nTGneg = 0.;
  foreach (reduction(min:TGmin_gas) reduction(+:nTGneg)) {
    if (f[] < F_ERR)
      TGmin_gas = min (TGmin_gas, TG[]);
    if (TG[] < 0.)
      nTGneg += 1.;
  }

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

  /**
  Columns of `OutputData-<maxlevel>`: t(1) Ms/Ms0(2) Tmax(3) TGmin_gas(4)
  nTGneg(5) dt(6). */

  fprintf (fp, "%g %g %g %g %g %g\n", t, solid_mass/solid_mass0, sT.max,
           TGmin_gas, nTGneg, dt);
  fflush (fp);
}

/**
The adapt criterion is `T` and the oxidiser. Do not adapt on `zmix - zsto`:
`flame.h` updates `zmix` and `zsto` only every `FLAME_PRINT_TIME`, so the
grid would follow a stale field between two updates. */

#if TREE
event adapt (i++) {
  scalar oxidiser = YGList_G[OpenSMOKE_IndexOfSpecies ("O2")];

  adapt_wavelet_leave_interface ({T, oxidiser}, {f},
    (double[]){5e0, 1e-2}, maxlevel, minlevel, 2);

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

  // OH. The dummy mechanism has no OH, so the column is then 0.
  scalar XOH[];
#if FATEHI_DUMMY
  foreach()
    XOH[] = 0.;
#else
  scalar XOH_G = XGList_G[OpenSMOKE_IndexOfSpecies ("OH")];
  scalar XOH_S = XGList_S[OpenSMOKE_IndexOfSpecies ("OH")];
  foreach()
    XOH[] = XOH_S[]*f[] + XOH_G[]*(1. - f[]);
#endif

  output_vtk ({f, T, XH2O, XCO2, XOH, u.x, u.y, zdiff}, n=(1<<maxlevel), fp=fvtk, linear=true);
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
Caution: `return 1` is what ends the run, not the time of the event.

`events()` in `$BASILISK/grid/events.h` keeps the loop alive while any event
still has a condition and is still alive. `event print_profile` carries
`t <= PROFILE_TEND`, so it holds the loop open, and a bare `event stop
(t = tend)` cannot end a run with `tend` below `PROFILE_TEND`. An action that
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
  return 1;
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
