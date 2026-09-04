/**
# The amplitude ladder: from the reduced case to `fatehi-combustion`

The slow oscillation of the release rate of the particle has the same
frequency in the reduced case and in the full case, but not the same
amplitude. At matched remaining mass the full case gives 3.5 to 4.1 per
cent and the reduced case gives 1.6 to 2.3 per cent. The probe temperature
differs by more: 23 K against 4.4 K. So the particle oscillates twice as
much, and the flame turns that oscillation into 2.5 times more temperature.

This case measures which ingredient carries the difference. The build with
no flags reproduces the reduced case `~/temp/expansion/test.c`. Each flag
adds one ingredient of `run/fatehi-combustion.c`.

| flag | reduced value | full value |
|---|---|---|
| `DA_VALUE` | 1e-12 | 1e-10 |
| `MOISTURE` | 0, pure biomass | 1, 6.1 % moisture and 0.4 % ash |
| `GRAVITY` | 0 | 1 |
| `SHAPE` | 0, sphere | 1, superquadric |
| `EMISSIVITY_DIBLASI` | 0, constant | 1, Di Blasi |
| the three transport flags | off | on |

`MAXLEVEL_VALUE` and `DT_VALUE` are separate: they test the grid and the
step, not an ingredient.

Compare the runs at matched remaining mass, never at matched time. The
cases burn at different rates, so the same instant is a different state.

Caution: `MOLAR_DIFFUSION`, `FICK_CORRECTED` and `MASS_DIFFUSION_ENTHALPY`
are existence-tested. Every site that this case reaches is `#ifdef`, 45 of
them in `src/`. A `#define MOLAR_DIFFUSION 0` therefore turns the flag ON.
The block below undefines them instead, so `-DMOLAR_DIFFUSION=0` means off,
and both `-DMOLAR_DIFFUSION` and `-DMOLAR_DIFFUSION=1` mean on. Do not
replace it with a `#ifndef ... 0` default. */

#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1

#if defined(MOLAR_DIFFUSION) && !MOLAR_DIFFUSION
# undef MOLAR_DIFFUSION
#endif

#if defined(FICK_CORRECTED) && !FICK_CORRECTED
# undef FICK_CORRECTED
#endif

#if defined(MASS_DIFFUSION_ENTHALPY) && !MASS_DIFFUSION_ENTHALPY
# undef MASS_DIFFUSION_ENTHALPY
#endif

/**
`TURN_OFF_HEAT_OF_REACTION` carries the same trap: `reactors.h` reads it with
`#ifdef` at lines 237 and 303, so `-DTURN_OFF_HEAT_OF_REACTION=0` would turn
it ON. Undefine it when the value is 0.

The flag keeps the kinetics of the species and zeroes the temperature source.
It therefore removes the endothermic arm of the loop and nothing else, which
is what separates the endothermic closure from the blowing closure. */

#if defined(TURN_OFF_HEAT_OF_REACTION) && !TURN_OFF_HEAT_OF_REACTION
# undef TURN_OFF_HEAT_OF_REACTION
#endif

/**
The log needs a value, and the four flags above carry none once they are
undefined. */

#ifdef TURN_OFF_HEAT_OF_REACTION
# define NOHEAT_ON 1
#else
# define NOHEAT_ON 0
#endif

#ifdef MOLAR_DIFFUSION
# define MOLAR_ON 1
#else
# define MOLAR_ON 0
#endif

#ifdef FICK_CORRECTED
# define FICK_ON 1
#else
# define FICK_ON 0
#endif

#ifdef MASS_DIFFUSION_ENTHALPY
# define MDE_ON 1
#else
# define MDE_ON 0
#endif

/**
The permeability of the reduced case. The full case leaves the default of
`darcy.h`, which is 1e-10, so the full particle is 100 times more
permeable. */

#ifndef DA_VALUE
# define DA_VALUE 1e-12
#endif

#ifndef MOISTURE
# define MOISTURE 0
#endif

#ifndef EMISSIVITY_DIBLASI
# define EMISSIVITY_DIBLASI 0
#endif

#ifndef GRAVITY
# define GRAVITY 0
#endif

#ifndef SHAPE
# define SHAPE 0 // 0 sphere, 1 superquadric
#endif

#ifndef MAXLEVEL_VALUE
# define MAXLEVEL_VALUE 10
#endif

/**
The reduced reference uses 5e-4. Keep this value, or the unflagged build
does not reproduce `~/temp/expansion/avg`. */

#ifndef DT_VALUE
# define DT_VALUE 5e-4
#endif

/**
## The timestep pair

`Tmax` is a single-valued function of `dt`: 16 to 29 K for each halving, with
R^2 = 0.89 in `test-shape`. `dt` itself is quantised, because `dtnext()` snaps
the step so that the run lands on the four `t += 0.01` events, so `dt` can
only take the values 0.01/n.

`CFLNUM` makes a constant step possible. With `CFLNUM = 2` the CFL never
binds, so `DT_VALUE` sets every step, and it still caps a runaway.

Caution: `DT_VALUE` must divide 0.01 exactly, or `dtnext()` subdivides the
step and `dt` is not constant. Use 2e-4 = 0.01/50, not 1.75e-4.

Caution: `CFL` is assigned in the `defaults` event of
`navier-stokes/centered.h`, which runs after `main()`. So this value is set in
`event init` below, never in `main()`. `TOLERANCE` has no such event and stays
in `main()`.

Caution: at ignition the peak velocity reaches about 2 m/s and the CFL binds
even at `CFLNUM = 2`. Branch the fixed-step runs from a plateau snapshot, and
check that column 2 of `expansion.dat` is constant before you quote them. */

#ifndef CFLNUM
# define CFLNUM 0.8
#endif

/**
## The blowing sweep

In the pyrolysis case the Stefan flow is not a perturbation. The blowing
ratio against the free stream is about 1.06, the Peclet number `v_w R/alpha`
is 2.5 to 4.4, and the blockage of the conductive flux moves by a factor of 4
over one cycle. `UIN_VALUE` changes the ratio, so it tests whether the loop
closes through the blowing.

`v_w` follows the rate of pyrolysis, not `Uin`, so `Uin` changes the ratio and
not the Peclet number. A larger `Uin` thins the layer and weakens the
blockage; a smaller `Uin` strengthens it.

Caution: `Uin` also changes the supply of heat, so it changes the rate of
burn. Compare at matched remaining mass, never at matched time. */

#ifndef UIN_VALUE
# define UIN_VALUE 0.13
#endif

/**
## Pyrolysis only

`PYROLYSIS_ONLY` removes the oxidiser and the reactions of the gas, and
reproduces `~/temp/10-fatehi/test.c`. The release still oscillates, 7 to 18
per cent peak-to-peak, so the flame does not close the loop. The case is much
cheaper and gives 6 cycles in 20 s, against 1.2 to 4.3 cycles in the runs
with a flame, so it is the correct case for the questions about the loop.

Both branches use `biomass/dummy-solid-gas`; see the note at `kinfolder`
below. The mechanism therefore holds O2 in either branch, but the
`PYROLYSIS_ONLY` build sets its mass fraction to zero everywhere. Keep the
lookups of O2 inside the `#else` branch anyway: they record which build wants
an oxidiser. */

#ifndef PYROLYSIS_ONLY
# define PYROLYSIS_ONLY 0
#endif

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h" 
#include "two-phase.h"

/**
`fatehi-combustion.c` puts `gravity.h` here, between `two-phase.h` and
`shrinking.h`. Same-name events run in reverse declaration order, so the
position of a module in this list changes when its events run. Keep the
order of the full case. */

#if GRAVITY
# include "gravity.h"
#endif

#include "superquadric.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"

//#include "flame.h"
#include "view.h"

const double Uin = UIN_VALUE; //inlet velocity
u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = neumann (0.);

double tend = 40;
int maxlevel = MAXLEVEL_VALUE, minlevel = 2;
double solid_mass0 = 0.;
double D0 = 8e-3, H0 = 8e-3;

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

int main() {

  /**
  Caution: under MPI every rank shares this stderr. Guard the message with
  `pid() == 0`, or the log carries one copy per rank. */

  if (pid() == 0)
    fprintf (stderr, "# ladder: MOLAR=%d FICK=%d MDE=%d MOISTURE=%d GRAVITY=%d"
                     " SHAPE=%d DIBLASI=%d Da=%g DT=%g maxlevel=%d"
                     " CFL=%g Uin=%g PYRO=%d NOHEAT=%d nranks=%d\n",
             MOLAR_ON, FICK_ON, MDE_ON, MOISTURE, GRAVITY, SHAPE,
             EMISSIVITY_DIBLASI, (double) DA_VALUE, (double) DT_VALUE,
             MAXLEVEL_VALUE, (double) CFLNUM, (double) UIN_VALUE,
             PYROLYSIS_ONLY, NOHEAT_ON, npe());

  lambdaSmodel = L_TENWOLDE;
  TS0 = 300.; TG0 = 1123.;
  rhoS = 1550; cpS = 1800;
  eps0 = 0.2;

  //rhoG    = 0.31;      // air ~1100 K, 1 atm
  //muG     = 4.5e-5;
  //lambdaG = 0.08;      cpG = 1200.;
  //lambdaS = 0.2;       cpS = 1500.;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_REACTION;

  DT = DT_VALUE;

  /**
  One mechanism for both branches.

  `biomass/dummy-solid` declares only `BIOMASS CHAR` as solid species. It has
  no `MOIST` and no `ASH`. `lambda_tenwolde` looks up both with the
  hard-error form (`solid-thermal-conductivity.h:171,173`), so a
  `PYROLYSIS_ONLY` build with that mechanism stopped at startup with "The
  requested species MOIST is not available".

  The oxidiser turns the gas phase off, not the mechanism. `PYROLYSIS_ONLY`
  sets `gas_start` to pure N2 and gives the boundaries the same value, so no
  reaction of the gas has an oxidiser. The case also defines no
  `GAS_PHASE_REACTIONS`, so `chemistry.h` integrates no gas kinetics at all.

  Two consequences. The mechanism carries 8 gas species against 3, so the
  species transport costs about 2.7 times more and a `PYROLYSIS_ONLY` run is
  slower than it was. And the solid mechanism adds the pair
  `MOIST => H2O` / `H2O => MOIST`, which stays active whatever `MOISTURE` is,
  so the pore water of pyrolysis can condense into `MOIST` even in a build
  that starts dry. Every flame case of the ladder already ran that way, so
  this makes the two branches consistent; but it means these runs are no
  longer comparable with `~/temp/{9,10,11}-fatehi`, which used
  `dummy-solid`. Compare new against new. */

  kinfolder = "biomass/dummy-solid-gas";
  shift_prod = true;

  L0 = 20*D0;
  origin (-L0/2, 0);

#if EMISSIVITY_DIBLASI
  emissivity = emissivity_diblasi;
#else
  emissivity = emissivity_constant;
#endif

  Da = (coord){DA_VALUE, DA_VALUE};

#if GRAVITY
  /**
  `gravity.h` declares `coord G = {0.,0.,0.}`. Without this line the header
  is present, the acceleration event runs, and the gravity is zero. */

  G.x = -9.81;
#endif

  init_grid(1 << min (maxlevel, 8));
  refine (circle (x, y, 4.*D0) > 0. && level < maxlevel);

  TOLERANCE = 1e-5;
  NITERMIN = 2;

  run();
}

double r0;
event init (i = 0) {
  scalar f0[];

  /**
  Caution: `navier-stokes/centered.h` assigns `CFL = 0.8` in its `defaults`
  event, which runs after `main()`. So the value belongs here. */

  CFL = CFLNUM;

#if SHAPE
  fraction (f0, superquadric (x, y, 20, 0.5*H0, 0.5*D0));
#else
  fraction (f0, circle(x, y, 0.5 * D0));
#endif

#if PYROLYSIS_ONLY
  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
#else
  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;
#endif

#if MOISTURE
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 0.935; // 93.5% biomass
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")]   = 0.061; // 6.1% moisture
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]     = 0.004; // 0.4% ash
#else
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;
#endif

  foreach()
    porosity[] = eps0*f0[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += f0[]*(1. - eps0)*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  TG[left] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);

  /**
  The mechanism carries O2 in either branch. The `PYROLYSIS_ONLY` build gives
  it a mass fraction of zero at the boundaries, so the particle sees none. */

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
#if PYROLYSIS_ONLY
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[left] = dirichlet (1.);
      YG[top] = dirichlet (1.);
    }
#else
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[left] = dirichlet (0.765);
      YG[top] = dirichlet (0.765);
    } else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      YG[left] = dirichlet (0.235);
      YG[top] = dirichlet (0.235);
    }
#endif
    else {
      YG[left] = dirichlet (0.);
      YG[top] = dirichlet (0.);
    }
  }

  if (restore (file = "last-snapshot", list = all)) {
    fprintf (stderr, "Restart file found!\n");
    restarted = true;
  } else {
    fprintf (stderr, "No restart file found, starting from scratch!\n");

    foreach() {
      f[] = f0[];
      porosity[] = eps0*f[];
    }
  }
}

/**
The H2O-weighted path-averaged temperature of the full case: the mole
fraction over a quarter of the domain.

Caution: the weight follows `MOLAR_DIFFUSION`, because `XGList_G` only
exists when that flag is on. The full case always uses the mole fraction.
So the change of weight is part of the transport rung, and only
`test-transport` and `test-full` weight the probes as the full case does. */

double T_H2O_weigthed_average (double x_interp, int n_samples = 1 << (maxlevel - 1),
                               const double length = L0/4.) {
#ifdef MOLAR_DIFFUSION
  scalar YH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
#else
  scalar YH2O = YGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
#endif

  double numerator = 0., denominator = 0.;
  coord pos, box[2] = {{x_interp, 0.}, {x_interp, length}}, nn = {1, n_samples};
  foreach_region (pos, box, nn, reduction(+:numerator) reduction(+:denominator)) {
    double yH2O_local = interpolate_linear (point, YH2O, pos.x, pos.y, pos.z);
    numerator   += yH2O_local;
    denominator += yH2O_local / interpolate_linear (point, T, pos.x, pos.y, pos.z);
  }

  if (denominator <= 0) // avoid division by 0
    return TG0;

  return numerator/denominator;
}

event output (t += 0.01) {

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, restarted ? "a" : "w");
  if (fp == NULL) {
    fprintf (stderr, "Error opening OutputData\n");
    exit(1);
  }

  /**
  These three points keep the layout of `~/temp/expansion/test.c`, so that
  `OutputData` stays comparable with the runs of the archive. The header of
  that case named them 1, 2 and 4 mm, but it sampled 2, 6 and 11 mm from the
  surface. The names below are the distances the case actually samples.

  The five points of the full case go to `TemperatureProfile.dat` instead.
  Do not add them here: `slow_flicker.py` reads the first six columns of
  this file by position. */

  double Tavg[3], sample_points[3] = {H0/2 + 2e-3, H0/2 + 6e-3, H0/2 + 11e-3};

  for (int ii = 0; ii < 3; ii++)
      Tavg[ii] = T_H2O_weigthed_average (sample_points[ii]);

  if (i == 0)
    fprintf (fp, "#t(1) Ms/Ms0(2) Tmax(3) Tavg_2mm(4) Tavg_6mm(5) Tavg_11mm(6)"
                 " dt(7) mgp_i(8) mgp_resa(9)\n");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[] - porosity[])*rhoS*dv();

  fprintf (fp, "%g %g %g %g %g %g %g %d %g\n", t, solid_mass/solid_mass0, statsf(T).max,
                                      Tavg[0], Tavg[1], Tavg[2], dt, mgp.i, mgp.resa);

  fflush(fp);
}

/**
## The probes of the full case

`run/fatehi-combustion.c` writes `TemperatureProfile.dat` with five
H2O-weighted path averages at 2, 4, 8, 11 and 15 mm from the surface. This
event writes the same five points in the same order, so that the amplitudes
compare directly with the runs under `~/temp/fatehi` and `slow_flicker.py`
reads them
with the branch it already has. The temperature gap between the two cases
is a factor 5, which is larger than the gap of the release rate, so this
file carries the quantity that this campaign must reduce.

`T_H2O_weigthed_average` is collective: it reduces over `foreach_region`.
Call it on every rank, and write on rank 0 only. */

event temperature_profile (t += 0.01) {

  double Tavg[5], sample_points[5] = {H0/2 + 2e-3, H0/2 + 4e-3, H0/2 + 8e-3,
                                      H0/2 + 11e-3, H0/2 + 15e-3};

  for (int ii = 0; ii < 5; ii++)
    Tavg[ii] = T_H2O_weigthed_average (sample_points[ii]);

  if (pid() == 0) {
    static FILE * fpT = NULL;
    if (!fpT) {
      fpT = fopen ("TemperatureProfile.dat", restarted ? "a" : "w");
      if (fpT == NULL) {
        fprintf (stderr, "Error opening TemperatureProfile.dat\n");
        exit (1);
      }
      fprintf (fpT, "#t(1) T2mm(2) T4mm(3) T8mm(4) T11mm(5) T15mm(6)\n");
    }
    fprintf (fpT, "%g %g %g %g %g %g\n",
             t, Tavg[0], Tavg[1], Tavg[2], Tavg[3], Tavg[4]);
    fflush (fpT);
  }
}

/**
Diagnostics for the phase-change expansion field, sampled every timestep.

The projection enforces div(uf) = -gas_source, so we log the two sides of that
identity separately:
  Qsrc = int gas_source dV   the source computed by the chemistry
  Qdiv = int div(uf) dV      the expansion the velocity field actually carries
and 'resmax', the largest pointwise violation. Note gas_source already carries
cm[], as does the discrete div(uf) below, so both integrate with sq(Delta).

'ur' probes the radial velocity on a 45 degree ray just outside the particle:
this is the quantity the velocity vectors show.

Columns 17 to 19 carry the instruments that the amplitude campaign needs.

  mdot    the mass release rate of the particle, integrated:
          `omega*(f - porosity)` is the rate per unit volume, because `omega`
          is the rate per cubic metre of solid material and `(f - porosity)`
          is `f(1 - eps)`, the solid volume fraction. This is the metric of
          the campaign. Log it, do not differentiate `Ms/Ms0`: a numerical
          derivative of a sampled signal changes the amplitude, and the
          campaign compares amplitudes across nine runs.

          Caution: `omega` is positive where the solid decomposes, so `mdot`
          must be positive and must agree with `-d(Ms)/dt`. Check the sign on
          the first run before you trust the column.

  ncells  the number of leaf cells. An adaptation event changes it, so this
          column separates a grid event from a physics event. The pyrolysis
          case loses its spectral peak at level 11, so the mesh is a suspect
          for the excitation of the mode, not only for its resolution.

  nsolid  the number of cells with `f > F_ERR`. It says whether the reacting
          volume switches cells. The gas-source campaign falsified that for
          the fast band. It is free here.
*/

event probe_expansion (t += 0.01) {
  double Qsrc = 0., Qdiv = 0., resmax = 0.;

  foreach (reduction(+:Qsrc) reduction(+:Qdiv) reduction(max:resmax)) {
    double d = 0.;
    foreach_dimension()
      d += uf.x[1] - uf.x[];
    d /= Delta;                        // = cm*div(u), same weighting as gas_source

    Qsrc += gas_source[]*sq(Delta);
    Qdiv += d*sq(Delta);
    resmax = max (resmax, fabs (d + gas_source[]));
  }

  stats so = statsf (omega);

  double ur[3];
  for (int k = 0; k < 3; k++) {
    double r = 0.5*D0 + (k + 1)*0.5e-3, c = cos(pi/4.), s = sin(pi/4.);
    ur[k] = interpolate (u.x, r*c, r*s)*c + interpolate (u.y, r*c, r*s)*s;
  }

  /**
  Particle-side temperature. statsf(T).max is useless here: it sits on TG0
  because the particle starts at TS0 < TG0 and is always the coldest region,
  so a maximum only ever reports the inlet boundary condition.

  We use T[] (assigned once per step as TS[]+TG[], both in tracer form, so it
  is the volume-weighted mixture temperature) together with f[]. Both are
  unambiguous wherever this event runs, unlike TS[]/TG[] whose normalisation
  depends on the position within the timestep.

    Tcore  coldest point, i.e. the centre of the particle
    Tbulk  f-weighted mean over the solid, the bulk particle temperature
    Tsurf  mean over interface cells: where the Arrhenius feedback bites first
  */

  double Tcore = statsf(T).min;
  double Tbulk = 0., fvol = 0., Tsurf = 0., nsurf = 0.;
  double mdot = 0., nsolid = 0.;

  foreach (reduction(+:Tbulk) reduction(+:fvol)
           reduction(+:Tsurf) reduction(+:nsurf)
           reduction(+:mdot) reduction(+:nsolid)) {

    /**
    `omega` is zero outside the solid, so this sum needs no test. */

    mdot += omega[]*(f[] - porosity[])*dv();

    if (f[] > F_ERR) {
      nsolid += 1.;
      Tbulk += T[]*f[]*dv();
      fvol  += f[]*dv();
    }
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      Tsurf += T[];
      nsurf += 1.;
    }
  }
  Tbulk = fvol   > 0. ? Tbulk/fvol  : 0.;
  Tsurf = nsurf  > 0. ? Tsurf/nsurf : 0.;

  /**
  Everything above is collective (reductions / interpolate), so every rank
  holds the same values. Only rank 0 writes: basilisk wraps popen for MPI but
  not fopen, so an unguarded fprintf would emit one copy of each line per rank.
  */

  if (pid() == 0) {
    static FILE * fe = fopen ("expansion.dat", restarted ? "a" : "w");
    if (fe == NULL) {
      fprintf (stderr, "Error opening expansion.dat\n");
      exit(1);
    }
    if (i == 0)
      fprintf (fe, "#t(1) dt(2) Qsrc(3) Qdiv(4) resmax(5) omega_min(6) omega_max(7)"
                   " ur_0.5mm(8) ur_1mm(9) ur_1.5mm(10) mgp_i(11) mgp_resa(12)"
                   " mgpsf_i(13) Tcore(14) Tbulk(15) Tsurf(16)"
                   " mdot(17) ncells(18) nsolid(19)\n");

    /**
    The new columns go on the end. `slow_flicker.py` reads this file with
    `load(path, 16)`, so it truncates to column 16 and keeps working. */

    fprintf (fe, "%g %g %g %g %g %g %g %g %g %g %d %g %d %g %g %g %g %ld %g\n",
             t, dt, Qsrc, Qdiv, resmax, so.min, so.max,
             ur[0], ur[1], ur[2], mgp.i, mgp.resa, mgpsf.i,
             Tcore, Tbulk, Tsurf,
             mdot, grid->tn, nsolid);
    fflush (fe);
  }
}

/**
## Angular profile along the particle surface

For a sphere blowing uniformly the *normal* velocity at the surface is the same
at every angle,

  un_pred = Q/(4 pi R^2),   Q = 2 pi Qdiv the total volumetric production,

which reduces to Qdiv/(2 R^2). |u| is *not* uniform even then, because the free
stream adds a tangential component that varies with angle. So 'un' answers "is
the blowing uniform?" while 'umag' is what the |u| screenshots show; logging
both separates the anomaly from the free-stream confound.

theta runs from the downstream pole (0 deg) through the equator (90) to the
upstream pole (180). Along each ray we locate the interface, sample the gas
just outside it, and take the *peak* of omega inside it -- peak rather than a
fixed depth, so the measurement follows the reaction front instead of sliding
off it as the front recedes. 'r_front' vs theta is then the front shape, and
'T_front' the temperature driving the local Arrhenius rate.

'T_gas' (column 10) is the temperature at the same point as 'un', just
outside the interface. The candidate mechanism of the slow oscillation is a
shield: a higher release rate raises the blowing, the blowing holds the hot
gas away from the surface, the surface cools, and the release falls. That
loop needs the temperature difference which drives the surface, and nothing
else in this case measures it. Read 'T_gas - T_front' against 'un'.

Every interpolate() here is collective and called an identical number of times
on every rank (the branch conditions depend only on reduced values), so the
arrays agree across ranks and only rank 0 writes.
*/

#define NANG 24        // angular samples over [0,pi]
#define NRAD 32        // radial samples along each ray

double profile_offset = 0.15;  // gas-side sampling offset, in units of R

event angular_profile (t += 0.01) {
  const double R = 0.5*D0, dr = profile_offset*R;

  double th[NANG], ri[NANG], un[NANG], um[NANG], om[NANG], rf[NANG], Tf[NANG];
  double Tg_out[NANG];

  for (int k = 0; k < NANG; k++) {
    double theta = (k + 0.5)*pi/NANG, c = cos(theta), s = sin(theta);

    /**
    Locate the interface along the ray. Marching outward to the first f < 0.5
    keeps this correct if the interface ever starts moving again.
    */

    double rint = R;
    for (int j = 0; j < NRAD; j++) {
      double rr = (j + 0.5)*1.5*R/NRAD;
      double ff = interpolate (f, rr*c, rr*s);
      if (ff != nodata && ff < 0.5) { rint = rr; break; }
    }

    /**
    Peak reaction rate along the ray, and the temperature where it peaks.
    */

    double ommax = 0., rfront = 0.;
    for (int j = 0; j < NRAD; j++) {
      double rr = (j + 0.5)*rint/NRAD;
      double o = interpolate (omega, rr*c, rr*s);
      if (o != nodata && fabs(o) > fabs(ommax)) {
        ommax = o; rfront = rr;
      }
    }

    /**
    Gas side, just outside the interface.
    */

    double ux = interpolate (u.x, (rint + dr)*c, (rint + dr)*s);
    double uy = interpolate (u.y, (rint + dr)*c, (rint + dr)*s);
    double Tg = interpolate (T, (rint + dr)*c, (rint + dr)*s);

    th[k] = theta*180./pi;
    ri[k] = rint;
    un[k] = ux*c + uy*s;                 // outward normal component = blowing
    um[k] = sqrt (sq(ux) + sq(uy));
    om[k] = ommax;
    rf[k] = rfront;
    Tf[k] = interpolate (T, rfront*c, rfront*s);
    Tg_out[k] = Tg;
  }

  /**
  Total production, for the uniform-blowing reference. Same weighting as in
  probe_expansion: the discrete divergence carries cm[], so it integrates
  with sq(Delta).
  */

  double Qdiv = 0.;
  foreach (reduction(+:Qdiv)) {
    double d = 0.;
    foreach_dimension()
      d += uf.x[1] - uf.x[];
    Qdiv += d*Delta;
  }
  double un_pred = Qdiv/(2.*sq(R));

  if (pid() == 0) {
    static FILE * fa = fopen ("angular.dat", restarted ? "a" : "w");
    if (fa == NULL) {
      fprintf (stderr, "Error opening angular.dat\n");
      exit(1);
    }
    if (i == 0)
      fprintf (fa, "#t(1) theta_deg(2) r_int(3) un(4) umag(5) omega_max(6)"
                   " r_front(7) T_front(8) un_pred(9) T_gas(10)\n");

    /**
    `T_gas` goes on the end, because `slow_flicker.py` reads this file with
    `load(path, 9)` and truncates to column 9. */

    for (int k = 0; k < NANG; k++)
      fprintf (fa, "%g %g %g %g %g %g %g %g %g %g\n",
               t, th[k], ri[k], un[k], um[k], om[k], rf[k], Tf[k], un_pred,
               Tg_out[k]);
    fflush (fa);
  }
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
  clear();
  view (theta=0, phi=0, psi=-pi/2., width = 1080, height = 1080);
  squares ("T", min = 300, max = 2000, spread = -1, linear = true);
  //isoline ("zmix - zsto", lw = 1.5, lc = {1., 1., 1.});
  draw_vof ("f", lw = 1.5);
  mirror ({0, 1}) {
    squares ("O2_G + O2_S", min = 0., max = 0.235, spread = -1, linear = true);
    //isoline ("zmix - zsto", lw = 1.5, lc = {1., 0., 0.});
    draw_vof ("f", lw = 1.5);
  }
  save ("movie.mp4");
}

event dump (t = 1; t += 1) {
  dump("last-snapshot");
}

event stop (t = tend);

/**
~~~gnuplot
~~~
**/
