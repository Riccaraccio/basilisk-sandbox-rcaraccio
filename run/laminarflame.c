/**
# Axisymmetric H2/N2-air coflow flame: the Toro case

A laminar axisymmetric hydrogen-air diffusion flame burns above a coaxial
burner. The central nozzle carries a H2/N2 mixture, the outer nozzle carries
air, and both streams enter at 0.5 m/s and 300 K. The case compares against
the measurements of [Toro et al. 2005](#toro2005combined) in
`data/toro2004/50cms/`:

| file | quantity |
|---|---|
| `axis-T.exp` | the temperature on the axis, CARS and Raman |
| `radial-<X>mm-T.exp` | the temperature at x = 3, 10, 20 and 30 mm |
| `radial-<X>mm-{H2,O2,H2O}.exp` | the mole fractions at the same stations |

This case is the port of `laminarflame.c` of the sandbox of Edoardo Cipriano.
That sandbox solves the gas with a low-Mach solver and a phase model. This
sandbox has no low-Mach solver, so the port uses the module stack of the
pyrolysis cases and runs it with no condensed phase.

## How the port works

The stack of this sandbox always carries two phases. This flame has one. The
port sets `f[] = 0` in every cell, and the whole solid branch then switches
off:

- `chemistry.h` integrates the solid reactor only where `f > F_ERR`, and the
  gas reactor everywhere `f < 1 - F_ERR`. With `f = 0` each cell is a
  homogeneous gas reactor, which is what a gas flame needs.
- `multicomponent-properties.h` fills the solid properties only where
  `f > F_ERR`, so no solid property function runs.
- `shrinking.h` builds `gas_source` from `omega*(f - porosity)`, which is 0.
  The velocity potential `psi` is then 0 and `ubf` equals `uf`.
- The interface solves of `int-temperature.h` and `int-concentration.h` act
  on interface cells. There are none.

The expansion of the gas does NOT come from `gas_source`. It comes from
`drhodt`, which `multicomponent-properties.h` builds from the heat release and
the change of the mixture molecular weight, and which
`navier-stokes/centered-phasechange.h` adds to the projection. So the flame
still pushes the flow.

## The kinetic scheme

Caution: do not point `kinfolder` at a gas-only scheme such as
`skeletal/hydrogen`. `memoryallocation-varprop.h` calls
`OpenSMOKE_NumberOfSolidSpecies()` in its `defaults` event, and
`int-temperature.h` calls `OpenSMOKE_IndexOfSolidSpeciesWithoutError()`. Both
read the solid map of OpenSMOKE++. A folder with no `kinetics.solid.xml`
leaves that map at NULL, and the run stops with a segmentation fault.

`biomass/Red-gas-2507` is the default here. It has a solid block, which
satisfies the constraint above, and its gas block carries the full H2/O2
chain (H, O, OH, HO2, H2O2), which the flame needs. The solid species stay at
0 because `f = 0`.

To use a different scheme, build with `-DKINFOLDER='"<folder>"'`.
`biomass/Solid-gas-88` also carries the chain, and it is more complete, but it
has 87 gas species against the 57 of the default, so it costs more.
`biomass/Solid-gas-2507` has 36 gas species and no OH and no HO2. Do not use
it for a flame. */

#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1

/**
The knobs below are overridable from the `Makefile`, so that an A/B keeps one
source.

Caution: `CORRECTIVE_CFL` is a floating constant. The preprocessor rejects a
floating constant in an `#if`, so test it at run time, never with `#if`. */

#ifndef KINFOLDER
# define KINFOLDER "biomass/Red-gas-2507"
#endif

#ifndef MAXLEVEL
# define MAXLEVEL 7
#endif

#ifndef TEND
# define TEND 0.3
#endif

#ifndef DT_VALUE
# define DT_VALUE 2e-4
#endif

#ifndef CFL_VALUE
# define CFL_VALUE 0.5
#endif

#ifndef CORRECTIVE_CFL
# define CORRECTIVE_CFL 0.8
#endif

/**
`FROZEN_CELL_GATE` spends one evaluation of the right-hand side to find that a
cell does not react, and then skips the stiff solve. Most of this domain holds
cold air, which is the state that the gate answers cheaply.

Caution: no run of THIS case has measured the gain. `fatehi-combustion.c`
measured 1.28 times more simulated time for the same wall clock, on a
different configuration. Repeat the measurement here before you quote it. */

#ifndef FROZEN_CELL_GATE
# define FROZEN_CELL_GATE 1
#endif

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "gravity.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "spark.h"
#include "view.h"

/**
## Simulation setup

The geometry, the inlet velocity, and the inlet composition follow the 50 cm/s
case of Toro et al. The mass fractions give a fuel stream of 0.0666 H2 in N2
and an oxidizer stream of air. */

#define FUEL_NOZZLE 9e-3
#define AIR_NOZZLE  95e-3
#define LENGTH      150e-3
#define U_AIR       0.5
#define H2_IN       0.066637
#define O2_IN       0.231442
#define T_IN        300.

/**
The inlet velocity of the fuel follows a parabolic profile, which is more
realistic than a flat profile for a laminar flow in a nozzle. The mean value
of the profile equals the velocity of the coflow. The central nozzle (y < R)
carries the fuel, the outer nozzle carries the air.

The domain is 150 mm wide and the outer nozzle is 95 mm wide, so the coflow
covers the whole inlet outside the fuel core. `AIR_NOZZLE` therefore only
records the burner of the experiment. */

double R = 0.5*FUEL_NOZZLE;

#define parabolic (2.*U_AIR*(1. - sq(y)/sq(R)))

u.n[left]  = dirichlet (y <= R ? parabolic : U_AIR);
u.t[left]  = dirichlet (0.);
p[left]    = neumann (0.);
pf[left]   = neumann (0.);
psi[left]  = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right]   = dirichlet (0.);
pf[right]  = dirichlet (0.);
psi[right] = neumann (0.);

u.n[top]   = neumann (0.);
u.t[top]   = neumann (0.);
p[top]     = dirichlet (0.);
pf[top]    = dirichlet (0.);
psi[top]   = dirichlet (0.);

/**
## The spark

The spark ignites the mixture at the lip of the fuel nozzle. `SPARK_HEATFLUX`
adds a heat source to `sGT`, which is the source of the gas energy equation.
Use this policy, not `SPARK_CONSTANT`.

`SPARK_CONSTANT` writes the temperature field directly. The density then
belongs to the old temperature, the divergence does not carry the jump, and
the projection fails within one step. A test at 1800 K broke the run at the
first spark step: 16 species solves and the pressure solve all stopped at 100
iterations, and the timestep became NaN.

The order of the events makes the heat source safe. `event reset_sources` of
`multicomponent-varprop.h` zeroes `sGT`, and it runs BEFORE `event chemistry`,
which is where `spark.h` adds the source. `event tracer_diffusion` reads `sGT`
after that. To confirm the order, read the registration order in the source
that `qcc -source` writes.

`spark.T` must be a scratch field, NOT `TG`. `SPARK_HEATFLUX` writes the
value of the source into that field before it accumulates it into `sGT`. */

scalar qspark[];

/**
`SPARK_Q` is a rate of temperature, in K/s, because `spark.h` multiplies it by
`rhoGv_G*cpGv_G`.

A sweep at level 7 gave this:

| `SPARK_Q` | result |
|---|---|
| 1e5 | no ignition. `Tmax` peaks at 810 K and decays to 553 K by t = 0.06 |
| 5e5 | ignition at t = 0.015. `Tmax` about 2500 K |
| 2e6 | ignition, but `Tmax` reaches 3120 K |

Caution: the three rows above ran at `SPARK_DURATION` of 0.02, which the
block below shows to be too long. Only the value that ignites carries over.
Repeat the sweep at 0.008 before you read the temperatures as settled.

5e5 is the default. Do not lower it below the value that ignites: a spark
that only heats the gas gives a plume which cools and leaves no flame.

## Do not hold the spark after the mixture ignites

`SPARK_DURATION` is 0.008, and the value matters more than `SPARK_Q`. The
kernel ignites 5 ms after the spark starts, so a spark which runs for 20 ms
keeps heating a gas that already burns.

An A/B at `SPARK_Q = 5e5` changed the duration alone. Both runs hold the same
state to t = 0.0196, because they only differ after the shorter spark ends.
The step then separates:

| t | `dt`, duration 0.020 | `dt`, duration 0.008 |
|---|---|---|
| 0.0180 | 6.1e-6 | 6.1e-6 |
| 0.0189 | 3.3e-5 | 3.5e-5 |
| 0.0196 | 3.4e-5 | 6.8e-5 |
| 0.0202 | **6.0e-9** | 2.3e-5 |
| 0.0275 | dead | 3.6e-5 |

The long spark collapses the step by four orders of magnitude and never
recovers. The short spark dips to 6.1e-6 at the instant the spark stops, then
recovers and holds 2.3e-5 to 3.6e-5.

The peak temperature follows the same order. At t = 0.02 the long spark gave
2494 K, above the adiabatic flame temperature of this diluted stream, which is
near 1960 K; the short spark gave 2284 K.

The long spark also gave 29 messages of `convergence ... not reached`, and the
short spark gave none. See the `TOLERANCE` block in `main()` for how to read
those messages.

So: if you raise `SPARK_Q`, do NOT also lengthen `SPARK_DURATION` to match.
Check the step in `tsolve.dat` before you trust a run that carries a new
spark.

## What the default gives

The run at the default spark holds the flame after the spark stops at
t = 0.018. One run reached t = 0.08 and it ended with status 0:

| t | `Tmax` | `dt` | messages |
|---|---|---|---|
| 0.015 | 2524 K | 1.7e-4 | 0 |
| 0.020 | 2284 K | 1.5e-5 | 0 |
| 0.025 | 2162 K | 2.0e-5 | 0 |
| 0.030 | **1988 K** | 5.5e-5 | 0 |
| 0.035 | 2013 K | 7.9e-5 | 0 |
| 0.040 | 2062 K | 9.8e-5 | 0 |
| 0.045 | 2174 K | 1.1e-4 | 0 |
| 0.050 | **2176 K** | 1.3e-4 | 0 |
| 0.060 | 2173 K | 1.4e-4 | 0 |
| 0.070 | 2160 K | 1.6e-4 | 0 |
| 0.080 | 2152 K | 1.7e-4 | 0 |

The step is healthy. It recovers from 1.5e-5 to 1.7e-4, which is the 1.9e-4
that the cold flow carried, to a factor 1.1. No solve of the whole run reports
a failure.

`Tmax` has two turns. It falls while the heat of the spark leaves and it
reaches a minimum of 1988 K at t = 0.030. It then rises to a maximum of
2176 K at t = 0.050. It then falls again, slowly: 2173, 2160, 2152, which is
about 2 K per 5 ms at the end, and the fall decelerates.

Caution: do NOT read the 1988 K of t = 0.030 as agreement with the adiabatic
flame temperature of this diluted stream, which is near 1960 K. The curve
crosses that value on the way up; it does not stop there. At t = 0.08 the peak
is 2152 K, which is 10 per cent above the adiabatic value.

Two readings of the 10 per cent, and this case does not separate them:

- Preferential diffusion. The Lewis number of H2 is near 0.3, and a hydrogen
  flame with differential diffusion carries a temperature above the adiabatic
  value where the flame curves. `MOLAR_DIFFUSION` and `FICK_CORRECTED` are on,
  so the solver carries that effect. Then the excess is physics.
- The three transport flags each carry a correction of the flux, and
  `run/test.c` measured defects in all three. Then part of the excess is the
  scheme.

To separate them, repeat the run with the three transport flags off and read
the peak again. A level ladder answers the refinement at the same time.

Caution: t = 0.08 is 62 ms after the spark, and `Tmax` still falls there. The
flame that Toro et al. measured is a steady flame, and `TEND` is 0.3. No run
of this case has yet reached that. Do not quote a comparison with the
measurement until one does. */

#ifndef SPARK_X
# define SPARK_X 5e-3
#endif

#ifndef SPARK_D
# define SPARK_D 4e-3
#endif

#ifndef SPARK_Q
# define SPARK_Q 5e5
#endif

#ifndef SPARK_START
# define SPARK_START 0.01
#endif

#ifndef SPARK_DURATION
# define SPARK_DURATION 0.008
#endif

const double tend = TEND;
int maxlevel = MAXLEVEL, minlevel = 2;

/**
`D0` is the length that scales the flame. This case does not read it, but
`flame.h` declares it `extern` and the other cases of this sandbox carry the
name, so keep it. Add `#include "flame.h"` to get the mixture fraction and the
flame diameter.

Caution: `flame.h` builds the fuel stream from `YGList_Int`, which is the
composition at the interface. This case has no interface, so `avg_interface`
returns 0 and the mixture fraction is meaningless. Give `flame.h` the
composition of the nozzle before you use it here. */

double D0 = FUEL_NOZZLE;

int main() {

  /**
  Caution: under MPI every rank shares this stderr. Guard the message with
  `pid() == 0`, or the log carries one copy per rank.

  Read this line before you quote a run. It prints what the build enabled, not
  what the directory name promises. */

  if (pid() == 0)
    fprintf (stderr, "# laminarflame: maxlevel=%d DT=%g CFL=%g Uair=%g"
                     " tend=%g corrCFL=%g frozen=%d kin=%s nranks=%d\n",
             MAXLEVEL, (double) DT_VALUE, (double) CFL_VALUE, U_AIR,
             (double) TEND, (double) CORRECTIVE_CFL, FROZEN_CELL_GATE,
             KINFOLDER, npe());

  /**
  The initial thermodynamic state of the gas. `TS0` is the solid temperature.
  No cell holds solid, so the value only keeps the fields finite. */

  Pref = 101325.;
  TG0 = T_IN, TS0 = T_IN;

  /**
  The properties of the solid. This case has no solid, so these values are
  placeholders that keep the header happy. */

  rhoS = 1500., eps0 = 0.;

  /**
  `two-phase.h` declares `rho1`, `rho2`, `mu1` and `mu2`.
  `variable-properties.h` overrides all four with the fields of the mixture,
  so the values below never reach the momentum equation. */

  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  /**
  `ZETA_SWELLING` holds `zeta` at 0. The partitioning has no meaning with no
  solid, and this policy costs the least. */

  zeta_policy = ZETA_SWELLING;

  L0 = LENGTH;
  G.x = -9.81;
  DT = DT_VALUE;

  kinfolder = KINFOLDER;

  /**
  The projection. `project_sf()` passes `TOLERANCE/sq(dt)` to `poisson()`, so
  the default 1e-3 with `dt = 2e-4` gives about 2.5e4 and the solve stops
  after one cycle at every step. Keep both lines together.

  Caution: no `defaults` event resets `TOLERANCE` or `NITERMIN`, so `main()`
  is the right place for them. `CFL` is the opposite case; see `event init`.

  Caution: `TOLERANCE` is ABSOLUTE, and 1e-5 comes from
  `fatehi-combustion.c`, which is a much gentler problem. The run that held
  the spark too long, and that then collapsed, gave 29 messages of
  `convergence ... not reached after 100 iterations` near t = 0.02. Read the
  two numbers of each message together before you call it a failure:

  | field | `res` | `sum` | res/sum |
  |---|---|---|---|
  | a species | 2.79 | -1.3e19 | 2e-19 |
  | a species | 0.041 | -3.5e16 | 1e-18 |
  | `p`, `pf`, `u.x` | 3.9e11 | 2.9e13 | 1e-2 |

  `sum` is the right-hand side of the solve. The species solves are converged
  to 1e-18 of it, and they cannot reach an absolute 1e-5 against a
  right-hand side of 1e19. The rows near 1e11 are the pressure solves, which
  carry the `1/dt^2` of `project_sf()`.

  So each species solve converged to 1e-18 of its own right-hand side. The
  message reports the absolute tolerance, which such a right-hand side puts
  out of reach. That is a scale of the tolerance.

  Caution: do NOT read that as "the messages are harmless". Every one of them
  came from the run which then collapsed, and the run at the default spark
  gave none. The messages did not cause the collapse, and they did not warn
  of it either, because they measure the wrong quantity.

  Read `tsolve.dat` instead. Column 2 is the step. A step which falls and does
  not recover is the failure; the messages are not. */

  TOLERANCE = 1e-5;
  NITERMIN = 2;

  init_grid (1 << min (maxlevel, 7));
  run();
}

event init (i = 0) {

  /**
  Caution: `navier-stokes/centered.h` assigns `CFL = 0.8` in its `defaults`
  event, which runs after `main()`, and `vof.h` then clamps it to 0.5. So the
  value belongs here. To confirm that it took effect, read the live variable,
  not the macro. */

  CFL = CFL_VALUE;

  /**
  No cell holds solid. `f` and `porosity` stay at 0 for the whole run, because
  `shrinking.h` advects them with `ubf`, and no source creates them. */

  foreach() {
    f[] = 0.;
    porosity[] = 0.;
    u.x[] = U_AIR;
  }

  /**
  The domain starts full of the oxidizer stream. `memoryallocation-varprop.h`
  reads `gas_start` in its own `init` event and fills every gas field with it.
  That event runs after this one, because this file includes the header. */

  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = O2_IN;
  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1. - O2_IN;

  /**
  The header refuses a solid composition that does not sum to 1. Give the
  whole mass to the ash, which is inert. No cell holds solid, so the value
  never enters a balance. */

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")] = 1.;

  /**
  The spark. See the block above the declaration of `qspark` for the reason
  behind the policy and behind the scratch field. */

  spark.T = qspark;
  spark.policy = SPARK_HEATFLUX;
  spark.position = (coord){SPARK_X, 0.5*FUEL_NOZZLE};
  spark.diameter = SPARK_D;
  spark.temperature = SPARK_Q;
  spark.time = SPARK_START;
  spark.duration = SPARK_DURATION;

  /**
  The boundary conditions of the temperature and of the species live here.
  Those fields do not exist in the global scope, because the header allocates
  them in its `defaults` event. */

  TG[left]   = dirichlet (T_IN);
  TG[top]    = dirichlet (T_IN);
  TG[right]  = neumann (0.);
  TG[bottom] = neumann (0.);

  for (int jj = 0; jj < NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("H2")) {
      YG[left] = dirichlet (y <= R ? H2_IN : 0.);
      YG[top]  = dirichlet (0.);
    }
    else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      YG[left] = dirichlet (y <= R ? 0. : O2_IN);
      YG[top]  = dirichlet (O2_IN);
    }
    else if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[left] = dirichlet (y <= R ? 1. - H2_IN : 1. - O2_IN);
      YG[top]  = dirichlet (1. - O2_IN);
    }
    else {
      YG[left] = dirichlet (0.);
      YG[top]  = dirichlet (0.);
    }
  }
}

/**
## Mesh adaptation

The grid follows the fuel, the temperature, and the velocity. The last cells
of the domain are coarsened, so that the outflow does not drive the cost. */

#if TREE
event adapt (i++) {
  scalar fuel = YGList_G[OpenSMOKE_IndexOfSpecies ("H2")];
  adapt_wavelet ({fuel, T, u.x, u.y},
      (double[]){1e-2, 1e1, 1e-1, 1e-1}, maxlevel, minlevel);
  unrefine (x >= 0.9*L0);
}
#endif

/**
## Log

One line per output step, so that a run can be followed while it goes.

Caution: `dt` here is the step of ONE instant, not a mean. Do not compare the
cost of two runs with it. A pair at t = 0.02 read 6.7e-5 against 1.5e-5, which
suggests a factor 4, while the step counts of the same pair were 186 against
197, which is 6 per cent. Compare the column `i`, which counts every step. */

event logfile (t += 0.005) {
  if (pid() == 0)
    fprintf (stderr, "%g %g %d %g\n", t, dt, i, statsf(T).max);
}

/**
## The temperature probes along the flame

`run/test.c` reads the temperature of the plume as a path average along a
line, because the measurement of Fatehi reads a beam and not a point. This
case uses the same two averages, at several heights of the flame, so that the
two cases compare with one reader.

The path of each probe is radial. It starts on the axis and it ends at
`PROBE_LENGTH`, which covers the flame and part of the coflow.

The H2O-weighted average is the value that a beam of an absorber sees:

$$
  T_{H_2O} = \dfrac{\sum x_{H_2O}}{\sum x_{H_2O}/T}
$$

Caution: the weight follows `MOLAR_DIFFUSION`, because `XGList_G` only exists
when that flag is on.

Caution: both functions reduce over `foreach_region`, thus both are
collective. Call them on every rank and write on rank 0 only. */

#ifndef PROBE_LENGTH
# define PROBE_LENGTH 15e-3
#endif

double T_H2O_weigthed_average (double x_interp,
                               int n_samples = 1 << (MAXLEVEL - 1),
                               const double length = PROBE_LENGTH) {
#ifdef MOLAR_DIFFUSION
  scalar YH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
#else
  scalar YH2O = YGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
#endif

  double numerator = 0., denominator = 0.;
  coord pos, box[2] = {{x_interp, 0.}, {x_interp, length}}, nn = {1, n_samples};
  foreach_region (pos, box, nn, reduction(+:numerator) reduction(+:denominator)) {
    double yH2O_local = interpolate_linear (point, YH2O, pos.x, pos.y, pos.z);
    double T_local = interpolate_linear (point, T, pos.x, pos.y, pos.z);
    if (T_local > 0. && T_local < nodata && yH2O_local > 0.) {
      numerator   += yH2O_local;
      denominator += yH2O_local/T_local;
    }
  }

  if (denominator <= 0.) // avoid division by 0
    return TG0;

  return numerator/denominator;
}

/**
The same path average with a uniform H2O concentration. The weight cancels,
and the average becomes the harmonic mean of the temperature along the path:

$$
  T_{uni} = \dfrac{N}{\sum 1/T}
$$

This form holds a value before the flame makes any H2O, where the weighted
probe returns `TG0`. Compare the two columns of one station against each
other; do not read one for the other.

Caution: every sample counts the same, so the coflow at 300 K dilutes this
average. The weighted form does not carry that dilution. */

double T_uniform_average (double x_interp,
                          int n_samples = 1 << (MAXLEVEL - 1),
                          const double length = PROBE_LENGTH) {
  double numerator = 0., denominator = 0.;
  coord pos, box[2] = {{x_interp, 0.}, {x_interp, length}}, nn = {1, n_samples};
  foreach_region (pos, box, nn, reduction(+:numerator) reduction(+:denominator)) {
    double T_local = interpolate_linear (point, T, pos.x, pos.y, pos.z);
    if (T_local > 0. && T_local < nodata) {
      numerator   += 1.;
      denominator += 1./T_local;
    }
  }

  if (denominator <= 0.) // avoid division by 0
    return TG0;

  return numerator/denominator;
}

/**
The stations are the four heights of the measurement of Toro et al., plus one
at 50 mm, which sits above the tip of the flame.

Caution: a station can read below `T_IN` while the front of the flame crosses
it. A run at `SPARK_Q = 5e5` gave 251 K at 20 mm at t = 0.015 s. This is an
undershoot of the advection of `TG` ahead of the front, not a fault of the
probe. Do not read the first samples of a station as a measurement. */

#define NPROBES 5
static const double probe_x[NPROBES] = {3e-3, 10e-3, 20e-3, 30e-3, 50e-3};

/**
`TemperatureProfile.dat` carries the weighted average of every station, then
the uniform average of the same stations in the same order. A reader which
stops at column 6 keeps the weighted set alone. */

event temperature_profile (t += 0.005) {

  double Tavg[NPROBES], Tuni[NPROBES];
  for (int k = 0; k < NPROBES; k++) {
    Tavg[k] = T_H2O_weigthed_average (probe_x[k]);
    Tuni[k] = T_uniform_average (probe_x[k]);
  }

  double Tmax = statsf(T).max;

  if (pid() == 0) {
    static FILE * fpT = NULL;
    if (!fpT) {
      fpT = fopen ("TemperatureProfile.dat", "w");
      if (fpT == NULL) {
        fprintf (stderr, "Error opening TemperatureProfile.dat\n");
        exit (1);
      }
      fprintf (fpT, "#t(1) Tavg_3mm(2) Tavg_10mm(3) Tavg_20mm(4) Tavg_30mm(5)"
                    " Tavg_50mm(6) Tuni_3mm(7) Tuni_10mm(8) Tuni_20mm(9)"
                    " Tuni_30mm(10) Tuni_50mm(11) Tmax(12) dt(13)\n");
    }
    fprintf (fpT, "%g", t);
    for (int k = 0; k < NPROBES; k++)
      fprintf (fpT, " %g", Tavg[k]);
    for (int k = 0; k < NPROBES; k++)
      fprintf (fpT, " %g", Tuni[k]);
    fprintf (fpT, " %g %g\n", Tmax, dt);
    fflush (fpT);
  }
}

/**
## Post-processing

### Profiles

`interpolate()` gives `nodata` when no rank holds the point. Guard every
value, or the file carries 1e30 and the plot reads it as a temperature. */

static double probe (scalar s, double x, double y)
{
  double v = interpolate (s, x, y);
  return (v == nodata) ? 0. : v;
}

/**
Write one radial profile of the temperature and of three mole fractions, from
the axis to `ymax`. */

static void radial_profile (const char * name, double xp, double ymax, int n)
{
  /**
  Caution: `XGList_G` only exists when `MOLAR_DIFFUSION` is on. The
  measurement gives mole fractions, so a build without that flag writes mass
  fractions in the same columns and does not compare with the data. */

#ifdef MOLAR_DIFFUSION
  scalar * list = XGList_G;
#else
  scalar * list = YGList_G;
#endif
  scalar T2 = T;
  scalar XH2  = list[OpenSMOKE_IndexOfSpecies ("H2")];
  scalar XO2  = list[OpenSMOKE_IndexOfSpecies ("O2")];
  scalar XH2O = list[OpenSMOKE_IndexOfSpecies ("H2O")];

  double * buf = (double *) malloc (5*n*sizeof (double));
  for (int j = 0; j < n; j++) {
    double yp = ymax*(double)j/(double)(n - 1);
    buf[5*j]     = yp;
    buf[5*j + 1] = probe (T2, xp, yp);
    buf[5*j + 2] = probe (XH2, xp, yp);
    buf[5*j + 3] = probe (XO2, xp, yp);
    buf[5*j + 4] = probe (XH2O, xp, yp);
  }

  if (pid() == 0) {
    FILE * fp = fopen (name, "w");
    fprintf (fp, "#y(1) T(2) xH2(3) xO2(4) xH2O(5)\n");
    for (int j = 0; j < n; j++)
      fprintf (fp, "%g %g %g %g %g\n", buf[5*j], buf[5*j + 1],
               buf[5*j + 2], buf[5*j + 3], buf[5*j + 4]);
    fclose (fp);
  }
  free (buf);
}

event profiles (t = tend) {
  char name[80];

  /**
  The axial profile, on the axis of symmetry. */

  scalar YH2  = YGList_G[OpenSMOKE_IndexOfSpecies ("H2")];
  scalar YO2  = YGList_G[OpenSMOKE_IndexOfSpecies ("O2")];
  scalar YH2O = YGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar YN2  = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];

  int n = 100;
  double * buf = (double *) malloc (6*n*sizeof (double));
  for (int j = 0; j < n; j++) {
    double xp = 140e-3*(double)j/(double)(n - 1);
    buf[6*j]     = xp;
    buf[6*j + 1] = probe (T, xp, 0.);
    buf[6*j + 2] = probe (YH2, xp, 0.);
    buf[6*j + 3] = probe (YH2O, xp, 0.);
    buf[6*j + 4] = probe (YO2, xp, 0.);
    buf[6*j + 5] = probe (YN2, xp, 0.);
  }

  if (pid() == 0) {
    sprintf (name, "AxialProfiles-%d", maxlevel);
    FILE * fp = fopen (name, "w");
    fprintf (fp, "#x(1) T(2) H2(3) H2O(4) O2(5) N2(6)\n");
    for (int j = 0; j < n; j++)
      fprintf (fp, "%g %g %g %g %g %g\n", buf[6*j], buf[6*j + 1],
               buf[6*j + 2], buf[6*j + 3], buf[6*j + 4], buf[6*j + 5]);
    fclose (fp);
  }
  free (buf);

  /**
  The radial profiles, at the four stations of the measurement. */

  double stations[4] = {3e-3, 10e-3, 20e-3, 30e-3};
  for (int k = 0; k < 4; k++) {
    sprintf (name, "RadialProfiles%dmm-%d",
             (int)(stations[k]*1e3 + 0.5), maxlevel);
    radial_profile (name, stations[k], 15e-3, 100);
  }

  /**
  The maps of the fields. */

  sprintf (name, "Maps-%d", maxlevel);
  FILE * fpm = pid() == 0 ? fopen (name, "w") : NULL;
  output_field ({YH2, YO2, YH2O, T, u.x, u.y}, fp = fpm, linear = true,
      box = {{0., 0.}, {0.85*LENGTH, 50e-3}});
  if (fpm)
    fclose (fpm);
}

/**
### Movie

The temperature, and the mass fractions of O2 and of H2O.

The ceiling is 2200 K, not the 1960 K of the adiabatic flame temperature.
`Tmax` reaches 2176 K by t = 0.05, so a ceiling at 2000 K saturates the whole
flame and hides the peak. */

event movie (t += 0.005) {
  clear();
  view (tx = -0.5);
  squares ("T", min = 300., max = 2200., linear = true);
  mirror ({0,-1}) {
    cells();
  }
  save ("temperature.mp4");

  clear();
  view (tx = -0.5);
  squares ("O2_G", min = 0., max = O2_IN, linear = true);
  mirror ({0,-1}) {
    squares ("H2O_G", min = 0., max = 0.173, linear = true);
  }
  save ("species.mp4");
}

event stop (t = tend) {
  return 1;
}

/**
## Results

The comparison with the measurement can show a displacement, which comes from
(i) a low level of refinement, (ii) the Soret effect, which the solver does
not carry, and (iii) a gas mechanism that is not tuned for hydrogen.

~~~gnuplot Temperature map
LEVEL = 7

set size ratio -1
unset key
unset xtics
unset ytics
unset colorbox
set pm3d
set pm3d map interpolate 3,3
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
                          0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, \
                      0.625 1 0.9333 0, 0.75 1 0.4392 0, \
                      0.875 0.9333 0 0, 1 0.498 0 0 )

set xlabel "x"
set ylabel "y"
set title "temperature [K]"
set colorbox

splot "Maps-".LEVEL u 1:2:6
~~~

~~~gnuplot H2 map
set title "mass fraction H2 [-]"
splot "Maps-".LEVEL u 1:2:3
~~~

~~~gnuplot O2 map
set title "mass fraction O2 [-]"
splot "Maps-".LEVEL u 1:2:4
~~~

~~~gnuplot H2O map
set title "mass fraction H2O [-]"
splot "Maps-".LEVEL u 1:2:5
~~~

~~~gnuplot u.x map
set title "u.x [m/s]"
splot "Maps-".LEVEL u 1:2:7
~~~

~~~gnuplot Path-averaged temperature against time
reset
set grid
set xlabel "time [s]"
set ylabel "path-averaged temperature [K]"
set key top left box width 1

plot "TemperatureProfile.dat" u 1:2 w l lw 2 t "H2O-weighted, x = 3 mm", \
     "TemperatureProfile.dat" u 1:3 w l lw 2 t "H2O-weighted, x = 10 mm", \
     "TemperatureProfile.dat" u 1:4 w l lw 2 t "H2O-weighted, x = 20 mm", \
     "TemperatureProfile.dat" u 1:5 w l lw 2 t "H2O-weighted, x = 30 mm", \
     "TemperatureProfile.dat" u 1:6 w l lw 2 t "H2O-weighted, x = 50 mm", \
     "TemperatureProfile.dat" u 1:12 w l lw 2 dt 2 lc -1 t "Tmax"
~~~

~~~gnuplot Uniform path average against time
reset
set grid
set xlabel "time [s]"
set ylabel "path-averaged temperature [K]"
set key top left box width 1

plot "TemperatureProfile.dat" u 1:7 w l lw 2 t "uniform, x = 3 mm", \
     "TemperatureProfile.dat" u 1:8 w l lw 2 t "uniform, x = 10 mm", \
     "TemperatureProfile.dat" u 1:9 w l lw 2 t "uniform, x = 20 mm", \
     "TemperatureProfile.dat" u 1:10 w l lw 2 t "uniform, x = 30 mm", \
     "TemperatureProfile.dat" u 1:11 w l lw 2 t "uniform, x = 50 mm"
~~~

~~~gnuplot Axial temperature profile
reset
LEVEL = 7
set grid
set xlabel "axial distance [mm]"
set ylabel "temperature [K]"

plot "AxialProfiles-".LEVEL u ($1*1e3):2 w l dt 1 lc -1 t "Temperature", \
     "../../data/toro2004/50cms/axis-T.exp" u 1:2 w p lc -1 t "Toro et al., 2005 - CARS", \
     "../../data/toro2004/50cms/axis-T.exp" u 3:4 w p lc -1 t "Toro et al., 2005 - Raman"
~~~

~~~gnuplot Radial profiles at x = 3 mm
reset
LEVEL = 7
set grid
set xlabel "radial distance [mm]"
set ylabel "temperature [K]"
set xr[-15:15]
set yr[200:2200]

set y2tics
set y2r[0:1]
set y2label "mole fractions [-]"

plot "RadialProfiles3mm-".LEVEL u ($1*1e3):2 w l dt 1 lc -1 t "Temperature", \
     "RadialProfiles3mm-".LEVEL u (-$1*1e3):2 w l dt 1 lc -1 notitle, \
     "RadialProfiles3mm-".LEVEL u ($1*1e3):3 w l dt 1 lc 1 t "H2" axis x1y2, \
     "RadialProfiles3mm-".LEVEL u (-$1*1e3):3 w l dt 1 lc 1 notitle axis x1y2, \
     "RadialProfiles3mm-".LEVEL u ($1*1e3):4 w l dt 1 lc 2 t "O2" axis x1y2, \
     "RadialProfiles3mm-".LEVEL u (-$1*1e3):4 w l dt 1 lc 2 notitle axis x1y2, \
     "RadialProfiles3mm-".LEVEL u ($1*1e3):5 w l dt 1 lc 3 t "H2O" axis x1y2, \
     "RadialProfiles3mm-".LEVEL u (-$1*1e3):5 w l dt 1 lc 3 notitle axis x1y2, \
     "../../data/toro2004/50cms/radial-3mm-T.exp" u 1:2 w p lc -1 t "Toro et al., 2005", \
     "../../data/toro2004/50cms/radial-3mm-H2.exp" w p lc 1 axis x1y2 notitle, \
     "../../data/toro2004/50cms/radial-3mm-O2.exp" w p lc 2 axis x1y2 notitle, \
     "../../data/toro2004/50cms/radial-3mm-H2O.exp" w p lc 3 axis x1y2 notitle
~~~

~~~gnuplot Radial profiles at x = 10 mm
reset
LEVEL = 7
set grid
set xlabel "radial distance [mm]"
set ylabel "temperature [K]"
set xr[-15:15]
set yr[200:2200]

set y2tics
set y2r[0:1]
set y2label "mole fractions [-]"

plot "RadialProfiles10mm-".LEVEL u ($1*1e3):2 w l dt 1 lc -1 t "Temperature", \
     "RadialProfiles10mm-".LEVEL u (-$1*1e3):2 w l dt 1 lc -1 notitle, \
     "RadialProfiles10mm-".LEVEL u ($1*1e3):3 w l dt 1 lc 1 t "H2" axis x1y2, \
     "RadialProfiles10mm-".LEVEL u (-$1*1e3):3 w l dt 1 lc 1 notitle axis x1y2, \
     "RadialProfiles10mm-".LEVEL u ($1*1e3):4 w l dt 1 lc 2 t "O2" axis x1y2, \
     "RadialProfiles10mm-".LEVEL u (-$1*1e3):4 w l dt 1 lc 2 notitle axis x1y2, \
     "RadialProfiles10mm-".LEVEL u ($1*1e3):5 w l dt 1 lc 3 t "H2O" axis x1y2, \
     "RadialProfiles10mm-".LEVEL u (-$1*1e3):5 w l dt 1 lc 3 notitle axis x1y2, \
     "../../data/toro2004/50cms/radial-10mm-T.exp" u 1:2 w p lc -1 t "Toro et al., 2005 - CARS", \
     "../../data/toro2004/50cms/radial-10mm-T.exp" u 3:4 w p lc -1 t "Toro et al., 2005 - Raman", \
     "../../data/toro2004/50cms/radial-10mm-H2.exp" w p lc 1 axis x1y2 notitle, \
     "../../data/toro2004/50cms/radial-10mm-O2.exp" w p lc 2 axis x1y2 notitle, \
     "../../data/toro2004/50cms/radial-10mm-H2O.exp" w p lc 3 axis x1y2 notitle
~~~

~~~gnuplot Radial profiles at x = 20 mm
reset
LEVEL = 7
set grid
set xlabel "radial distance [mm]"
set ylabel "temperature [K]"
set xr[-15:15]
set yr[200:2200]

set y2tics
set y2r[0:1]
set y2label "mole fractions [-]"

plot "RadialProfiles20mm-".LEVEL u ($1*1e3):2 w l dt 1 lc -1 t "Temperature", \
     "RadialProfiles20mm-".LEVEL u (-$1*1e3):2 w l dt 1 lc -1 notitle, \
     "RadialProfiles20mm-".LEVEL u ($1*1e3):3 w l dt 1 lc 1 t "H2" axis x1y2, \
     "RadialProfiles20mm-".LEVEL u (-$1*1e3):3 w l dt 1 lc 1 notitle axis x1y2, \
     "RadialProfiles20mm-".LEVEL u ($1*1e3):4 w l dt 1 lc 2 t "O2" axis x1y2, \
     "RadialProfiles20mm-".LEVEL u (-$1*1e3):4 w l dt 1 lc 2 notitle axis x1y2, \
     "RadialProfiles20mm-".LEVEL u ($1*1e3):5 w l dt 1 lc 3 t "H2O" axis x1y2, \
     "RadialProfiles20mm-".LEVEL u (-$1*1e3):5 w l dt 1 lc 3 notitle axis x1y2, \
     "../../data/toro2004/50cms/radial-20mm-T.exp" u 1:2 w p lc -1 t "Toro et al., 2005 - CARS", \
     "../../data/toro2004/50cms/radial-20mm-T.exp" u 3:4 w p lc -1 t "Toro et al., 2005 - Raman", \
     "../../data/toro2004/50cms/radial-20mm-H2.exp" w p lc 1 axis x1y2 notitle, \
     "../../data/toro2004/50cms/radial-20mm-O2.exp" w p lc 2 axis x1y2 notitle, \
     "../../data/toro2004/50cms/radial-20mm-H2O.exp" w p lc 3 axis x1y2 notitle
~~~

~~~gnuplot Radial profiles at x = 30 mm
reset
LEVEL = 7
set grid
set xlabel "radial distance [mm]"
set ylabel "temperature [K]"
set xr[-15:15]
set yr[200:2200]

set y2tics
set y2r[0:1]
set y2label "mole fractions [-]"

plot "RadialProfiles30mm-".LEVEL u ($1*1e3):2 w l dt 1 lc -1 t "Temperature", \
     "RadialProfiles30mm-".LEVEL u (-$1*1e3):2 w l dt 1 lc -1 notitle, \
     "RadialProfiles30mm-".LEVEL u ($1*1e3):3 w l dt 1 lc 1 t "H2" axis x1y2, \
     "RadialProfiles30mm-".LEVEL u (-$1*1e3):3 w l dt 1 lc 1 notitle axis x1y2, \
     "RadialProfiles30mm-".LEVEL u ($1*1e3):4 w l dt 1 lc 2 t "O2" axis x1y2, \
     "RadialProfiles30mm-".LEVEL u (-$1*1e3):4 w l dt 1 lc 2 notitle axis x1y2, \
     "RadialProfiles30mm-".LEVEL u ($1*1e3):5 w l dt 1 lc 3 t "H2O" axis x1y2, \
     "RadialProfiles30mm-".LEVEL u (-$1*1e3):5 w l dt 1 lc 3 notitle axis x1y2, \
     "../../data/toro2004/50cms/radial-30mm-T.exp" u 1:2 w p lc -1 t "Toro et al., 2005 - CARS", \
     "../../data/toro2004/50cms/radial-30mm-T.exp" u 3:4 w p lc -1 t "Toro et al., 2005 - Raman", \
     "../../data/toro2004/50cms/radial-30mm-H2.exp" w p lc 1 axis x1y2 notitle, \
     "../../data/toro2004/50cms/radial-30mm-O2.exp" w p lc 2 axis x1y2 notitle, \
     "../../data/toro2004/50cms/radial-30mm-H2O.exp" w p lc 3 axis x1y2 notitle
~~~

## References

~~~bib
@article{toro2005combined,
  title={Combined experimental and computational study of laminar, axisymmetric hydrogen--air diffusion flames},
  author={Toro, VV and Mokhov, AV and Levinsky, HB and Smooke, MD},
  journal={Proceedings of the Combustion Institute},
  volume={30},
  number={1},
  pages={485--492},
  year={2005},
  publisher={Elsevier}
}
~~~
*/
