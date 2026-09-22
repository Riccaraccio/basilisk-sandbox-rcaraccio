/**
# Chemistry solver
This header contains the implementation of the chemistry solver.
The general form of the equation is we solve an equation of the form:
$$\frac{dMs_{i}}{dt} = \sum^{NR}_j R_{ji} \nu_{ji}(1-\epsilon)$$
where y is the mass (or mass fractions) of species i, R is the reaction rate and nu is the 
stoichiometric coefficient of species i in reaction j.

We use the OpenSMOKE++ library to solve the ODE system as this offers
a wide variety of stiff solvers optimized for chemical kinetics.

We also have the option to use an explicit solver for non-stiff problems.
*/

#include "reactors.h"

extern scalar zeta;
extern scalar T;
extern scalar porosity;

#ifdef STORE_SOURCES
/**
Optional diagnostic: when `STORE_SOURCES` is defined, the full per-cell
`data.sources` vector (length NEQ) computed by the reactor is copied into the
user-provided field list `sourcesList`. */
extern scalar * sourcesList;
#endif

/**
`GAS_CHEMISTRY_STRANG` selects the Strang split of the gas chemistry. The
default is 0, the Lie split. See "The Strang split of the gas chemistry"
below. The default stays outside `TURN_OFF_REACTIONS`, so a case can print
the flag in every build. */

#ifndef GAS_CHEMISTRY_STRANG
# define GAS_CHEMISTRY_STRANG 0
#endif

#if GAS_CHEMISTRY_STRANG
/**
The link to `chem-split-probe.h`. The probe sets `strang_probe_armed` on the
step that it measures. The second half then records the largest intrinsic
gas temperature before it runs (`strang_Tmax_tr`) and its heat release
(`strang_Q2`), in the same form as the column `Qchem` of the probe. When no
second half runs, the two values stay 0. */

bool strang_probe_armed = false;
double strang_Tmax_tr = 0., strang_Q2 = 0.;
#endif

#ifndef TURN_OFF_REACTIONS

/**
## Explicit ODE solvers
These are simple implementations of explicit ODE solvers: Euler and
Runge-Kutta 4(5). They can be used for non-stiff problems.
*/
void ODESolverEXP (odefunction ode, unsigned int neq, double dt, double* y, void* args) {

  double dy[neq];
  ode(y, dt, dy, args);

  for (int jj=0; jj<neq; jj++)
    y[jj] += dt*dy[jj];
}

void RungeKutta45EXP(odefunction ode, unsigned int neq, double dt, double *y, void *args) {

  // Allocate arrays for the k values and temporary y values
  double k1[neq], k2[neq], k3[neq], k4[neq], k5[neq], k6[neq];
  double ytmp[neq];

  // Coefficients for the RK45 method
  const double a2 = 1.0 / 4.0;
  const double a3 = 3.0 / 8.0;
  const double a4 = 12.0 / 13.0;
  const double a6 = 1.0 / 2.0;

  const double b31 = 3.0 / 32.0;
  const double b32 = 9.0 / 32.0;

  const double b41 = 1932.0 / 2197.0;
  const double b42 = -7200.0 / 2197.0;
  const double b43 = 7296.0 / 2197.0;

  const double b51 = 439.0 / 216.0;
  const double b52 = -8.0;
  const double b53 = 3680.0 / 513.0;
  const double b54 = -845.0 / 4104.0;

  const double b61 = -8.0 / 27.0;
  const double b62 = 2.0;
  const double b63 = -3544.0 / 2565.0;
  const double b64 = 1859.0 / 4104.0;
  const double b65 = -11.0 / 40.0;

  // Coefficients for the 5th order solution
  const double c1 = 16.0 / 135.0;
  const double c3 = 6656.0 / 12825.0;
  const double c4 = 28561.0 / 56430.0;
  const double c5 = -9.0 / 50.0;
  const double c6 = 2.0 / 55.0;

  // Step 1: Calculate k1 = f(t, y)
  ode(y, 0, k1, args);

  // Step 2: Calculate k2 = f(t + a2*dt, y + a2*k1*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * a2 * k1[j];
  ode(ytmp, a2 * dt, k2, args);

  // Step 3: Calculate k3 = f(t + a3*dt, y + b31*k1*dt + b32*k2*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * (b31 * k1[j] + b32 * k2[j]);
  ode(ytmp, a3 * dt, k3, args);

  // Step 4: Calculate k4 = f(t + a4*dt, y + b41*k1*dt + b42*k2*dt + b43*k3*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * (b41 * k1[j] + b42 * k2[j] + b43 * k3[j]);
  ode(ytmp, a4 * dt, k4, args);

  // Step 5: Calculate k5 = f(t + a5*dt, y + b51*k1*dt + b52*k2*dt + b53*k3*dt + b54*k4*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * (b51 * k1[j] + b52 * k2[j] + b53 * k3[j] + b54 * k4[j]);
  ode(ytmp, dt, k5, args);

  // Step 6: Calculate k6 = f(t + a6*dt, y + b61*k1*dt + b62*k2*dt + b63*k3*dt + b64*k4*dt + b65*k5*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * (b61 * k1[j] + b62 * k2[j] + b63 * k3[j] + b64 * k4[j] + b65 * k5[j]);
  ode(ytmp, a6 * dt, k6, args);

  // Update y using the 5th order solution
  for (int j = 0; j < neq; j++) {
    y[j] += dt * (c1 * k1[j] + c3 * k3[j] + c4 * k4[j] + c5 * k5[j] + c6 * k6[j]);
    y[j] = y[j] < 0 ? 0 : y[j]; // Ensure non-negativity
  }

  y[NGS+NSS] = clamp (y[NGS+NSS], 0., 1.); // Ensure boundness for porosity
}

event init (i = 0) {
  OpenSMOKE_InitODESolver ();
}

event cleanup (t = end) {
  OpenSMOKE_CleanODESolver ();
}

/**
## Diagnostics of the solid source

`SOLID_SOURCE_DIAG` adds two fields that explain a flicker of `gas_source`.

`solid_diag[]` records what the solid branch did with the cell: 0 no solid,
1 the reactor ran, 2 a guard skipped the cell, 3 the solve returned a
non-finite state. `omega[]` starts every step at zero, and only the case 1
writes it. The cases 2 and 3 therefore remove the source of that cell for one
step, and the cell returns the next step. That is a switch, not a rate.

`dTS_step[]` records `|TS_end - TS_start|` over the step, per unit of `f`. It
decides whether the end-state sampling of `omega` matters. `omega` is read at
the converged end state, and the Arrhenius factor is exponential in `TS`. With
`Ea/R = 15000` K at `TS = 800` K, a rise of 10 K changes that factor by about
26 percent, and a rise of 1 K by about 2 percent. Below 1 K the sampling of
`omega` cannot explain a visible flicker, because the solid conversion is much
slower than the step that the gas phase imposes. */

#ifndef SOLID_SOURCE_DIAG
# define SOLID_SOURCE_DIAG 0
#endif

#if SOLID_SOURCE_DIAG
scalar solid_diag[], dTS_step[];
#endif

event reset_sources (i++) {
  foreach() {
    omega[] = 0.;
#if SOLID_SOURCE_DIAG
    solid_diag[] = 0.;
    dTS_step[] = 0.;
#endif
  }
}

#ifdef CHEMISTRY_LOG
scalar t_solid[], t_gas[];
#endif

#ifdef BINNING
/**
Scale the gas tracers (species + temperature) of the current cell by a common
factor. `1/(1-f)` maps the VOF-tracer form (`Y*(1-f)`) to the actual mass
fractions the reactor expects; `(1-f)` is the inverse. */

static void scale_gas_tracers (Point point, double factor) {
  for (int jj = 0; jj < NGS; jj++) {
    scalar YG = YGList_G[jj];
    YG[] *= factor;
  }
  TG[] *= factor;
}
#endif

#ifdef VARPROP
/**
## Gas-phase reaction source for the low-Mach divergence

`DYDtG_G` [kg/m3/s] and `DTDtG` [W/m3] carry the gas-phase reaction
contribution into `drhodt`, and therefore into the right-hand side of the
pressure Poisson equation. Two forms are available.

`gas_source_averaged = false` re-evaluates the reactor right-hand side at the
converged end-of-step state. That is the rate at one single state. Near a stiff
flame a cell alternates between "reacting" and "burnt out" from one step to the
next, so this rate flickers, and the flicker goes directly into the velocity
field. The velocity field then moves the flame, which changes the rate again.

`gas_source_averaged = true` (the default) uses the step-averaged rate,
`(state_end - state_start)/dt`. The integrator already produced both states, so
this costs one subtraction and no extra call to the reactor. It is the exact
mean of the same quantity over the step, it is conservative, and it removes the
end-state sensitivity.

`rhoGv_G` and `cpGv_G` do not change during the chemistry event. The same
values therefore weight the start state and the end state, and the source stays
an exact `rhoGv_G*dY/dt` and `rhoGv_G*cpGv_G*dT/dt`. These weights are the
start values: the last `update_properties()` before this event is the one of
the `adapt` event of the previous step.

Caution: `divu2` in `multicomponent-properties.h` does not divide by the same
values. `update_divergence()` runs after the second `update_properties()` of
the step (the `tracer_diffusion` event of `multicomponent-varprop.h`), which
reads the state after the chemistry. So the numerator has the start weights
and the denominators have the end values. See "The exact expansion" below.

## The weight of the increment

The exact quantity is not `rho*(Y_end - Y_start)/dt`. It is

    (1/dt) * integral of rho(tau)*dY/dtau dtau

so the weight must represent `rho` over the whole step, not at one end of it.
`gas_source_rho_mean = false` (the default) uses `rhoGv_G`, the value at the
step start. In a burning cell the gas expands and `rho` falls by about 30% in
one step, so the start value weights the increment too much.

`gas_source_rho_mean = true` uses the mean of the start and the end values,
`0.5*(rho_start + rho_end)`. `test/gas-source-cell.c` measures both against a
sub-stepped reference over 13 states, with `T` from 1200 to 2100 K, the fuel
mass fraction from 0.02 to 0.20, and `dt` from 2e-6 to 2e-4 s:

    weight       mean ratio to the reference    worst
    rho_start              1.174                1.234
    mean                   1.012                1.025

`cp` stays at the start value. The mean of `cp` changes the result by 0.2%,
which does not pay for the extra call to the property library. The mean of
`rho` needs no call at all: `1/MW = sum_j Y_j/MW_j` gives `rho` from the ideal
gas law with pure arithmetic.

Caution: do not read the end-state density from `data.rhog` after the solve.
`reactors.h` starts the reactor with `UserDataODE data = *(UserDataODE *)args`,
so the reactor writes `rhog` and `cpg` in a local copy and the caller keeps the
old values. `data.sources` behaves differently because it is a pointer. Even
with a pointer, the last evaluation of the right-hand side is a trial point of
the stiff solver, not the converged end state.

Caution: with the averaged form, `TURN_OFF_HEAT_OF_REACTION` also removes the
heat release from the expansion source. The instantaneous form keeps it, because
it fills `sources[NGS]` before it zeroes `dy[NGS]`.

Compile with `-DGAS_SOURCE_AVERAGED=0` to select the instantaneous form and
with `-DGAS_SOURCE_RHO_MEAN=1` to select the mean weight, or assign
`gas_source_averaged` and `gas_source_rho_mean` in `main()` to override the
compiled defaults. The mean weight applies to the averaged form only. The
instantaneous form ignores it.

## The exact expansion

Both forms above feed `divu2` in `update_divergence()`, which divides them by
`TG`, `rhoGv_G` and `cpGv_G`, and multiplies the species part by `MWmixG_G`.
All four hold values after the chemistry. `rhoGv_G`, `cpGv_G` and `MWmixG_G`
come from the `update_properties()` call in the `tracer_diffusion` event of
`multicomponent-varprop.h`, which reads the state after the chemistry. `TG`
also holds the advection of the step. The numerator has the start weights
`rhoGv_G` and `cpGv_G` of this event. So the numerator and the denominators
come from two time levels, and the default path gives too much expansion:

    code  = rho_0*cp_0*(T_end - T_0)/(dt*T_end*rho_end*cp_end)
          + (rho_0/rho_end)*MW_end*sum_j (Y_end,j - Y_0,j)/(MW_j*dt)
    exact = ln(rho_0/rho_end)/dt

`test/gas-source-cell.c` measures `code/exact` in one cell with no flow. With
the dummy kinetics, 13 burning states and `dt` from 2e-6 to 2e-4 s, the
default path gives 1.16 to 1.39. `gas_source_rho_mean` does not repair that,
because it changes the weight and not the time level of the denominators.
An earlier version of this comment and of the test used the start values as
denominators. That gave an expansion 8 to 27 percent under the exact value,
which was not correct.

The exact step mean of the expansion at constant pressure needs no
denominators. The expansion rate is `-d(ln rho)/dt`, so its mean over the step
is the closed form

    ln(rho_start/rho_end)/dt

and both densities follow from the ideal gas law with `1/MW = sum_j Y_j/MW_j`.
The pressure cancels in the ratio. `GAS_SOURCE_EXACT` selects this form, and
it has no denominators, so it removes the error of the two time levels above.
The default of `GAS_SOURCE_EXACT` stays 0. The
chemistry event then writes `cm[]*ln(rho_start/rho_end)/dt` to `drhodt_chem`,
per unit volume of gas, and `update_divergence()` adds it to `divu2` with the
same `(1-f)` weight as the other terms. The reaction part no longer passes
through `DYDtG_G` and `DTDtG`. With this flag `gas_source_averaged` and
`gas_source_rho_mean` have no effect.

The same flag switches on the filter of the divergence source in
`navier-stokes/centered-phasechange.h`. Set `gas_source_filter_passes = 0`
there to keep the exact source and remove the filter.

The exact form covers the gas-phase reactions of the external gas only. The
pore gas inside the solid keeps the source vector of the reactor, because its
temperature equation carries the heat capacity of the solid and of the gas
together, and the closed form does not apply there.

Caution: `TURN_OFF_HEAT_OF_REACTION` zeroes the temperature increment of the
reactor, so with the exact form the heat release also leaves the expansion.
This matches the averaged form.
*/

#ifndef GAS_SOURCE_EXACT
# define GAS_SOURCE_EXACT 0
#endif

#if defined(BINNING) && GAS_SOURCE_EXACT
# error "GAS_SOURCE_EXACT is not available with BINNING. The binning path\
 does not keep the start state of each cell."
#endif

#ifndef GAS_SOURCE_AVERAGED
# define GAS_SOURCE_AVERAGED 1
#endif

#ifndef GAS_SOURCE_RHO_MEAN
# define GAS_SOURCE_RHO_MEAN 0
#endif

#if defined(BINNING) && GAS_SOURCE_RHO_MEAN
# error "GAS_SOURCE_RHO_MEAN needs the end state at the time of the start-state\
 subtraction. The binning path subtracts the start state before the solve and\
 does not keep it, so the mean weight is not available there yet."
#endif

bool gas_source_averaged = GAS_SOURCE_AVERAGED;
bool gas_source_rho_mean = GAS_SOURCE_RHO_MEAN;

/**
The gas density at the end state, from the ideal gas law. `1/MW` is the sum of
`Y_j/MW_j`, so this needs no call to the property library. Returns 0 if the
state is not usable, and the caller then keeps the start value. */

static double gas_end_state_density (Point point, const double * yend) {
  double invMW = 0.;
  for (int jj = 0; jj < NGS; jj++)
    invMW += (yend[jj] > 0. ? yend[jj] : 0.)/gas_MWs[jj];
  double T = yend[NGS];
  if (!(invMW > 0.) || !(T > 0.))
    return 0.;
  return (Pref + p[])/(R_GAS*1000.*T*invMW);
}

#if GAS_SOURCE_EXACT
/**
The exact step mean of the expansion, `ln(rho_start/rho_end)`, from the two
states of the reactor. Returns 0 if one of the states is not usable. */

static double gas_log_expansion (Point point, const double * ystart,
                                 const double * yend) {
  double rho_start = gas_end_state_density (point, ystart);
  double rho_end = gas_end_state_density (point, yend);
  if (!(rho_start > 0.) || !(rho_end > 0.))
    return 0.;
  return log (rho_start/rho_end);
}
#else // !GAS_SOURCE_EXACT
/**
Instantaneous form: one extra evaluation of the reactor at the state `ys`. */

static void gas_sources_instantaneous (Point point, const double * ys) {
  UserDataODE data;
  data.P = Pref + p[];
  data.T = ys[NGS];
  double sources[NGS + 1];
  data.sources = sources;

  double dy_tmp[NGS + 1];
  gas_batch_nonisothermal_constantpressure (ys, dt, dy_tmp, &data);

  for (int jj = 0; jj < NGS; jj++) {
    scalar DYDtGjj = DYDtG_G[jj];
    DYDtGjj[] += sources[jj]*cm[];
  }
  DTDtG[] += sources[NGS]*cm[];
}

/**
One half of the averaged form. Call it with `sgn = -1` on the pre-reaction
state and with `sgn = +1` on the post-reaction state. The split lets the
binning path avoid a per-cell copy of the whole composition vector.

Caution: `rho` and `cp` are arguments, not reads of `rhoGv_G`/`cpGv_G`. The
two calls must use the same values, otherwise the result is
`rho_end*Y_end - rho_start*Y_start`, which folds the density change into the
species source. In the binning path `binning_remap()` updates `rhoGv_G` and
`cpGv_G` between the two calls, so the caller keeps the pre-reaction values. */

static void gas_sources_accumulate_state (Point point, const double * ys,
                                          double rho, double cp, double sgn) {
  if (!(dt > 0.) || !(rho > 0.) || !(cp > 0.))
    return;

  double w = sgn*cm[]/dt;
  for (int jj = 0; jj < NGS; jj++) {
    scalar DYDtGjj = DYDtG_G[jj];
    DYDtGjj[] += rho*ys[jj]*w;
  }
  DTDtG[] += rho*cp*ys[NGS]*w;
}
#endif // GAS_SOURCE_EXACT

/**
Convenience wrapper for the per-cell path, which holds both states and where
`rhoGv_G`/`cpGv_G` do not change over the chemistry event. */

static void accumulate_gas_sources (Point point, const double * ystart,
                                    const double * yend) {
#if GAS_SOURCE_EXACT
  if (dt > 0.)
    drhodt_chem[] += cm[]*gas_log_expansion (point, ystart, yend)/dt;
#else
  if (gas_source_averaged) {
    double rho = rhoGv_G[], cp = cpGv_G[];

    /**
    The mean weight. Both calls below still use one common `rho`, which is
    what keeps the result an increment of `Y` and not an increment of
    `rho*Y`. */

    if (gas_source_rho_mean) {
      double rho_end = gas_end_state_density (point, yend);
      if (rho_end > 0.)
        rho = 0.5*(rho + rho_end);
    }

    gas_sources_accumulate_state (point, ystart, rho, cp, -1.);
    gas_sources_accumulate_state (point, yend,   rho, cp, +1.);
  }
  else
    gas_sources_instantaneous (point, yend);
#endif
}
#endif

/**
## The gate that skips the cells the step cannot change

`FROZEN_CELL_GATE` spends **one** evaluation of the reactor right-hand side to
decide whether the stiff solve of a gas cell can change the state over the
step. The measurements below come from `test/bench88.c` and
`test/gate-ignition.c` with `biomass/Solid-gas-88` (87 gas species, 33 solid
species) at `dt = 2.4e-4` s.

The gas branch solves 88 equations in every cell that holds gas. In the free
stream the mixture is air at the inlet temperature and it does not react, but
the stiff solve still costs 1.8 to 3.8 ms, which is about 100 evaluations. The
gate applies the explicit update instead, which is exact at that size, and
skips the solve. It replaces those 100 evaluations with one.

The margin is large. Air at 1123 K gives `max|dY|` of 2e-17 and `|dT|` of
3e-13 K over the step. Hot air with 0.5 per cent of CO and 0.5 per cent of tar
gives 2e-4 and 3e-2 K.

## An ignition never closes the gate

`test/gate-ignition.c` marches a batch reactor at the step of the production
case and compares the gated trajectory with the ungated one. An ignition is
the case that a rate test can miss: the species move slowly through the
induction period, and the mixture then runs away inside one step.

With the pyrolysis gas of the biomass in air, from 800 K to 1300 K, the gate
stays open on **every one of the 3000 steps**, at every temperature. The
ignition delay and the whole temperature trajectory are equal to the last bit,
even at 800 K where the induction lasts 0.48 s, which is 2000 steps.

At the tolerances below the gate closes on hot air only under about 1e-9 of
fuel by mass. A mixture with 1e-8 of fuel still gets the solve. At the level
where the gate does close, the mixture moves by 2e-9 K over 0.72 s, and the
gated trajectory follows the solved one to 1.4e-9 K.

## Why the tolerances are as tight as they are

The tolerance is the size of the state that the gate throws away on one step,
so it is also the size of the perturbation that the gate feeds to the rest of
the solver. A 2-D A/B on `run/test.c` at `maxlevel` 8 measured it. With
`FROZEN_CELL_YTOL` at 1e-10 the mass, `Tmax`, `dt` and the count of the
pressure iterations stay equal to the printed precision over 0.6 s, but the
projection residual moves by 5e-5 in relative terms from t = 0.33, and the two
`Tavg` probes nearest the surface then differ, because they are a ratio of two
integrals that are both near zero when the plume arrives.

At 1e-15 the same run is equal to the ungated run in every column of every
row. The tolerances are therefore the tight values: they still sit 40 times
above the rate of the free stream, so the gate keeps its work, and they leave
the trajectory reproducible. Raise them only if a run needs the speed more
than it needs a run-to-run comparison.

## What it is worth

A 2-D `fatehi-combustion` at `maxlevel` 7 with the 88-species scheme, over a
fixed wall budget, reaches **1.28 times** the simulated time of the same case
with the gate off. The gate takes the whole domain on the first step, about
42 per cent of the gas cells by step 10, and 9 per cent once the plume is
established. At `FROZEN_CELL_YTOL` of 1e-10 the same case reaches 1.43 times,
with the drift above.

Caution: compare two builds only inside one batch of runs. The absolute time
of this case swings 15 per cent from one run to the next on the same binary,
while the ratio inside a batch repeats to 1 per cent.

`FROZEN_CELL_GATE` is **off by default**. Set it to 1 to switch the gate on.
*/

#ifndef FROZEN_CELL_GATE
# define FROZEN_CELL_GATE 0
#endif

#ifndef FROZEN_CELL_YTOL
# define FROZEN_CELL_YTOL 1e-15
#endif

#ifndef FROZEN_CELL_TTOL
# define FROZEN_CELL_TTOL 1e-11
#endif

#if FROZEN_CELL_GATE
/**
The number of cells that the gate skipped over the last step. It makes the
gain of the gate visible in a production log. The `foreach` loop below carries
a `reduction` clause for it, so the value is the total over every rank and
every thread. */

int frozen_cell_gate_n = 0;
#endif

/**
## The Strang split of the gas chemistry

Item TL-1 of `~/discretization-report/time-level-review.md`. The default
step is a Lie split. The `chemistry` event integrates the gas reactor over
the full `dt` at the start of the step, and the advection and the implicit
diffusion follow. A constant-dt ladder from the plateau at level 10 showed
that this split carries the whole `Tmax(dt)` law: the jump of `Tmax` over
the chemistry is 50.9, 28.8, 15.7 and 8.4 K at `dt` 4e-4, 2e-4, 1e-4 and
5e-5 s, thus first order in `dt`.

`GAS_CHEMISTRY_STRANG` 1 makes the split symmetric:

    R(dt/2)   the `chemistry` event, first in the step
    T(dt)     VOF, advection, interface, species and temperature solves
    R(dt/2)   the `tracer_diffusion` event below, after the solves

The second half runs before `shrinking.h` puts `uf` back, before the
momentum and the projection, and before `end_timestep` and `adapt`. So the
output, the probes of `end_timestep` and the refinement criterion all read
the state after the second half.

What the split covers:

* The gas-phase reactions of the external gas (`YGList_G`, `TG`). This is
  where the flame is.
* NOT the solid reactor. The solid, the pore gas (`YGList_S`, `TS`) and
  `porosity` stay on the full `dt` in the first call, as before. The pore gas
  reacts inside the same stiff system as the solid, with the heat capacity of
  the solid in the temperature equation, so a split of the pore gas needs a
  split of the whole solid reactor. The solid changes slowly: `dTS_step` is
  0.05 to 0.11 K per step, and the split error of the solid is 0.1 to 0.3 %
  of `omega` (review, section 3.1). So `omega`, `zeta`, `prod`, `ubf` and
  `gas_source` come from the same full-dt solve as at 0, at the same point of
  the step. The projection of the step reads the same `gas_source`.

The expansion of the gas reactions (the chemistry part of `drhodt`):

* The first half adds its increment to `DTDtG` and `DYDtG_G` (or to
  `drhodt_chem` under `GAS_SOURCE_EXACT`) exactly as the full step does. The
  weights divide by the step `dt`, not by `dt/2`. So the increment of the
  first half becomes its share of the mean rate of the step.
  `update_divergence()` then puts it into `drhodt`, with no change.
* The second half runs after `update_divergence()`. It therefore adds its
  share directly to `drhodt`, with the weight `(1-f)` and the factor `cm`
  that `update_divergence()` gives the chemistry part. The projection of the
  same step reads `drhodt` in `advection_term` and in `projection`, and both
  run after this event. So the projection of step n receives the sum of the
  two half increments divided by `dt`, which is what the review asks for.
  No increment enters two projections, and no increment is lost.
* Under `GAS_SOURCE_EXACT` both halves use `ln(rho_start/rho_end)/dt`.
* In the default path the second half uses the linear form
  `[(T_1 - T_0)/T_0 + MW_0*sum_j (Y_1,j - Y_0,j)/MW_j]/dt`, with the start
  state of the half as divisor. This is the form of `update_divergence()`
  when the weight and the divisors come from the same level: `rho*cp` of the
  weight cancels. It does not carry the two-level error of TL-3, which the
  first half keeps. The difference to the log form is second order in
  `(T_1 - T_0)/T_0`, which is below 0.03 in one half step on the plateau.
* `gas_source_rho_mean` acts on the first half only.
* The instantaneous form (`gas_source_averaged = false`) evaluates a rate,
  not an increment. The split does not support it, and the run stops at
  `init`.

The heat of reaction is not counted two times. The gas reactor puts its heat
into `TG` only. `data.sources` of the gas branch stays `NULL`, and the
expansion source of each half comes from the state change of that half only.

The state that the second half reads:

* `TG` and `YGList_G` in tracer form, `* (1 - f)`, with `f` after the VOF
  sweep. The `tracer_diffusion` event of `multicomponent-varprop.h` has put
  them back into tracer form and has set `T = TS + TG`. The reactor divides
  by the same `1 - f`, and this event sets `T = TS + TG` again at its end.
* `p` of the start of the step, because the projection has not run.
* `rhoGv_G` and `cpGv_G` of the `update_properties()` call before the solves.
  The reactor recomputes `rho` and `cp` from the state at each evaluation
  under `VARPROP`, so these are start values only. The expansion of the
  second half does not use them.
* The properties are not recomputed after the second half. The momentum of
  this step reads the properties of the start of the step, as at 0. `adapt`
  then computes the properties from the state after the second half, so the
  next step starts with consistent properties.

Caution: with `PROPS_AFTER_SOLVES` the properties of the momentum come from
the state before the second half.

Caution: `drhodt-budget.h` reads `drhodt` before the second half, so its
reference `|drhodt|` does not hold the expansion of the second half.

What the split can change, and what it cannot:

* At a constant `dt` the Strang sequence is the Lie sequence with one more
  half step at the end: `R_h T R_h R_h T R_h = R_h T R T R_h`. So the state
  at the end of a Strang step is `R(dt/2)` applied to the state after the
  transport of a Lie run. The split changes the time level at which the
  output, `adapt` and the probes read the state. It also moves the
  expansion of each half into the projection of its own step. It does not
  change the sequence of the reactor and the transport. So expect a
  smaller `Tmax(dt)` law from a ladder with the split, not zero.
* The transport is first order in time: the diffusion solves are backward
  Euler. The whole step therefore stays first order. The split removes only
  the first-order error of the splitting.
* `test/strang-cell.c` (one stirred cell, exact transport) gives, against a
  fine reference: Lie -42.8, -21.2, -10.5, -5.2 K at `dt` 8e-4, 4e-4,
  2e-4, 1e-4 (first order), and Strang +1.76, +0.50, +0.19, +0.11 K. Below
  about 0.1 K the tolerance of the stiff solver sets the error.

Cost: the gas reactor runs two times in each step, each time over `dt/2`.
Each call of the Gear solver has a fixed start cost, so a cell that does
not react costs about 2 times. A burning cell costs about 1.1 times (0.4 to
1.8, `test/strang-cell.c`), because the solver takes fewer internal steps
over a shorter interval. The smoke run of `run/test.c` at level 8 from
t = 0 to 0.3 s (no flame yet) gave 42.4 s of chemistry at 0 and 69.6 s at
1 (the second half 33.7 s), thus 1.64 times the chemistry and about 1.7
times the gas part. The solid reactor does not change. With
`CHEMISTRY_LOG` the second half prints its time on a line that starts with
`S2`.

Restart: the split adds no field. A snapshot of a run at 0 restarts with the
split and the reverse.

`FROZEN_CELL_GATE` tests each half with its own `dt/2`. `BINNING` is not
available with the split. */

#if GAS_CHEMISTRY_STRANG && defined(BINNING)
# error "GAS_CHEMISTRY_STRANG is not available with BINNING."
#endif

#if GAS_CHEMISTRY_STRANG && TURN_OFF_GAS_REACTIONS
# warning "GAS_CHEMISTRY_STRANG does nothing with TURN_OFF_GAS_REACTIONS."
#endif

#if GAS_CHEMISTRY_STRANG && !TURN_OFF_GAS_REACTIONS
# ifdef VARPROP
/**
The expansion of the second half, per unit volume of the cell. See the list
above for the form. */

static void strang_second_half_expansion (Point point, const double * ystart,
                                          const double * yend)
{
#  ifndef NO_EXPANSION
  if (!(dt > 0.))
    return;
  double rate = 0.;
#   if GAS_SOURCE_EXACT
  rate = gas_log_expansion (point, ystart, yend);
#   else
  double invMW0 = 0., dinvMW = 0.;
  for (int jj = 0; jj < NGS; jj++) {
    invMW0 += (ystart[jj] > 0. ? ystart[jj] : 0.)/gas_MWs[jj];
    dinvMW += (yend[jj] - ystart[jj])/gas_MWs[jj];
  }
  if (!(invMW0 > 0.) || !(ystart[NGS] > 0.))
    return;
  rate = (yend[NGS] - ystart[NGS])/ystart[NGS] + dinvMW/invMW0;
#   endif
  drhodt[] -= (1. - f[])*cm[]*rate/dt;
#  endif // !NO_EXPANSION
}
# endif // VARPROP
#endif // GAS_CHEMISTRY_STRANG && !TURN_OFF_GAS_REACTIONS

#if !defined(BINNING) && !TURN_OFF_GAS_REACTIONS
/**
## The sweep of the gas-phase reactions

This function holds the loop of the gas-phase reactions of the external gas.
`dtc` is the time over which the reactor integrates. The `chemistry` event
gives `dt`, or `dt/2` under `GAS_CHEMISTRY_STRANG`. The expansion source
always divides by the step `dt`, so each half gives its part of the mean of
the step. `second_half` is true only for the second half of the Strang split.
That half writes its expansion directly into `drhodt`, because
`update_divergence()` has already run. See "The Strang split of the gas
chemistry" above. */

static void gas_phase_reactions (double dtc, bool second_half)
{
#if FROZEN_CELL_GATE
  foreach (reduction(+:frozen_cell_gate_n)) {
#else
  foreach() {
#endif
    if (f[] < 1. - F_ERR) {
      double temperature = TG[]/(1. - f[]);
      if (!(temperature > 273.) || !(temperature < 3500.))
        continue;

      // Freshly-uncovered cells can carry an all-zero composition: the RHS
      // clamps each species to >= 0, so an empty vector reaches the
      // mole-fraction conversion as MW = 1/0 (mirrors the solid-branch gate).
      double ygsum_seed = 0.;
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        ygsum_seed += YG[];
      }
      if (!(ygsum_seed > 0.))
        continue;

      double y0ode[NGS + 1]; // NGS + T
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        y0ode[jj] = YG[]/(1. - f[]);
      }
      y0ode[NGS] = temperature;

      /**
      Keep the pre-reaction state: the step-averaged expansion source needs
      both ends of the step. The copy is local and is discarded below. */

#ifdef VARPROP
      double ystart[NGS + 1];
      for (int jj = 0; jj < NGS + 1; jj++)
        ystart[jj] = y0ode[jj];
#endif

      UserDataODE data;
      data.P = Pref + p[];
      data.T = y0ode[NGS];
      data.sources = NULL; // do not fill sources during integration; predict after the solve
# ifdef VARPROP
      data.rhog = rhoGv_G[];
      data.cpg = cpGv_G[];
# else
      data.rhog = rhoG;
      data.cpg = cpG;
# endif

#if FROZEN_CELL_GATE
      /**
      One evaluation decides whether this cell reacts at all. When it does not,
      the explicit update carries the whole change of the step, and the stiff
      solve has nothing to add. The gate then continues to the write-back
      below, so the expansion source and the fields stay on the same path as a
      solved cell. */

      bool frozen = false;
      if (dtc > 0.) {
        double dy_gate[NGS + 1];
        gas_batch_nonisothermal_constantpressure (y0ode, dtc, dy_gate, &data);

        double dYmax = 0.;
        for (int jj = 0; jj < NGS; jj++)
          dYmax = fmax (dYmax, fabs (dy_gate[jj]));

        if (dYmax*dtc < FROZEN_CELL_YTOL &&
            fabs (dy_gate[NGS])*dtc < FROZEN_CELL_TTOL) {
          frozen = true;
          frozen_cell_gate_n++;
          for (int jj = 0; jj < NGS + 1; jj++)
            y0ode[jj] += dtc*dy_gate[jj];
        }
      }

      if (!frozen)
#endif
      /**
        Using an explicit solver for gas-phase reactions is not
        recommended as they are usually stiff.
        */
      OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure, NGS + 1, dtc, y0ode, &data);

      /**
      Keep the pre-integration state if the solve diverged (see the solid
      branch above). */

      bool valid = true;
      for (int jj = 0; jj < NGS + 1; jj++)
        if (!isfinite (y0ode[jj]))
          valid = false;

      if (!valid)
        continue;

      /**
        The expansion source is taken over the whole step, as
        `(state_end - state_start)/dt`, which is conservative. Set
        `gas_source_averaged = false` to recover the older instantaneous form,
        evaluated at the converged end-of-step state. */

# ifdef VARPROP
#  if GAS_CHEMISTRY_STRANG
      if (second_half)
        strang_second_half_expansion (point, ystart, y0ode);
      else
#  endif
        accumulate_gas_sources (point, ystart, y0ode);
# endif

      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        YG[] = (y0ode[jj] > 0.) ? y0ode[jj]*(1. - f[]) : 0.;
      }
      TG[] = y0ode[NGS]*(1. - f[]);
    }
  }
}
#endif // !BINNING && !TURN_OFF_GAS_REACTIONS

event chemistry (i++) {

#ifdef CHEMISTRY_LOG
  reset ({t_solid, t_gas}, 0.);
  double time_mpi[npe()];

  for (int pe = 0; pe < npe(); pe++)
    time_mpi[pe] = 0.;

  struct timespec start, end;
  clock_gettime (CLOCK_MONOTONIC, &start);
#endif

#if FROZEN_CELL_GATE
  frozen_cell_gate_n = 0;
#endif

#ifdef SOLVE_TEMPERATURE
  odefunction batch = &solid_batch_nonisothermal_constantpressure;
  unsigned int NEQ = NGS + NSS + 1 + 1; //NGS + NSS + porosity + T
#else
  odefunction batch = &solid_batch_isothermal_constantpressure;
  unsigned int NEQ = NGS + NSS + 1;
#endif
  /**
  ## Solid-gas reactions
  We solve the solid-gas reaction system in each cell where there is
  solid present (i.e. f > F_ERR). The system is solved in terms of mass
  because the volume of the solid phase is variable due to porosity changes.
  */
  foreach ()
    if (f[] > F_ERR) {
      double temperature = TS[]/f[];
#if SOLID_SOURCE_DIAG
      solid_diag[] = 2.;   // a guard below can still skip this cell
#endif
      // Reject two FPE triggers before mutating state, both of which make the
      // gas-species mole-fraction conversion in the RHS divide by sum(y/MW)==0:
      //  - sliver-garbage temperature (TS/f outside a physical window);
      //  - an empty reactor seed (gasmass = YG/f * rhoGvh * porosity, built below).
      double rhoGvh_seed;
      #ifdef VARPROP
      rhoGvh_seed = rhoGv_S[];
      #else
      rhoGvh_seed = rhoG;
      #endif
      double ygsum_seed = 0.;
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_S[jj];
        ygsum_seed += YG[];
      }
      if (!(temperature > 273.) || !(temperature < 3500.) ||
          !(ygsum_seed > 0.) || !(porosity[] > 0.) || !(rhoGvh_seed > 0.))
        continue;

      porosity[] /= f[];
#if SOLID_SOURCE_DIAG
      solid_diag[] = 1.;
#endif

      double y0ode[NEQ];
      UserDataODE data;
      data.P = Pref + p[];
#ifdef VARPROP
      data.rhos = rhoSv[];
      data.rhog = rhoGv_S[];
#else
      data.rhos = rhoS;
      data.rhog = rhoG;
#endif
      data.zeta = zeta[];
#ifdef SOLVE_TEMPERATURE
# ifdef VARPROP
      data.cps = cpSv[];
      data.cpg = cpGv_S[];
# else
      data.cps = cpS;
      data.cpg = cpG;
# endif
#endif
      double sources[NEQ];
#ifdef STORE_SOURCES
      for (int jj = 0; jj < NEQ; jj++) 
        sources[jj] = 0.; // solid-species slots are never written
#endif
      data.sources = NULL; // do not fill sources during integration; predict after the solve

      double gasmass[NGS];
      double rhoGvh;
      #ifdef VARPROP
      rhoGvh = rhoGv_S[];
      #else
      rhoGvh = rhoG;
      #endif

      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_S[jj];
        gasmass[jj] = YG[]/f[]*rhoGvh*porosity[];
        y0ode[jj] = gasmass[jj];
      }

      double solidmass[NSS];
      for (int jj = 0; jj < NSS; jj++) {
        scalar YS = YSList[jj];
        solidmass[jj] = YS[]/f[]*rhoS*(1. - porosity[]);
        y0ode[jj+NGS] = solidmass[jj];
      }

      y0ode[NGS+NSS] = porosity[];

#ifdef SOLVE_TEMPERATURE
      y0ode[NGS+NSS+1] = TS[]/f[];
#endif

#ifdef EXPLICIT_REACTIONS
    // ODESolverEXP (batch, NEQ, dt, y0ode, &data);
      RungeKutta45EXP (batch, NEQ, dt, y0ode, &data);
#else //default
      OpenSMOKE_ODESolver (batch, NEQ, dt, y0ode, &data);
#endif

      /**
      A diverged solve returns a non-finite state; writing it back poisons the
      fields (and the next step's source prediction) far worse than losing one
      cell's reaction step, so keep the pre-integration state instead. */

      bool valid = true;
      for (int jj = 0; jj < NEQ; jj++)
        if (!isfinite (y0ode[jj]))
          valid = false;

      if (!valid) {
        porosity[] *= f[]; // undo the tracer-form conversion above
#if SOLID_SOURCE_DIAG
        solid_diag[] = 3.;
#endif
        continue;
      }

#if SOLID_SOURCE_DIAG && defined(SOLVE_TEMPERATURE)
      dTS_step[] = fabs (y0ode[NGS+NSS+1] - temperature);
#endif

      /**
      The source term is predicted once, at the converged end-of-step state
      (exact as dt -> 0), rather than being accumulated from the solver's
      internal RHS evaluations. */

      data.sources = sources;
      double dy_tmp[NEQ];
      batch (y0ode, dt, dy_tmp, &data);

      double totgasmass = 0;
      for (int jj = 0; jj < NGS; jj++)
        totgasmass += fmax (0., y0ode[jj]);

      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_S[jj];
        YG[] = (totgasmass < 1e-8) ? 0. : fmax (0., y0ode[jj])/totgasmass*f[];
      }

      double totsolidmass = 0;
      for (int jj = 0; jj < NSS; jj++)
        totsolidmass += fmax (0., y0ode[jj+NGS]);

      for (int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        YS[] = (totsolidmass < 1e-8) ? 0. : fmax (0., y0ode[jj+NGS])/totsolidmass*f[];
      }

      porosity[] = y0ode[NGS+NSS]*f[];

#ifdef VARPROP
      for (int jj = 0; jj < NGS; jj++) {
        scalar DYDtGjj = DYDtG_S[jj];
        DYDtGjj[] += sources[jj]*cm[];
      }
#endif

#ifdef SOLVE_TEMPERATURE
      TS[] = y0ode[NGS+NSS+1]*f[];
# ifdef VARPROP
      DTDtS[] += sources[NGS+NSS+1]*cm[];
# endif
#endif
      omega[] = sources[NGS+NSS];
#ifdef STORE_SOURCES
      for (int jj = 0; jj < NEQ; jj++) {
        scalar src = sourcesList[jj];
        src[] = sources[jj];
      }
#endif
    }

  /**
  ## Gas-phase reactions
  We solve the gas-phase reaction system in every cell that contains gas
  (i.e. f < 1 - F_ERR). The system is solved in terms of mass fraction.
  */

#ifdef BINNING
# ifndef VARPROP
#   error "BINNING requires VARPROP (it uses rhoGv_G/cpGv_G and the DYDtG_G/DTDtG sources)"
# endif

  /**
  The bin partitioning is driven by a case-provided list of thermochemical
  `targets` and a per-target tolerance `eps` (mixed-radix bin id, see
  binning.h). */

  extern scalar * targets;
  extern double * eps;

  /**
  Flag the pure-gas cells (the same set integrated by the non-binning branch)
  and convert their gas fields from VOF-tracer form (`Y*(1-f)`) to the actual
  mass fractions the reactor expects. */

  scalar gasmask[], rho0[], cp0[];
  foreach() {
    gasmask[] = (f[] < 1. - F_ERR && TG[] > 0.) ? 1. : 0.;
    rho0[] = rhoGv_G[], cp0[] = cpGv_G[];
    if (gasmask[]) {
      scale_gas_tracers (point, 1./(1. - f[]));

      /**
      First half of the step-averaged source: subtract the pre-reaction state.
      The loop below adds the post-reaction state over the same `gasmask`, so
      every cell that gets the subtraction also gets the addition. `rho0` and
      `cp0` keep the weights that `binning_remap()` is about to overwrite. */

      if (gas_source_averaged) {
        double ystart[NGS + 1];
        for (int jj = 0; jj < NGS; jj++) {
          scalar YG = YGList_G[jj];
          ystart[jj] = YG[];
        }
        ystart[NGS] = TG[];
        gas_sources_accumulate_state (point, ystart, rho0[], cp0[], -1.);
      }
    }
  }

  /**
  Agglomerate the flagged cells into bins of similar thermochemical state and
  integrate the stiff chemistry ODE once per bin. `bin->phi[j]` holds the
  mass-averaged value of `fields[j]`: entries `[0..NGS-1]` are the gas species
  and entry `[NGS]` is the temperature. */

  scalar * fields = list_concat (YGList_G, {TG});

  BinTable * table = binning (fields, targets, eps, rhoGv_G, cpGv_G, gasmask);

#ifdef CHEMISTRY_LOG
  static FILE * fp = NULL;
  if (!fp) {
    char name[20];
    sprintf (name, "bin-%d", pid());
    fp = fopen (name, "w");
  }
  fprintf (fp, "%g %ld %ld\n", t, grid->n, binning_stats(table).nactive);
  fflush (fp);
#endif

  foreach_bin (table) {
    double y0ode[NGS + 1];
    for (size_t j = 0; j < bin->nfields; j++)
      y0ode[j] = bin->phi[j];

    UserDataODE data;
    data.P = Pref + bin_average (bin, p);
    data.sources = NULL;
    data.rhog = bin->rho;
    data.cpg = bin->cp;

    OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
        NGS + 1, dt, y0ode, &data);

    for (size_t j = 0; j < bin->nfields; j++)
      bin->phi[j] = (j < (size_t)NGS) ? fmax (0., y0ode[j]) : y0ode[j];

    bin->rho = data.rhog;
    bin->cp = data.cpg;
  }

  binning_remap (table, fields, rhoGv_G, cpGv_G);
  binning_cleanup (table);
  free (fields), fields = NULL;

  /**
  Second half of the step-averaged source: add the post-reaction state. Then
  restore the VOF-tracer form of the gas fields. */

  foreach() {
    if (gasmask[]) {
      double y0ode[NGS + 1]; // NGS + T
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        y0ode[jj] = YG[];
      }
      y0ode[NGS] = TG[];

      if (gas_source_averaged)
        gas_sources_accumulate_state (point, y0ode, rho0[], cp0[], +1.);
      else
        gas_sources_instantaneous (point, y0ode);

      scale_gas_tracers (point, 1. - f[]); // restore VOF-tracer form
    }
  }
#else // !BINNING

  /**
  `TURN_OFF_GAS_REACTIONS` removes this whole loop, so no gas cell reacts.
  The gas sources that the loop feeds are reset in every step, so they stay
  at 0. The pore gas of the solid branch is switched in `reactors.h`. */

# if !TURN_OFF_GAS_REACTIONS
  gas_phase_reactions (GAS_CHEMISTRY_STRANG ? 0.5*dt : dt, false);
# endif // !TURN_OFF_GAS_REACTIONS
#endif // BINNING

#ifdef CHEMISTRY_LOG
  clock_gettime (CLOCK_MONOTONIC, &end);
  time_mpi[pid()] = (end.tv_sec - start.tv_sec) +
                    (end.tv_nsec - start.tv_nsec)*1e-9;
@if _MPI
  if (pid() == 0) {
    MPI_Reduce(MPI_IN_PLACE, time_mpi, npe(), MPI_DOUBLE,
        MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(time_mpi, NULL, npe(), MPI_DOUBLE,
        MPI_SUM, 0, MPI_COMM_WORLD);
  }
@endif

  fprintf (stderr, "%g ", t);

  for (int pe = 0; pe < npe(); pe++)
    fprintf (stderr, "%g ", time_mpi[pe]);

  fprintf (stderr, "\n");
#endif
}

#if GAS_CHEMISTRY_STRANG && !defined(BINNING) && !TURN_OFF_GAS_REACTIONS
/**
## The second half of the Strang split

The split needs the averaged form of the expansion, see above. Stop the run
at the start if a case has changed it in `main()`. */

# if defined(VARPROP) && !GAS_SOURCE_EXACT
event init (i = 0) {
  if (!gas_source_averaged) {
    fprintf (stderr, "GAS_CHEMISTRY_STRANG needs gas_source_averaged ="
                     " true. Stop.\n");
    exit (1);
  }
}
# endif

/**
This event runs after the `tracer_diffusion` events of
`multicomponent-varprop.h`, because `chemistry.h` comes before them and
same-name events run in reverse order of the declaration. It runs before the
`tracer_diffusion` event of `shrinking.h`, which the case includes before
`multicomponent-varprop.h`. Check the order with `qcc -events`. */

event tracer_diffusion (i++) {

  if (!(dt > 0.))
    return 0;

# ifdef CHEMISTRY_LOG
  struct timespec s2start, s2end;
  clock_gettime (CLOCK_MONOTONIC, &s2start);
# endif

  /**
  The probe records the state before this half. `TG` is in tracer form, so
  the intrinsic value is `TG/(1-f)`. */

  if (strang_probe_armed) {
    scalar s2TG[];
    double Tmax = -HUGE;
    foreach (reduction(max:Tmax)) {
      s2TG[] = TG[];
      double fG = 1. - f[];
      if (fG > F_ERR)
        Tmax = max (Tmax, TG[]/fG);
    }

    gas_phase_reactions (0.5*dt, true);

    double Q2 = 0.;
    foreach (reduction(+:Q2)) {
      double fG = 1. - f[];
      if (fG > F_ERR) {
# ifdef VARPROP
        double rc = rhoGv_G[]*cpGv_G[];
# else
        double rc = rhoG*cpG;
# endif
        Q2 += rc*(TG[] - s2TG[])/dt*dv();
      }
    }
    strang_Tmax_tr = (Tmax > -HUGE) ? Tmax : 0.;
    strang_Q2 = Q2;
  }
  else
    gas_phase_reactions (0.5*dt, true);

  /**
  The output and `adapt` read `T`. The gas reactor changes `TG` only. */

# ifdef SOLVE_TEMPERATURE
  foreach()
    T[] = TS[] + TG[];
# endif

# ifdef CHEMISTRY_LOG
  clock_gettime (CLOCK_MONOTONIC, &s2end);
  double s2time = (s2end.tv_sec - s2start.tv_sec) +
                  (s2end.tv_nsec - s2start.tv_nsec)*1e-9;
  mpi_all_reduce (s2time, MPI_DOUBLE, MPI_SUM);
  if (pid() == 0)
    fprintf (stderr, "S2 %g %g\n", t, s2time);
# endif
}
#endif // GAS_CHEMISTRY_STRANG && !BINNING && !TURN_OFF_GAS_REACTIONS

#endif // TURN_OFF_REACTIONS
