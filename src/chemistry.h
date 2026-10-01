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
`GAS_CHEMISTRY_STRANG` selects the Strang split of the gas chemistry, and
`FROZEN_CELL_GATE` selects the frozen-cell gate of the gas branch. Both are 1
by default. Set one to 0 to get the reference path for a comparison. See "The
Strang split of the gas chemistry" and "The gate that skips the cells the
step cannot change" below. The defaults stay outside `TURN_OFF_REACTIONS`, so
a case can print the flags in every build.

The `BINNING` path uses the same Strang split. It does not use the gate. */

#ifndef GAS_CHEMISTRY_STRANG
# define GAS_CHEMISTRY_STRANG 1
#endif

#ifndef FROZEN_CELL_GATE
# define FROZEN_CELL_GATE 1
#endif

/**
`binning.h` stops the build if a case includes it after this file. */

#define CHEMISTRY_H 1

#if defined(BINNING) && !defined(VARPROP)
# error "BINNING requires VARPROP (it uses rhoGv_G/cpGv_G and drhodt_chem)"
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

event reset_sources (i++) {
  foreach() {
    omega[] = 0.;
  }
}

#ifdef CHEMISTRY_LOG
scalar t_solid[], t_gas[];
#endif

#ifdef VARPROP
/**
## The gas-phase reaction source of the low-Mach divergence

At constant pressure the expansion rate of the gas is `-d(ln rho)/dt`, so
its exact mean over the step is

    ln(rho_start/rho_end)/dt

Both densities follow from the ideal gas law with `1/MW = sum_j Y_j/MW_j`.
The form thus needs no call to the property library, and it has no
denominator from another time level. The chemistry event writes
`cm[]*ln(rho_start/rho_end)/dt` to `drhodt_chem`, per unit volume of gas.
`update_divergence()` adds it to `divu2` with the same `(1-f)` weight as the
other terms. `test/gas-source-cell.c` compares it with a sub-stepped
reference.

The filter of the divergence source in `navier-stokes/centered-phasechange.h`
goes with this form. Set `gas_source_filter_passes = 0` to remove the filter.

The exact form covers the external gas only. The pore gas inside the solid
keeps the source vector of the reactor, because its temperature equation
carries the heat capacity of the solid and of the gas together.

Caution: `TURN_OFF_HEAT_OF_REACTION` zeroes the temperature increment of the
reactor, so the heat release also leaves the expansion.

Caution: do not read the end-state density from `data.rhog` after the solve.
`reactors.h` starts the reactor with `UserDataODE data = *(UserDataODE *)args`,
so the reactor writes `rhog` in a local copy and the caller keeps the old
value.

The `BINNING` path uses the same form. It keeps the start density of each
cell in a field, because the reactor of a bin does not give the state of a
cell.

The functions below exist only when the gas reactor runs
(`TURN_OFF_GAS_REACTIONS` 0). */

#if !TURN_OFF_GAS_REACTIONS

/**
The gas density from the ideal gas law, with `invMW = 1/MW = sum_j Y_j/MW_j`.
This needs no call to the property library. Returns 0 if the state is not
usable, and the caller then adds no expansion. */

static double gas_ideal_density (double P, double invMW, double T) {
  if (!(invMW > 0.) || !(T > 0.))
    return 0.;
  return P/(R_GAS*1000.*T*invMW);
}

#ifndef BINNING
/**
The gas density of the state `ys` (`NGS` mass fractions, then the
temperature). */

static double gas_state_density (Point point, const double * ys) {
  double invMW = 0.;
  for (int jj = 0; jj < NGS; jj++)
    invMW += (ys[jj] > 0. ? ys[jj] : 0.)/gas_MWs[jj];
  return gas_ideal_density (Pref + p[], invMW, ys[NGS]);
}
#endif // !BINNING

/**
`ln(rho_start/rho_end)` from the two densities. Returns 0 if one of them is
not usable. */

static double gas_log_density_ratio (double rho_start, double rho_end) {
  if (!(rho_start > 0.) || !(rho_end > 0.))
    return 0.;
  return log (rho_start/rho_end);
}

/**
The exact step mean of the expansion, `ln(rho_start/rho_end)`, from the two
states of the reactor. Returns 0 if one of the states is not usable. */

#ifndef BINNING
static double gas_log_expansion (Point point, const double * ystart,
                                 const double * yend) {
  return gas_log_density_ratio (gas_state_density (point, ystart),
                                gas_state_density (point, yend));
}
#endif // !BINNING

/**
Add the increment `ln(rho_start/rho_end)` of one call of the gas reactor to
the expansion of the step. The increment divides by the step `dt`, also when
the reactor runs over `dt/2`. The first call of the step writes to
`drhodt_chem`. The second half of the Strang split writes to `drhodt`,
because `update_divergence()` has already run. See "The Strang split of the
gas chemistry" below for the weights. */

static void add_gas_expansion (Point point, double lnratio, bool second_half)
{
  if (!(dt > 0.))
    return;
#if GAS_CHEMISTRY_STRANG
  if (second_half) {
# ifndef NO_EXPANSION
    drhodt[] -= (1. - f[])*cm[]*lnratio/dt;
# endif
    return;
  }
#endif
  drhodt_chem[] += cm[]*lnratio/dt;
}
#endif // !TURN_OFF_GAS_REACTIONS
#endif // VARPROP

/**
## The gate that skips the cells the step cannot change

The gate spends **one** evaluation of the reactor right-hand side to decide
whether the stiff solve of a gas cell can change the state over the step.
With `biomass/Solid-gas-88` the gas branch solves 88 equations in every cell
that holds gas. In the free stream the mixture is air that does not react,
but the stiff solve still costs 1.8 to 3.8 ms, about 100 evaluations. When
the change of the step is below the tolerances, the gate applies the explicit
update, which is exact at that size, and skips the solve.

The margin is large. Air at 1123 K gives `max|dY|` of 2e-17 and `|dT|` of
3e-13 K over a step of 2.4e-4 s. Hot air with 0.5 per cent of CO and 0.5 per
cent of tar gives 2e-4 and 3e-2 K.

An ignition does not close the gate. `test/gate-ignition.c` marches a batch
reactor with the pyrolysis gas in air from 800 K to 1300 K. The gate stays
open on every one of the 3000 steps, and the gated trajectory is equal to the
ungated one to the last bit, also through an induction of 0.48 s at 800 K.

The tolerances are the size of the state that the gate throws away on one
step, thus of the perturbation that it feeds to the solver. At
`FROZEN_CELL_YTOL` = 1e-10 a 2-D run at level 8 kept the mass, `Tmax` and
`dt`, but the projection residual moved by 5e-5 from t = 0.33 s. At 1e-15 the
run is equal to the ungated run in every column. The tolerances below are
40 times above the rate of the free stream, so the gate keeps its work and
the trajectory stays reproducible.

`FROZEN_CELL_GATE` 0 removes the gate. Every gas cell then goes to the stiff
solve.

A 2-D `fatehi-combustion` at level 7 with 88 species reaches 1.28 times the
simulated time of the ungated case over a fixed wall budget (1.43 times at
`FROZEN_CELL_YTOL` = 1e-10). Caution: compare two builds only inside one batch
of runs. The absolute time swings 15 per cent from one run to the next. */

#ifndef FROZEN_CELL_YTOL
# define FROZEN_CELL_YTOL 1e-15
#endif

#ifndef FROZEN_CELL_TTOL
# define FROZEN_CELL_TTOL 1e-11
#endif

/**
The number of cells that the gate skipped over the last step, over both
halves of the Strang split. The value is the total over every rank and every
thread. With `CHEMISTRY_LOG` the last sweep of the step prints it on a line
that starts with `GATE`. A case can also print it in its own log. */

int frozen_cell_gate_n = 0;

/**
## The Strang split of the gas chemistry

A Lie split integrates the gas reactor over the full `dt` at the start of
the step, and the advection and the implicit diffusion follow. At level 10
this split carries the whole `Tmax(dt)` law: the jump of `Tmax` over the
chemistry is 50.9, 28.8, 15.7 and 8.4 K at `dt` 4e-4, 2e-4, 1e-4 and 5e-5 s,
thus first order in `dt`.

The gas chemistry therefore uses a symmetric (Strang) split when
`GAS_CHEMISTRY_STRANG` is 1 (the default):

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
  `porosity` stay on the full `dt` in the first call. The pore gas
  reacts inside the same stiff system as the solid, with the heat capacity of
  the solid in the temperature equation, so a split of the pore gas needs a
  split of the whole solid reactor. The solid changes slowly: `TS` changes
  by 0.05 to 0.11 K per step, and the split error of the solid is 0.1 to
  0.3 % of `omega`. So `omega`, `zeta`, `prod`, `ubf` and `gas_source` come
  from one full-dt solve at the start of the step. The projection of the
  step reads the same `gas_source`.

The expansion of the gas reactions (the chemistry part of `drhodt`):

* Both halves use the exact form `ln(rho_start/rho_end)/dt` of the half.
  The form divides by the step `dt`, not by `dt/2`, so the increment of
  each half becomes its share of the mean rate of the step.
* The first half adds its increment to `drhodt_chem`.
  `update_divergence()` then puts it into `drhodt`, with no change.
* The second half runs after `update_divergence()`. It therefore adds its
  share directly to `drhodt`, with the weight `(1-f)` and the factor `cm`
  that `update_divergence()` gives the chemistry part. The projection of the
  same step reads `drhodt` in `advection_term` and in `projection`, and both
  run after this event. So the projection of step n receives the sum of the
  two half increments divided by `dt`. No increment enters two
  projections, and no increment is lost.

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
  this step reads the properties of the start of the step. `adapt`
  then computes the properties from the state after the second half, so the
  next step starts with consistent properties.

What the split can change, and what it cannot:

* At a constant `dt` the Strang sequence is the Lie sequence with one more
  half step at the end: `R_h T R_h R_h T R_h = R_h T R T R_h`. So the state
  at the end of a Strang step is `R(dt/2)` applied to the state after the
  transport of a Lie run. The split changes the time level at which the
  output, `adapt` and the probes read the state. It also moves the
  expansion of each half into the projection of its own step. It does not
  change the sequence of the reactor and the transport. So expect a
  smaller `Tmax(dt)` law with the split, not zero.
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
over a shorter interval. At level 8 before the ignition, the chemistry cost
1.64 times the Lie split. The solid reactor does not change. With
`CHEMISTRY_LOG` the second half prints its time on a line that starts with
`S2`.

Restart: the split adds no field.

The gate tests each half with its own `dt/2`. The `BINNING` path uses the
same split, with no gate: it bins the cells again in each half.

`GAS_CHEMISTRY_STRANG` 0 gives the Lie split: the `chemistry` event
integrates the gas reactor over the full `dt`, and the build has no second
half. */

#if !TURN_OFF_GAS_REACTIONS
/**
## The cells where the gas reactor runs

The gas reactor runs in a cell that holds gas (`f < 1 - F_ERR`), when the
temperature of the gas is in the window of the thermodynamic data and the
gas composition is not empty. A cell that the interface has just uncovered
can carry an all-zero composition. The right-hand side clamps each species
to 0 or more, so an empty vector reaches the conversion to mole fractions as
`MW = 1/0`. The solid branch has the same check. The per-cell path and the
`BINNING` path use this function, so they select the same cells. */

static bool gas_reactor_cell (Point point)
{
  if (!(f[] < 1. - F_ERR))
    return false;

  double temperature = TG[]/(1. - f[]);
  if (!(temperature > 273.) || !(temperature < 3500.))
    return false;

  double ygsum = 0.;
  for (int jj = 0; jj < NGS; jj++) {
    scalar YG = YGList_G[jj];
    ygsum += YG[];
  }
  return ygsum > 0.;
}
#endif // !TURN_OFF_GAS_REACTIONS

#if !defined(BINNING) && !TURN_OFF_GAS_REACTIONS
/**
## The sweep of the gas-phase reactions

This function holds the loop of the gas-phase reactions of the external gas.
`dtc` is the time over which the reactor integrates: `dt/2` for each half
of the Strang split, `dt` for the Lie split. The expansion source
always divides by the step `dt`, so each half gives its part of the mean of
the step. `second_half` is true only for the second half of the Strang split.
That half writes its expansion directly into `drhodt`, because
`update_divergence()` has already run. See "The Strang split of the gas
chemistry" above. */

static void gas_phase_reactions (double dtc, bool second_half)
{
  foreach (reduction(+:frozen_cell_gate_n)) {
    if (gas_reactor_cell (point)) {
      double temperature = TG[]/(1. - f[]);

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

      /**
      One evaluation decides whether this cell reacts at all. When it does not,
      the explicit update carries the whole change of the step, and the stiff
      solve has nothing to add. The gate then continues to the write-back
      below, so the expansion source and the fields stay on the same path as a
      solved cell. */

      bool frozen = false;
#if FROZEN_CELL_GATE
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
#endif // FROZEN_CELL_GATE

      if (!frozen)
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
        The expansion source comes from the two ends of the half,
        `ln(rho_start/rho_end)/dt`. */

# ifdef VARPROP
      add_gas_expansion (point, gas_log_expansion (point, ystart, y0ode),
                         second_half);
# endif

      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        YG[] = (y0ode[jj] > 0.) ? y0ode[jj]*(1. - f[]) : 0.;
      }
      TG[] = y0ode[NGS]*(1. - f[]);
    }
  }

# ifdef CHEMISTRY_LOG
  if ((second_half || !GAS_CHEMISTRY_STRANG) && pid() == 0)
    fprintf (stderr, "GATE %g %d\n", t, frozen_cell_gate_n);
# endif
}
#endif // !BINNING && !TURN_OFF_GAS_REACTIONS

#if defined(BINNING) && !TURN_OFF_GAS_REACTIONS
/**
## The sweep of the gas-phase reactions with binning

This function replaces the per-cell sweep above when the case includes
`binning.h`. It has the same arguments, the same cells
(`gas_reactor_cell()`) and the same expansion source
(`add_gas_expansion()`), so the Strang split and `TURN_OFF_GAS_REACTIONS`
work in the same way.

The reaction step follows section 2 of Goldin, Ren and Zahirovic, Combust.
Theory Model. 13 (2009) 721 (`papers/Goldin2009.md`):

1. `binning()` puts the cells with a similar state into one bin. The case
   gives the `targets` and a tolerance `eps` for each target (see
   `binning.h`).
2. The state of a bin is the mean of the cells, with the weight of the gas
   mass `rho*(1-f)*dV` (Eq. 3). The mean covers the mass fractions and the
   total enthalpy `h`, not the temperature. The temperature of the bin comes
   from `h(T, Y) = h` with a Newton solve.
3. The stiff solver integrates the bin over `dtc` (Eq. 4).
4. Each cell gets the change of its bin (Eq. 5), scaled by one factor
   `alpha` (0 to 1) for the whole bin. `alpha` is the largest value that
   keeps every mass fraction of every cell of the bin in [0, 1]. One factor
   for all species keeps the elements and the sum of the mass fractions.
5. Each cell gets the change of the total enthalpy of its bin, also scaled
   by `alpha`, and its temperature comes from a Newton solve. For an
   adiabatic reactor the change is zero, and each cell keeps its own total
   enthalpy: this is Eq. 6. The gas reactor of `reactors.h` also has the
   optically thin radiation, so the change of a bin is its radiation loss.
   With `TURN_OFF_HEAT_OF_REACTION` the reactor does not change the
   temperature, so the cells keep their temperature.

When the bins hold one cell each, the step is equal to the per-cell step,
apart from the gate and the tolerance of the Newton solves.

The gate of the per-cell sweep does not apply here. Thus
`frozen_cell_gate_n` stays 0. */

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

/**
The gas density of the current cell, from its fields. Call it only between
the two calls of `scale_gas_tracers()`, where the fields hold the actual
values.

Caution: do not replace this function with a point function that copies the
state of the cell into a local array, followed by a call of
`gas_state_density()` on that array in the same `foreach`. That pattern
makes `qcc` crash with a segmentation fault (in `replace_break()`). */

static double gas_cell_density (Point point) {
  double invMW = 0.;
  for (int jj = 0; jj < NGS; jj++) {
    scalar YG = YGList_G[jj];
    invMW += (YG[] > 0. ? YG[] : 0.)/gas_MWs[jj];
  }
  return gas_ideal_density (Pref + p[], invMW, TG[]);
}

/**
The total (sensible plus formation) enthalpy in J/kg of the gas with the mass
fractions `y` at the temperature `T`. A negative mass fraction counts as 0. */

static double gas_mass_enthalpy (const double * y, double T)
{
  double yc[NGS], x[NGS], MW;
  for (int jj = 0; jj < NGS; jj++)
    yc[jj] = y[jj] > 0. ? y[jj] : 0.;
  OpenSMOKE_MoleFractions_From_MassFractions (x, &MW, yc);
  OpenSMOKE_GasProp_SetTemperature (T);
  return OpenSMOKE_GasProp_Enthalpy (x);
}

/**
The temperature at which the gas with the mass fractions `y` has the total
enthalpy `h`, from a Newton solve that starts at `Tguess`. `dh/dT` is the
`cp` of the mixture. The solve stays in [`TMIN_ODE`, `TMAX_ODE`], the range
of the reactor. It returns false if it does not converge, and then `*T` is
`Tguess`. */

#define GAS_ENTHALPY_NEWTON_MAXIT 50
#define GAS_ENTHALPY_NEWTON_RTOL 1e-10

static bool gas_temperature_from_enthalpy (double h, const double * y,
                                           double Tguess, double * T)
{
  double yc[NGS], x[NGS], MW;
  for (int jj = 0; jj < NGS; jj++)
    yc[jj] = y[jj] > 0. ? y[jj] : 0.;
  OpenSMOKE_MoleFractions_From_MassFractions (x, &MW, yc);

  double Tn = clamp_temperature (Tguess);
  for (int it = 0; it < GAS_ENTHALPY_NEWTON_MAXIT; it++) {
    OpenSMOKE_GasProp_SetTemperature (Tn);
    double res = OpenSMOKE_GasProp_Enthalpy (x) - h;
    double cp = OpenSMOKE_GasProp_HeatCapacity (x);
    if (!(cp > 0.) || !isfinite (res))
      break;
    double dT = - res/cp;
    Tn = clamp_temperature (Tn + dT);
    if (fabs (dT) < GAS_ENTHALPY_NEWTON_RTOL*Tn) {
      *T = Tn;
      return true;
    }
  }
  *T = Tguess;
  return false;
}

/**
The statistics of the last call of the binning sweep, over all ranks:
`binning_alpha_min` is the smallest `alpha` of the bins,
`binning_nlimited` is the number of bins with `alpha < 1` (see
`BINNING_ALPHA_ROUNDOFF` below), `binning_nbadbin`
is the number of bins that the sweep did not change (the solve diverged or
the Newton solve of the bin did not converge), and `binning_nbadcell` is the
number of cells where the Newton solve did not converge. With
`CHEMISTRY_LOG` they go to `bin-<pid>` (values of that rank). */

double binning_alpha_min = 1.;

/**
A bin with one cell, or with cells of equal state, can get `alpha` below 1
by a rounding error only (the mean of the bin differs from the cell in the
last bit). `binning_nlimited` does not count a bin with `alpha` above
`1 - BINNING_ALPHA_ROUNDOFF`. */

#define BINNING_ALPHA_ROUNDOFF 1e-10

/**
The margin of the realizability condition of `alpha`. Section 2 of Goldin et
al. keeps every mass fraction of every cell in [0, 1]. With that strict
condition, a change far below the tolerance of the stiff solver sets
`alpha`:

* In the free stream of `fatehi-combustion` (dummy mechanism, level 8),
  one cell of the air bin has `H2O = 0` and the bin mean is 6e-14. The bin
  takes about 1e-20 of H2O (the reverse reaction of `H2 + 0.5 O2 = H2O`), so
  `alpha = 0` and the whole air bin does not react.
* In a bin near the particle, H2 at 1.5e-20 in one cell gives `alpha = 0.70`,
  so the heat release of the whole bin is 30 per cent too slow.

The condition therefore lets a cell go out of [0, 1] by `BINNING_YTOL`, and
the clamp removes the excess. The default is the absolute tolerance of the
OpenSMOKE ODE solver (1e-10, `ODE_Parameters`). The per-cell sweep clamps
the output of the solver to 0 in the same way, so both paths have an error of
the same size in the conservation. A real limit (a bin that takes more of a
species than one of its cells has) still gives `alpha < 1`. Set
`BINNING_YTOL` to 0 for the strict condition of the paper. */

#ifndef BINNING_YTOL
# define BINNING_YTOL 1e-10
#endif
long binning_nlimited = 0, binning_nbadbin = 0, binning_nbadcell = 0;

static void gas_phase_reactions (double dtc, bool second_half)
{
  extern scalar * targets;
  extern double * eps;

  /**
  Flag the cells where the reactor runs. Change their gas fields from the
  tracer form to the actual values. Keep the start density of each cell for
  the expansion source, the gas mass `wG` for the weight of the means, and
  the total enthalpy `hG`. */

  scalar gasmask[], rho0[], wG[], hG[];
  foreach() {
    gasmask[] = gas_reactor_cell (point) ? 1. : 0.;
    rho0[] = 0., wG[] = 0., hG[] = 0.;
    if (gasmask[]) {
      scale_gas_tracers (point, 1./(1. - f[]));
      rho0[] = gas_cell_density (point);
      wG[] = rho0[]*(1. - f[]);
      double ys[NGS];
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        ys[jj] = YG[];
      }
      hG[] = gas_mass_enthalpy (ys, TG[]);
    }
  }

  /**
  The entries of `bin->phi`: `[0..NGS-1]` the mass fractions, `[NGS]` the
  total enthalpy, `[NGS+1]` the temperature. The mean temperature is only the
  first guess of the Newton solve of the bin. The sweep then sets
  `bin->phi0[NGS+1]` to the temperature of the bin. */

  const int ih = NGS, iT = NGS + 1;
  scalar * fields = list_concat (YGList_G, {hG, TG});

  BinTable * table = binning (fields, targets, eps, wG, unity, gasmask);

  /**
  `alpha` of each bin. A bin with `alpha < 0` is a bin that the sweep does
  not change. */

  double * alpha = malloc ((table->nbins > 0 ? table->nbins : 1)*sizeof (double));

  foreach_bin (table) {
    alpha[i] = -1.;

    double ybar[NGS];
    for (int jj = 0; jj < NGS; jj++)
      ybar[jj] = bin->phi[jj];
    double Tbar;
    if (!gas_temperature_from_enthalpy (bin->phi[ih], ybar, bin->phi[iT],
                                        &Tbar))
      continue;
    bin->phi0[iT] = Tbar;

    double y0ode[NGS + 1];
    for (int jj = 0; jj < NGS; jj++)
      y0ode[jj] = ybar[jj];
    y0ode[NGS] = Tbar;

    /**
    `rho` and `cp` are start values only. Under `VARPROP` the reactor
    computes them again from the state at each evaluation, in a local copy
    of `data`. So do not read them back from `data` after the solve. */

    UserDataODE data;
    data.P = Pref + bin_average (bin, p);
    data.T = Tbar;
    data.sources = NULL;
    double invMW = 0.;
    for (int jj = 0; jj < NGS; jj++)
      invMW += (ybar[jj] > 0. ? ybar[jj] : 0.)/gas_MWs[jj];
    data.rhog = gas_ideal_density (data.P, invMW, Tbar);
    double xbar[NGS], MWbar;
    OpenSMOKE_MoleFractions_From_MassFractions (xbar, &MWbar, y0ode);
    OpenSMOKE_GasProp_SetTemperature (Tbar);
    data.cpg = OpenSMOKE_GasProp_HeatCapacity (xbar);

    OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
        NGS + 1, dtc, y0ode, &data);

    /**
    If the solve diverged, the sweep does not change the cells of this bin,
    as in the per-cell sweep.

    Clamp the mass fractions of the bin to 0 or more, as the per-cell sweep
    does. The stiff solver can return a negative value of the size of its
    tolerance for a species that is absent. Without the clamp, a cell where
    that species is exactly 0 gives `alpha = 0`, and the whole bin does not
    react. The clamp changes the sum of the mass fractions by the tolerance
    of the solver only. */

    bool valid = true;
    for (int jj = 0; jj < NGS + 1; jj++)
      if (!isfinite (y0ode[jj]))
        valid = false;
    if (!valid)
      continue;

    for (int jj = 0; jj < NGS; jj++) {
      y0ode[jj] = fmax (0., y0ode[jj]);
      bin->phi[jj] = y0ode[jj];
    }
    bin->phi[ih] = gas_mass_enthalpy (y0ode, y0ode[NGS]);
    bin->phi[iT] = y0ode[NGS];
    alpha[i] = 1.;
  }

  /**
  ## The remap (Eq. 5 and 6)

  First the largest `alpha` that keeps every mass fraction of every cell of
  the bin in [`-BINNING_YTOL`, `1 + BINNING_YTOL`], then the change of each
  cell, with a clamp to [0, 1]. See `BINNING_YTOL` for the margin. */

  double amin = 1.;
  long nlimited = 0, nbadbin = 0, nbadcell = 0;
  foreach_bin (table) {
    if (alpha[i] < 0.) {
      nbadbin++;
      continue;
    }

    double a = 1.;
    foreach_bin_cell (bin)
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        double dY = bin->phi[jj] - bin->phi0[jj];
        if (dY < 0. && YG[] + a*dY < - BINNING_YTOL)
          a = (max (YG[], 0.) + BINNING_YTOL)/(- dY);
        else if (dY > 0. && YG[] + a*dY > 1. + BINNING_YTOL)
          a = (1. + BINNING_YTOL - min (YG[], 1.))/dY;
      }
    alpha[i] = a;
    amin = min (amin, a);
    if (a < 1. - BINNING_ALPHA_ROUNDOFF)
      nlimited++;

    foreach_bin_cell (bin) {
      double ycell[NGS];
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        YG[] += a*(bin->phi[jj] - bin->phi0[jj]);
        YG[] = clamp (YG[], 0., 1.);
        ycell[jj] = YG[];
      }
#ifndef TURN_OFF_HEAT_OF_REACTION
      double hnew = hG[] + a*(bin->phi[ih] - bin->phi0[ih]);
      double Tguess = TG[] + a*(bin->phi[iT] - bin->phi0[iT]);
      double Tnew;
      if (!gas_temperature_from_enthalpy (hnew, ycell, Tguess, &Tnew))
        nbadcell++;
      TG[] = Tnew;
#endif
    }
  }
  free (alpha);

# ifdef CHEMISTRY_LOG
  /**
  The file `bin-<pid>` gets one line for each call: the time, the number of
  cells, the number of bins, 1 for the second half of the Strang split (else
  0), the smallest `alpha`, the number of bins with `alpha < 1`, the number
  of bins that the sweep did not change, and the number of cells where the
  Newton solve did not converge. All values are those of this rank. */

  static FILE * fp = NULL;
  if (!fp) {
    char name[20];
    sprintf (name, "bin-%d", pid());
    fp = fopen (name, "w");
  }
  fprintf (fp, "%g %ld %zu %d %g %ld %ld %ld\n", t, grid->n,
           binning_stats (table).nactive, second_half ? 1 : 0,
           amin, nlimited, nbadbin, nbadcell);
  fflush (fp);
# endif

  mpi_all_reduce (amin, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (nlimited, MPI_LONG, MPI_SUM);
  mpi_all_reduce (nbadbin, MPI_LONG, MPI_SUM);
  mpi_all_reduce (nbadcell, MPI_LONG, MPI_SUM);
  binning_alpha_min = amin;
  binning_nlimited = nlimited;
  binning_nbadbin = nbadbin;
  binning_nbadcell = nbadcell;

  binning_cleanup (table);
  free (fields), fields = NULL;

  /**
  Add the expansion of each cell from its two densities. Then change the gas
  fields back to the tracer form. */

  foreach() {
    if (gasmask[]) {
      double lnratio = gas_log_density_ratio (rho0[], gas_cell_density (point));
      add_gas_expansion (point, lnratio, second_half);

      scale_gas_tracers (point, 1. - f[]);
    }
  }
}
#endif // BINNING && !TURN_OFF_GAS_REACTIONS

event chemistry (i++) {

#ifdef CHEMISTRY_LOG
  reset ({t_solid, t_gas}, 0.);
  double time_mpi[npe()];

  for (int pe = 0; pe < npe(); pe++)
    time_mpi[pe] = 0.;

  struct timespec start, end;
  clock_gettime (CLOCK_MONOTONIC, &start);
#endif

  frozen_cell_gate_n = 0;

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
        continue;
      }

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
  With `BINNING` the solver integrates bins of cells, not single cells.

  `TURN_OFF_GAS_REACTIONS` removes this call, so no gas cell reacts.
  The gas sources that the call feeds are reset in every step, so they stay
  at 0. The pore gas of the solid branch is switched in `reactors.h`. */

#if !TURN_OFF_GAS_REACTIONS
  gas_phase_reactions (GAS_CHEMISTRY_STRANG ? 0.5*dt : dt, false);
#endif // !TURN_OFF_GAS_REACTIONS

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

#if GAS_CHEMISTRY_STRANG && !TURN_OFF_GAS_REACTIONS
/**
## The second half of the Strang split

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
#endif // GAS_CHEMISTRY_STRANG && !TURN_OFF_GAS_REACTIONS

#endif // TURN_OFF_REACTIONS
