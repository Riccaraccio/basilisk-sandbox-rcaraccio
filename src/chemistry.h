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

The `BINNING` path does not keep the start state of each cell, so it cannot
use the exact form. It keeps the step-averaged form
`rho*(state_end - state_start)/dt` (the default) or, with
`gas_source_averaged = false`, the rate at the end state. */

#ifdef BINNING
bool gas_source_averaged = true;

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
#endif // BINNING

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

/**
Convenience wrapper for the per-cell path, which holds both states and where
`rhoGv_G`/`cpGv_G` do not change over the chemistry event. */

static void accumulate_gas_sources (Point point, const double * ystart,
                                    const double * yend) {
  if (dt > 0.)
    drhodt_chem[] += cm[]*gas_log_expansion (point, ystart, yend)/dt;
}
#endif

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
The number of cells that the gate skipped over the last step. It makes the
gain of the gate visible in a production log. The `foreach` loop below carries
a `reduction` clause for it, so the value is the total over every rank and
every thread. */

int frozen_cell_gate_n = 0;

/**
## The Strang split of the gas chemistry

A Lie split integrates the gas reactor over the full `dt` at the start of
the step, and the advection and the implicit diffusion follow. At level 10
this split carries the whole `Tmax(dt)` law: the jump of `Tmax` over the
chemistry is 50.9, 28.8, 15.7 and 8.4 K at `dt` 4e-4, 2e-4, 1e-4 and 5e-5 s,
thus first order in `dt`.

The gas chemistry therefore uses a symmetric (Strang) split:

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

The gate tests each half with its own `dt/2`. The `BINNING` path keeps the
Lie split. */

#if !defined(BINNING) && !TURN_OFF_GAS_REACTIONS
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
  rate = gas_log_expansion (point, ystart, yend);
  drhodt[] -= (1. - f[])*cm[]*rate/dt;
#  endif // !NO_EXPANSION
}
# endif // VARPROP
#endif // !BINNING && !TURN_OFF_GAS_REACTIONS

#if !defined(BINNING) && !TURN_OFF_GAS_REACTIONS
/**
## The sweep of the gas-phase reactions

This function holds the loop of the gas-phase reactions of the external gas.
`dtc` is the time over which the reactor integrates, `dt/2` for each half
of the Strang split. The expansion source
always divides by the step `dt`, so each half gives its part of the mean of
the step. `second_half` is true only for the second half of the Strang split.
That half writes its expansion directly into `drhodt`, because
`update_divergence()` has already run. See "The Strang split of the gas
chemistry" above. */

static void gas_phase_reactions (double dtc, bool second_half)
{
  foreach (reduction(+:frozen_cell_gate_n)) {
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
      if (second_half)
        strang_second_half_expansion (point, ystart, y0ode);
      else
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
  gas_phase_reactions (0.5*dt, false);
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

#if !defined(BINNING) && !TURN_OFF_GAS_REACTIONS
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
#endif // !BINNING && !TURN_OFF_GAS_REACTIONS

#endif // TURN_OFF_REACTIONS
