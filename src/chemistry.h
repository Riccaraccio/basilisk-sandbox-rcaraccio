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
  foreach()
    omega[] = 0.;
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
an exact `rhoGv_G*dY/dt` and `rhoGv_G*cpGv_G*dT/dt`. This is what `divu2` in
`multicomponent-properties.h` needs: it divides by the same two fields, so the
result is the mole-number change rate and the thermal expansion rate.

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
*/

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

/**
Convenience wrapper for the per-cell path, which holds both states and where
`rhoGv_G`/`cpGv_G` do not change over the chemistry event. */

static void accumulate_gas_sources (Point point, const double * ystart,
                                    const double * yend) {
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
}
#endif

event chemistry (i++) {

#ifdef CHEMISTRY_LOG
  reset ({t_solid, t_gas}, 0.);
  double time_mpi[npe()];

  for (int pe = 0; pe < npe(); pe++)
    time_mpi[pe] = 0.;

  struct timespec start, end;
  clock_gettime (CLOCK_MONOTONIC, &start);
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

  foreach() {
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
        Using an explicit solver for gas-phase reactions is not
        recommended as they are usually stiff.
        */
      OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure, NGS + 1, dt, y0ode, &data);

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
      accumulate_gas_sources (point, ystart, y0ode);
# endif

      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        YG[] = (y0ode[jj] > 0.) ? y0ode[jj]*(1. - f[]) : 0.;
      }
      TG[] = y0ode[NGS]*(1. - f[]);
    }
  }
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
#endif // TURN_OFF_REACTIONS
