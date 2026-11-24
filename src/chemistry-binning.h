/**
# The sweep of the gas-phase reactions with binning

`chemistry.h` includes this file when the case includes `binning.h` and
`TURN_OFF_GAS_REACTIONS` is 0. Do not include it directly. It uses the
helpers of `chemistry.h`: `gas_reactor_cell()`, `add_gas_expansion()`,
`gas_ideal_density()`, `gas_log_density_ratio()` and the field `binid`. */

/**
## The sweep of the gas-phase reactions with binning

This function replaces the per-cell sweep of `chemistry-cell.h` when the
case includes `binning.h`. It has the same arguments, the same cells
(`gas_reactor_cell()`) and the same expansion source
(`add_gas_expansion()`), so the Strang split and `TURN_OFF_GAS_REACTIONS`
work in the same way.

The reaction step follows section 2 of Goldin, Ren and Zahirovic, Combust.
Theory Model. 13 (2009) 721 (`papers/Goldin2009.md`):

1. `binning()` puts the cells with a similar state into one bin (section 3
   of the paper, see `binning.h`). The coordinates of the bins are
   `binning_M` targets: T, then the species of `binning_forced`, then the
   species with the largest gas mass. `binning_tol` is the tolerance,
   relative to the current range of each target.
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

`BINNING_GATE` 1 applies the frozen-cell gate of the per-cell sweep before
step 1, with the same test and the same tolerances. A frozen cell gets the
explicit update of the gate, and it does not go into a bin. The range of
each coordinate then covers the cells that react only. With
`BINNING_GATE` 0 (the default) every reactor cell goes into a bin, and
`frozen_cell_gate_n` stays 0. */

#if !defined(VARPROP)
# error "BINNING requires VARPROP (it uses rhoGv_G/cpGv_G and drhodt_chem)"
#endif

#ifndef BINNING_GATE
# define BINNING_GATE 0
#endif

/**
## The bin of each cell: the field `binid`

`binid` gives the bin of each gas cell after the last call of the gas
chemistry. In a step with the Strang split, these are the bins of the second
half. A cell that is not in a bin (no gas, or a frozen cell of
`BINNING_GATE`) holds `nodata`. The field is a global field, thus `dump()`
writes it, and jview can show it.

The numbers are local to each MPI rank, because each rank makes its own
bins. Two ranks can give the same number to two different bins.

The default refine and restriction functions average the values. A mean of
two bin numbers is not a bin number. Thus a new cell of `adapt` gets the
number of its parent, and a coarse cell gets the number of one of its
children. The next call of the chemistry writes the field again in all
cells. */

scalar binid[];

static void restriction_binid (Point point, scalar s)
{
  double v = nodata;
  foreach_child()
    if (s[] != nodata)
      v = s[];
  s[] = v;
}

event init (i = 0) {
  binid.refine = binid.prolongation = refine_injection;
  binid.restriction = restriction_binid;
  foreach()
    binid[] = nodata;
}

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
`BINNING_ALPHA_ROUNDOFF` below), `binning_ncells` and `binning_nbins` are
the number of cells of the mask (`N_c` of the paper) and of bins (`N_CA`),
`binning_nbadbin`
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
long binning_ncells = 0, binning_nbins = 0;

/**
## The parameters of the binning

Set them in an `init` event of the case (or `defaults`):

* `binning_tol`: the tolerance of the paper. Each coordinate is cut into
  `1/binning_tol` equal intervals of its current range (Eq. 7). The paper
  uses 0.01, 0.05 and 0.1.
* `binning_M`: the number of coordinates. The paper uses 4 and 11. A value
  above `NGS + 1` counts as `NGS + 1`.
* `binning_temperature`: T is the first coordinate (the paper does this by
  default; exclude it only for an isothermal system).
* `binning_forced`: the names of the minor species that the case forces into
  the coordinates, before the most abundant species. A NULL entry ends the
  list. The paper forces NO and CO in its cases A to C and no species in D
  and E. The default is no forced species. A forced species counts only when
  its gas mass in the domain is above 0 at that step.

The pressure is not a coordinate (the paper excludes it by default). */

double binning_tol = 0.05;
int binning_M = 4;
bool binning_temperature = true;
#define BINNING_MAXFORCED 16
const char * binning_forced[BINNING_MAXFORCED + 1] = {NULL};

/**
## The coordinates of the bins (section 3 of the paper)

At each call: T, then the forced species, then the species with the largest
total gas mass `sum Y w dV` over the cells of the mask, until there are
`binning_M` coordinates. `w` is the gas mass of the cell per unit volume.
With `CHEMISTRY_LOG`, rank 0 writes the coordinates to `bin-coords` when
they change. */

static scalar * binning_select_targets (scalar mask, scalar w)
{
  double gasmass[NGS];
  for (int jj = 0; jj < NGS; jj++)
    gasmass[jj] = 0.;
  foreach (serial)
    if (mask[]) {
      double m = w[]*dv();
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        gasmass[jj] += YG[]*m;
      }
    }
  mpi_all_reduce_array (gasmass, MPI_DOUBLE, MPI_SUM, NGS);

  int M = min (binning_M, NGS + (binning_temperature ? 1 : 0));
  bool used[NGS];
  for (int jj = 0; jj < NGS; jj++)
    used[jj] = false;
  int chosen[NGS + 1], n = 0;
  NOT_UNUSED (chosen);

  scalar * list = NULL;
  if (binning_temperature && n < M) {
    list = list_append (list, TG);
    chosen[n++] = -1;
  }

  for (int k = 0; k < BINNING_MAXFORCED && binning_forced[k] && n < M; k++) {
    int jj = OpenSMOKE_IndexOfSpecies (binning_forced[k]);
    if (jj < 0 || jj >= NGS) {
      static bool warned = false;
      if (!warned && pid() == 0)
        fprintf (stderr, "WARNING: binning_forced: no gas species %s\n",
                 binning_forced[k]);
      warned = true;
      continue;
    }
    if (!used[jj] && gasmass[jj] > 0.) {
      list = list_append (list, YGList_G[jj]);
      used[jj] = true;
      chosen[n++] = jj;
    }
  }

  while (n < M) {
    int best = -1;
    for (int jj = 0; jj < NGS; jj++)
      if (!used[jj] && gasmass[jj] > 0. &&
          (best < 0 || gasmass[jj] > gasmass[best]))
        best = jj;
    if (best < 0)
      break;
    list = list_append (list, YGList_G[best]);
    used[best] = true;
    chosen[n++] = best;
  }

# ifdef CHEMISTRY_LOG
  static int * last = NULL, nlast = -1;
  if (!last)
    last = malloc ((NGS + 1)*sizeof (int));
  bool changed = (n != nlast);
  for (int k = 0; k < n && !changed; k++)
    if (chosen[k] != last[k])
      changed = true;
  if (changed) {
    if (pid() == 0) {
      static FILE * fc = NULL;
      if (!fc)
        fc = fopen ("bin-coords", "w");
      fprintf (fc, "%g", t);
      for (int k = 0; k < n; k++)
        fprintf (fc, " %s", chosen[k] < 0 ? "T" :
                 OpenSMOKE_NamesOfSpecies (chosen[k]));
      fprintf (fc, "\n");
      fflush (fc);
    }
    for (int k = 0; k < n; k++)
      last[k] = chosen[k];
    nlast = n;
  }
# endif

  return list;
}

static void gas_phase_reactions (double dtc, bool second_half)
{
  /**
  Flag the cells where the reactor runs. Change their gas fields from the
  tracer form to the actual values. Keep the start density of each cell for
  the expansion source, the gas mass `wG` for the weight of the means, and
  the total enthalpy `hG`. */

  scalar gasmask[], rho0[], wG[], hG[], frozen[];
  foreach (reduction(+:frozen_cell_gate_n)) {
    gasmask[] = gas_reactor_cell (point) ? 1. : 0.;
    rho0[] = 0., wG[] = 0., hG[] = 0., frozen[] = 0.;
    if (gasmask[]) {
      scale_gas_tracers (point, 1./(1. - f[]));
      rho0[] = gas_cell_density (point);
      wG[] = rho0[]*(1. - f[]);

#if BINNING_GATE
      /**
      The gate of the per-cell sweep (see `FROZEN_CELL_GATE`). A frozen cell
      leaves the mask of the bins, and `frozen` keeps it in the loop of the
      expansion source at the end. The fields hold the actual values
      here. */

      if (dtc > 0.) {
        double yg[NGS + 1], dy_gate[NGS + 1];
        for (int jj = 0; jj < NGS; jj++) {
          scalar YG = YGList_G[jj];
          yg[jj] = YG[];
        }
        yg[NGS] = TG[];

        UserDataODE data;
        data.P = Pref + p[];
        data.T = TG[];
        data.sources = NULL;
        data.rhog = rhoGv_G[];
        data.cpg = cpGv_G[];
        gas_batch_nonisothermal_constantpressure (yg, dtc, dy_gate, &data);

        double dYmax = 0.;
        for (int jj = 0; jj < NGS; jj++)
          dYmax = fmax (dYmax, fabs (dy_gate[jj]));

        if (dYmax*dtc < FROZEN_CELL_YTOL &&
            fabs (dy_gate[NGS])*dtc < FROZEN_CELL_TTOL) {
          for (int jj = 0; jj < NGS; jj++) {
            scalar YG = YGList_G[jj];
            double y = yg[jj] + dtc*dy_gate[jj];
            YG[] = y > 0. ? y : 0.;
          }
          TG[] = yg[NGS] + dtc*dy_gate[NGS];
          gasmask[] = 0.;
          frozen[] = 1.;
          frozen_cell_gate_n++;
          continue;
        }
      }
#endif // BINNING_GATE

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
  total enthalpy, `[NGS+1]` the temperature, `[NGS+2]` the pressure `p`.
  The mean temperature is only the first guess of the Newton solve of the
  bin. The sweep then sets `bin->phi0[NGS+1]` to the temperature of the
  bin. */

  const int ih = NGS, iT = NGS + 1, ip = NGS + 2;
  scalar * fields = list_concat (YGList_G, {hG, TG, p});
  scalar * targets = binning_select_targets (gasmask, wG);

  BinTable * table = binning (fields, targets, binning_tol, wG, gasmask);
  scalar cellbin = table->id;

  /**
  `alpha` of each bin. A bin with `alpha < 0` is a bin that the sweep does
  not change. */

  int nb = table->nbins;
  double * alpha = malloc ((nb > 0 ? nb : 1)*sizeof (double));

# ifdef BINNING_MPI_DIAG
  double * cost = calloc (nb > 0 ? nb : 1, sizeof (double));
# endif

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
    data.P = Pref + bin->phi[ip];
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

# ifdef BINNING_MPI_DIAG
    struct timespec c0, c1;
    clock_gettime (CLOCK_MONOTONIC, &c0);
# endif
    OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
        NGS + 1, dtc, y0ode, &data);
# ifdef BINNING_MPI_DIAG
    clock_gettime (CLOCK_MONOTONIC, &c1);
    cost[i] = (c1.tv_sec - c0.tv_sec) + (c1.tv_nsec - c0.tv_nsec)*1e-9;
# endif

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
  cell, with a clamp to [0, 1]. See `BINNING_YTOL` for the margin. Each
  cell lowers the `alpha` of its bin to its own limit, so the result is the
  minimum over the cells and does not depend on their order.

  A cell can arrive with a mass fraction out of [0, 1] by more than the
  margin (the transport can make it slightly negative). Its bounds are then
  `min (Y, 0) - BINNING_YTOL` and `max (Y, 1) + BINNING_YTOL`: the reaction
  must not move it further out, but it does not have to bring it back. The
  clamp to [0, 1] then removes the excess, as the per-cell sweep does. Without
  this, the condition is false from the start and gave `alpha = 1.5e16` in
  the air bin of `fatehi-combustion`. */

# ifdef CHEMISTRY_LOG
  /**
  `alim[b]` is the species that sets `alpha` of the bin `b` (-1: none). */

  int * alim = malloc ((nb > 0 ? nb : 1)*sizeof (int));
  for (int b = 0; b < nb; b++)
    alim[b] = -1;
# endif

  foreach (serial) {
    int b = cellbin[];
    if (b >= 0 && alpha[b] > 0.) {
      Bin * bin = table->bins + b;
      double a = alpha[b];
      int jlim = -1;
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        double dY = bin->phi[jj] - bin->phi0[jj];
        double lower = min (YG[], 0.) - BINNING_YTOL;
        double upper = max (YG[], 1.) + BINNING_YTOL;
        if (dY < 0. && YG[] + a*dY < lower)
          a = min (a, (YG[] - lower)/(- dY)), jlim = jj;
        else if (dY > 0. && YG[] + a*dY > upper)
          a = min (a, (upper - YG[])/dY), jlim = jj;
      }
# ifdef CHEMISTRY_LOG
      if (jlim >= 0)
        alim[b] = jlim;
# endif
      alpha[b] = max (a, 0.);
    }
  }

  long nbadcell = 0;
  foreach (reduction(+:nbadcell)) {
    int b = cellbin[];
    if (b >= 0 && alpha[b] >= 0.) {
      Bin * bin = table->bins + b;
      double a = alpha[b];
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

  double amin = 1.;
  long nlimited = 0, nbadbin = 0;
  for (int b = 0; b < nb; b++) {
    if (alpha[b] < 0.)
      nbadbin++;
    else {
      amin = min (amin, alpha[b]);
      if (alpha[b] < 1. - BINNING_ALPHA_ROUNDOFF)
        nlimited++;
    }
  }

# ifdef CHEMISTRY_LOG
  /**
  Rank 0 writes one line to `bin-limit` for each call, with the values of
  all ranks:

      t half nbins nlimited n1e-6 n1e-3 n0.1 n0.5 mlim mlost name:count ...

  * `n<x>`: the bins with `alpha` below `x`;
  * `mlim`: the gas mass of the limited bins over the gas mass of all bins;
  * `mlost`: the sum of `mass*(1 - alpha)` over the gas mass of all bins.
    It is the part of the reaction that the remap does not give to the
    cells;
  * `name:count`: the 5 species that set `alpha` in the most bins.

  The bins that the sweep did not change (`alpha < 0`) do not count. */

  const double athr[4] = {1e-6, 1e-3, 0.1, 0.5};
  long lcount[NGS + 5];
  for (int k = 0; k < NGS + 5; k++)
    lcount[k] = 0;
  double lmass[3] = {0., 0., 0.}; // all bins, limited bins, mass*(1 - alpha)
  for (int b = 0; b < nb; b++) {
    if (alpha[b] < 0.)
      continue;
    double m = table->bins[b].mass;
    lmass[0] += m;
    lmass[2] += m*(1. - alpha[b]);
    if (alpha[b] < 1. - BINNING_ALPHA_ROUNDOFF) {
      lmass[1] += m;
      if (alim[b] >= 0)
        lcount[alim[b]]++;
      for (int k = 0; k < 4; k++)
        if (alpha[b] < athr[k])
          lcount[NGS + k]++;
    }
  }
  lcount[NGS + 4] = nb - nbadbin;
  long nlim_all = nlimited;
  mpi_all_reduce_array (lcount, MPI_LONG, MPI_SUM, NGS + 5);
  mpi_all_reduce_array (lmass, MPI_DOUBLE, MPI_SUM, 3);
  mpi_all_reduce (nlim_all, MPI_LONG, MPI_SUM);
  free (alim);

  if (pid() == 0) {
    static FILE * flim = NULL;
    if (!flim)
      flim = fopen ("bin-limit", "w");
    double mall = lmass[0] > 0. ? lmass[0] : 1.;
    fprintf (flim, "%g %d %ld %ld %ld %ld %ld %ld %.4g %.4g", t,
             second_half ? 1 : 0, lcount[NGS + 4], nlim_all,
             lcount[NGS], lcount[NGS + 1], lcount[NGS + 2], lcount[NGS + 3],
             lmass[1]/mall, lmass[2]/mall);
    bool used[NGS];
    for (int jj = 0; jj < NGS; jj++)
      used[jj] = false;
    for (int r = 0; r < 5; r++) {
      int best = -1;
      for (int jj = 0; jj < NGS; jj++)
        if (!used[jj] && lcount[jj] > 0 &&
            (best < 0 || lcount[jj] > lcount[best]))
          best = jj;
      if (best < 0)
        break;
      used[best] = true;
      fprintf (flim, " %s:%ld", OpenSMOKE_NamesOfSpecies (best), lcount[best]);
    }
    fprintf (flim, "\n");
    fflush (flim);
  }
# endif

  free (alpha);

# ifdef CHEMISTRY_LOG
  /**
  The file `bin-<pid>` gets one line for each call: the time, the number of
  cells of the mask (`N_c` of the paper), the number of bins (`N_CA`), 1 for
  the second half of the Strang split (else 0), the smallest `alpha`, the
  number of bins with `alpha < 1`, the number of bins that the sweep did not
  change, and the number of cells where the Newton solve did not converge.
  All values are those of this rank, apart from the last one (all ranks). */

  static FILE * fp = NULL;
  if (!fp) {
    char name[20];
    sprintf (name, "bin-%d", pid());
    fp = fopen (name, "w");
  }
  fprintf (fp, "%g %ld %d %d %g %ld %ld %ld\n", t, table->nmask, nb,
           second_half ? 1 : 0, amin, nlimited, nbadbin, nbadcell);
  fflush (fp);
# endif

# ifdef BINNING_MPI_DIAG
  /**
  With `BINNING_MPI_DIAG` (and MPI), rank 0 writes the comparison of local
  and global bins to `binning-mpi` (see `binning_mpi_report()` in
  `binning.h`). */

@if _MPI
  static FILE * fmpi = NULL;
  if (!fmpi && pid() == 0)
    fmpi = fopen ("binning-mpi", "w");
  binning_mpi_report (table, cost, t, second_half ? 1 : 0, fmpi);
@endif
  free (cost);
# endif

  long ncells = table->nmask, nbins = nb;
  mpi_all_reduce (amin, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (nlimited, MPI_LONG, MPI_SUM);
  mpi_all_reduce (nbadbin, MPI_LONG, MPI_SUM);
  mpi_all_reduce (ncells, MPI_LONG, MPI_SUM);
  mpi_all_reduce (nbins, MPI_LONG, MPI_SUM);
  binning_ncells = ncells;
  binning_nbins = nbins;
  binning_alpha_min = amin;
  binning_nlimited = nlimited;
  binning_nbadbin = nbadbin;
  binning_nbadcell = nbadcell;

# if defined(CHEMISTRY_LOG) && BINNING_GATE
  if ((second_half || !GAS_CHEMISTRY_STRANG) && pid() == 0)
    fprintf (stderr, "GATE %g %d\n", t, frozen_cell_gate_n);
# endif

  binning_ids (table, binid);
  binning_cleanup (table);
  free (fields), fields = NULL;
  free (targets), targets = NULL;

  /**
  Add the expansion of each cell from its two densities. Then change the gas
  fields back to the tracer form. The frozen cells of `BINNING_GATE` are not
  in `gasmask`, but they need both steps. */

  foreach() {
    if (gasmask[] || frozen[]) {
      double lnratio = gas_log_density_ratio (rho0[], gas_cell_density (point));
      add_gas_expansion (point, lnratio, second_half);

      scale_gas_tracers (point, 1. - f[]);
    }
  }
}
