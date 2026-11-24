/**
# The stiff solve of the gas reactor, with the DRG reduction

`chemistry.h` includes this file only with `GAS_DRG` 1, before the two
sweeps of the gas-phase reactions. Do not include it directly. Both sweeps
call `gas_stiff_solve()` for each stiff solve of the gas reactor: one cell
in `chemistry-cell.h`, one bin in `chemistry-binning.h`. Without `GAS_DRG`,
`chemistry.h` gives a `gas_stiff_solve()` that integrates the full
mechanism.

Here, the Directed Relation Graph (DRG) of the interface
(`OpenSMOKE_DRG_*`) reduces the mechanism for the state at the start of
each solve:

1. The analysis finds the important species. A reaction is important only
   when all its species are important.
2. The stiff solver integrates only the important species and the
   temperature. The Jacobian thus has `n + 1` columns, not `NGS + 1`. This
   gives most of the speed-up. The interface computes the kinetic constants
   of the important reactions only, but it computes the rates of all the
   reactions.
3. The unimportant species stay constant. In the reduced mechanism their
   formation rate is exactly zero, because no important reaction contains
   them. Their concentrations still go into the third-body terms and into
   the mixture properties.

Caution: a species with a zero mass fraction at the start of the solve is
not important. It stays at zero for the full solve, also when the full
mechanism forms it. With a long step through an ignition, this changes the
temperature by tens of kelvin, also with `eps = 0` (`test/drg-cell.c`). In
the cells of a flame, the step is short and the radical pool exists, and
the error stays at the tolerance of the solver (`test/drg-dump.c`).

The speed-up depends on the size of the mechanism. With 55 or 87 species,
the reduction gives no speed-up. With 402 species and `eps = 1e-3`, it gives
1.8 times (`test/drg-dump.c` on the cells of a flame).

The frozen-cell gate and the solid reactor call the right-hand side outside
of this solve. They always use the full mechanism (see `gas_drg_active` in
`reactors.h`). */

#if _OPENMP
# error "GAS_DRG uses one DRG object for each process: do not use it with OpenMP"
#endif

/**
## The settings

Set the settings in an `init` event of the case. The `init` event of the
case runs before the `init` event of this file, because Basilisk runs
same-name events in reverse order of declaration.

* `drg_key_species`, `drg_nkey`: the names of the key species, for example
  the fuel and O2. These species are always important. You must set them.
* `drg_eps`, `drg_T`, `drg_nT`: the thresholds of the DRG and the upper
  temperature bounds [K], in ascending order. The analysis uses
  `drg_eps[i]` for the first `i` with `T < drg_T[i]`. Above the last bound,
  it uses 1e-2. The default is one interval with `eps = 1e-2`.
* `drg_T_threshold` [K]: below this temperature, only the key species are
  important and no reaction is important. The default is 300 K.
* `drg_R_threshold` [kmol/m3/s]: a species with its production and
  destruction rates below this value does not add species to the graph. The
  default is 0.

Example:

~~~literatec
event init (i = 0) {
  static const char * keys[] = {"IC4H9OH", "O2"};
  drg_key_species = keys, drg_nkey = 2;
}
~~~
*/

const char ** drg_key_species = NULL;
int drg_nkey = 0;

static double drg_eps_default[] = {1e-2};
static double drg_T_default[] = {1e5};
double * drg_eps = drg_eps_default;
double * drg_T = drg_T_default;
int drg_nT = 1;

double drg_T_threshold = 300.;
double drg_R_threshold = 0.;

/**
## The statistics

`gas_stiff_solve()` adds to these counters on each rank. The mean number of
important species is `drg_nspecies/drg_nsolve`. The case can print and reset
them. */

long drg_nsolve = 0;
double drg_nspecies = 0., drg_nreactions = 0.;

event init (i = 0) {
  if (drg_key_species == NULL || drg_nkey < 1) {
    fprintf (stderr, "drg.h: set drg_key_species and drg_nkey "
             "in an init event of the case\n");
    exit (1);
  }
  OpenSMOKE_DRG_Init (drg_key_species, drg_nkey, drg_eps, drg_T, drg_nT,
                      drg_T_threshold, drg_R_threshold);
}

/**
## The reduced system

The reduced state is the mass fractions of the important species, then the
temperature. `drg_reduced_rhs()` puts the reduced state into a copy of the
full state, calls the gas right-hand side of `reactors.h`, and keeps the
derivatives of the reduced state.

The three variables below are valid only during a reduced solve.
`gas_stiff_solve()` points them to its local arrays. */

static int drg_n = 0;
static const int * drg_index = NULL;
static const double * drg_yfull = NULL;

static void drg_reduced_rhs (const double * yr, const double dt, double * dyr,
                             void * args)
{
  double y[NGS + 1], dy[NGS + 1];
  for (int jj = 0; jj < NGS + 1; jj++)
    y[jj] = drg_yfull[jj];
  for (int k = 0; k < drg_n; k++)
    y[drg_index[k]] = yr[k];
  y[NGS] = yr[drg_n];

  gas_batch_nonisothermal_constantpressure (y, dt, dy, args);

  for (int k = 0; k < drg_n; k++)
    dyr[k] = dy[drg_index[k]];
  dyr[drg_n] = dy[NGS];
}

/**
The analysis uses the concentrations of the full state, with the same
conversion and the same clamps as the gas right-hand side. */

static void drg_analysis (const double * y, double P)
{
  double T = clamp_temperature (y[NGS]);
  P = clamp_pressure (P);

  double Y[NGS], x[NGS], c[NGS], MW;
  for (int jj = 0; jj < NGS; jj++)
    Y[jj] = y[jj] > 0. ? y[jj] : 0.;
  OpenSMOKE_MoleFractions_From_MassFractions (x, &MW, Y);

  double ctot = P/(R_GAS*1000.*T); // kmol/m3
  for (int jj = 0; jj < NGS; jj++)
    c[jj] = ctot*x[jj];

  OpenSMOKE_DRG_Analysis (T, P, c);
}

/**
## The solve

`y` holds the `NGS` mass fractions and the temperature. On return, it holds
the state at the end of `dtc`. When the solve diverges, `y` can hold values
that are not finite, as with `OpenSMOKE_ODESolver()`. The sweeps check this
after the call. */

static void gas_stiff_solve (double * y, double dtc, UserDataODE * data)
{
  drg_analysis (y, data->P);

  /**
  Copy the indices: the next analysis can make the pointer of the
  interface invalid. */

  int n = OpenSMOKE_DRG_NumberImportantSpecies ();
  const int * important = OpenSMOKE_DRG_IndicesImportantSpecies ();
  int index[NGS];
  for (int k = 0; k < n; k++)
    index[k] = important[k];

  double yfull[NGS + 1], yr[NGS + 1];
  for (int jj = 0; jj < NGS + 1; jj++)
    yfull[jj] = y[jj];
  for (int k = 0; k < n; k++)
    yr[k] = y[index[k]];
  yr[n] = y[NGS];

  drg_nsolve++;
  drg_nspecies += n;
  drg_nreactions += OpenSMOKE_DRG_NumberImportantReactions ();

  drg_n = n, drg_index = index, drg_yfull = yfull;
  gas_drg_active = true;
  OpenSMOKE_ODESolver (&drg_reduced_rhs, n + 1, dtc, yr, data);
  gas_drg_active = false;
  drg_index = NULL, drg_yfull = NULL;

  for (int k = 0; k < n; k++)
    y[index[k]] = yr[k];
  y[NGS] = yr[n];
}
