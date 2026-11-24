/**
# Does the DRG reduction keep the cells of a real flame?

`drg-cell.c` tests an ignition. This case tests the cells that the gas
reactor gets in a flow. It reads the gas cells of one call of the gas
reactor, from a dump of `CHEMISTRY_DUMP` (`chemistry-cell.h`), and makes one
step of the reactor for each cell:

* `full` calls `OpenSMOKE_ODESolver` with the full mechanism.
* `drg` calls `gas_stiff_solve()` of `drg.h`, for each threshold `eps`.

The default dump is the second Strang half of step 134065 of the Fatehi case
(`run/fatehi-combustion.c`, `biomass/Solid-gas-88`, steady flame at
t = 15 s, 8071 gas cells), in `~/cuSMOKE/bench/results/fatehi_20261113`.
Set these environment variables to change the input:

* `DRG_DUMP`: the merged `.cells` file. The `.meta` file must be next to it.
* `DRG_DTC`: the interval of the reactor [s]. The default is the half step
  of that call, 2.769e-5 s.
* `DRG_STRIDE`: use one cell out of `DRG_STRIDE`. The default is 1.
* `DRG_KIN`: the kinetic scheme of the test. The default is
  `biomass/Solid-gas-88`. It can have more gas species than the dump: the
  species that are not in the dump start at zero.
* `DRG_KEYS`: the key species, separated by commas. The default is the list
  below, for the Fatehi flame.
* `DRG_WARMUP`: the number of full steps of `dtc` that each cell makes before
  the comparison. The comparison starts from the state after them. Use it
  when the mechanism of the test differs from the mechanism of the run: the
  species that are not in the dump start at zero, and the warm-up forms them.
  The default is 0.
* `DRG_REPEAT`: the number of times the case repeats each solve of the
  cells. The CPU time is the shortest of the repetitions. The default is 1.

The case skips the cells that the frozen-cell gate skipped in the run.

The error of a cell is the largest error of its components (Y and T) in
units of the tolerance of the run's solver, `|y_drg - y_full|/(rtol |y_full|
+ atol)`, with rtol = 1.19e-5 and atol = 1e-10. A value near 1 or less is at
the level of the error of the solver itself.

The error of the main components, `errM`, takes only T and the species
with `Y > 1e-4` in the full result. A species that the full mechanism forms
at a small Y and the DRG drops gives a large error, but no `errM`.

This case does not stop with an error. It only prints the results.
*/

#define F_ERR 1e-10
#define SOLVE_TEMPERATURE 1
#define GAS_DRG 1

scalar f[];

#include "run.h"

attribute {
  scalar * tracers, c;
  bool inverse;
}

scalar p[];
scalar porosity[];
scalar zeta[];

double rhoG = 0.3;
double rhoS = 1550.;

#include "memoryallocation-varprop.h"
#include "chemistry.h"
#include <time.h>
#include <stdint.h>

/**
The default key species: O2 and the main fuels of the Fatehi flame gas, by
their largest mass fraction in the dump. `DRG_KEYS` replaces them. */

#define MAXKEY 64
static const char * keys[MAXKEY] = {"O2", "CO", "H2", "CH4", "C6H10O5",
                                    "FURFURAL", "C5H8O4", "MLINO", "C6H6O3",
                                    "CH2OHCHO", "C2H2", "C2H4"};
static int NKEY = 12;

static void read_keys (void) {
  const char * env = getenv ("DRG_KEYS");
  if (!env)
    return;
  char * list = strdup (env);
  NKEY = 0;
  for (char * k = strtok (list, ","); k && NKEY < MAXKEY; k = strtok (NULL, ","))
    keys[NKEY++] = k;
}

#define RTOL 1.19e-5
#define ATOL 1e-10

static void gas_props (const double * y, UserDataODE * data) {
  double x[NGS], mf[NGS], MW;
  for (int jj = 0; jj < NGS; jj++)
    mf[jj] = fmax (0., y[jj]);
  OpenSMOKE_GasProp_SetTemperature (clamp_temperature (y[NGS]));
  OpenSMOKE_GasProp_SetPressure (clamp_pressure (data->P));
  OpenSMOKE_MoleFractions_From_MassFractions (x, &MW, mf);
  data->rhog = data->P/(R_GAS*1000.*y[NGS])*MW;
  data->cpg = OpenSMOKE_GasProp_HeatCapacity (x);
}

static double cpu_time (void) {
  struct timespec ts;
  clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &ts);
  return ts.tv_sec + 1e-9*ts.tv_nsec;
}

/**
Stop the run when a read fails. Do not put a read into `assert()`: with
`NDEBUG`, the read does not run. */

static void must (bool ok, const char * what) {
  if (!ok) {
    fprintf (stderr, "drg-dump: %s\n", what);
    exit (1);
  }
}

static int cmp_double (const void * a, const void * b) {
  double x = *(const double *) a, y = *(const double *) b;
  return (x > y) - (x < y);
}

/**
The percentile `q` (0 to 1) of `n` values. The function sorts `v`. */

static double percentile (double * v, int n, double q) {
  if (n == 0)
    return 0.;
  qsort (v, n, sizeof (double), cmp_double);
  int i = (int) (q*(n - 1) + 0.5);
  return v[i];
}

/**
One step of the reactor for one cell, with the full mechanism or with
`gas_stiff_solve()`. */

static void step (double * y, double P, double dtc, bool reduced) {
  UserDataODE data;
  data.P = P;
  data.zeta = 0.;
  data.sources = NULL;
  gas_props (y, &data);
  data.T = y[NGS];
  if (reduced)
    gas_stiff_solve (y, dtc, &data);
  else
    OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
                         NGS + 1, dtc, y, &data);
}

int main() {
  TS0 = 800.; TG0 = 1123.;
  cpS = 1500.; cpG = 1100.;
  kinfolder = getenv ("DRG_KIN") ? getenv ("DRG_KIN") : "biomass/Solid-gas-88";
  read_keys ();
  DT = 1.;
  init_grid (1);
  run();
}

event init (i = 0) {
  drg_key_species = keys, drg_nkey = NKEY;

  restarted = true;
  foreach() {
    f[] = 1.; porosity[] = 0.2; zeta[] = 0.; p[] = 0.;
    TS[] = TS0; TG[] = 0.; T[] = TS0;
  }
}

event check (i = 0) {

  char dumpfile[1024];
  const char * env = getenv ("DRG_DUMP");
  if (env)
    snprintf (dumpfile, sizeof (dumpfile), "%s", env);
  else
    snprintf (dumpfile, sizeof (dumpfile), "%s/cuSMOKE/bench/results/"
              "fatehi_20261113/merged/chemdump-134065-2.cells",
              getenv ("HOME"));
  double dtc = getenv ("DRG_DTC") ? atof (getenv ("DRG_DTC")) : 2.769e-5;
  int stride = getenv ("DRG_STRIDE") ? atoi (getenv ("DRG_STRIDE")) : 1;
  int repeat = getenv ("DRG_REPEAT") ? atoi (getenv ("DRG_REPEAT")) : 1;
  int warmup = getenv ("DRG_WARMUP") ? atoi (getenv ("DRG_WARMUP")) : 0;
  if (repeat < 1)
    repeat = 1;

  /**
  Read the header and match the names of the species with the mechanism. */

  FILE * fp = fopen (dumpfile, "rb");
  if (!fp) {
    fprintf (stderr, "drg-dump: cannot open %s\n", dumpfile);
    exit (1);
  }
  int64_t ncells;
  int32_t nsp;
  must (fread (&ncells, sizeof (int64_t), 1, fp) == 1, "cannot read ncells");
  must (fread (&nsp, sizeof (int32_t), 1, fp) == 1, "cannot read nsp");
  must (nsp <= NGS, "the dump has more species than the mechanism");
  int map[nsp];
  for (int jj = 0; jj < nsp; jj++) {
    char name[33] = {0};
    must (fread (name, 1, 32, fp) == 32, "cannot read a species name");
    map[jj] = OpenSMOKE_IndexOfSpecies (name);
    must (map[jj] >= 0 && map[jj] < NGS, "a species of the dump is not in the mechanism");
  }

  double * cells = malloc (ncells*(4 + nsp)*sizeof(double));
  must (fread (cells, sizeof(double), ncells*(4 + nsp), fp) ==
        ncells*(4 + nsp), "cannot read the cells");
  fclose (fp);

  char metafile[1024];
  snprintf (metafile, sizeof (metafile), "%.*s.meta",
            (int) (strlen (dumpfile) - 6), dumpfile);
  fp = fopen (metafile, "rb");
  if (!fp) {
    fprintf (stderr, "drg-dump: cannot open %s\n", metafile);
    exit (1);
  }
  double * meta = malloc (ncells*4*sizeof(double));
  must (fread (meta, sizeof(double), ncells*4, fp) == ncells*4,
        "cannot read the meta file");
  fclose (fp);

  /**
  The cells of the test: all the cells that the run solved, one out of
  `stride`. */

  int n = 0;
  int * list = malloc (ncells*sizeof(int));
  for (int64_t c = 0; c < ncells; c += stride)
    if (meta[4*c + 1] == 0.)
      list[n++] = c;

  fprintf (stderr, "# key species:");
  for (int k = 0; k < NKEY; k++)
    fprintf (stderr, " %s", keys[k]);
  fprintf (stderr, "\n");
  fprintf (stderr, "# %s\n# %ld cells in the dump, %d in the test, "
           "dtc = %g s, %d gas species (%d in the dump), %d key species, "
           "%d repetitions, %d warm-up steps\n",
           dumpfile, (long) ncells, n, dtc, NGS, (int) nsp, NKEY, repeat,
           warmup);

  /**
  The full solve of each cell. */

  double * y0 = malloc (n*(NGS + 1)*sizeof(double));
  double * yfull = malloc (n*(NGS + 1)*sizeof(double));
  double * P = malloc (n*sizeof(double));
  bool * hot = malloc (n*sizeof(bool));
  int nhot = 0;
  double dTstep = 0.;
  for (int k = 0; k < n; k++) {
    double * row = cells + list[k]*(4 + nsp);
    double * y = y0 + k*(NGS + 1);
    for (int jj = 0; jj < NGS; jj++)
      y[jj] = 0.;
    for (int jj = 0; jj < nsp; jj++)
      y[map[jj]] = row[4 + jj];
    y[NGS] = row[3];
    P[k] = meta[4*list[k]];
    hot[k] = y[NGS] > 1500.;
    nhot += hot[k];
  }

  /**
  The warm-up: full steps that change the start state of each cell. */

  for (int k = 0; k < n; k++)
    for (int w = 0; w < warmup; w++)
      step (y0 + k*(NGS + 1), P[k], dtc, false);
  if (warmup > 0) {
    nhot = 0;
    for (int k = 0; k < n; k++) {
      hot[k] = y0[k*(NGS + 1) + NGS] > 1500.;
      nhot += hot[k];
    }
  }

  double cpu_full = HUGE, c0;
  for (int r = 0; r < repeat; r++) {
    c0 = cpu_time ();
    for (int k = 0; k < n; k++) {
      double * y = yfull + k*(NGS + 1);
      for (int jj = 0; jj < NGS + 1; jj++)
        y[jj] = y0[k*(NGS + 1) + jj];
      step (y, P[k], dtc, false);
      dTstep = fmax (dTstep, fabs (y[NGS] - y0[k*(NGS + 1) + NGS]));
    }
    cpu_full = fmin (cpu_full, cpu_time () - c0);
  }

  fprintf (stderr, "# full: cpu %.3g s (%.3g ms per cell), largest |T_end - T_start| %.3g K, "
           "%d cells above 1500 K\n", cpu_full, 1e3*cpu_full/n, dTstep, nhot);
  fprintf (stderr, "# eps  ns_mean  ns_max  nr_mean  "
           "err_median  err_p99  err_max  err_max_hot  "
           "dT_max[K]  cells_err>10  cpu[s]  speedup  "
           "errM_median  errM_p99  errM_max\n");

  /**
  The DRG solve of each cell, for each threshold. */

  double epss[] = {0., 1e-3, 1e-2, 5e-2, 1e-1};
  double * err = malloc (n*sizeof(double));
  double * errM = malloc (n*sizeof(double));
  for (int e = 0; e < 5; e++) {
    double Tmax = 1e5;
    OpenSMOKE_DRG_Init (keys, NKEY, &epss[e], &Tmax, 1,
                        drg_T_threshold, drg_R_threshold);
    drg_nsolve = 0, drg_nspecies = 0., drg_nreactions = 0.;

    /**
    The repetitions before the last one measure only the CPU time. The last
    one also gives the errors and the statistics. */

    double cpu = HUGE;
    for (int r = 0; r < repeat - 1; r++) {
      c0 = cpu_time ();
      for (int k = 0; k < n; k++) {
        double y[NGS + 1];
        for (int jj = 0; jj < NGS + 1; jj++)
          y[jj] = y0[k*(NGS + 1) + jj];
        step (y, P[k], dtc, true);
      }
      cpu = fmin (cpu, cpu_time () - c0);
    }
    drg_nsolve = 0, drg_nspecies = 0., drg_nreactions = 0.;

    int nsmax = 0, nbad = 0;
    double dTmax = 0., errhot = 0., cpu_last = 0.;
    for (int k = 0; k < n; k++) {
      double y[NGS + 1];
      for (int jj = 0; jj < NGS + 1; jj++)
        y[jj] = y0[k*(NGS + 1) + jj];

      double ns0 = drg_nspecies;
      c0 = cpu_time ();
      step (y, P[k], dtc, true);
      cpu_last += cpu_time () - c0;
      nsmax = max (nsmax, (int) (drg_nspecies - ns0));

      const double * yf = yfull + k*(NGS + 1);
      double ek = 0., ekM = 0.;
      for (int jj = 0; jj < NGS + 1; jj++) {
        double e = fabs (y[jj] - yf[jj])/(RTOL*fabs (yf[jj]) + ATOL);
        ek = fmax (ek, e);
        if (jj == NGS || fabs (yf[jj]) > 1e-4)
          ekM = fmax (ekM, e);
      }
      err[k] = ek;
      errM[k] = ekM;
      dTmax = fmax (dTmax, fabs (y[NGS] - yf[NGS]));
      if (hot[k])
        errhot = fmax (errhot, ek);
      nbad += ek > 10.;
    }

    cpu = fmin (cpu, cpu_last);
    double emax = percentile (err, n, 1.);
    double eM50 = percentile (errM, n, 0.5), eM99 = percentile (errM, n, 0.99),
      eMmax = percentile (errM, n, 1.);
    fprintf (stderr, "%g %.1f %d %.1f %.3g %.3g %.3g %.3g %.3g %d %.3g %.2f"
             " %.3g %.3g %.3g\n",
             epss[e], drg_nspecies/drg_nsolve, nsmax,
             drg_nreactions/drg_nsolve,
             percentile (err, n, 0.5), percentile (err, n, 0.99), emax,
             errhot, dTmax, nbad, cpu, cpu_full/cpu, eM50, eM99, eMmax);
  }

  free (cells), free (meta), free (list), free (y0), free (yfull);
  free (P), free (hot), free (err), free (errM);
}

event stop (i = 1);
