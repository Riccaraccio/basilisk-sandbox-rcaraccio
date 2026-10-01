/**
# Binning against the per-cell solve through an ignition

This test checks the `BINNING` path of `chemistry.h` (Goldin, Ren and
Zahirovic 2009, `papers/Goldin2009.md`) where the binning is most difficult:
an ignition, where a bin can mix a cell that runs away with cells that are
still in their induction period.

The domain is a grid of independent batch reactors at constant pressure,
with no transport and with gas only (`f = 0`). Each cell holds a mixture of
the pyrolysis gas of `gate-ignition.c` and air. The fraction `s` of
pyrolysis gas grows along x, and the temperature grows along y from 900 K
to 1050 K. At the production step (2.4e-4 s), a cell at 950 K ignites after
about 13 ms, a cell at 1000 K after 3 ms, and a cell at 900 K after 0.11 s.
So over 80 steps the ignition passes through the grid, and some cells stay
in induction.

The reference is one stiff solve for each cell and each step: the calls of
the per-cell sweep of `chemistry.h`, without the gate. The test computes it
one time and keeps the whole trajectory. Each case of the binning then calls
`gas_phase_reactions()` of `chemistry.h` itself, in two ways:

* local: one step from the reference state, at each step. This is the
  error of the reaction step of the paper (Eq. 10),
  `E_i = |T_i(binning) - T_i(reference)|/T_i(reference)`;
* global: the binning marches from the initial state, as in a run. The
  error then includes the errors of all earlier steps, and the ignition
  time of a cell can move.

The cases follow the paper. A, B and C have `tol` = 0.01, 0.05 and 0.1
with M = 4, and E has `tol` = 0.05 with M = 11, with the coordinates of
section 3 (T, then the most abundant species). B is also the case D of the
paper, because no species is forced. F is B with the forced species H2O2
and OH, which mark the induction and the ignition (not in the paper). The
last case has one state in each bin (`tol` = 1e-9, all species and T as
coordinates).

What the cases show (2026-10-01):

* The error of an ignition is the spurious mixing of the bins, which the
  paper describes for the unburnt cells next to a flame. With M = 4, a bin
  holds cells with the same T and major species but different radicals.
  Two experiments, not kept in the code, gave the same local errors: the
  remap of Eq. 5 with no `alpha` (case C then ignites 249 cells, also cells
  that never ignite in the reference) and one `alpha` for each cell (210).
  With one `alpha` for each bin, as in the paper, case C ignites 103 of the
  184 cells. `alpha` thus moves the error, it does not make it.
* The case with one state in each bin does not give the reference to the
  rounding error: the relative difference is 5e-5 (local) and 3e-4 (global)
  in this ignition, the size of the tolerance of the stiff solver amplified
  by the ignition. The log therefore tests it against 1e-3 and 1e-2, and it
  tests that every cell ignites at the same step as in the reference.

The log gives no time, and the errors with two digits only, so that it can
be compared with `binning-ignition.ref`. A machine with another compiler or
another OpenSMOKE build can change the last digit; regenerate the reference
there only after the case with one state in each bin says `yes`. */

#define F_ERR 1e-10
#define SOLVE_TEMPERATURE 1

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
double rhoS = 1500.;

/**
The fields of `variable-properties.h` and `multicomponent-properties.h` that
`chemistry.h` reads under `VARPROP`. This test has no flow solver and no
transport, so it declares them here, as `gas-source-cell.c` does. */

#define VARPROP
scalar rhoGv_G[], rhoGv_S[], rhoSv[];
scalar muGv_G[], muGv_S[];
scalar lambdaGv_G[], lambdaGv_S[], lambdaSv[];
scalar cpGv_G[], cpGv_S[], cpSv[];
scalar drhodt[], drhodt_chem[];
scalar DTDtS[];

/**
Caution: `DYDtG_S` must hold `NGS` fields, also when no cell is solid.
Basilisk runs the body of each `foreach` one time in a stencil pass, through
all branches, so the solid branch of the `chemistry` event reads the list. A
`NULL` list makes the test crash there. */

scalar * DYDtG_S = NULL;

#include "binning.h"
#include "memoryallocation-varprop.h"

/**
All fields are 0, so the `chemistry` event of the first step finds no gas
cell, and only the event `ignition` below changes the state. */

event init (i = 0) {
  restarted = true;
  for (int jj = 0; jj < NGS; jj++) {
    scalar s = new scalar;
    DYDtG_S = list_append (DYDtG_S, s);
  }
  reset (DYDtG_S, 0.);
  foreach() {
    f[] = 0.; porosity[] = 0.; zeta[] = 0.; p[] = 0.;
    TS[] = 0.; TG[] = 0.; T[] = 0.;
    for (scalar s in YGList_G) s[] = 0.;
    for (scalar s in YGList_S) s[] = 0.;
    for (scalar s in YSList) s[] = 0.;
    rhoGv_G[] = 1., cpGv_G[] = 1000.;
  }
}

#include "chemistry.h"

#define NX 16
#define NSTEPS 80
#define DTB 2.4e-4
#define RISE 200.

int main() {
  TS0 = 800.; TG0 = 1123.;
  cpS = 1500.; cpG = 1100.;
  kinfolder = "biomass/Solid-gas-88";
  DT = DTB;
  init_grid (NX);
  run();
}


/**
## The state of the cells

A state is `NGS` mass fractions and the temperature. A trajectory keeps the
state of every cell at every step, in the order of `foreach (serial)`. */

#define NV (NGS + 1)

static double * state (double * S, long ncells, int k, long c) {
  return S + ((long) k*ncells + c)*NV;
}

static void initial_state (double * y, double xc, double yc)
{
  for (int jj = 0; jj < NV; jj++)
    y[jj] = 0.;
  double s = 0.25 + 0.75*xc;
  y[OpenSMOKE_IndexOfSpecies ("N2")]  = 0.767*(1. - s) + 0.700*s;
  y[OpenSMOKE_IndexOfSpecies ("O2")]  = 0.233*(1. - s) + 0.212*s;
  y[OpenSMOKE_IndexOfSpecies ("CO")]  = 0.050*s;
  y[OpenSMOKE_IndexOfSpecies ("H2")]  = 0.004*s;
  y[OpenSMOKE_IndexOfSpecies ("CH4")] = 0.014*s;
  y[OpenSMOKE_IndexOfSpecies ("CO2")] = 0.010*s;
  y[OpenSMOKE_IndexOfSpecies ("H2O")] = 0.010*s;
  y[NGS] = 900. + 150.*yc;
}

static void load_state (const double * S, long ncells, int k)
{
  long c = 0;
  foreach (serial) {
    const double * ys = S + ((long) k*ncells + c)*NV;
    for (int jj = 0; jj < NGS; jj++) {
      scalar YG = YGList_G[jj];
      YG[] = ys[jj];
    }
    TG[] = ys[NGS];
    c++;
  }
}

static void save_temperature (double * Tc)
{
  long c = 0;
  foreach (serial)
    Tc[c++] = TG[];
}

/**
The reference step of one cell: the calls of the per-cell sweep of
`chemistry.h` (with `f = 0` and no gate). */

static void reference_step (double * y)
{
  double y0ode[NV];
  for (int jj = 0; jj < NV; jj++)
    y0ode[jj] = y[jj];

  UserDataODE data;
  data.P = Pref;
  data.T = y0ode[NGS];
  data.sources = NULL;
  data.rhog = 1.;
  data.cpg = 1000.;
  OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure, NV, DTB,
                       y0ode, &data);

  for (int jj = 0; jj < NV; jj++)
    if (!isfinite (y0ode[jj]))
      return;
  for (int jj = 0; jj < NGS; jj++)
    y[jj] = y0ode[jj] > 0. ? y0ode[jj] : 0.;
  y[NGS] = y0ode[NGS];
}

/**
The first step at which the temperature of a cell is `RISE` above its
initial value, or -1. */

static int ignition_step (const double * Tk, int k, int previous,
                          double T0)
{
  if (previous >= 0)
    return previous;
  return Tk[0] > T0 + RISE ? k : -1;
}

event ignition (i = 0) {

  long ncells = 0;
  foreach (serial)
    ncells++;

  /**
  ## The reference trajectory */

  double * S = malloc ((long) (NSTEPS + 1)*ncells*NV*sizeof (double));
  long c = 0;
  foreach (serial) {
    initial_state (state (S, ncells, 0, c), x - X0, y - Y0);
    c++;
  }
  for (int k = 0; k < NSTEPS; k++)
    for (c = 0; c < ncells; c++) {
      double * y1 = state (S, ncells, k + 1, c);
      const double * y0 = state (S, ncells, k, c);
      for (int jj = 0; jj < NV; jj++)
        y1[jj] = y0[jj];
      reference_step (y1);
    }

  int * kref = malloc (ncells*sizeof (int));
  int nign = 0, kmin = NSTEPS, kmax = -1;
  for (c = 0; c < ncells; c++) {
    double T0 = state (S, ncells, 0, c)[NGS];
    kref[c] = -1;
    for (int k = 1; k <= NSTEPS; k++)
      kref[c] = ignition_step (state (S, ncells, k, c) + NGS, k, kref[c], T0);
    if (kref[c] >= 0) {
      nign++;
      kmin = min (kmin, kref[c]);
      kmax = max (kmax, kref[c]);
    }
  }

  fprintf (stderr, "# binning-ignition: %ld cells, %d steps of %g s, "
           "Solid-gas-88 (%d species)\n", ncells, NSTEPS, DTB, NGS);
  fprintf (stderr, "# reference: %d cells ignite (T rises by %g K), "
           "first at step %d, last at step %d\n", nign, RISE, kmin, kmax);
  fprintf (stderr, "# case tol M | NCA/Nc alpha_min | Eloc_max Eloc_mean |"
           " Eglob_max Eglob_end_mean | nign dk_max dk_mean\n");
  fprintf (stderr, "# E = |T - Tref|/Tref (Eq. 10); nign: cells that"
           " ignite; dk: shift of the ignition in steps (%d: one of the"
           " two does not ignite)\n", NSTEPS);

  /**
  ## The cases of the binning */

  struct { const char * name; double tol; int M; bool forced; } cases[] = {
    {"A", 0.01, 4, false}, {"B", 0.05, 4, false}, {"C", 0.10, 4, false},
    {"E", 0.05, 11, false}, {"F", 0.05, 4, true}, {"one", 1e-9, 1000, false}
  };
  int ncases = sizeof (cases)/sizeof (cases[0]);

  double * Tc = malloc (ncells*sizeof (double));
  int * kbin = malloc (ncells*sizeof (int));

  for (int n = 0; n < ncases; n++) {
    binning_tol = cases[n].tol;
    binning_M = cases[n].M;
    binning_forced[0] = cases[n].forced ? "H2O2" : NULL;
    binning_forced[1] = cases[n].forced ? "OH" : NULL;
    binning_forced[2] = NULL;

    /**
    Local error (Eq. 10): one step from the reference state. */

    double elmax = 0., elsum = 0., amin = 1.;
    long nca = 0, nc = 0;
    for (int k = 0; k < NSTEPS; k++) {
      load_state (S, ncells, k);
      gas_phase_reactions (DTB, false);
      nca += binning_nbins, nc += binning_ncells;
      amin = min (amin, binning_alpha_min);
      save_temperature (Tc);
      for (c = 0; c < ncells; c++) {
        double Tr = state (S, ncells, k + 1, c)[NGS];
        double e = fabs (Tc[c] - Tr)/Tr;
        elmax = max (elmax, e);
        elsum += e;
      }
    }

    /**
    Global error: the binning marches from the initial state. */

    double egmax = 0., egend = 0.;
    load_state (S, ncells, 0);
    for (c = 0; c < ncells; c++)
      kbin[c] = -1;
    for (int k = 1; k <= NSTEPS; k++) {
      gas_phase_reactions (DTB, false);
      save_temperature (Tc);
      for (c = 0; c < ncells; c++) {
        double Tr = state (S, ncells, k, c)[NGS];
        double e = fabs (Tc[c] - Tr)/Tr;
        egmax = max (egmax, e);
        if (k == NSTEPS)
          egend += e/ncells;
        kbin[c] = ignition_step (Tc + c, k, kbin[c],
                                 state (S, ncells, 0, c)[NGS]);
      }
    }

    int nignb = 0, dkmax = 0;
    double dksum = 0.;
    long ndk = 0;
    for (c = 0; c < ncells; c++) {
      if (kbin[c] >= 0)
        nignb++;
      if (kbin[c] >= 0 && kref[c] >= 0) {
        int dk = abs (kbin[c] - kref[c]);
        dkmax = max (dkmax, dk);
        dksum += dk, ndk++;
      }
      else if (kbin[c] != kref[c])
        dkmax = NSTEPS;
    }

    /**
    The case with one state in each bin gives a test, not the errors (see
    above for the thresholds). */

    if (cases[n].tol < 1e-6)
      fprintf (stderr, "%s %g all | %.2f %.2g | local error < 1e-3: %s |"
               " global error < 1e-2: %s | same ignition steps: %s\n",
               cases[n].name, cases[n].tol, (double) nca/nc, amin,
               elmax < 1e-3 ? "yes" : "no", egmax < 1e-2 ? "yes" : "no",
               nignb == nign && dkmax == 0 ? "yes" : "no");
    else
      fprintf (stderr, "%s %g %d%s | %.2f %.2g | %.2g %.2g | %.2g %.2g |"
               " %d %d %.2g\n",
               cases[n].name, cases[n].tol, cases[n].M,
               cases[n].forced ? "+H2O2,OH" : "", (double) nca/nc,
               amin, elmax, elsum/(NSTEPS*ncells), egmax, egend,
               nignb, dkmax, ndk ? dksum/ndk : 0.);
  }

  /**
  Clear the gas, so that the `chemistry` event of the next step finds no
  cell. */

  foreach()
    TG[] = 0.;

  free (S), free (kref), free (Tc), free (kbin);
}

/**
The `chemistry` event of the second step still reads `DYDtG_S` (stencil
pass), so free it only here. */

event stop (i = 1) {
  delete (DYDtG_S), free (DYDtG_S), DYDtG_S = NULL;
  return 1;
}
