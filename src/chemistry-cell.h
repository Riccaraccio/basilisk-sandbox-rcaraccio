/**
# The per-cell sweep of the gas-phase reactions

`chemistry.h` includes this file when the case does not include
`binning.h` and `TURN_OFF_GAS_REACTIONS` is 0. Do not include it directly.
It uses the helpers of `chemistry.h`: `gas_reactor_cell()`,
`add_gas_expansion()` and the tolerances of the frozen-cell gate. */

/**
## The dump of the gas cells

`CHEMISTRY_DUMP` N > 0 writes the gas cells of the first N steps after the
start of the run, for cuSMOKE (M4.1, `make_reference --cells`). Each call of
the gas reactor writes three files on each rank:

    chemdump-<iter>-<half>-<pid>.cells
    chemdump-<iter>-<half>-<pid>.end
    chemdump-<iter>-<half>-<pid>.meta

`<half>` is 1 for the `chemistry` event and 2 for the second half of the
Strang split. It is 0 for the Lie split. The `.cells` file has the format of
`export.py`: int64 ncells, int32 nsp, nsp names of 32 bytes, then for each
cell the doubles x, y, w, T and Y[nsp]. The state is the state that goes into
the reactor, before the gate. `w` is the gas volume of the cell,
`(1 - f)*Delta^2*cm`. The `.end` file has the same format and the same cells
in the same order, with the state that comes out of the reactor: the
solver's output, or the explicit update of a cell that the gate skipped.
It is the raw state, before the write-back clips `Y < 0` to 0. A solve that
diverged leaves its non-finite values there, and the run keeps the state
that went in. The `.meta` file holds four doubles for each cell, in the
same order: P, the gate flag (1 if the gate skipped the solve), the CPU time
of the solve in seconds, and the number of evaluations of the right-hand
side. Rank 0 prints the interval of each dumped call on a line
`CHEMDUMP <iter> <half> dtc <dtc>` of the log: the replay takes it as `--dt`.

Caution: a case that runs on `run.h` alone and calls `dtnext (DT)` without
assigning the result keeps `dt = 1` (the loop of `run.h` never sets `dt`),
so its reactor integrates over 1 s, or 0.5 s per half. Read `dtc` from the
log, not from `DT`.

The dump does not change the run. The solver calls the right-hand side
through `chemistry_dump_rhs()`, which counts the call and calls the same
function. `misc/chemdump-merge.py` joins the files of all the ranks.

Caution: do not use the dump with OpenMP. The cells of one rank go to one
file in the order of the loop. */

#ifndef CHEMISTRY_DUMP
# define CHEMISTRY_DUMP 0
#endif

#if CHEMISTRY_DUMP
# if _OPENMP
#  error "CHEMISTRY_DUMP needs a build without OpenMP"
# endif
#include <stdint.h>
#include <time.h>

static int chemistry_dump_calls = 0;
static long chemistry_dump_nrhs = 0;

static void chemistry_dump_rhs (const double * y, const double dt,
                                double * dy, void * args)
{
  chemistry_dump_nrhs++;
  gas_batch_nonisothermal_constantpressure (y, dt, dy, args);
}

static double chemistry_dump_cpu (void)
{
  struct timespec ts;
  clock_gettime (CLOCK_THREAD_CPUTIME_ID, &ts);
  return ts.tv_sec + 1e-9*ts.tv_nsec;
}

static FILE * chemistry_dump_open (const char * ext, bool second_half)
{
  char name[96];
  snprintf (name, sizeof (name), "chemdump-%d-%d-%d.%s", iter,
            GAS_CHEMISTRY_STRANG ? (second_half ? 2 : 1) : 0, pid(), ext);
  FILE * fp = fopen (name, "wb");
  if (!fp)
    fprintf (stderr, "chemistry: cannot open %s for the dump\n", name);
  return fp;
}

/**
The header of a `.cells` or `.end` file. The number of cells is 0 here:
`chemistry_dump_close()` writes it at the end, from the size of the file. */

static void chemistry_dump_header (FILE * fp)
{
  int64_t n0 = 0;
  int32_t nsp = NGS;
  fwrite (&n0, sizeof (n0), 1, fp);
  fwrite (&nsp, sizeof (nsp), 1, fp);
  for (int jj = 0; jj < NGS; jj++) {
    char name[32] = {0};
    strncpy (name, OpenSMOKE_NamesOfSpecies (jj), sizeof (name) - 1);
    fwrite (name, 1, sizeof (name), fp);
  }
}

static void chemistry_dump_close (FILE * fp)
{
  long header = sizeof (int64_t) + sizeof (int32_t) + 32*NGS;
  int64_t n = (ftell (fp) - header)/(long) (sizeof (double)*(4 + NGS));
  fseek (fp, 0, SEEK_SET);
  fwrite (&n, sizeof (n), 1, fp);
  fclose (fp);
}
#endif // CHEMISTRY_DUMP

/**
## The sweep of the gas-phase reactions

This function holds the loop of the gas-phase reactions of the external gas.
`dtc` is the time over which the reactor integrates: `dt/2` for each half
of the Strang split, `dt` for the Lie split. The expansion source
always divides by the step `dt`, so each half gives its part of the mean of
the step. `second_half` is true only for the second half of the Strang split.
That half writes its expansion directly into `drhodt`, because
`update_divergence()` has already run. See "The Strang split of the gas
chemistry" in `chemistry.h`. */

static void gas_phase_reactions (double dtc, bool second_half)
{
#if CHEMISTRY_DUMP
  FILE * dump_cells = NULL, * dump_end = NULL, * dump_meta = NULL;
  if (chemistry_dump_calls < (GAS_CHEMISTRY_STRANG ? 2 : 1)*CHEMISTRY_DUMP) {
    chemistry_dump_calls++;
    dump_cells = chemistry_dump_open ("cells", second_half);
    dump_end = chemistry_dump_open ("end", second_half);
    dump_meta = chemistry_dump_open ("meta", second_half);
    if (dump_cells && dump_end && dump_meta) {
      chemistry_dump_header (dump_cells);
      chemistry_dump_header (dump_end);
      if (pid() == 0)
        fprintf (stderr, "CHEMDUMP %d %d dtc %.17g\n", iter,
                 GAS_CHEMISTRY_STRANG ? (second_half ? 2 : 1) : 0, dtc);
    }
    else {
      if (dump_cells)
        fclose (dump_cells), dump_cells = NULL;
      if (dump_end)
        fclose (dump_end), dump_end = NULL;
      if (dump_meta)
        fclose (dump_meta), dump_meta = NULL;
    }
  }
#endif // CHEMISTRY_DUMP

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

#if CHEMISTRY_DUMP
      /**
      Keep the state that goes into the reactor for the dump. */

      double dump_y[NGS + 1], dump_cpu = 0., dump_nrhs = 0.;
      if (dump_cells)
        for (int jj = 0; jj < NGS + 1; jj++)
          dump_y[jj] = y0ode[jj];
#endif

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
      {
#if CHEMISTRY_DUMP
        if (dump_cells) {
          long nrhs0 = chemistry_dump_nrhs;
          double cpu0 = chemistry_dump_cpu ();
          OpenSMOKE_ODESolver (&chemistry_dump_rhs, NGS + 1, dtc, y0ode, &data);
          dump_cpu = chemistry_dump_cpu () - cpu0;
          dump_nrhs = chemistry_dump_nrhs - nrhs0;
        }
        else
#endif
        OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure, NGS + 1, dtc, y0ode, &data);
      }

#if CHEMISTRY_DUMP
      if (dump_cells) {
        double row[4 + NGS];
        row[0] = x, row[1] = y;
        row[2] = (1. - f[])*sq(Delta)*cm[];
        row[3] = dump_y[NGS];
        for (int jj = 0; jj < NGS; jj++)
          row[4 + jj] = dump_y[jj];
        fwrite (row, sizeof (double), 4 + NGS, dump_cells);
        row[3] = y0ode[NGS];
        for (int jj = 0; jj < NGS; jj++)
          row[4 + jj] = y0ode[jj];
        fwrite (row, sizeof (double), 4 + NGS, dump_end);
        double meta[4] = {data.P, frozen ? 1. : 0., dump_cpu, dump_nrhs};
        fwrite (meta, sizeof (double), 4, dump_meta);
      }
#endif

      /**
      Keep the pre-integration state if the solve diverged (see the solid
      branch in the `chemistry` event of `chemistry.h`). */

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

#if CHEMISTRY_DUMP
  if (dump_cells) {
    chemistry_dump_close (dump_cells);
    chemistry_dump_close (dump_end);
    fclose (dump_meta);
  }
#endif // CHEMISTRY_DUMP

# ifdef CHEMISTRY_LOG
  if ((second_half || !GAS_CHEMISTRY_STRANG) && pid() == 0)
    fprintf (stderr, "GATE %g %d\n", t, frozen_cell_gate_n);
# endif
}
