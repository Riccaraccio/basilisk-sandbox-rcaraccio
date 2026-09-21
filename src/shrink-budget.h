/**
# Probe NEW-1 and 9d: the solid volume that the VOF sweep removes

Items NEW-1, NEW-4 and KNOWN 9d of
`~/discretization-report/formulation-coherence.md`, section 3. This header
changes nothing in the run. It writes one file.

## The question

The chemistry assigns a part `zeta` of the decomposition to the shrinkage of
the solid. `project_sv()` (`velocity-potential.h`) turns that part into the
source

    prod = omega*f*zeta*cm/rhoS

of the Poisson equation of the velocity potential `psi`, and then into the
solid velocity `ubf`. The VOF sweep of `$BASILISK/vof.h` removes the solid
material with the discrete divergence of `ubf`. Three known defects make the
volume that the sweep removes different from the volume that the chemistry
assigned:

* NEW-1. With `shift_prod = true`, `shift_field()` moves `prod` from each cut
  cell to the pure solid cells of its 3x3 stencil. The receiver removes the
  donated volume with ITS porosity, not with the porosity of the donor.
* KNOWN 9d. `vof.h:353` sets `cc = (c > 0.5)`. A cut cell with no pure
  receiver keeps its `prod`, and it removes no volume at all when `f <= 0.5`.
* NEW-4. The stop test of the `psi` solve is absolute, `TOLERANCE_SOLID`
  = 1e-5 in the units of `prod`. The residual of that solve enters the
  removed volume of every cell with `f > 0.5`.

The measured drift of the solid mass is -1.5 per cent at t = 10 s. This probe
says which of the three defects carries it.

## The discrete identity

For the tracer `porosity` of the volume fraction `f`, the sweep gives

    c[]  += dt*(flux[] - flux[1] + cc[]*(ubf.x[1] - ubf.x[]))/(cm[]*Delta)
    t[]  += dt*(tflux[] - tflux[1] + tc[]*(ubf.x[1] - ubf.x[]))/(cm[]*Delta)

with `cc = (c > 0.5)` and `tc = (c > 0.5) ? t/c : 0`. The fluxes cancel in
the sum over the domain, and `cc - tc = (c > 0.5)*(1 - eps)` with
`eps = porosity/f`. So the solid material of the domain changes by

    d/dt sum (f - porosity)*cm*sq(Delta)
      = sum_{f > 0.5} (1 - eps)*D*sq(Delta),   D = sum_d (ubf.d[1] - ubf.d[])/Delta

The Poisson solve gives `D = -(prod - r)`, where `r` is its residual. The
probe therefore measures

    Mtgt = dt*rhoS*sum_all (1 - eps)*prod0*sq(Delta)
    Mact = dt*rhoS*sum_{f > 0.5} (1 - eps)*(-D)*sq(Delta)

`prod0` is the value of `prod` BEFORE the shift. `Mtgt` is the solid mass
that the chemistry assigned to the shrinkage in this step. `Mact` is the mass
that the sweep really removes. Both are in kg per radian, as
`solid_mass0` of `run/test.c` is.

## The three buckets

`Mdiff = Mtgt - Mact` splits into three parts and one closure error:

* `Mshift` the porosity of the receiver against the porosity of the donor
  (NEW-1). It covers every cut cell that has at least one pure receiver:

      Mshift = dt*rhoS*sum_donors prod0*(s_donor - mean s_receiver)*sq(Delta)

  with `s = 1 - eps`. It is positive when the donors hold less solid material
  than the ring of full cells that receives from them.
* `Mnorecv` the cut cells with no pure receiver and `f <= 0.5` (KNOWN 9d).
  Those cells keep their `prod` and remove nothing. The part is always
  positive: the solid stays too heavy. With `shift_prod = false` every cut
  cell keeps its `prod`, so this part then covers every cut cell with
  `f <= 0.5` and `Mshift` is zero.
* `Mres` the residual of the `psi` solve (NEW-4), summed over the cells with
  `f > 0.5`. Either sign.
* `Mclose = Mdiff - Mshift - Mnorecv - Mres` the closure error. It holds the
  parts that the three formulas above do not cover, for example a shift over
  a jump of the level, where the donor and the receiver do not have the same
  size. Read it as a check: a large `Mclose` says that the split is not
  complete for that state.

The file holds the value of each part for the step, and the running integral
of each part from the start of the run.

## How to read `shrinkbudget.dat`

    #t(1) i(2) dt(3) Mtgt(4) Mact(5) Mdiff(6) Mshift(7) Mnorecv(8) Mres(9)
     Mclose(10) Ctgt(11) Cact(12) Cdiff(13) Cshift(14) Cnorecv(15) Cres(16)
     Cclose(17) ncut(18) nnorecv(19) nlow(20)

The file holds one line for each output time. Columns 4 to 10 are the
values of the step that starts at that time, in kg per radian. Columns 11 to
17 are their running integrals over the run. The running integrals cover
every step, also the steps that the file does not show. Column 18 gives the number of
cut cells, column 19 the number of cut cells with no pure receiver, and
column 20 the number of those that also have `f <= 0.5`.

`Cdiff` is the quantity to read first. Divide it by `solid_mass0` and compare
it with the drift of the mass balance. `Cdiff > 0` means that the sweep
removes less than the chemistry assigned, so the solid of the field is too
heavy. Then read which of `Cshift`, `Cnorecv` and `Cres` carries `Cdiff`.

Caution: a restart sets every running integral back to zero. Compare the
integrals of one continuous run only.

Caution: `zeta` lags `omega` by one step, and `set_zeta()` runs at the end of
the `phasechange` event of `shrinking.h`. The probe reads `zeta` before that
call, so it reads the same value that `project_sv()` uses.

## Flags

* `SHRINK_BUDGET` 1 turns the probe on. The default is 0.
* `SHRINK_BUDGET_EVERY_STEP` 1 writes one line for every step instead of one
  line for each output time. The default is 0.
* `SHRINK_BUDGET_FILE` the name of the file.

## Cost

The running integrals need every step, so the probe works at every step: two
sweeps of the grid, and one loop over the 3x3 stencil of each cut cell. The
probe adds one field to each cell. It writes one line for each output
time. */

#ifndef SHRINK_BUDGET
# define SHRINK_BUDGET 0
#endif

#if SHRINK_BUDGET

#include <time.h>

#ifndef SHRINK_BUDGET_EVERY_STEP
# define SHRINK_BUDGET_EVERY_STEP 0
#endif

#ifndef SHRINK_BUDGET_FILE
# define SHRINK_BUDGET_FILE "shrinkbudget.dat"
#endif

/**
`sb_prod0` holds `prod` before the shift. The solid fraction of the matrix,
`1 - eps`, needs no field: `shrinking.h` clamps `f` and `porosity` in place
before `project_sv()`, so the `vof` event below reads the clamped values
that the loop of `phasechange` computes in local variables. */

scalar sb_prod0[];

static bool sb_want = false;
static double sb_cpu = 0.;
static int sb_ncall = 0;

/**
The arm of the writer. It runs at the start of the step that starts at an
output time. The interval is the one of the output events of `run/test.c`,
so the event adds no new event time and does not change `dt`. */

event shrink_budget_arm (t += 0.01) {
  sb_want = true;
}

double sb_Mtgt = 0., sb_Ctgt = 0.;
double sb_Cact = 0., sb_Cshift = 0., sb_Cnorecv = 0., sb_Cres = 0.;
double sb_Cdiff = 0., sb_Cclose = 0.;

event defaults (i = 0) {
  sb_prod0.nodump = true;
}

/**
This instance of `phasechange` runs after the `chemistry` event and before
the `phasechange` event of `shrinking.h`, because same-name events run in
reverse declaration order. So it reads the `omega` and the `zeta` that
`project_sv()` reads.

`shrinking.h` clamps `f` and `porosity` before it builds `prod`. The loop
below repeats that clamp in local variables. It does not write `f` and it
does not write `porosity`. */

event phasechange (i++) {
  clock_t c0 = clock();
  double Mtgt = 0.;

  foreach (reduction(+:Mtgt)) {
    double fc = clamp (f[], 0., 1.);
    fc = (fc > F_ERR) ? fc : 0.;
    fc = (fc > 1. - F_ERR) ? 1. : fc;
    double pc = clamp (porosity[], 0., 1.);
    pc = (fc > F_ERR) ? pc : 0.;

    double sc = (fc > F_ERR) ? 1. - pc/fc : 0.;
    sb_prod0[] = omega[]*fc*zeta[]*cm[]/rhoS;
    Mtgt += sc*sb_prod0[]*sq(Delta);
  }

  sb_Mtgt = dt*rhoS*Mtgt;
  sb_Ctgt += sb_Mtgt;
  sb_cpu += (double)(clock() - c0)/CLOCKS_PER_SEC;
}

/**
This instance of `vof` runs first in the `vof` event, thus before the
exchange of `uf` and `ubf` of `shrinking.h` and before the sweep of
`vof.h`. So `f` and `porosity` still hold the state of the start of the
step, and `prod`, `psi` and `ubf` hold the values that the sweep uses.

The probe reads `ubf` and not `uf`, so the order against the exchange of
`shrinking.h` does not matter. */

#if TREE
# define SB_CELL (is_leaf(cell) && !is_boundary(cell))
#else
# define SB_CELL (!is_boundary(point))
#endif

event vof (i++) {
  clock_t c0 = clock();
  double Mact = 0., Mshift = 0., Mnorecv = 0., Mres = 0.;
  double ncut = 0., nnorecv = 0., nlow = 0.;

  foreach (reduction(+:Mact) reduction(+:Mshift) reduction(+:Mnorecv)
           reduction(+:Mres) reduction(+:ncut) reduction(+:nnorecv)
           reduction(+:nlow)) {
    double d = 0.;
    foreach_dimension()
      d += ubf.x[1] - ubf.x[];
    d /= Delta;                     // = cm*div(ubf), as gas_source is weighted

    double s = (f[] > F_ERR) ? 1. - porosity[]/f[] : 0.;

    if (f[] > 0.5) {
      Mact += -d*s*sq(Delta);
      Mres += (d + prod[])*s*sq(Delta);
    }

    /**
    The buckets of the cut cells. The count of the receivers repeats the rule
    of `shift_field()`: a pure solid leaf of the 3x3 stencil, inside the
    domain. */

    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      ncut += 1.;
      int count = 0;
      double ssum = 0.;
      if (shift_prod) {
        foreach_neighbor (1)
          if (SB_CELL && f[] > 1. - F_ERR) {
            count++;
            ssum += 1. - porosity[]/f[];
          }
      }
      if (count > 0)
        Mshift += sb_prod0[]*(s - ssum/count)*sq(Delta);
      else {
        nnorecv += 1.;
        if (f[] <= 0.5) {
          nlow += 1.;
          Mnorecv += sb_prod0[]*s*sq(Delta);
        }
      }
    }
  }

  Mact *= dt*rhoS;
  Mshift *= dt*rhoS;
  Mnorecv *= dt*rhoS;
  Mres *= dt*rhoS;

  double Mdiff = sb_Mtgt - Mact;
  double Mclose = Mdiff - Mshift - Mnorecv - Mres;

  sb_Cact += Mact;
  sb_Cshift += Mshift;
  sb_Cnorecv += Mnorecv;
  sb_Cres += Mres;
  sb_Cdiff += Mdiff;
  sb_Cclose += Mclose;

  bool write = sb_want || SHRINK_BUDGET_EVERY_STEP;
  sb_want = false;
  if (pid() == 0 && write) {
    static FILE * fp = NULL;
    if (!fp) {
      fp = fopen (SHRINK_BUDGET_FILE, restarted ? "a" : "w");
      if (fp == NULL) {
        fprintf (stderr, "Error opening %s\n", SHRINK_BUDGET_FILE);
        exit (1);
      }
      if (!restarted)
        fprintf (fp, "#t(1) i(2) dt(3) Mtgt(4) Mact(5) Mdiff(6) Mshift(7)"
                     " Mnorecv(8) Mres(9) Mclose(10) Ctgt(11) Cact(12)"
                     " Cdiff(13) Cshift(14) Cnorecv(15) Cres(16) Cclose(17)"
                     " ncut(18) nnorecv(19) nlow(20)\n");
    }
    fprintf (fp, "%g %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
                 " %g %g %g\n",
             t, i, dt, sb_Mtgt, Mact, Mdiff, Mshift, Mnorecv, Mres, Mclose,
             sb_Ctgt, sb_Cact, sb_Cdiff, sb_Cshift, sb_Cnorecv, sb_Cres,
             sb_Cclose, ncut, nnorecv, nlow);
    fflush (fp);
  }
  sb_cpu += (double)(clock() - c0)/CLOCKS_PER_SEC;
  sb_ncall++;
}

/**
The cost of the probe, in the log at the end of the run. */

event shrink_budget_cost (t = end) {
  if (pid() == 0)
    fprintf (stderr, "# shrink-budget: %g s CPU in %d steps\n",
             sb_cpu, sb_ncall);
}

#undef SB_CELL

#endif // SHRINK_BUDGET
