/**
# The ghost cells of phase values

The test checks `phase-prolongation.h`. A narrow band of refinement follows
a circle, so the parents next to the band are cut by the interface, as in
the phase-change cases. Three species hold the value of the gas phase
outside the circle, with $\sum Y = 1$ there and $Y = 0$ inside.

For each choice of restriction and prolongation, the test gives the largest
$|\sum Y - 1|$ that an open face sees on its two sides, when both sides
hold the phase. A ghost cell at a change of level is one of those sides.

1. The functions of the VOF tracers (`restriction_volume_average`,
   `vof_concentration_refine`).
2. The plain average and the bilinear prolongation.
3. `restriction_phase` and `refine_phase_bilinear`: the error must be at
   round-off.

Then the coupled solve of `coupled-diffusion.h` runs with choice 3, and it
must keep $\sum Y = 1$ in every gas cell to the solver tolerance. */

#include "grid/quadtree.h"
#include "run.h"
#include "fractions.h"

scalar * interfaces = NULL;
face vector uf[];
#include "vof.h"

#define F_ERR 1e-10
#include "phase-prolongation.h"
#include "coupled-diffusion.h"
#include "utils.h"

#define NS 3

scalar f[], Ya[], Yb[], Yc[], Da[], Db[], Dc[];
scalar thetac[];
face vector fs[], alpha[];
double MWs[NS] = {2., 28., 162.};

static void initial (scalar * YL)
{
  scalar a = YL[0], b = YL[1], c = YL[2];
  foreach() {
    if (1. - f[] > F_ERR) {
      double s = 0.5*(1. + tanh ((x - 0.5)/0.05));
      a[] = 0.05 + 0.6*s;
      b[] = 0.85 - 0.8*s;
      c[] = 1. - a[] - b[];
    }
    else
      a[] = b[] = c[] = 0.;
  }
}

static double face_error (scalar * YL, int * nbad)
{
  double e = 0.;
  int n = 0;
  foreach_face (reduction(max:e) reduction(+:n)) {
    if (fm.x[] > 0. && 1. - fs.x[] > 0. &&
        1. - f[-1] > F_ERR && 1. - f[] > F_ERR) {
      double sm = 0., sp = 0.;
      for (scalar Y in YL)
        sm += Y[-1], sp += Y[];
      double ef = max (fabs (sm - 1.), fabs (sp - 1.));
      e = max (e, ef), n += (ef > 1e-12);
    }
  }
  *nbad = n;
  return e;
}

static void set_functions (scalar * YL, int choice)
{
  for (scalar Y in YL) {
    if (choice == 1) {
      set_restriction (Y, restriction_volume_average);
      set_prolongation (Y, vof_concentration_refine);
    }
    else if (choice == 2) {
      set_restriction (Y, restriction_average);
      set_prolongation (Y, refine_bilinear);
    }
    else {
      set_restriction (Y, restriction_phase);
      set_prolongation (Y, refine_phase_bilinear);
    }
  }
}

int main()
{
  init_grid (16);
  double R = 0.2 + 1e-3;
  refine (fabs (sqrt (sq(x - 0.5) + sq(y - 0.5)) - R) < Delta && level < 7);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(R) - sq(x - 0.5) - sq(y - 0.5);
  fractions (phi, f, fs);
  restriction ({f});

  int nl[8] = {0};
  foreach()
    nl[level]++;
  fprintf (stderr, "cells per level: 4:%d 5:%d 6:%d 7:%d\n",
           nl[4], nl[5], nl[6], nl[7]);

  scalar * YL = {Ya, Yb, Yc}, * DL = {Da, Db, Dc};
  for (scalar Y in YL) {
    Y.c = f;
    Y.inverse = true;
  }

  const char * name[4] = {"", "VOF tracer", "average + bilinear",
                          "phase functions"};
  double err[4];
  for (int choice = 1; choice <= 3; choice++) {
    set_functions (YL, choice);
    initial (YL);
    int nbad;
    err[choice] = face_error (YL, &nbad);
    fprintf (stderr, "%-20s: max |sum Y - 1| seen by a face %.3e, faces %d\n",
             name[choice], err[choice], nbad);
  }
  assert (err[3] < 1e-12);
  assert (err[2] > 1e-6);

  /**
  The coupled solve with the phase functions. */

  set_functions (YL, 3);
  foreach() {
    thetac[] = cm[]*max (1. - f[], F_ERR);
    Da[] = 1., Db[] = 0.2, Dc[] = 0.05;
  }
  foreach_face()
    alpha.x[] = (1. - fs.x[])*fm.x[];

  for (int molar = 1; molar >= 0; molar--) {
    initial (YL);
    double emax = 0., ymin = HUGE;
    int cycles = 0;
    for (int s = 0; s < 3; s++) {
      mgstats st = coupled_diffusion (YL, 1e-3, DL, alpha, thetac,
                                      MWs = molar ? MWs : NULL,
                                      tolerance = 1e-8);
      assert (st.resa <= 1e-8);
      cycles = max (cycles, st.i);
      foreach (reduction(max:emax) reduction(min:ymin))
        if (1. - f[] > F_ERR) {
          double sum = 0.;
          for (scalar Y in YL)
            sum += Y[], ymin = min (ymin, Y[]);
          emax = max (emax, fabs (sum - 1.));
        }
    }
    fprintf (stderr, "coupled solve, %s, phase functions: cycles %d, "
             "max |sum Y - 1| %.3e, min Y %.3e\n",
             molar ? "molar" : "Fick ", cycles, emax, ymin);
    assert (emax < 1e-9 && ymin >= 0.);
  }
}
