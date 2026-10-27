/**
# The coupled diffusion on a tree

The test checks `coupled_diffusion()` of `coupled-diffusion.h` on a quadtree
with three levels of refinement. The cases of `test/coupled-diffusion.c`
test the accuracy in 1-D. This test checks the parts that only a tree has:
the restriction of the coefficients and the fluxes between a coarse and a
fine cell.

1. With equal diffusivities, the result must be equal to `diffusion()`.
2. With the diffusivities of `test/coupled-diffusion.c`, upwind face values
   and a large step, both forms must converge, keep $\sum Y = 1$, keep the
   mass of each species, and keep the species positive. */

#include "grid/quadtree.h"
#include "coupled-diffusion.h"
#include "diffusion.h"
#include "utils.h"

#define NS 3

scalar Ya[], Yb[], Yc[], Da[], Db[], Dc[];
scalar thetac[]; // not theta: utils.h has a global theta
face vector alpha[];
double MWs[NS] = {2., 28., 100.};

static void initial (scalar * YL)
{
  scalar a = YL[0], b = YL[1], c = YL[2];
  foreach() {
    double r = sqrt (sq(x - 0.5) + sq(y - 0.5));
    double s = 0.5*(1. + tanh ((r - 0.25)/0.03));
    a[] = 0.05 + 0.8*s;
    b[] = 0.85 - 0.8*s;
    c[] = 0.1;
  }
}

static void set_diffusivities (scalar * DL, double d0, double d1, double d2)
{
  scalar a = DL[0], b = DL[1], c = DL[2];
  foreach()
    a[] = d0, b[] = d1, c[] = d2;
}

int main()
{
  init_grid (32);
  refine (fabs (sqrt (sq(x - 0.5) + sq(y - 0.5)) - 0.25) < 0.08 && level < 7);
  foreach()
    thetac[] = cm[];
  foreach_face()
    alpha.x[] = fm.x[];

  int ncells = 0, lmin = 100, lmax = 0;
  foreach (reduction(+:ncells) reduction(min:lmin) reduction(max:lmax))
    ncells++, lmin = min (lmin, level), lmax = max (lmax, level);
  fprintf (stderr, "cells %d, levels %d to %d\n", ncells, lmin, lmax);

  scalar * YL = {Ya, Yb, Yc}, * DL = {Da, Db, Dc};
  scalar * RL = list_clone (YL);
  double dt = 0.002;

  /**
  ## 1. Equal diffusivities */

  set_diffusivities (DL, 0.3, 0.3, 0.3);
  initial (YL);
  initial (RL);
  for (int s = 0; s < 5; s++) {
    for (scalar R in RL) {
      scalar th[];
      face vector Df[];
      foreach()
        th[] = thetac[];
      foreach_face()
        Df.x[] = alpha.x[]*0.3;
      diffusion (R, dt, D = Df, theta = th, tolerance = 1e-9);
    }
    coupled_diffusion (YL, dt, DL, alpha, thetac, tolerance = 1e-9);
  }
  double d = 0.;
  for (int i = 0; i < NS; i++) {
    scalar Y = YL[i], R = RL[i];
    foreach (reduction(max:d))
      d = max (d, fabs (Y[] - R[]));
  }
  fprintf (stderr, "equal D: max |coupled - diffusion()| = %.3e\n", d);
  assert (d < 1e-9);

  /**
  ## 2. Different diffusivities, large steps */

  set_diffusivities (DL, 1., 0.2, 0.05);
  for (int molar = 0; molar <= 1; molar++) {
    initial (YL);
    double mass0[NS];
    for (int i = 0; i < NS; i++) {
      scalar Y = YL[i];
      mass0[i] = statsf (Y).sum;
    }
    int imax = 0;
    double resmax = 0.;
    for (int s = 0; s < 5; s++) {
      mgstats st = coupled_diffusion (YL, 5.*dt, DL, alpha, thetac,
                                      MWs = molar ? MWs : NULL,
                                      tolerance = 1e-9);
      imax = max (imax, st.i), resmax = max (resmax, st.resa);
    }
    double se = 0., ym = HUGE, dm = 0.;
    foreach (reduction(max:se) reduction(min:ym)) {
      double s = 0.;
      for (scalar Y in YL)
        s += Y[], ym = min (ym, Y[]);
      se = max (se, fabs (s - 1.));
    }
    for (int i = 0; i < NS; i++) {
      scalar Y = YL[i];
      dm = max (dm, fabs (statsf (Y).sum - mass0[i]));
    }
    fprintf (stderr, "%s, upwind, 5 steps of %g: cycles %d, res %.2e, "
             "|sum Y - 1| %.2e, d(mass) %.2e, min Y %.4e\n",
             molar ? "molar" : "Fick ", 5.*dt, imax, resmax, se, dm, ym);
    assert (resmax <= 1e-9 && se < 1e-9 && dm < 1e-10 && ym >= 0.);
  }

  delete (RL), free (RL);
}
