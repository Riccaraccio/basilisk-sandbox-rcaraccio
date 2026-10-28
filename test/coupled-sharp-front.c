/**
# The coupled solve at a sharp front of a major species

A step of 0.98 of the light species `Ya` against the carrier `Yb`, with
$D_a/D_b = 5$, at an angle to the grid. The cells on the two sides of the
front have very different $Y^n$, thus the face values $[Y^n]_f$ of the
corrective term are far from the cell values.

The first version of the relaxation of `coupled-diffusion.h` took the cell
value in the diagonal block. On this problem it diverged for $Fo \geq 4$
and stopped with a SIGFPE. The test checks that each solve converges and
keeps $\sum Y = 1$, on a uniform grid at level 7 and on a grid with a change
of level, in the Fick and the molar forms. The front crosses the change of
level at $(0.5, 0.5)$.

$Fo = D_a \Delta t/\Delta^2$ with the finest $\Delta$.

## The result (2026-10-28)

* With the exact diagonal block, every solve on the uniform grid converges
  without a split (at most 64 cycles).
* On the tree, the Fick solves at $Fo = 16$ and $64$ and the molar solves at
  $Fo = 16$ and $64$ fail once and split (`cd_nsplit`), then converge.
* The test does not assert $Y \geq 0$. The converged solution of the Fick
  form on the uniform grid has $Y_c$ down to -2.8e-5 at $Fo = 16$ and to
  -3.2e-4 at $Fo = 64$, thus 2.8 and 32 per cent of the step of $Y_c$. The
  cause is the linearized term $[Y_i^n]_f (U - U^n)$: it is not monotone
  when $U$ changes much in one step. A step of 0.8, or a front of one cell,
  stays positive (`test/coupled-positivity.c` gives the smooth cases). */

#include "grid/quadtree.h"
#include "run.h"
#include "coupled-diffusion.h"

scalar Ya[], Yb[], Yc[], Da[], Db[], Dc[], thetac[];
face vector alpha[];
double MWs[3] = {2., 28., 162.};

int main()
{
  scalar * YL = {Ya, Yb, Yc}, * DL = {Da, Db, Dc};
  double Fos[] = {1., 4., 16., 64.};
  int maxcycles = 0;
  for (int adaptive = 0; adaptive <= 1; adaptive++) {
    init_grid (64);
    if (adaptive)
      refine (x < 0.5 && level < 7);
    else
      refine (level < 7);
    for (int molar = 0; molar <= 1; molar++) {
      fprintf (stderr, "%s %s:", adaptive ? "tree   " : "uniform",
               molar ? "molar" : "Fick ");
      for (int iF = 0; iF < 4; iF++) {
        double dt = Fos[iF]*sq(1./128.), tol = 1e-13/dt;
        foreach() {
          double s = x + 0.25*(y - 0.5) < 0.5;
          Ya[] = 0.01 + 0.98*s;
          Yc[] = 1e-3*s;
          Yb[] = 1. - Ya[] - Yc[];
          Da[] = 1., Db[] = 0.2, Dc[] = 0.1;
          thetac[] = cm[];
        }
        foreach_face()
          alpha.x[] = fm.x[];
        mgstats st = coupled_diffusion (YL, dt, DL, alpha, thetac,
                                        MWs = molar ? MWs : NULL,
                                        tolerance = tol);
        double ymin = HUGE, es = 0.;
        foreach (reduction(min:ymin) reduction(max:es)) {
          double sum = 0.;
          for (scalar Y in YL)
            sum += Y[], ymin = min (ymin, Y[]);
          es = max (es, fabs (sum - 1.));
        }
        fprintf (stderr, "  Fo %-2g cycles %2d min Y %+.1e |sum-1| %.1e",
                 Fos[iF], st.i, ymin, es);
        assert (st.resa <= tol);
        assert (es < 1e-10);
        maxcycles = max (maxcycles, st.i);
      }
      fprintf (stderr, "\n");
    }
  }
  fprintf (stderr, "largest number of cycles: %d, splits: %d\n",
           maxcycles, cd_nsplit);
  assert (cd_nsplit > 0);
}
