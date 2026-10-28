/**
# Negative values from the coarse-fine fluxes of an implicit diffusion

The name does not start with `coarse-`: `runtest` removes `coarse-*` and
`fine-*` in the directory of the test, thus also the binary.

An implicit diffusion step keeps $Y \geq 0$ when its matrix is an M-matrix.
On a tree, the fine cell next to a coarser cell C computes its flux against
a ghost value. With the bilinear prolongation, the ghost value is
$$
Y_g = (9 Y_C + 3 Y_P + 3 Y_T + Y_D)/16
$$
$P$ is the coarse parent of the fine cell, $T$ and $D$ are the tangential
and the diagonal neighbours of C. C receives the sum of the fine fluxes, so
its row has the entries $-(3/16)\,\alpha/\Delta_f$ on $Y_P$, $Y_T$ and
$-(1/16)\,\alpha/\Delta_f$ on $Y_D$, with the wrong sign for an M-matrix.

The test puts a step of one species next to a change of level, and gives
the smallest $Y$ after one implicit step for a range of cell Fourier numbers
$Fo = D\,\Delta t/\Delta_f^2$. The left half of the box is at level 7, the
right half at level 6. A uniform grid at level 7 is the control.

* `fine`: $Y = 1$ for $x < 0.5 - k\Delta_f$, thus the step is $k$ fine
  cells away from the change of level, on the fine side.
* `coarse`: $Y = 1$ for $x > 0.5 + k\Delta_c$, on the coarse side.
* `tilt`: as `fine`, but the step makes an angle with the change of level.
* `peak`: $Y = 1$ only in the column of fine cells $k$ cells away from the
  change of level, as a radical in a thin flame.
* `spot`: $Y = 1$ only in one fine cell of that column.
* `ystep`: $Y = 1$ for $y > 0.5$ on both sides, thus a front that crosses
  the change of level at a right angle. (`k` has no effect.)

The solvers:

* `bilinear`: `diffusion()` with the default prolongation.
* `injection`: `diffusion()` with `refine_injection`, thus $Y_g = Y_C$.
* `coupled`: `coupled_diffusion()` with three species, the step on the
  trace species, the bilinear prolongation. The carrier `Yb` takes the
  complement, so the corrective flux is $(D_c - D_b)\nabla Y_c$.

## The result (2026-10-28)

No solver gives a negative value in any configuration, for $Fo$ from 0.25
to 64. The smallest value is at the level of the solver tolerance (6e-14 of
the height of the step). The result is the same with `-DAXI=1`, where the
change of level and the profiles are along the radius. An M-matrix is a
sufficient condition, not a necessary one: the wrong-sign entries are too
small to make the inverse negative.

Thus the coarse-fine fluxes do not explain the negative gas mass fractions
of `run/fatehi-combustion.c`, although the smallest one is in a coarse cell
(level 8 of 11) there. The test stays as a guard: it fails if a change of
the solvers makes a value negative.

Caution: give `diffusion()` the argument `theta = cm` in an axisymmetric
case. Without it, $\lambda = -1/\Delta t$ has no metric. The first coarse
cell next to the change of level then goes to -1.4e-3 of the step at
$Fo = 0.25$. */

#include "grid/quadtree.h"
#if AXI
# include "axi.h"
/* the change of level and the profiles are along the radius y */
# define XC y
# define YC x
#else
# define XC x
# define YC y
#endif
#include "run.h"
#include "diffusion.h"
#include "coupled-diffusion.h"

scalar Ya[], Yb[], Yc[], Da[], Db[], Dc[], thetac[];
face vector alpha[];

enum { FINE, COARSE, TILT, PEAK, SPOT, YSTEP, NSHAPE };
enum { BILINEAR, INJECTION, COUPLED };

static double profile (int shape, double x, double y, int k)
{
  double df = 1./128., dc = 1./64.;
  switch (shape) {
  case FINE:   return x < 0.5 - k*df;
  case COARSE: return x > 0.5 + k*dc;
  case TILT:   return x + 0.25*(y - 0.5) < 0.5 - k*df;
  case PEAK:   return x > 0.5 - (k + 1)*df && x < 0.5 - k*df;
  case SPOT:   return (x > 0.5 - (k + 1)*df && x < 0.5 - k*df &&
                       y > 0.5 - df/2. && y < 0.5 + df/2.);
  default:     return y > 0.5 + 0.3*df;
  }
}

static void make_grid (bool adaptive)
{
  init_grid (64);
  if (adaptive)
    refine (XC < 0.5 && level < 7);
  else
    refine (level < 7);
#if AXI
  /* the test calls no run(), so the metric event of axi.h must be called */
  event ("metric");
#endif
}

/**
One implicit step. The function returns the smallest $Y$ of the species
with the step, and its level and position. */

static double one_step (int solver, int shape, int k, double Fo,
                        int * lmin, double * xmin)
{
  double dt = Fo*sq(1./128.);
  scalar * YL = {Ya, Yb, Yc}, * DL = {Da, Db, Dc};
  for (scalar Y in YL)
    set_prolongation (Y, solver == INJECTION ? refine_injection :
                      refine_bilinear);
  foreach() {
    /* the trace species Yc: the step goes to exactly zero */
    Yc[] = profile (shape, XC, YC, k)*1e-3;
    Ya[] = 0.2;
    Yb[] = 0.8 - Yc[];
    Da[] = 1., Db[] = 0.2, Dc[] = 1.;
    thetac[] = cm[];
  }
  foreach_face()
    alpha.x[] = fm.x[];

  scalar Y = Yc;
  if (solver == COUPLED)
    coupled_diffusion (YL, dt, DL, alpha, thetac, tolerance = 1e-13/dt);
  else {
    face vector Dface[];
    foreach_face()
      Dface.x[] = fm.x[];
    /* diffusion() multiplies theta by 1/dt in place: give it a copy */
    scalar th[];
    foreach()
      th[] = cm[];
    diffusion (Y, dt, Dface, theta = th, tolerance = 1e-13/dt);
  }

  double ymin = HUGE;
  foreach (reduction(min:ymin))
    ymin = min (ymin, Y[]);
  *lmin = 0, *xmin = 0.;
  foreach (serial)
    if (Y[] == ymin)
      *lmin = level, *xmin = XC;
  return ymin/1e-3;
}

int main()
{
  const char * sname[] = {"fine", "coarse", "tilt", "peak", "spot", "ystep"};
  const char * vname[] = {"bilinear", "injection", "coupled"};
  double Fos[] = {0.25, 1., 4., 16., 64.};
  int ks[] = {0, 1, 2, 4};

  /**
  The control: a uniform grid. No value can be negative. */

  make_grid (false);
  double umin[3] = {HUGE, HUGE, HUGE};
  for (int solver = BILINEAR; solver <= COUPLED; solver += 2)
    for (int shape = FINE; shape < NSHAPE; shape++)
      for (int ik = 0; ik < 4; ik++)
        for (int iF = 0; iF < 5; iF++) {
          int l; double xm;
          umin[solver] = min (umin[solver],
                              one_step (solver, shape, ks[ik], Fos[iF],
                                        &l, &xm));
        }
  fprintf (stderr, "uniform grid: smallest Y/Ymax, bilinear %.3e, "
           "coupled %.3e\n", umin[BILINEAR], umin[COUPLED]);

  /**
  The tree. The table gives the smallest $Y$, relative to the height of the
  step, and the level and $x$ of that cell. */

  make_grid (true);
  double worst[3][NSHAPE] = {{0}};
  for (int solver = BILINEAR; solver <= COUPLED; solver++)
    for (int shape = FINE; shape < NSHAPE; shape++)
      for (int ik = 0; ik < 4; ik++) {
        fprintf (stderr, "%-9s %-6s k=%d:", vname[solver], sname[shape],
                 ks[ik]);
        for (int iF = 0; iF < 5; iF++) {
          int l; double xm;
          double m = one_step (solver, shape, ks[ik], Fos[iF], &l, &xm);
          worst[solver][shape] = min (worst[solver][shape], m);
          if (m < -1e-10)
            fprintf (stderr, "  Fo %-4g %+.2e (L%d x %.4f)", Fos[iF], m, l,
                     xm);
          else
            fprintf (stderr, "  Fo %-4g %+.2e            ", Fos[iF],
                     max (m, 0.));
        }
        fprintf (stderr, "\n");
      }

  for (int solver = BILINEAR; solver <= COUPLED; solver++)
  {
    fprintf (stderr, "%-9s worst Y/Ymax:", vname[solver]);
    for (int shape = FINE; shape < NSHAPE; shape++)
      fprintf (stderr, " %s %+.2e", sname[shape], worst[solver][shape]);
    fprintf (stderr, "\n");
  }

  assert (umin[BILINEAR] > -1e-10 && umin[COUPLED] > -1e-10);
  for (int solver = BILINEAR; solver <= COUPLED; solver++)
    for (int shape = FINE; shape < NSHAPE; shape++)
      assert (worst[solver][shape] > -1e-10);
}
