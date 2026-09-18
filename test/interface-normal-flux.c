/**
# The interface flux of a sphere, with one plane or with two

This is the check of item 8 of the report on the consistency of the
discretization. It needs no chemistry, no flow and no OpenSMOKE.

## The issue under test

The interface heat and species sources of `multicomponent-varprop.h` take
the gradient from `ebmgrad`, which uses the `mycs` plane
(`interface_normal`). The old code took the area, the centroid and the
weights of the conductivity from `facet_normal (point, fS, fsS)`, which
uses the face fractions. So one cut cell used two planes.
`interface_source_normal()` in `intgrad.h` now gives the normal of the
sources, and `INTERFACE_NORMAL_MYCS` selects it (1 `mycs`, 0 the old
`facet_normal`).

## The set-up

A sphere of radius `R` in axisymmetric coordinates, as in `run/test.c`. The
temperature is linear in the radius on each side of the interface:

    TS = Ti + gS (r - R)     (solid, r < R)
    TG = Ti + gG (r - R)     (gas,   r > R)

The interface value is `Ti`. The exact heat that leaves the solid through
the interface is `lambdaS gS 4 pi R^2`. The exact heat that enters the gas
is `- lambdaG gG 4 pi R^2`, because `ebmgrad` gives the gradient of the gas
along the normal that points into the solid.

The test builds the two interface sources the way
`interface_temperature_sources()` builds them, and it compares them with the
exact values. It does this for three normals:

- `old`: `facet_normal (point, fS, fsS)`, the code before the fix.
- `new`: `interface_normal (point, fS)`, the plane of `ebmgrad`.
- `src`: `interface_source_normal (point, fS, fsS)`, the function that the
  solver calls. It must be equal to `new` when `INTERFACE_NORMAL_MYCS` is 1
  and equal to `old` when it is 0.

## The errors

All errors are relative.

- `Q`: the error of the heat, integrated over the whole interface. This is
  the check that the report asks for. The larger of the two sides is given.
- `A cell`: the error of the interface area, cell by cell, in the L1 norm:
  `sum |A_num - A_exact| / sum A_exact`. The exact area of the sphere in one
  cell comes from a fine quadrature of the circle. This error does not
  contain the gradient, so it isolates the plane.
- `Q cell`: the error of the heat, cell by cell, in the same L1 norm, with
  the exact heat `lambda g A_exact` of each cell.
- `angle`: the angle between the two normals, mean and largest over the cut
  cells, in degrees.

## The levels

The level `L` of a row gives the cell size of `run/test.c` at `maxlevel = L`:
`Delta = 20 D0/2^L`, with `D0 = 8 mm`. The domain of this test is smaller
(`L0 = 1.25 D0`), so the grid has `2^(L-4)` cells on a side. The run is
therefore short at every level.

The position of the sphere relative to the grid changes the error. Each
level uses `NOFF` positions of the centre along the axis. The table gives
the mean over the positions, and the largest value for the angle.

The column `fb` counts the cut cells where `facet_normal()` has no normal
and returns the diagonal `(1/2, 1/2)`. The column `drop` counts the cut
cells where `ebmgrad` has no stencil and returns a zero gradient.

## What it checks

1. `src` is equal to `new` (flag 1) or to `old` (flag 0) in every cell.
2. The error of the integrated heat of `src` falls by `RATE_MIN` or more
   from one level to the next.
3. The errors of `src` (integrated heat, area of each cell) are not more
   than `REGRESSION_MAX` times the errors of `old`.

The test does not ask the `mycs` plane to be MORE accurate than the old
plane. On this smooth sphere it is not: read the results below.

## What it measured, 2026-09-18

Flag 1, 8 positions per level. The build with flag 0 gives the same `old`
and `new` columns, and `Q src` is then equal to `Q old`.

| L | R/dx | Q old | Q new | Acell old | Acell new | Qcell old | Qcell new | angle mean | angle max | fb | drop |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 9 | 12.8 | 1.26e-2 | 1.28e-2 | 1.48e-2 | 1.53e-2 | 2.24e-2 | 2.27e-2 | 0.28 | 4.50 | 0 | 0 |
| 10 | 25.6 | 5.80e-3 | 5.83e-3 | 8.19e-3 | 8.49e-3 | 1.18e-2 | 1.21e-2 | 0.12 | 2.98 | 0 | 0 |
| 11 | 51.2 | 2.67e-3 | 2.63e-3 | 4.62e-3 | 4.78e-3 | 6.29e-3 | 6.49e-3 | 0.12 | 87.1 | 2 | 1 |
| 12 | 102.4 | 1.34e-3 | 1.39e-3 | 2.19e-3 | 2.26e-3 | 3.04e-3 | 3.13e-3 | 0.17 | 94.1 | 8 | 3 |
| 13 | 204.8 | 6.86e-4 | 7.54e-4 | 1.09e-3 | 1.12e-3 | 1.48e-3 | 1.52e-3 | 0.07 | 87.9 | 4 | 1 |

(Levels 12 and 13: `-DLEVEL_MIN=12 -DLEVEL_MAX=13`.)

1. The two normals agree to 0.1 to 0.3 degrees on average. The largest
   angle in a well-cut cell is 3 to 4.5 degrees. The angles near 90 degrees
   are the `fb` cells, with `f < 1e-6`, where the old normal is the
   diagonal.
2. The error of the integrated heat is first order for both normals, and
   the two differ by 1 to 10 per cent of the error. The one-sided gradient
   of `ebmgrad` sets this error, not the plane: along a straight line that
   does not go through the centre, the radius is not linear.
3. The area of the old normal is about 3 per cent more accurate, cell by
   cell. The face fractions average the planes of two cells, and on a
   smooth interface that is a good normal.

So the fix gives one plane per cell. It does not give a more accurate flux
on a smooth interface.

The test prints the table on standard error and returns 1 if a check
fails. */

#include "grid/multigrid.h"
#include "axi.h"
#include "run.h"

/**
`intgrad.h` reads the `inverse` attribute of the field. The full solver sets
it in `memoryallocation-varprop.h`. Here we declare it and set it by hand. */

attribute {
  bool inverse;
}

#include "fractions.h"
#include "intgrad.h"

#ifndef F_ERR
# define F_ERR 1.e-10
#endif

#ifndef LEVEL_MIN
# define LEVEL_MIN 9
#endif
#ifndef LEVEL_MAX
# define LEVEL_MAX 11
#endif
#ifndef NOFF
# define NOFF 8
#endif
#ifndef NQUAD
# define NQUAD 20000   // quadrature points of the circle, per cell
#endif

scalar f[], fS[], fG[], TS[], TG[], Aexact[];
face vector fsS[], fsG[];

const double D0 = 8e-3;
const double R = 4e-3;
const double Ti = 800., gS = 2e4, gG = -5e4;  // K and K/m
const double lambdaS = 0.2, lambdaG = 0.05;    // W/m/K

int level, nfail = 0;
double xoff;

/**
The three normals, in the order of the columns. */

enum { OLD, NEW, SRC, NNORM };

coord pick_normal (Point point, int k)
{
  if (k == OLD)
    return facet_normal (point, fS, fsS);
  if (k == NEW)
    return interface_normal (point, fS);
  return interface_source_normal (point, fS, fsS);
}

/**
The exact area of the sphere in one cell, for the whole revolution. The
circle is `(xoff + R cos(th), R sin(th))` and the area element is
`2 pi R^2 sin(th) dth`. The quadrature covers a window of angle around the
cell, which is larger than the cell. */

double exact_cell_area (double xc, double yc, double h)
{
  double r = sqrt (sq(xc - xoff) + sq(yc));
  if (fabs (r - R) > h)
    return 0.;
  double thc = atan2 (yc, xc - xoff), w = 2.*h/R;
  double th0 = max (0., thc - w), th1 = min (pi, thc + w);
  double dth = (th1 - th0)/NQUAD, A = 0.;
  for (int q = 0; q < NQUAD; q++) {
    double th = th0 + (q + 0.5)*dth;
    double px = xoff + R*cos(th), py = R*sin(th);
    if (fabs (px - xc) <= h/2. && fabs (py - yc) <= h/2.)
      A += 2.*pi*sq(R)*sin(th)*dth;
  }
  return A;
}

/**
The errors of one normal. The loop is the loop of
`interface_temperature_sources()`, with the isotropic conductivity. The
weights of the anisotropic conductivity reduce to `lambda` when the two
components are equal. */

typedef struct {
  double Q, Acell, Qcell;
  int ndrop;
} Errors;

Errors interface_errors (int k)
{
  double QS = 0., QG = 0.;
  double dA = 0., sA = 0., dQS = 0., dQG = 0.;
  int nd = 0;
  foreach (reduction(+:QS) reduction(+:QG) reduction(+:dA) reduction(+:sA)
           reduction(+:dQS) reduction(+:dQG) reduction(+:nd)) {
    double AS = 0., qS = 0., qG = 0.;
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      coord n = pick_normal (point, k), p;
      double alpha = plane_alpha (fS[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      bool success = false;
      double Strgrad = ebmgrad (point, TS, fS, fG, fsS, fsG, false, Ti,
                                &success);
      double Gtrgrad = ebmgrad (point, TG, fS, fG, fsS, fsG, true, Ti,
                                &success);
      if (Strgrad == 0. || Gtrgrad == 0.)
        nd++;

      n.x = fabs(n.x); n.y = fabs(n.y);
      double lambda1vh = n.x/(n.x + n.y)*lambdaS + n.y/(n.x + n.y)*lambdaS;
      double lambda2vh = n.x/(n.x + n.y)*lambdaG + n.y/(n.x + n.y)*lambdaG;

      double aov = area*(y + p.y*Delta)/(Delta*y)*cm[];

      /**
      The sources are per unit volume and carry the metric `cm`, so
      `source*Delta^2` is the heat of the cell per radian. */

      AS = 2.*pi*aov*sq(Delta);
      qS = lambda1vh*Strgrad*AS;
      qG = lambda2vh*Gtrgrad*AS;
    }
    QS += qS;
    QG += qG;
    dA += fabs (AS - Aexact[]);
    sA += Aexact[];
    dQS += fabs (qS - lambdaS*gS*Aexact[]);
    dQG += fabs (qG + lambdaG*gG*Aexact[]);
  }

  double Aex = 4.*pi*sq(R);
  double QSex = lambdaS*gS*Aex, QGex = - lambdaG*gG*Aex;
  Errors e;
  e.Q = max (fabs (QS - QSex)/fabs (QSex), fabs (QG - QGex)/fabs (QGex));
  e.Acell = dA/sA;
  e.Qcell = max (dQS/(fabs (lambdaS*gS)*sA), dQG/(fabs (lambdaG*gG)*sA));
  e.ndrop = nd;
  return e;
}

/**
Check 1: the normal of the solver is the one that the flag selects. The
same loop measures the angle between the old and the new normal. */

double angle_mean, angle_max;
int nfallback;

int check_source_normal (void)
{
  int nbad = 0, ncut = 0, nfb = 0;
  double amean = 0., amax = 0.;
  foreach (reduction(+:nbad) reduction(+:ncut) reduction(+:amean)
           reduction(max:amax) reduction(+:nfb))
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      coord a = interface_source_normal (point, fS, fsS);
#if INTERFACE_NORMAL_MYCS
      coord b = interface_normal (point, fS);
#else
      coord b = facet_normal (point, fS, fsS);
#endif
      foreach_dimension()
        if (a.x != b.x)
          nbad++;

      /**
      `facet_normal()` has no normal when the face fractions of the two
      sides of the cell are equal in each direction. This is the case when
      `face_fraction()` sets all four to 0, that is when `f < 1e-6`. It
      then returns the diagonal `(1/2, 1/2)`, whatever the geometry is.
      Count these cells. */

      double nn = 0.;
      foreach_dimension()
        nn += fabs (fsS.x[] - fsS.x[1]);
      if (nn == 0.)
        nfb++;

      coord n1 = facet_normal (point, fS, fsS);
      coord n2 = interface_normal (point, fS);
      normalize (&n1);
      normalize (&n2);
      double c = 0.;
      foreach_dimension()
        c += n1.x*n2.x;
      double ang = acos (clamp (c, -1., 1.))*180./pi;
      amean += ang;
      amax = max (amax, ang);
      ncut++;
    }
  angle_mean += ncut ? amean/ncut/NOFF : 0.;
  angle_max = max (angle_max, amax);
  nfallback += nfb;
  return nbad;
}

Errors err[NNORM];
int ndrop;

/**
The limits of checks 2 and 3. */

#ifndef RATE_MIN
# define RATE_MIN 1.5   // smallest drop of the error from one level to the next
#endif
#ifndef REGRESSION_MAX
# define REGRESSION_MAX 1.1   // largest ratio of the src error to the old error
#endif

int main()
{
  TG.inverse = true;
  TS.inverse = false;

  fprintf (stderr, "# INTERFACE_NORMAL_MYCS %d, %d positions per level\n",
           INTERFACE_NORMAL_MYCS, NOFF);
  fprintf (stderr, "# relative errors, mean over the positions; "
           "angle between the two normals in degrees\n");
  fprintf (stderr, "#%3s %5s | %9s %9s %9s | %9s %9s | %9s %9s |"
           " %5s %6s | %4s %4s\n",
           "L", "R/dx", "Q old", "Q new", "Q src", "Acell old", "Acell new",
           "Qcell old", "Qcell new", "ang", "angmax", "fb", "drop");

  double Qprev = HUGE;
  for (level = LEVEL_MIN; level <= LEVEL_MAX; level++) {
    L0 = 1.25*D0;
    origin (-L0/2., 0.);
    N = 1 << (level - 4);
    for (int k = 0; k < NNORM; k++)
      err[k] = (Errors){0};
    angle_mean = angle_max = 0.;
    nfallback = ndrop = 0;
    for (int o = 0; o < NOFF; o++) {
      xoff = (o + 0.5)/NOFF*L0/N;
      run();
    }

    fprintf (stderr, " %3d %5.1f | %9.2e %9.2e %9.2e | %9.2e %9.2e |"
             " %9.2e %9.2e | %5.2f %6.2f | %4d %4d\n",
             level, R/(L0/N), err[OLD].Q, err[NEW].Q, err[SRC].Q,
             err[OLD].Acell, err[NEW].Acell, err[OLD].Qcell, err[NEW].Qcell,
             angle_mean, angle_max, nfallback, ndrop);

    /**
    Check 2. The integrated heat of the solver converges. */

    if (err[SRC].Q*RATE_MIN > Qprev) {
      fprintf (stderr, "FAIL: at level %d the heat error does not fall by "
               "%g\n", level, (double) RATE_MIN);
      nfail++;
    }
    Qprev = err[SRC].Q;

    /**
    Check 3. The normal of the solver is not worse than the old normal, for
    the integrated heat and for the area of each cell. */

    if (err[SRC].Q > REGRESSION_MAX*err[OLD].Q ||
        err[SRC].Acell > REGRESSION_MAX*err[OLD].Acell) {
      fprintf (stderr, "FAIL: at level %d the normal of the solver is worse "
               "than the old normal\n", level);
      nfail++;
    }
  }

  fprintf (stderr, "\n%s: %d checks failed\n", nfail ? "FAIL" : "PASS",
           nfail);
  return nfail ? 1 : 0;
}

/**
The sphere and the two linear profiles. The value of a cell is the value of
the linear profile at the centre of the cell, in the cells of each phase
only. The solver holds the value of one phase there too: it divides the
tracer form by the phase fraction before it calls `ebmgrad`. */

event init (i = 0)
{
  fraction (f, sq(R) - sq(x - xoff) - sq(y));
  foreach() {
    fS[] = f[];
    fG[] = 1. - f[];
    double r = sqrt (sq(x - xoff) + sq(y));
    TS[] = fS[] > F_ERR ? Ti + gS*(r - R) : 0.;
    TG[] = fG[] > F_ERR ? Ti + gG*(r - R) : 0.;
    Aexact[] = exact_cell_area (x, y, Delta);
  }
  face_fraction (fS, fsS);
  face_fraction (fG, fsG);
}

/**
`run()` frees the grid when it returns, so the errors must be taken in an
event. This event runs after `init`, once per position. */

event measure (i = 0)
{
  if (check_source_normal()) {
    fprintf (stderr, "FAIL: interface_source_normal() does not follow "
             "INTERFACE_NORMAL_MYCS\n");
    nfail++;
  }

  for (int k = 0; k < NNORM; k++) {
    Errors e = interface_errors (k);
    err[k].Q     += e.Q/NOFF;
    err[k].Acell += e.Acell/NOFF;
    err[k].Qcell += e.Qcell/NOFF;
    if (k == NEW)
      ndrop = max (ndrop, e.ndrop);
  }
}

event stop (i = 0) {
  return 1;
}
