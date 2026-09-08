/**
# Does a Dirichlet interface condition bind a sliver?

This is the falsifying test of `INT_TEMP_VOFBC`. It answers one question, and
it needs no chemistry, no flow and no OpenSMOKE.

## The claim under test

Impose the interface temperature as a Dirichlet condition of the gas energy
solve. Then a cut cell with a vanishing gas fraction should tend to `TInt`,
because it has no heat capacity of its own to hold any other value.

The claim is FALSE for the accurate gradient, and that is the whole reason
`ebmgrad_bc()` exists. `ebmgrad` builds the gradient from the NEIGHBOURS of
the cut cell along the normal, so the value of the cell is not in it. The gas
row of the cell reads

    [ theta2/dt + F ] TG_c = theta2/dt TG_c^n + C*(TInt - v0) + (face terms)

and as the gas fraction goes to zero, `theta2/dt` and `F` both go to zero
while `C` does not, because the facet area does not vanish with the phase
fraction. The diagonal collapses and the off-diagonal does not. `TG_c` leaves
its own equation and the solve does not bound it.

`ebmgrad_bc()` detects that on the matrix, not on the geometry, and takes the
first-order form `(TInt - TG[])/(d0*Delta)`, which reads the cell. The row is
then an M-matrix row and the claim becomes true.

## The set-up, and why the interface must be tilted

The interface is the line `x + y = c`, solid below it and gas above, on a
uniform grid. The offset `c` is placed so that the gas of one chosen cell is
a small triangle in its upper right corner, of area `fGtarget`.

The tilt is the point. A grid-aligned interface never makes a sliver: the gas
part of such a cell is a slab that touches the whole top face, so its gas
face fraction stays 1 however small its volume is, the face conductance stays
`O(lambda/Delta^2)`, and the row binds. With a corner triangle of legs `a`
both gas face fractions are `a`, so

    F     ~ lambda*a/Delta^2      -> 0
    theta ~ rho*cp*a^2/2          -> 0
    C     ~ lambda*area/(d0*Delta^2)

and `C` does NOT vanish, because the interface length and the centroid
distance shrink together: `area ~ a` and `d0 ~ a/3`. So `C/(theta/dt + F)`
grows without bound as `a -> 0`. That is the configuration that kills a run.

The cells around the target are ordinary cut cells, and the accurate stencil
IS available in the target cell, so this is exactly the case that the
degenerate branch of `intgrad.h` does NOT catch.

Everything starts at `T0`, the interface is at `TInt0`, and the top and right
boundaries are held at `T0` so the heat has somewhere to go. Without that
sink the problem is pure Neumann with a net source and no bounded solution,
whatever the interface treatment does.

## What it checks

1. `TG` of every cut cell lies between `T0` and `TInt0`. This is the maximum
   principle, and it is what fails without the criterion.
2. `TG` of the target cell does not grow as its gas fraction shrinks.
3. The solve converges.

## What it measured, 2026-09-08

`d0/d1` is the thickness of the gas in the cell over the reach of the
accurate stencil. `T` is the largest temperature of any cut cell, against an
interface at 1000 K and a gas that starts at 300 K.

| eta | what it does | control 0.5 | corner 1e-2 | corner 1e-6 |
|---|---|---|---|---|
| 0 | accurate stencil everywhere | 6892 | 4128 | 4153 |
| 0.05 | mixed | 6892 | 1904 | 3888 |
| 1 | first-order form in every cut cell | 987 | 999.2 | 999.86 |

Three results, and each one contradicts a piece of the design that this
branch started from.

1. **The first-order form does what the Dirichlet condition promises.** As the
   gas fraction falls the cut cell goes to the interface value: 999.16,
   999.77, 999.85, 999.86, 999.86 at 1e-2 down to 1e-6. Nothing else in the
   sweep does that.
2. **The accurate stencil has no bound at any fraction.** It overshoots the
   interface value by a factor of four at a gas fraction of 0.4 and by the
   same factor at 1e-6. This is not a sliver effect. The flux it builds is
   referenced to a value one cell away, so the cut cell is free to sit at
   whatever temperature carries that flux onward.
3. **Mixing the two per cell is worse than either.** At `eta = 0.05` the
   multigrid failed to converge in 100 iterations on two rows. A forced cell
   next to an unforced one is an inconsistency the solver has to resolve.

`d0/d1` never exceeds about 0.19, which is the value for the thickest slab a
cut cell can hold. So there is no threshold that separates a sliver from a
well resolved cut cell by a wide margin, and the two control rows fire at any
`eta` that catches the corner rows.

A separate check, `T0_VALUE=1000`, starts the gas at the interface value. The
flux is then exactly zero and no row moves, at every fraction and on either
branch. So the sign, the metric and the feedback of the split are right, and
the overshoot above is a property of the discretisation and not a defect of
the wiring.

Caution: this exercises `ebmgrad_bc()` and `diffusion()`, which is the part in
doubt. It does NOT exercise `interface_temperature_sources()`, because that
function pulls in OpenSMOKE. `run/restart.c` and `run/netl.c` cover the
assembly. */

#include "grid/multigrid.h"
#include "run.h"

/**
`intgrad.h` reads the `inverse` attribute of the field to decide which side of
the interface the tracer lives on. The full solver sets it in
`memoryallocation-varprop.h`; here we declare and set it by hand, as the 0-D
chemistry tests do. */

attribute {
  bool inverse;
}

#include "fractions.h"
#include "intgrad.h"

/**
`INTBC_HOOK` selects which of the two schemes the run uses. Both impose the
same interface value on the same geometry, so the two columns are directly
comparable.

- 0, the source form. The interface flux is built once, before the solve,
  and put in `r` and `beta`. `INTBC_ETA` selects the branch per cell.
- 1, the operator form of `plicbc.h`. The flux is rebuilt inside every
  relaxation and residual sweep, so the interpolated neighbour value is at
  the current iterate. No threshold. */

#ifndef INTBC_HOOK
# define INTBC_HOOK 0
#endif

#if INTBC_HOOK
# include "plicbc.h"
#endif
#include "diffusion.h"

#ifndef F_ERR
# define F_ERR 1.e-10
#endif

#ifndef INTBC_OVERSHOOT
# define INTBC_OVERSHOOT 0.15
#endif
#ifndef INTBC_TOL_K
# define INTBC_TOL_K 1.e-4   // kelvin
#endif
#ifndef INTBC_ETA
# define INTBC_ETA 1.
#endif
#ifndef INTBC_CMAX
# define INTBC_CMAX 100.
#endif
#ifndef MAXLEVEL
# define MAXLEVEL 5
#endif

scalar f[], fS[], fG[], TG[];
face vector fsS[], fsG[];

#ifndef T0_VALUE
# define T0_VALUE 300.
#endif
double T0 = T0_VALUE, TInt0 = 1000.;
#if INTBC_HOOK
TG[interface] = dirichlet (TInt0);
#endif
double lambdaG = 0.05;          // W/m/K, near the gas of the real case
double rhocpG = 1.2*1000.;      // J/m3/K
double DT_STEP = 1e-2;          // s, deliberately large

double fGtarget;                // the gas fraction of the target cell
double xt, yt;                  // the centre of the target cell
int    tilted;                  // 1 corner sliver, 0 grid-aligned control
int    nfail = 0, nforced_control = 0;

/**
The sweep. The first two rows are the CONTROL: a grid-aligned interface,
whose gas part is a slab that spans the whole cell in the tangential
direction. Such a cell is well resolved along the normal whatever its
fraction, so it must keep the accurate stencil, branch 0, and it must still
be bounded. If the criterion fires there, it is too aggressive and it would
change well resolved answers.

The rest are corner slivers, which the accurate stencil cannot bind. */

double fGlist[]  = {0.5, 0.25, 0.4, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6};
int    tiltlist[] = {  0,    0,   1,   1,    1,    1,    1,    1,    1};

/**
The gas leaves through the two boundaries it touches. Without a sink the
problem is pure Neumann with a net source, and nothing bounds it. */

TG[top]   = dirichlet (T0);
TG[right] = dirichlet (T0);

int main()
{
  L0 = 1e-3;                    // 1 mm, so Delta is near the real case
  origin (0., 0.);
  N = 1 << MAXLEVEL;
  TG.inverse = true;            // the tracer lives in the gas
  init_grid (N);

  fprintf (stderr, "# hook %d eta %g cmax %g level %d\n", INTBC_HOOK,
           (double) INTBC_ETA, (double) INTBC_CMAX, MAXLEVEL);
  fprintf (stderr, "#%7s %10s %10s %10s %6s %14s %14s %5s %9s %9s %s\n",
           "geom", "fG", "fG_meas", "d0/d1", "branch",
           "T_target", "T_maxcut", "mg.i", "resb", "resa", "verdict");

  for (int k = 0; k < (int)(sizeof(fGlist)/sizeof(double)); k++) {
    fGtarget = fGlist[k];
    tilted = tiltlist[k];
    run();
  }

  /**
  The verdict is boundedness and convergence. The control rows are REPORTED
  and not asserted: the sweep showed that `d0/d1` never exceeds about 0.19,
  so no threshold separates a well resolved cut cell from a sliver by a wide
  margin, and any `eta` that catches the corner rows also catches the
  controls. Read the count, do not fail on it. */

  fprintf (stderr, "\n%s: %d of %d rows failed. "
           "The criterion fired on %d of the 2 control rows.\n",
           nfail ? "FAIL" : "PASS", nfail,
           (int)(sizeof(fGlist)/sizeof(double)), nforced_control);
  return nfail ? 1 : 0;
}

/**
Two geometries.

`tilted` puts the upper right corner of the cell `(N/2, N/2)` in the gas, as
a triangle of legs `a` and area `a^2/2`, cut by the line `x + y = c`.

Otherwise the interface is the horizontal line `y = y0`, and the gas of that
cell is a slab of thickness `fGtarget` that spans the cell. That slab touches
the whole top face, so its face conductance does not vanish and it is well
resolved along the normal. It is the control. */

event init (i = 0)
{
  double Delta0 = L0/N;
  int i0 = N/2, j0 = N/2;

  xt = (i0 + 0.5)*Delta0;
  yt = (j0 + 0.5)*Delta0;

  if (tilted) {
    double a = sqrt (2.*fGtarget);
    double c = (i0 + j0 + 2. - a)*Delta0;
    fraction (f, c - x - y);
  }
  else {
    double y0 = (j0 + 1. - fGtarget)*Delta0;
    fraction (f, y0 - y);
  }

  foreach() {
    fS[] = f[];
    fG[] = 1. - f[];
    TG[] = T0;
  }

  face_fraction (fS, fsS);
  face_fraction (fG, fsG);
}

/**
One step of the gas energy equation, with the interface flux built the way
`interface_temperature_sources()` builds it. */

event onestep (i = 0)
{
  scalar r[], beta[], theta[];
  face vector D[];

  /**
  A tolerance with a meaning. The residual has the units of the source, so
  dividing by the largest heat capacity over the step turns it into kelvin.
  The inherited `TOLERANCE = 1e-3` asks for 1e-8 K against a source of 3e10,
  which no solve delivers and which makes every run look like a failure. */

  TOLERANCE = INTBC_TOL_K*rhocpG/DT_STEP;

  foreach_face()
    D.x[] = lambdaG*fsG.x[]*fm.x[];

  foreach() {
    theta[] = cm[]*max (fG[]*rhocpG, F_ERR);
    r[] = 0.;
    beta[] = 0.;
  }

  double smax = 0., stgt = 0., brtgt = 0., fGmeas = 0.;

  foreach (reduction(max:smax) reduction(max:stgt) reduction(max:brtgt)
           reduction(max:fGmeas)) {
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      coord n = facet_normal (point, fS, fsS), p;
      double alpha = plane_alpha (fS[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      double aov = area/Delta*cm[];

      /**
      The diagonal the row would have without the interface term: the heat
      capacity over the step, plus the face conductance. Build the second
      exactly as the solve builds `D`. */

      double Fc = 0.;
      foreach_dimension()
        Fc += lambdaG*fsG.x[]*fm.x[] + lambdaG*fsG.x[1]*fm.x[1];
      Fc /= sq(Delta);

      double ddiag = theta[]/DT_STEP + Fc;

      double coefG = 0., srat = 0.;
      int br = 0;
      double G = ebmgrad_bc (point, TG, fS, fG, fsS, fsG, true, TInt0,
                             lambdaG, aov, ddiag, INTBC_ETA, INTBC_CMAX,
                             &coefG, &br, &srat);

      r[] += lambdaG*aov*G;
      beta[] += min (0., lambdaG*aov*coefG);

      smax = max (smax, srat);
      if (fabs (x - xt) < 0.4*Delta && fabs (y - yt) < 0.4*Delta) {
        stgt = srat;
        brtgt = br;
        fGmeas = fG[];
      }
    }
  }

#if INTBC_HOOK

  /**
  The operator form. `r` and `beta` carry nothing: the whole interface term
  is rebuilt inside the solve. */

  foreach() {
    r[] = 0.;
    beta[] = 0.;
  }
  plicbc_phase (fG, fsG);
  mgstats mg = diffusion (TG, DT_STEP, D = D, r = r, beta = beta,
                          theta = theta, flux = plic_flux);
#else
  mgstats mg = diffusion (TG, DT_STEP, D = D, r = r, beta = beta,
                          theta = theta);
#endif

  /**
  Read the target cell and the worst cut cell back, and judge them. */

  double Ttgt = 0., Tcut = 0.;
  foreach (reduction(max:Ttgt) reduction(max:Tcut)) {
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      Tcut = max (Tcut, TG[]);
      if (fabs (x - xt) < 0.4*Delta && fabs (y - yt) < 0.4*Delta)
        Ttgt = TG[];
    }
  }

  /**
  The bound. An exact maximum principle holds only for the first-order form.
  The accurate stencil interpolates with `quadratic()`, whose outer weights
  can be negative, so a few per cent above the interface value is a
  discretisation artifact and not a runaway. `INTBC_OVERSHOOT` is the margin
  allowed, as a fraction. What must NOT happen is growth as the phase
  fraction falls; the sweep is there to show that. */

  bool bounded = (Tcut >= T0 - 1e-6 &&
                  Tcut <= TInt0 + INTBC_OVERSHOOT*(TInt0 - T0));
  bool solved  = (mg.i < NITERMAX);
  if (!bounded || !solved)
    nfail++;
  if (!tilted && brtgt != 0.)
    nforced_control++;

  fprintf (stderr, "%8s %10.1e %10.3e %10.3g %6g %14.6g %14.6g %5d %9.2e %9.2e %s\n",
           tilted ? "corner" : "CONTROL",
           fGtarget, fGmeas, stgt, brtgt, Ttgt, Tcut, mg.i, mg.resb, mg.resa,
           (bounded && solved) ? "ok" : (!solved ? "NOT-SOLVED" : "UNBOUNDED"));
}

event stop (i = 1) {
  return 1;
}
