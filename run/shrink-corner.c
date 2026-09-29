/**
# Smear of the pellet corner under shrinkage: an isolated test

At level 11, `fatehi-combustion.c` shows solid in gas cells near the corner of
the pellet by t = 4 s, with f from 0.1 to 0.6. At that time the pellet lost
only 0.4 % of its mass, so the physical shrinkage is much less than one cell.
Level 10 shows almost no smear. Shrinkage only pulls the interface inward, so
this solid is a numerical defect.

In that case only two operations change `f`: the VOF advection with the solid
velocity `ubf`, and the adaptation of the grid. This case keeps both, with the
geometry, the grid, the `psi` boundary conditions and the adaptation call of
`fatehi-combustion.c`. It removes everything else:

- no chemistry: `omega` is a fixed function of a model temperature;
- no gas expansion: `rhoG = rhoS`, so `gas_source` is 0 and `u` stays 0;
- no temperature, species or radiation.

## The model source

The model temperature `theta` is 1 at the gas and 0 deep in the solid:

  theta = 1 - g(dx) g(dy),   g(d) = G0 + (1 - G0) erf (d/DELTA)

where `dx` and `dy` are the depths below the flat faces. On a face, theta is
1 - G0. At the corner, theta is 1 - G0^2, because two faces heat it. The rate
follows an Arrhenius-like law, `shape = exp (SHARP (theta - 1))`, so the corner
reacts `exp (SHARP G0 (1 - G0))` times faster than a face: 7.4 times at the
default values.

The case scales `omega` at each step so that the solid loses mass at the
constant rate `MDOT` (1/s):

  sum (omega (f - porosity) dv) = MDOT m0

The default 1e-3 1/s is the rate of `fatehi-11` near t = 5 s, so t = 4 s gives
a mass loss of 0.4 %, as in the production run.

## The switches

| flag | default | effect |
|---|---|---|
| `MAXLEVEL` | 11 | the finest level of the production domain |
| `ZOOM` | 2 | the domain is `20 D0/2^ZOOM`, and each level is `ZOOM` less |
| `SHIFT_PROD` | 1 | `shift_prod` of `velocity-potential.h` |
| `ZETA_POLICY` | `ZETA_REACTION` | the shrinkage policy |
| `OMEGA_SHAPE` | 1 | 1: the corner model above; 0: uniform `omega` |
| `ADAPT` | 1 | 0: the grid stays as `init` builds it; set `INIT_FIX = 1` too |
| `INIT_FIX` | 0 | 0: the level 8 start of production; 1: refine before `fraction()` |
| `MDOT` | 1e-3 | 0: no shrinkage, so only the adaptation acts on `f` |
| `TEND` | 4 | the end time |
| `DT_VALUE` | 5e-4 at level 10 | the fixed step; it halves at each level |

`DT_VALUE` matches the step of the production runs: 5e-4 at level 10, and
about 2.5e-4 at level 11. The flow CFL set it there. Here `u` is 0, so the
value is the step.

Caution: the diagnostics below read `f` of the neighbours and locate points
on one process. Run this case serial.

## The output

`debris.dat` holds one line each 0.1 s. `corner-<t>.txt` holds a map of `f`
near the corner each second: `#` is f > 0.99, `.` is f < 0.01, and a digit is
10 f. Read the columns in the header of `debris.dat`. */

#define NO_ADVECTION_DIV 1

#ifndef MAXLEVEL
# define MAXLEVEL 11
#endif

/**
`ZOOM` makes the domain `2^ZOOM` times smaller than the domain of
`fatehi-combustion.c`, and removes `ZOOM` levels. The cell size of each level
stays the same, so `MAXLEVEL` keeps its production meaning: `MAXLEVEL = 11`
gives the cell of 78 um at every zoom. This case solves no temperature, so the
far field is not necessary.

Caution: `psi` is 0 on the left and top boundaries. At `ZOOM = 2` the left
boundary is 16 mm from the pellet face, which is 4 half-widths. A larger zoom
puts that boundary near the pellet and changes `ubf` there. Check the result
against `ZOOM = 1` before you use a larger zoom. */

#ifndef ZOOM
# define ZOOM 2
#endif

#ifndef SHIFT_PROD
# define SHIFT_PROD 1
#endif

#ifndef ZETA_POLICY
# define ZETA_POLICY ZETA_REACTION
#endif

#ifndef OMEGA_SHAPE
# define OMEGA_SHAPE 1
#endif

#ifndef ADAPT
# define ADAPT 1
#endif

#ifndef MDOT
# define MDOT 1e-3
#endif

#ifndef INIT_FIX
# define INIT_FIX 0
#endif

#ifndef TEND
# define TEND 4.
#endif

#ifndef DT_VALUE
# define DT_VALUE (5e-4*pow (2., 10 - MAXLEVEL))
#endif

/**
The model temperature. `DELTA` is the depth of the heated layer, `G0` sets
the temperature of a face, and `SHARP` sets how much the corner dominates. */

#ifndef DELTA
# define DELTA 5e-4
#endif

#ifndef G0
# define G0 0.5
#endif

#ifndef SHARP
# define SHARP 8.
#endif

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke.h"
#include "constant-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "superquadric.h"

/**
The boundary conditions of `fatehi-combustion.c`, without the inflow. The
`psi` lines are the same as in that case. */

u.n[left]  = dirichlet (0.);
u.t[left]  = dirichlet (0.);
p[left]    = neumann (0.);
psi[left]  = dirichlet (0.);

psi[top]   = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right]   = dirichlet (0.);
psi[right] = neumann (0.);

int maxlevel = MAXLEVEL - ZOOM, minlevel = max (2 - ZOOM, 1);
double D0 = 8e-3, H0 = 8e-3;

scalar omega[];
scalar fS[];
face vector fsS[];

/**
`shape[]` is the normalized rate. The adaptation reads it in the place of the
temperature of the production case, so the heated layer stays refined. */

scalar shape[];

double m0 = 0.;          // the initial solid mass
double V0 = 0.;          // the initial solid volume
double Vshrink = 0.;     // the volume that the source removed, sum of dt Q

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

int main()
{
  if (pid() == 0)
    fprintf (stderr, "# shrink-corner: MAXLEVEL=%d zoom=%d level=%d L0=%g"
             " shift=%d zeta=%d shape=%d"
             " adapt=%d MDOT=%g DT=%g TEND=%g DELTA=%g G0=%g SHARP=%g\n",
             MAXLEVEL, ZOOM, MAXLEVEL - ZOOM, 20*D0/(1 << ZOOM),
             SHIFT_PROD, (int) ZETA_POLICY, OMEGA_SHAPE, ADAPT,
             (double) MDOT, (double) DT_VALUE, (double) TEND,
             (double) DELTA, (double) G0, (double) SHARP);

  eps0 = 0.2;
  rhoS = 1550.;
  rhoG = rhoS;            // gas_source = 0
  muG = 1e-5;
  rho1 = rho2 = 1.;
  mu1 = mu2 = 1.;

  zeta_policy = ZETA_POLICY;
  shift_prod = SHIFT_PROD;

  DT = DT_VALUE;

  L0 = 20*D0/(1 << ZOOM);
  origin (-L0/2, 0.);

  /**
  The same grid as `fatehi-combustion.c`.

  Caution: this `refine()` does nothing. `run()` calls `init_grid (N)` again
  (`$BASILISK/run.h:17`), and `init_grid` of the tree frees the grid, so the
  run starts on a uniform grid at level 8. `fatehi-combustion.c` has the same
  lines, so the case keeps them with `INIT_FIX = 0` to copy that state.
  Level `8 - ZOOM` of this domain has the cell size of level 8 there. With
  `INIT_FIX = 1`, `event init` refines before `fraction()`. */

  init_grid (1 << min (maxlevel, 8 - ZOOM));
#if !INIT_FIX
  refine (circle (x, y, 4.*D0) > 0. && level < maxlevel);
#endif

  TOLERANCE = 1e-5;
  NITERMIN = 2;

  run();
}

event init (i = 0)
{
  /**
  A disc of radius 0.75 D0 holds the pellet: the corner is at 0.71 D0. At
  level 11 the disc holds about 9 000 cells. */

#if INIT_FIX
  refine (circle (x, y, 0.75*D0) > 0. && level < maxlevel);
#endif

  fraction (f, superquadric (x, y, 20, 0.5*H0, 0.5*D0));

  foreach()
    porosity[] = eps0*f[];

  m0 = V0 = 0.;
  foreach (reduction(+:m0) reduction(+:V0)) {
    m0 += (f[] - porosity[])*rhoS*dv();
    V0 += f[]*dv();
  }
}

/**
## The model source

The depths `dx` and `dy` are negative outside the pellet. The clamp at 0 gives
a cell outside the pellet the rate of the nearest face or corner, as a hot
fragment of solid in the gas would have. */

static double model_shape (double x, double y)
{
#if OMEGA_SHAPE
  double dx = max (0.5*H0 - fabs (x), 0.);
  double dy = max (0.5*D0 - y, 0.);
  double gx = G0 + (1. - G0)*erf (dx/DELTA);
  double gy = G0 + (1. - G0)*erf (dy/DELTA);
  double theta = 1. - gx*gy;
  return exp (SHARP*(theta - 1.));
#else
  return 1.;
#endif
}

/**
This event takes the place of the chemistry. It also adds the volume that the
previous step removed: `prod[]` of `velocity-potential.h` still holds the
shifted source of that step, and the Poisson equation gives
`div (ubf) = -prod`, so the rate of volume change is `-sum (prod Delta^2)`.

The porosity update is the one of `run/temp.c`: the part `1 - zeta` of the
reaction opens pores and does not shrink the pellet. */

event chemistry (i++)
{
  if (i > 0) {
    double Q = 0.;
    foreach (reduction(+:Q))
      Q += prod[]*sq(Delta);
    Vshrink += dt*Q;
  }

  /**
  `shape[]` is 0 in the gas. Otherwise the adaptation refines the gas along
  the lines of the faces, where the clamped depths change the value. */

  double norm = 0.;
  foreach (reduction(+:norm)) {
    shape[] = f[] > F_ERR ? model_shape (x, y) : 0.;
    norm += shape[]*(f[] - porosity[])*dv();
  }

  foreach() {
    omega[] = (f[] > F_ERR && norm > 0.) ? MDOT*m0*shape[]/norm : 0.;
    if (f[] > F_ERR) {
      porosity[] = porosity[]/f[];
      porosity[] += omega[]*(1. - porosity[])*(1. - zeta[])/rhoS*dt;
      porosity[] *= f[];
    }
  }
}

#if TREE && ADAPT
event adapt (i++)
{
  adapt_wavelet_leave_interface ({shape, porosity}, {f},
    (double[]){1e-2, 1e-2}, maxlevel, minlevel, 2);
}
#endif

/**
## The diagnostics

`err` is the L1 distance between `f` and the exact initial shape on the
present grid, `sum (|f - fex| dv)/V0`. In a case without shrinkage it measures
the defect directly. With shrinkage it also holds the real change of shape,
which is at most `Vshrink/V0`.

Caution: `Vout` misses the defect of the level 8 start. That start moves both
faces about one cell inward, so the fragments stay inside the initial square.
Read `err` and the maps, not `Vout`, for that defect.

A cell is *outside* when its centre is more than one cell outside the initial
square, `max (|x| - H0/2, y - D0/2) > Delta`. Shrinkage cannot put solid
there, so every bit of `f` in such a cell is a defect. The superquadric of
order 20 lies inside the square (0.14 mm inside it on the diagonal), so the
test does not count the solid of the initial corner.

A mixed cell is *detached* when f > 0.01 and no cell of its 3x3 block holds
f >= 0.5. `shift_field()` can lose source, so the case compares
`sum (prod Delta^2)` after the shift with the same sum before the shift. */

event logfile (t = 0; t += 0.1)
{
  double V = 0., Vout = 0., fmax_out = 0., Qpost = 0., Qpre = 0.;
  double prodmax = 0., areasrc = 0.;
  int nout = 0, nmixed = 0, ndetached = 0, ncells = 0;

  foreach (reduction(+:V) reduction(+:Vout) reduction(max:fmax_out)
           reduction(+:nout) reduction(+:nmixed) reduction(+:ndetached)
           reduction(+:Qpost) reduction(+:Qpre) reduction(max:prodmax)
           reduction(+:areasrc) reduction(+:ncells)) {
    ncells++;
    V += f[]*dv();

    double out = max (fabs (x) - 0.5*H0, y - 0.5*D0);
    if (out > Delta && f[] > F_ERR) {
      Vout += f[]*dv();
      if (f[] > 0.01)
        nout++;
      if (f[] > fmax_out)
        fmax_out = f[];
    }

    if (f[] > 0.01 && f[] < 1. - F_ERR) {
      nmixed++;
      bool anchored = false;
      foreach_neighbor (1)
        if (f[] >= 0.5)
          anchored = true;
      if (!anchored)
        ndetached++;
    }

    Qpost += prod[]*sq(Delta);
    Qpre += omega[]*f[]*zeta[]*cm[]/rhoS*sq(Delta);
    if (prod[] > prodmax)
      prodmax = prod[];
    if (prod[] > 0.)
      areasrc += sq(Delta);
  }

  /**
  `lost` is the part of the source that the shift removed. `peak` is the
  largest source against the mean source over the cells that hold a source:
  a large value marks a point sink. */

  scalar fex[];
  fraction (fex, superquadric (x, y, 20, 0.5*H0, 0.5*D0));
  double err = 0.;
  foreach (reduction(+:err))
    err += fabs (f[] - fex[])*dv();

  double lost = Qpre > 0. ? 1. - Qpost/Qpre : 0.;
  double peak = Qpost > 0. ? prodmax*areasrc/Qpost : 0.;

  static FILE * fp = NULL;
  if (fp == NULL) {
    fp = fopen ("debris.dat", "w");
    fprintf (fp, "#t(1) i(2) ncells(3) V/V0-1(4) Vshrink/V0(5) Vout/V0(6)"
             " nout(7) fmax_out(8) nmixed(9) ndetached(10) lost(11)"
             " peak(12) err(13)\n");
  }
  fprintf (fp, "%g %d %d %.6e %.6e %.6e %d %.4f %d %d %.4f %.1f %.6e\n",
           t, i, ncells, V/V0 - 1., Vshrink/V0, Vout/V0, nout, fmax_out,
           nmixed, ndetached, lost, peak, err/V0);
  fflush (fp);
}

/**
The map of `f` near the corner, on a lattice of the finest cell size. Each
character is the leaf cell that holds the lattice point. The window covers
x and y from 1.5 mm to 5.8 mm; the corner is at (4 mm, 4 mm). */

static void corner_map (const char * name)
{
  FILE * fp = fopen (name, "w");

  double h = L0/(1 << maxlevel), lo = 1.5e-3, hi = 5.8e-3;
  int n = (hi - lo)/h;
  for (int jy = n - 1; jy >= 0; jy--) {
    for (int jx = 0; jx < n; jx++) {
      double xp = lo + (jx + 0.5)*h, yp = lo + (jy + 0.5)*h;
      Point point = locate (xp, yp);
      char c = '?';
      if (point.level >= 0) {
        double v = f[];
        c = v > 0.99 ? '#' : v < 0.01 ? '.' : '0' + (int) (10.*v);
      }
      fputc (c, fp);
    }
    fputc ('\n', fp);
  }
  fclose (fp);
}

event cornermap (t = 0; t += 1.)
{
  char name[80];
  sprintf (name, "corner-%g.txt", t);
  corner_map (name);
}

/**
The map after the first step shows what the first adaptation did to the
corner. With `INIT_FIX = 0` that adaptation builds the finest levels from the
level 8 pellet. */

event cornermap1 (i = 1)
{
  corner_map ("corner-i1.txt");
}

/**
Caution: `return 1` ends the run. The `logfile` event has no end time, so a
bare `event stop (t = TEND)` does not stop the loop. */

event stop (t = TEND)
{
  return 1;
}
