/**
# The interface value as a boundary condition of the solve

This imposes a Dirichlet condition on the PLIC interface of a cut cell, from
inside the Poisson operator. It is a port of `plic_flux`
(`basilisk-sandbox-ecipriano/src/gradients.h`) onto the fields of this
sandbox, with one change, which the last section describes.

It needs the two patches in `basilisk-patches/`. Without them the callback
does not exist and nothing here compiles.

## Why the operator, and not a source term

`diffusion()` offers `r` and `beta`, and both are constant during the solve.
So neither can carry a term that depends on the neighbours at step `n+1`.
Freeze the interface flux into `r` and the cut-cell row reads

    [ theta/dt + F ] T_c = theta/dt T_c^n + C*(TInt - v0^n) + (face couplings)

where `v0` is an interpolated value about one cell away. Because `v0` is held
at step `n`, the source does not shut off as the cell heats: it keeps
delivering the flux the OLD field asked for. The cell then settles wherever it
must to pass that flux onward, which can be far above `TInt`, and no Dirichlet
condition should allow that. `test/intbc-sliver.c` measures the overshoot at a
factor of four, at every gas fraction from 0.4 down to 1e-6.

`poisson.h` calls the flux callback inside every relaxation and every residual
sweep, so `v0` is at the current iterate. The feedback then closes inside the
step: as the cell and its neighbours heat, `v0` rises toward `TInt` and the
flux goes to zero. That is the whole reason to patch Basilisk rather than to
use `beta`.

## Usage

A case attaches the condition to the `interface` boundary id:

    TS[interface] = dirichlet (TInt[]);

Then, before each solve, name the phase that the field lives in and pass the
callback:

    plicbc_phase (fS, fsS);
    diffusion (TS, dt, D = lambda1f, r = sST, theta = theta1, flux = plic_flux);

`plicbc_phase()` copies the fraction pair into the globals the callback reads.
The two phases use the same boundary id and differ only in that copy, so the
order matters: set the phase immediately before its own solve.

Caution: `face_fraction()` has no refinement procedure
(`$BASILISK/fractions.h:621`), and `relax` runs on every level. `cs1` gets a
prolongation and a restriction here. `fs1` gets neither, because none exists,
so on a tree the coarse levels see a face fraction that was restricted as a
plain field. That degrades the preconditioner, not the answer: the residual
is evaluated on the leaves. Watch `mgS.i` in `tsolve.dat`.

## The one change from the original ends with

    grad = (grad == nodata || coef != 0.) ? 0. : grad;
    ...
    return 0.;

so it discards the diagonal coefficient and, in a degenerate cell, drops the
interface flux entirely. That reproduces the behaviour of `intgrad.h`, where
the same cell silently receives no surface heat.

This version returns the coefficient, as `embed.h` does. Two consequences.

1. A degenerate cell keeps its flux, in the first-order form
   `(TInt - T_c)/(d0*Delta)`, which reads the cell and therefore lands on the
   diagonal. `poisson.h` adds it there with `d += e*sq(Delta)`.
2. That cell obeys a maximum principle: it lands between `TInt`, its old
   value and its neighbours, for any fraction and any `dt`.

`PLICBC_DROP_DEGENERATE` restores the original behaviour for a comparison. Do
not run physics with it: a cut cell that receives no surface heat is a hole in
the energy balance that nothing reports. */

#ifndef PLICBC_H
#define PLICBC_H

#include "intgrad.h"

/**
The boundary id that a case attaches the interface condition to. */

bid interface;

/**
The phase that the field being solved lives in. `plicbc_phase()` writes them
and the callback reads them. */

scalar cs1[];
face vector fs1[];

/**
`dirichlet` and `dirichlet_homogeneous` have to accept the `data` pointer that
`s.boundary[]` passes, so that the same expression serves a domain boundary
and this one. With `data` set the value is the boundary value itself and the
flag says the condition is Dirichlet; without it, the usual ghost formula.

Caution: the homogeneous form must report Dirichlet with a value of ZERO, and
not `-s[]`. `mg_cycle` relaxes the correction on every level, and the part
proportional to `s[]` is what the callback returns for the diagonal. Return it
twice and the diagonal is doubled. */

macro2
double dirichlet (double expr, Point point = point,
                  scalar s = _s, bool * data = data)
{
  return data ? *((bool *)data) = true, expr : 2.*expr - s[];
}

macro2
double dirichlet_homogeneous (double expr, Point point = point,
                              scalar s = _s, bool * data = data)
{
  return data ? *((bool *)data) = true, 0 : - s[];
}

/**
Name the phase of the next solve. */

void plicbc_phase (scalar c, face vector fc)
{
  foreach()
    cs1[] = c[];
  foreach_face()
    fs1.x[] = fc.x[];

#if TREE
  cs1.prolongation = fraction_refine;
  cs1.refine = fraction_refine;
#endif

  /**
  BOTH fields, and the face vector is not optional. `relax` runs on every
  level of the multigrid and reads `fs1` there, but `face_fraction()` fills
  the leaves only. Restrict just `cs1`, as the original does, and the coarse
  levels differentiate against a face fraction that was never written.

  Measured on `test/intbc-sliver.c`: without this line the solve stalls, at a
  relative residual of 5e-3 on a grid-aligned interface. `poisson()` restricts
  its own face vector `alpha` for the same reason. */

  restriction ({cs1, fs1});
}

/**
The smallest fraction that carries an interface condition. Below it the cell
is treated as pure and gets none.

Caution: this is a hole, not a bound. A cut cell below the threshold receives
no surface heat at all, which is what `TG_FGMIN_MODE 1` does. Keep it at the
value that only excludes cells the geometry cannot describe. */

#ifndef PLICBC_TOL
# define PLICBC_TOL 1.e-10
#endif

/**
The callback. `poisson.h` uses the two returned numbers as

    flux = *val + (return value)*s[]

putting the first in the numerator and the second on the diagonal. */

double plic_flux (Point point, scalar s, face vector D, double * val)
{
  *val = 0.;
  if (cs1[] < PLICBC_TOL || cs1[] > 1. - PLICBC_TOL)
    return 0.;

  bool dirichlet = false;
  double bc = s.boundary[interface] (point, point, s, &dirichlet);
  if (!bc && !dirichlet)
    return 0.;

  coord m = interface_normal (point, cs1), p;
  double alpha = plane_alpha (cs1[], m);
  double area = plane_area_center (m, alpha, &p);
#if AXI
  area *= (y + p.y*Delta);
#endif

  /**
  The distance from the centroid of this phase to the interface, for the
  first-order form. Compute it before `normalize`, because `plane_alpha`
  returns the alpha of the box normal. */

  double d0 = interface_phase_distance (point, cs1, m, alpha, false);

  normalize (&m);

  bool third = false;
#ifdef INTGRAD_3rd
  third = true;
#endif

  double coef = 0.;
  /**
  `force` is false: use the accurate stencil wherever it exists, exactly as
  the original does. The first-order form is taken only where the stencil is
  missing, and there it arrives through `coef`, not through a threshold. */

  double grad = concentration_gradient (point, s, cs1, fs1, m, p, bc, third,
                                        &coef, d0, false);

#if PLICBC_DROP_DEGENERATE
  if (grad == nodata || coef != 0.)
    grad = 0., coef = 0.;
#endif

  /**
  Recover the conductivity without the face fractions and without the metric.
  `D` carries both, exactly as `embed_flux` assumes. */

  double Da = 0., fa = 0.;
  foreach_dimension() {
    Da += D.x[] + D.x[1];
    fa += fm.x[]*fs1.x[] + fm.x[1]*fs1.x[1];
  }
  double mua = Da/(fa + 1.e-30);

  *val = - mua*grad*area/Delta;
  return - mua*coef*area/Delta;
}

#endif // PLICBC_H
