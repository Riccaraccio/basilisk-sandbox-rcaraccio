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

This port keeps that behaviour. It was changed once, to return the
coefficient as `embed.h` does, so that a degenerate cell kept its flux on the
diagonal. That was reverted, for two measured reasons.

1. It is not needed. `PLICBC_TOL` excludes the cells that misbehave, and the
   heat it costs is 0.1 per cent of `QS`. See the note at `PLICBC_TOL`.
2. It does not work. `relax` and `residual` of `poisson.h` disagree on the
   sign of `e*a`, which is invisible while `e` is non-zero only in the rare
   degenerate branch and fatal when a branch returns it on many cells at once.
   `run/restart.c` stopped at the FIRST step. See
   `basilisk-patches/README.md`. */

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
The smallest phase fraction that carries an interface condition. Below it the
cell is treated as pure and gets none. This is the small-cell treatment, and
it is the same one `embed.h` and the original use: EXCLUDE the cell, do not
try to bound it.

The value 1e-3 is Alexandre's. Measured on `test-vofbcm`, MOISTURE, level 10,
8 ranks, over t = 5.0 to 5.97 against a build with 1e-10:

| PLICBC_TOL | QS | QG | outcome |
|---|---|---|---|
| 1e-10 | 0.008521 | -0.001395 | stops at t = 5.97, TG_min -2.3e6 |
| 1e-3 | 0.008530 | -0.001404 | passes t = 6.12, TG_min 525 |

So the exclusion changes the interfacial heat by +0.1 per cent on the solid
side and +0.6 per cent on the gas side, and it removes the failure. The cells
it drops hold under 0.1 per cent of the interfacial gas volume: they carried
the instability, not the heat.

In that run `SG_max` reached 287 with 93 cells above an exchange number of 1,
and nothing happened. That is the point of the whole scheme: with the
interface condition inside the operator, `S` stops being a stability number.

Caution: this IS an energy hole, in the way `TG_FGMIN_MODE 1` is. It is small
here because the excluded volume is small. Check `QS` against a run with a
lower value before you raise it on a new case. */

#ifndef PLICBC_TOL
# define PLICBC_TOL 1.e-3
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

  /**
  The interface value, from the condition the case attached to the `interface`
  boundary id.

  Caution: a field with no such condition keeps the default of a new `bid`,
  which is `symmetry`. That returns `s[]` and leaves the flag false, so
  without this test the flux is built with the CELL's own value as the
  interface value, the gradient collapses to an internal one, and the
  interface exchange silently goes to nothing. It reads as a converged,
  well behaved run with no heat transfer. Fail here instead. */

  bool dirichlet = false;
  double bc = s.boundary[interface] (point, point, s, &dirichlet);
  if (!dirichlet) {
    fprintf (stderr, "plic_flux: no Dirichlet condition on the interface for "
             "this field. Set s[interface] = dirichlet(...) where the field "
             "exists, not in a defaults event.\n");
    assert (dirichlet);
  }

  coord m = interface_normal (point, cs1), p;
  double alpha = plane_alpha (cs1[], m);
  double area = plane_area_center (m, alpha, &p);
#if AXI
  area *= (y + p.y*Delta);
#endif

  normalize (&m);

  bool third = false;
#ifdef INTGRAD_3rd
  third = true;
#endif

  double coef = 0.;
  double grad = concentration_gradient (point, s, cs1, fs1, m, p, bc, third,
                                        &coef, 0., false);


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

  /**
  No diagonal, ever. `concentration_gradient` is called with `force = false`,
  so `coef` is zero in every branch: the accurate estimates do not read the
  cell, and the degenerate branch reports no flux. The assertion states the
  invariant rather than trusting it.

  Caution: do NOT return `-mua*coef*area/Delta` here, even though it looks
  more general. `relax` and `residual` of `poisson.h` disagree on the sign of
  that term, so a non-zero value makes the multigrid relax toward one operator
  and measure the residual of another. It is invisible until a degenerate cell
  appears, and then the solve diverges in ONE step. A measured `test-vofbcm`
  ran to t = 10.03 with the solid solve converging in 2 cycles, then reported
  a residual of 1.5e35 and stopped inside the GSL Jacobian of
  `EqTemperature`. Read `basilisk-patches/README.md` before you change this. */

  assert (coef == 0.);
  return 0.;
}

#endif // PLICBC_H
