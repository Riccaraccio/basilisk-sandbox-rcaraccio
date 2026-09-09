# Patches to the Basilisk install

`$BASILISK` is not a git repository on this machine, so a change to it leaves
no record. These two patches are that record. Apply them after any
reinstallation or update of Basilisk, or the temperature solve of
`INT_TEMP_VOFBC` loses its interface condition and falls back silently to the
frozen source.

    cd $BASILISK
    patch -p0 < .../src/basilisk-patches/poisson-flux-hook.patch
    patch -p0 < .../src/basilisk-patches/diffusion-flux-hook.patch

The originals are kept at `~/basilisk/basilisk-patched-backup/`.

Source: Alexandre's tree, `implicit/alexandre/basilisk/src`. The base version
matches the install here, so the patches apply cleanly.

## What they do

`poisson.h` calls a flux callback inside every relaxation and every residual
sweep. That callback is how `embed.h` imposes a Dirichlet condition on a cut
cell. It was compiled out unless `EMBED` was defined, and this sandbox must
not define `EMBED`, because `embed.h` owns the fields `cs` and `fs` and the
metric `cm` and `fm`.

The patch removes the four `#if EMBED` guards around the callback. The hook
then exists for any caller, and `poisson()` takes the callback in its existing
`flux` argument. `diffusion.h` gains the same argument and forwards it.

## Why they are safe

The behaviour of every existing case is unchanged.

- With no `flux` argument the pointer is `NULL` and `if (p->embed_flux)` is
  false, so the operator is exactly the one before the patch. The cost is one
  test of a null pointer per cell per sweep.
- With `EMBED` defined and no `flux`, `poisson()` still installs `embed_flux`,
  because the patch sets `p.embed_flux = flux` BEFORE the `#if EMBED` block
  and the block still overwrites it. Read that ordering carefully if you ever
  re-apply these by hand: swapping the two lines would break every embed case
  in the install, and nothing would report it.

## Why the interface condition cannot be built without them

The interface flux has to be evaluated on the CURRENT iterate of the solve,
not on the field of step `n`. `diffusion()` offers `r` and `beta`, and both
are constant during the solve, so neither can carry a term that depends on the
neighbours at step `n+1`. Only the callback can. Measured consequence of the
frozen form: `test/intbc-sliver.c` overshoots the interface value fourfold.

A copy of `poisson.h` inside `src/` was considered and rejected. A quoted
include resolves against the directory of the including file, so
`$BASILISK/viscosity.h` would still take the original while our `diffusion.h`
took the copy. Both would land in one translation unit and every symbol of
the file would be defined twice.

## A latent inconsistency in `poisson.h`, recorded but NOT patched

`relax` and `residual` do not solve the same equation when the flux callback
returns a non-zero diagonal coefficient `e`.

    relax     n -= c*Delta^2, d += e*Delta^2   =>  lambda*a + div - (c + e*a) = b
    residual  res[] += c - e*a[]               =>  lambda*a + div + e*a - c = b

They agree on `c` and disagree on the sign of `e*a`. `viscosity-embed.h` uses
the same callback and is self-consistent: its relax puts `-dt*c` in the
numerator and `+dt*d` in the denominator, and its residual does
`res.x[] -= dt*(c + d*u.x[])`.

This is invisible in normal use, because `embed_flux` returns a non-zero `e`
only in its degenerate branch, which is rare and local. It becomes fatal when
a callback returns `e` on many cells at once: the multigrid relaxes toward one
operator and measures the residual of another. A measured `run/restart.c`
stopped at the FIRST step.

The patches here do NOT change it. This sandbox no longer needs a diagonal
from the callback: `plic_flux` returns 0 and the small cells are excluded by
`PLICBC_TOL`, as `embed.h` and the original do. The note is kept because the
next person who returns a non-zero `e` will lose a day to it.

If you ever do need it, the fix is `res[] += c + e*a[];` in both branches of
`residual`, and it must be measured against every embed case in the install.
