/**
# The interface temperature on the diagonal

`INT_TEMP_VOFBC` makes `TInt` a Dirichlet condition of each temperature solve,
instead of a frozen source. The code is in `multicomponent-varprop.h`, in
`interface_temperature_sources()`. This file holds the design note, the
counters and the event that writes them. It changes no field.

## What the flag changes

`ebmgrad` is affine in the interface value, so the interface gradient splits
exactly into a part that does not involve this cell and a part that does:

    grad = G + coef*T[]

Without the flag the whole product `lambda*aov*grad` goes to the source `r`.
With it, the first part goes to `r` and the second to the diagonal, through
the `beta` argument of `diffusion()`:

    sST[]    += lambda1vh*aov*G
    betaST[] += lambda1vh*aov*coef

The split is exact. It is NOT a deferred correction: nothing is added and
subtracted, so it lags nothing and it withholds no heat. This is the
difference from `INT_TEMP_ROBIN`, which added an artificial conductance `K`
and had to carry the withheld heat in a debt field to stay conservative.

## Why the split alone is not enough

`coef` is zero whenever the accurate stencil is available, because that
stencil reads the NEIGHBOURS of the cut cell along the normal and never the
cell itself. So for most cut cells the flag changes nothing, bit for bit.

That is also the defect, and the reason is not the one this branch started
from. Write the gas row of a cut cell:

    [ theta2/dt + F ] T_c = theta2/dt T_c^n + C*(TInt - v0) + (face couplings)

`v0` is an interpolated value about one cell away, and `T_c` is not in the
interface term at all. The row is a well posed statement about the
neighbours, and it says nothing about how far `T_c` may go. The cell settles
wherever it must to pass the imposed flux onward, and that can be far above
`TInt`, which no Dirichlet condition should allow.

Caution: the row does NOT become singular, and the solve does NOT fail.
`test/intbc-sliver.c` measured the interface conductance over the diagonal of
the row at 0.33, unchanged from a gas fraction of 0.4 down to 1e-6, while the
answer overshot the interface value fourfold at every one of them. That ratio
is bounded for a reason: for a convex cut the interface area is at most the
sum of the face fractions of the phase, so the interface conductance is at
most the face conductance. Do not build a criterion on it.

So a Dirichlet condition imposed at the facet does not, on its own, bind a
cut cell. `embed.h` has the same hole and plugs it twice: the `v[0] == nodata`
branch of `dirichlet_gradient_x`, which puts the cell back in its row, and
`fractions_cleanup()`, which deletes such cells by setting `cs = 0`. The
second is closed to us, because `f` is a physical volume fraction and not a
mask.

## The criterion

Take the first branch when the accurate stencil cannot describe this phase.
Per cut cell and per phase:

    d0 = the thickness of the phase, centroid to interface, in units of Delta
    d1 = the reach of the accurate stencil, 1/(h*Delta), about one cell

`d0` comes from `plane_center`, and `d1` from the exact affine slope of
`ebmgrad`, so both are free of any extra geometry. When `d0/d1` is below
`INT_TEMP_ETA` the stencil reads a value from well outside the layer it
claims to differentiate, and the cell takes the first-order form that reads
itself:

    grad = (TInt - T[])/(d0*Delta)

The row is then

    (theta/dt + F + C) T_c = theta/dt T_c^n + C*TInt + (face couplings)

an M-matrix row: `T_c` lands between `TInt`, `T_c^n` and its neighbours, for
any phase fraction and any `dt`. As the capacity and the face conductance
vanish it gives `T_c = TInt` exactly, which is what the Dirichlet condition
promised. `test/intbc-sliver.c` measures 999.16, 999.77, 999.85, 999.86 K
against an interface at 1000 K, as the gas fraction falls from 1e-2 to 1e-6.

## Why `INT_TEMP_ETA` is not `INT_TEMP_ROBIN_SMAX`

`SMAX` scaled a physical flux: it put `A/SMAX` on the diagonal, which is less
than the conductance the geometry asks for, and the cell then kept only
`SMAX/S` of the heat the interface gave it. A measured run lost 20 to 25 K of
surface temperature and 25 per cent of the mass loss rate, and no bookkeeping
of the withheld energy restored a response the cell was not allowed to make.
Refining the grid did not retire it.

`INT_TEMP_ETA` selects between two consistent discretisations of the same
flux. Neither is scaled, and no heat is withheld. A cell above the threshold
keeps the accurate stencil, bit for bit; a cell below it takes a first-order
form of the same Dirichlet condition. Refinement moves cells out of the
second branch, so the influence of the number goes away with `Delta`.

Caution: that argument bounds the ERROR of the threshold, not its effect on
a given run. The table below shows the effect is a cliff, not a slope, so
the number still has to be measured on a case.

## Measured: there is a plateau, and there is a cliff

`run/restart.c`, 200 steps at level 7, against the build without the flag, at
matched time:

| eta | forced | mass | radius | timestep |
|---|---|---|---|---|
| 0 | 0 of 26 | 0.0000 % | 0.0000 % | 1.00 |
| 0.02 | 8 of 26 | -0.096 % | -0.016 % | 0.94 |
| 0.05 | 8 of 26 | -0.096 % | -0.016 % | 0.94 |
| 0.10 | 16 of 26 | -2.06 % | -0.400 % | 0.36 |
| 1 | 26 of 26 | -3.99 % | -0.726 % | 0.36 |

`eta = 0` reproduces the build without the flag exactly, on every column and
on the timestep. That is not a coincidence: with no cell forced, `coef` is
zero, `beta` stays zero and `r` is the same expression as before.

0.02 and 0.05 force the same eight cells, so the default sits on a plateau.
The cliff is at 0.10, where the forced set doubles, the mass loss moves by
two per cent and the timestep falls to a third. Both of those are the
signature of a scheme that has begun to rewrite the answer rather than to
bound it, and it is the same signature `INT_TEMP_ROBIN_SMAX = 1` gave.

The largest `d0/d1` of that case is 0.373. So a real surface holds cut cells
far better resolved than any geometry in `test/intbc-sliver.c`, whose ratios
all sit below 0.19. Take the threshold from a case, never from that test.

`INT_TEMP_VOFBC_CMAX` caps `C` at that multiple of `Ddiag`, default 100. The
maximum principle holds for any positive `C`, so the cap costs nothing in
boundedness; it stops a very thin cell from wrecking the condition number of
the multigrid. `G` and `coef` are scaled together, so the capped pair is
still the flux `C*(TInt - T[])` and the interface value it relaxes to is
unchanged.

## Picard is not optional here

Without the flag, `sST` and `sGT` are both built from one `TInt` and from the
same step-`n` gradients, so the heat that leaves the solid equals the heat
that enters the gas, cell by cell, to the tolerance of the root find. That
scheme is exactly conservative and merely lagged.

With `TInt` frozen and the two phases solved apart, that identity is gone.
After the solves the imbalance of a cut cell is

    aov*[ -q_rad(TInt^n) + lambda1vh*gS(T^{n+1}) + lambda2vh*gG(T^{n+1}) ]

which is proportional to the change of the near-interface temperatures over
one step, so it is first order in `dt`. Only the Picard fixed point removes
it. `INT_TEMP_PICARD` is therefore required, and `INT_TEMP_PICARD_MAXITER = 0`
is an error under this flag.

The outer map now contracts, which it did not before. Each phase solve obeys
a discrete maximum principle with respect to `TInt`, so the gain of one pass
is at most

    (lambda1vh*hS + lambda2vh*hG)/D  =  1 - 4*sigma*eps*TInt^3/D

with `D = dF/dTInt`. The radiation term alone keeps it below one. Without the
flag the same map has gain of order `S`, the exchange number, which a
measured run put at 100.

Caution: the bound needs the M-matrix property, which is a theorem in the
first-order branch and a guideline in the third-order and vof-averaged
branches, whose `quadratic()` weights can be negative.

## What the flag does NOT change

`ijc_CoupledTemperature()`, `EqTemperature` and the radiation term
`sigma*eps*(Tbulk^4 - TInt^4)` are untouched, and the radiation is not
linearised. `TInt` is frozen during each linear solve, so the operator stays
linear whatever produced it.

`TG_FGMIN` is left alone. Its mode 2 puts the full conductance on the
diagonal, but that line lives inside `#if INT_TEMP_ROBIN`, so with this flag
and no Robin it is inert. The criterion above covers the same cells and does
not need a fraction threshold.

## The open question

Under the criterion the operator applies `C*(TInt - T_c)` while
`EqTemperature` still measures the balance with the accurate stencil. The two
are then different functions of `TInt`, and the Picard fixed point satisfies
neither exactly. `INT_TEMP_VOFBC_CONSISTENT`, off by default, gives
`EqTemperature` the same criterion and the same form. Read `rel_max` of
`intres.dat` on a build with it and a build without it before you choose.

## How to read `vofbc.dat`

    #t(1) dt(2) nint(3) nforce(4) ndegen(5) ncap(6) Smax(7) betamin(8)

- `nforce` the cut cells that took the first-order form because the criterion
  fired. When this is zero the flag changes nothing at all, and the run must
  reproduce a build without it to every printed digit.
- `ndegen` the cut cells whose accurate stencil was missing. These take the
  first-order form whatever the criterion says. A build without the flag
  DROPS the interface flux of these cells and reports nothing, so a non-zero
  count here also measures a defect of the old path.
- `ncap` the cells whose `C` hit `INT_TEMP_VOFBC_CMAX`.
- `Smax` the largest `C_acc/Ddiag` over the cut cells, BEFORE the criterion
  is applied. This is the quantity that used to decide whether the run
  survived. It may stay large, and that is the point: the scheme must be
  indifferent to it.
- `betamin` the most negative diagonal contribution, in the units of the
  source divided by a temperature.

Caution: this file shows which branch each cell took. It does not show that
the answer is right. The proof is `rel_max` of `intres.dat`, and `Tsurf` and
`mdot` against a build without the flag.
*/

double ITV_nint, ITV_nforce, ITV_ndegen, ITV_ncap, ITV_Smax, ITV_betamin;

event vofbc_output (i++, last) {
  static FILE * fp = NULL;
  if (!fp) {
    fp = fopen ("vofbc.dat", restarted ? "a" : "w");
    if (!fp) {
      fprintf (stderr, "Error opening vofbc.dat\n");
      return 0;
    }
    if (!restarted)
      fprintf (fp, "#t(1) dt(2) nint(3) nforce(4) ndegen(5) ncap(6)"
                   " Smax(7) betamin(8)\n");
  }

  fprintf (fp, "%g %g %g %g %g %g %g %g\n",
           t, dt, ITV_nint, ITV_nforce, ITV_ndegen, ITV_ncap,
           ITV_Smax, ITV_betamin);
  fflush (fp);

  return 0;
}
