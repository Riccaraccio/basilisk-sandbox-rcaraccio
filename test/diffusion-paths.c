/**
# The two diffusion paths on a grid ladder, gas only

`run/fatehi-combustion.c` at level 10 gives a different answer with
`COUPLED_DIFFUSION` 0 and 1: Tmax differs by 32 K, and the gap is not a time
error. This test checks if the two paths of `multicomponent-varprop.h`
differ in the gas alone, with no solid and no interface.

The case is 1-D. The grid is a 2-D channel 2 times longer than its height,
with the default symmetry conditions at the top and the bottom, thus the
solution does not change along y. The stack does not compile on
`multigrid1D.h`. With a channel 4 or more times longer than its height, the
solve of `pf` diverges at the first step (not investigated). A mixture of `H2`, `TAR` (MW 162) and `CO2` in `N2` fills
the left half, and `N2` fills the right half. The front has the width
`W`. The temperature is uniform, and the reactions are off. The left end is
a wall, and the right end is an outlet. The stack and the transport flags
are the ones of fatehi: `MOLAR_DIFFUSION`, `FICK_CORRECTED`,
`MASS_DIFFUSION_ENTHALPY` and the gas expansion.

The step `DT_VALUE` is fixed and the same at each level. It must stay below
the limit of the corrective flux of the old path, `corrective_dtmax`. The run
writes `dtmin` and `corrective_dtmax` at the end, so that a reader can check
this.

At `tend` the run writes `prof`: x, `YG` of each species, `TG` and `u.x`.
`misc/diffusion-paths.py` compares the profiles of the levels and of the two
paths.

Build flags: `-DLEVEL=<n>`, `-DCOUPLED_DIFFUSION=0|1`, `-DDT_VALUE=<s>`,
`-DTEND_VALUE=<s>`.

## The result (2026-11-02)

Levels 5 to 10, `dt = 1e-5`, `tend = 2e-3` (211 steps). The limit of the
corrective flux was at least 1.3e-4, thus it did not bind. No solve gave a
warning. Level 7 has the cell size of fatehi at level 10 (1.5625e-4 m).

The mean of |error| against level 10 of the same path, for `TAR` (0.45 at
the start):

| level | old path | coupled path |
|---|---|---|
| 7 | 9.0e-5 | 9.1e-5 |
| 8 | 2.2e-5 | 2.7e-5 |
| 9 | 4.4e-6 | 8.5e-6 |

* The old path converges with second order (1.9 to 2.3 for every species).
  The coupled path gives 1.7 to 1.9 for `TAR`, `CO2` and `N2`: the upwind
  face value of the corrective flux adds a first-order part.
* The difference of the two paths for `TAR`, `CO2` and `N2` halves at each
  level and does not change when `dt` is halved. It is the first-order part
  above, and the two paths go to the same limit.
* The difference for `H2` and `u` stays near 2.7e-6 and 2.2e-5 above level
  8, and it falls when `dt` is halved: it is a difference of the time
  discretizations.
* At level 7, the largest difference of the two paths is 5e-4 for `TAR`,
  2.6e-5 for `H2` and 1.6e-4 m/s for `u` (the largest `u` is 0.09 m/s).

Thus in the gas alone the two paths agree to about 0.1 % at the cell size of
fatehi. This does not explain the 32 K of fatehi. Both paths are accurate
here. */

#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1
#define TURN_OFF_REACTIONS 1

#ifndef LEVEL
# define LEVEL 7
#endif
#ifndef DT_VALUE
# define DT_VALUE 1e-5
#endif
#ifndef TEND_VALUE
# define TEND_VALUE 2e-3
#endif

#include "grid/multigrid.h"
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"

const double tend = TEND_VALUE;
const double W = 1e-3;      // width of the front [m]

u.n[left]  = dirichlet (0.);
p[left]    = neumann (0.);
pf[left]   = neumann (0.);
psi[left]  = dirichlet (0.);

u.n[right] = neumann (0.);
p[right]   = dirichlet (0.);
pf[right]  = dirichlet (0.);
psi[right] = dirichlet (0.);

int main() {
  TS0 = 1123.; TG0 = 1123.;
  rhoS = 1550.;
  eps0 = 0.2;

  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_REACTION;
  kinfolder = "biomass/dummy-solid-gas";

  L0 = 0.02;
  DT = DT_VALUE;
  TOLERANCE = 1e-6;
  NITERMIN = 2;

  dimensions (nx = 2, ny = 1);
  init_grid (1 << LEVEL);
  run();
}

event init (i = 0) {
  restarted = true;

  int iN2 = OpenSMOKE_IndexOfSpecies ("N2");
  double yA[NGS];
  for (int jj = 0; jj < NGS; jj++)
    yA[jj] = 0.;
  yA[OpenSMOKE_IndexOfSpecies ("H2")]  = 0.05;
  yA[OpenSMOKE_IndexOfSpecies ("TAR")] = 0.45;
  yA[OpenSMOKE_IndexOfSpecies ("CO2")] = 0.20;
  yA[iN2] = 0.30;

  foreach() {
    f[] = 0., porosity[] = 0.;
    double s = 0.5*(1. - tanh ((x - L0/2.)/W)); // 1 on the left
    for (int jj = 0; jj < NGS; jj++) {
      scalar YS = YGList_S[jj], YG = YGList_G[jj], YI = YGList_Int[jj];
      YS[] = 0.;
      YG[] = s*yA[jj] + (1. - s)*(jj == iN2);
      YI[] = (jj == iN2);
    }
    for (int jj = 0; jj < NSS; jj++) {
      scalar YSol = YSList[jj];
      YSol[] = 0.;
    }
    TS[] = 0.;
    TG[] = TG0;
    T[]  = TG0;
    TInt[] = 0.;
    u.x[] = 0.;
  }
}

double dtmin = HUGE, cdtmin = HUGE;

event logfile (i++) {
  dtmin = min (dtmin, dt);
#ifdef FICK_CORRECTED
  cdtmin = min (cdtmin, corrective_uodx > 0. ? 0.5/corrective_uodx : HUGE);
#endif
}

event stop (t = tend) {
  FILE * fp = fopen ("prof", "w");
  fprintf (fp, "# x");
  for (int jj = 0; jj < NGS; jj++)
    fprintf (fp, " %s", OpenSMOKE_NamesOfSpecies (jj));
  fprintf (fp, " TG u\n");
  foreach()
    if (y < L0/(1 << LEVEL)) { // the first row: the solution is 1-D
    fprintf (fp, "%.10g", x);
    for (int jj = 0; jj < NGS; jj++) {
      scalar YG = YGList_G[jj];
      fprintf (fp, " %.12g", YG[]);
    }
    fprintf (fp, " %.12g %.12g\n", TG[], u.x[]);
    }
  fclose (fp);
  fprintf (stderr, "level %d coupled %d: steps %d, t %g, dt min %g,"
           " corrective limit min %g\n", LEVEL, COUPLED_DIFFUSION, i, t,
           dtmin, cdtmin);
  return 1;
}
