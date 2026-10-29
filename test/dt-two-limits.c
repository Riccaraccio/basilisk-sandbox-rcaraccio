/**
# The timestep with the two stability events

`centered-phasechange.h` and `shrinking.h` each have a `stability` event. The
event of `shrinking.h` runs first. Before the fix, each event called
`timestep()` and `dtnext()`. The velocity `ubf` is zero outside the solid, so
the first call returned only the growth value of `timestep()` and set `tnext`
from it. The second call got the limit of `uf`. When this limit was between
half of `tnext - t` and `tnext - t`, `dtnext()` cut the step in half. A
restart of fatehi at level 10 with `COUPLED_DIFFUSION` (2026-10-28) had
`dt = 0.55 dt_uf` in 94 % of the steps.

This test has no solid (`f = 0`), thus `ubf = 0`, and a uniform flow `U`.
`DT` is 2.4 times the limit of `uf`, as in fatehi. After the growth from the
first step, `dt` must stay at the limit of `uf`, `CFL*Delta/U`. The test
counts the steps to `tend` and gives the smallest `dt` between t = 0.1 and
0.25. It does not use the last `dt`: `dtnext()` changes the last steps so
that the run stops at `tend`.

Result (2026-10-29), `CFL*Delta/U = 2.5e-3`, `tend = 0.3`:

| build | steps | smallest dt/limit |
|---|---|---|
| before the fix | 212 | 0.56 |
| fix | 131 | 0.99 |
*/

#define NO_ADVECTION_DIV 1
#define TURN_OFF_REACTIONS 1
#define SOLVE_TEMPERATURE 1

#ifndef LEVEL
# define LEVEL 7
#endif

#include "grid/multigrid.h"
#include "navier-stokes/centered-phasechange.h"

#define MULTICOMPONENT
#include "opensmoke.h"
#include "constant-properties.h"

#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"

const double U = 0.05;      // velocity [m/s]
const double tend = 0.3;

u.n[left]  = dirichlet (U);
u.t[left]  = dirichlet (0.);
p[left]    = neumann (0.);
psi[left]  = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right]   = dirichlet (0.);
psi[right] = neumann (0.);

double dtlim;

int main() {
  rhoG = 0.31;
  muG  = 4.5e-5;
  rhoS = 1550.;
  eps0 = 0.2;

  TS0 = 300.; TG0 = 300.;
  lambdaG = 0.08; cpG = 1200.;
  lambdaS = 0.2;  cpS = 1500.;

  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_REACTION;
  kinfolder = "biomass/dummy-solid";

  L0 = 0.04;
  TOLERANCE = 1e-6;

  init_grid (1 << LEVEL);
  run();
}

/**
`CFL` is set in `init`: `centered.h` resets it in `defaults`. */

event init (i = 0) {
  restarted = true;

  CFL = 0.4;
  dtlim = CFL*(L0/(1 << LEVEL))/U;
  DT = 2.4*dtlim;

  int iN2 = OpenSMOKE_IndexOfSpecies ("N2");
  foreach() {
    f[] = 0., porosity[] = 0.;
    for (int jj = 0; jj < NGS; jj++) {
      scalar YS = YGList_S[jj], YG = YGList_G[jj], YI = YGList_Int[jj];
      YS[] = 0.;
      YG[] = (jj == iN2);
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
    u.x[] = U;
  }
  TG[left] = dirichlet (TG0);
  for (int jj = 0; jj < NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == iN2)
      YG[left] = dirichlet (1.);
    else
      YG[left] = dirichlet (0.);
  }
}

double ratiomin = HUGE;

event logfile (i++) {
  static FILE * fp = fopen ("dt.dat", "w");
  fprintf (fp, "%d %g %g %g\n", i, t, dt, dt/dtlim);
  if (t > 0.1 && t < 0.25)
    ratiomin = min (ratiomin, dt/dtlim);
}

event stop (t = tend) {
  fprintf (stderr, "limit %g, steps %d, smallest dt/limit %.3f\n",
           dtlim, i, ratiomin);
  assert (ratiomin > 0.9);
  assert (i < 1.1*tend/dtlim + 40);
  return 1;
}
