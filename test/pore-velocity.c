/**
# The velocity of the pore species

This test measures the velocity at which `multicomponent-varprop.h` moves the
gas species in the pores (`YGList_S`). It is the check of item 4 of the
discretization report.

A porous slab fills the full height of a 2D channel. A uniform flow `U`
enters from the left. The flow is 1D, so the superficial velocity is `U` in
the slab and in the gas. The pores start full of `TAR`, and the gas outside
is `N2`. The reactions are off, so `TAR` is a passive species.

The pore gas moves at the interstitial velocity `U/eps`. The plug flow
therefore flushes the slab of length `Ls` in `eps*Ls/U`. A code that moves the
pore species with `U` flushes it in `Ls/U`, which is `1/eps` times longer.

The test writes `pore.dat`: the time, the mass of `TAR` in the pores over its
initial value, the timestep, and the CFL limit of the pore velocity. At the
end it writes the times at which the mass falls to 50 % and to 10 %.

Build the two versions of the transport:

    -DPORE_SPECIES_INTERSTITIAL=1   (the default) u/eps
    -DPORE_SPECIES_INTERSTITIAL=0   the previous code, u

`-DDT_VALUE=<dt>` caps the timestep. Use it to run the two versions with the
same timestep.

Result (2026-09-19, level 7, `eps0 = 0.2`, `Ls/U = 0.2 s`):

| build | dt | t50 [s] | t10 [s] |
|---|---|---|---|
| interstitial (1) | 6.25e-4, pore limit | 0.0221 | 0.0386 |
| previous (0) | 6.25e-4, `DT` cap | 0.1020 | 0.1959 |
| previous (0) | 3.05e-3, free | 0.1010 | 0.1946 |
| interstitial, `SLAB_SHIFT=0.5` | 6.25e-4 | 0.0242 | 0.0411 |
| previous, `SLAB_SHIFT=0.5` | 6.25e-4 | 0.1017 | 0.1960 |

The ratio of the flush times is 0.22 (t50) and 0.20 (t10), which is `eps0`.
The plug flow gives t50 = 0.02 s and t10 = 0.036 s. The build 0 gives the
dump of commit 93ed0a0 bit for bit. */

#define NO_ADVECTION_DIV 1
#define TURN_OFF_REACTIONS 1

/**
The stack does not build without `SOLVE_TEMPERATURE`. The temperature is
uniform and stays uniform, because the reactions are off. */

#define SOLVE_TEMPERATURE 1

#ifndef DT_VALUE
# define DT_VALUE 1.
#endif

#ifndef LEVEL
# define LEVEL 7
#endif

#ifndef SLAB_SHIFT
# define SLAB_SHIFT 0.   // shift of the slab in cells: 0.5 gives cut cells
#endif

#include "grid/multigrid.h"
#include "navier-stokes/centered-phasechange.h"

#define MULTICOMPONENT
#include "opensmoke.h"
#include "constant-properties.h"

#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"

const double U = 0.05;      // superficial velocity [m/s]
const double Ls = 0.01;     // length of the slab [m]
const double xs = 0.01;     // left face of the slab [m]
const double tend = 0.3;

u.n[left]  = dirichlet (U);
u.t[left]  = dirichlet (0.);
p[left]    = neumann (0.);
psi[left]  = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right]   = dirichlet (0.);
psi[right] = neumann (0.);

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
  shift_prod = true;

  Da = (coord){1e-9, 1e-9};

  L0 = 0.04;
  DT = DT_VALUE;
  TOLERANCE = 1e-6;

  init_grid (1 << LEVEL);
  run();
}

double mass0 = 0.;

double pore_tar_mass (void) {
  scalar YT = YGList_S[OpenSMOKE_IndexOfSpecies ("TAR")];
  double m = 0.;
  foreach (reduction(+:m))
    if (f[] > F_ERR)
      m += YT[]/f[]*porosity[]*rhoG*dv();   // YT and porosity: tracer form
  return m;
}

event init (i = 0) {

  /**
  `restarted` stops the uniform fill of `memoryallocation-varprop.h`, which
  gives the pores and the gas the same composition. So this event fills
  every field. */

  restarted = true;

  double x0 = xs + SLAB_SHIFT*L0/(1 << LEVEL);
  fraction (f, min (x - x0, x0 + Ls - x));
  foreach()
    porosity[] = eps0*f[];

  int iN2 = OpenSMOKE_IndexOfSpecies ("N2");
  int iTAR = OpenSMOKE_IndexOfSpecies ("TAR");
  foreach() {
    for (int jj = 0; jj < NGS; jj++) {
      scalar YS = YGList_S[jj], YG = YGList_G[jj], YI = YGList_Int[jj];
      YS[] = (jj == iTAR)*f[];
      YG[] = (jj == iN2)*(1. - f[]);
      YI[] = (jj == iN2);
    }
    for (int jj = 0; jj < NSS; jj++) {
      scalar YSol = YSList[jj];
      YSol[] = (jj == 0)*f[];
    }
  }

  foreach() {
    TS[] = TS0*f[];
    TG[] = TG0*(1. - f[]);
    T[]  = TS[] + TG[];
    TInt[] = (f[] > F_ERR && f[] < 1. - F_ERR) ? TS0 : 0.;
  }
  TG[left] = dirichlet (TG0);

  for (int jj = 0; jj < NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == iN2)
      YG[left] = dirichlet (1.);
    else
      YG[left] = dirichlet (0.);
  }

  foreach()
    u.x[] = U;

  mass0 = pore_tar_mass();
}

double t50 = -1., t10 = -1., mprev = 1., tprev = 0.;

event logfile (i++) {
  static FILE * fp = fopen ("pore.dat", "w");
  if (i == 0)
    fprintf (fp, "#t m/m0 dt pore_dtmax\n");

  double m = pore_tar_mass()/mass0;
  fprintf (fp, "%g %.6g %g %g\n", t, m, dt, pore_dtmax);
  fflush (fp);

  // linear interpolation of the crossing times
  if (t50 < 0. && m <= 0.5)
    t50 = tprev + (t - tprev)*(mprev - 0.5)/(mprev - m);
  if (t10 < 0. && m <= 0.1)
    t10 = tprev + (t - tprev)*(mprev - 0.1)/(mprev - m);
  mprev = m, tprev = t;
}

event stop (t = tend) {
  double umean = 0., n = 0.;
  foreach_face (x, reduction(+:umean) reduction(+:n))
    if (x > xs + 0.25*Ls && x < xs + 0.75*Ls)
      umean += uf.x[], n += 1.;
  fprintf (stderr, "PORE_SPECIES_INTERSTITIAL=%d eps0=%g U=%g Ls=%g level=%d"
           " DT=%g\n", PORE_SPECIES_INTERSTITIAL, eps0, U, Ls, LEVEL,
           (double) DT_VALUE);
  fprintf (stderr, "uf in the slab %g, t50 %g, t10 %g,"
           " plug flow eps*Ls/U %g, Ls/U %g, steps %d\n",
           umean/n, t50, t10, eps0*Ls/U, Ls/U, i);
#if DUMP_END
  dump (file = "end-dump"); // for a bit for bit comparison of two builds
#endif
  return 1;
}
