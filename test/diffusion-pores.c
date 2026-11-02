/**
# The two diffusion paths in the pore gas

`test/diffusion-paths.c` showed that the two paths of
`multicomponent-varprop.h` agree to about 0.1 % in the gas alone. This test
is the same case with the domain full of porous solid (`f = 1`, porosity
`EPS`). The mixture of `H2`, `TAR` and `CO2` in `N2` fills the pores of the
left half, and `N2` fills the pores of the right half. There is no
interface and no external gas.

In the pore gas the solve of each species divides by
`theta1 = cm*rho*porosity`. The corrective step of the old path
(`fick_correct_species()`) divides by `cm*rho` only, unless
`CORRECTIVE_CAPACITY` is 1. The test runs three variants:

* `-DCOUPLED_DIFFUSION=0`: the old path as in production;
* `-DCOUPLED_DIFFUSION=0 -DCORRECTIVE_CAPACITY=1`: the old path with the
  capacity of the solve;
* `-DCOUPLED_DIFFUSION=1`: the coupled path.

At `tend` the run writes `prof`: x, `YGList_S` of each species (the value of
the pore gas), `TS` and `u.x`. `misc/diffusion-paths.py` compares them.

Build flags: `-DLEVEL=<n>`, `-DEPS=<porosity>`, `-DDT_VALUE=<s>`,
`-DTEND_VALUE=<s>`, `-DNO_EXPANSION=1`, and the two flags above.

## The result (2026-11-02)

The test found two defects of the old path in the pore gas.

1. The `tracer_diffusion` event called
   `uodx = max (uodx, fick_correct_species (YGList_S, ...))`. `max()` is a
   macro of qcc (`common.h`), and it evaluates its second argument again when
   that argument is the larger one. The corrective step of the pore gas thus
   ran two times in each step where the pore gas had the larger `|u_c|/Delta`.
   Commit 2c6a61e (2026-09-30) brought this in. The event now calls each
   function once.
2. The corrective step divides by `cm*rho`, and the solve divides by
   `cm*rho*porosity`. `CORRECTIVE_CAPACITY=1` gives the step the capacity of
   the solve.

Before the fix of 1, the weight of the corrective step of the pore gas was
2 `eps` against the solve. The paths agreed at `eps = 0.5` by chance, and
the gap changed sign between `eps = 0.2` and 1. With `NO_EXPANSION` the gap
stayed: it is in the species equations.

Checks after the fix of 1 (level 6, `dt = 1e-5`, `tend = 2e-3`):

* `eps = 1`, `NO_EXPANSION`: the old path in the pores gives the profiles of
  the old path in the gas (`test/diffusion-paths.c`) bit for bit. Before the
  fix, the largest difference of `TAR` was 1.6e-2. The coupled path gives
  the two bit for bit before and after the fix.
* `eps = 0.2`: the mean of |old - coupled| for `TAR`:

| level | `CORRECTIVE_CAPACITY=0` | `CORRECTIVE_CAPACITY=1` |
|---|---|---|
| 6 | 7.0e-4 | 1.1e-4 |
| 7 | 7.0e-4 | 5.5e-5 |
| 8 | 7.1e-4 | 2.7e-5 |

With both fixes, the old path converges to the coupled path with first
order, as in the gas (the upwind face value of the coupled path). Without
the capacity, the gap does not change with the level. The largest local
difference at level 8 before the two fixes was 4.6e-3 for `TAR` (0.45 at
the start). */

#ifndef EPS
# define EPS 0.2
#endif

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

#ifndef CORRECTIVE_CAPACITY
# define CORRECTIVE_CAPACITY -1 // the coupled path does not use it
#endif

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
  eps0 = EPS;

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
    f[] = 1., porosity[] = eps0;   // tracer form: eps*f
    double s = 0.5*(1. - tanh ((x - L0/2.)/W)); // 1 on the left
    for (int jj = 0; jj < NGS; jj++) {
      scalar YS = YGList_S[jj], YG = YGList_G[jj], YI = YGList_Int[jj];
      YS[] = s*yA[jj] + (1. - s)*(jj == iN2);
      YG[] = 0.;
      YI[] = (jj == iN2);
    }
    for (int jj = 0; jj < NSS; jj++) {
      scalar YSol = YSList[jj];
      YSol[] = (jj == OpenSMOKE_IndexOfSolidSpecies ("CHAR"));
    }
    TS[] = TS0;
    TG[] = 0.;
    T[]  = TS0;
    TInt[] = 0.;
    u.x[] = 0.;
  }
}

double dtmin = HUGE, cdtmin = HUGE;

/**
`mass` gives, at each step, the mass of `H2`, `TAR` and `CO2` in the pores
per unit depth, and the largest |sum Y - 1|. These species do not reach the
outlet before `tend`, so a consistent scheme keeps their mass. */

event logfile (i++) {
  static FILE * fm = fopen ("mass", "w");
  static int jH2 = -1, jTAR, jCO2;
  if (jH2 < 0) {
    jH2 = OpenSMOKE_IndexOfSpecies ("H2");
    jTAR = OpenSMOKE_IndexOfSpecies ("TAR");
    jCO2 = OpenSMOKE_IndexOfSpecies ("CO2");
  }
  scalar Y1 = YGList_S[jH2], Y2 = YGList_S[jTAR], Y3 = YGList_S[jCO2];
  double m1 = 0., m2 = 0., m3 = 0., se = 0.;
  foreach (reduction(+:m1) reduction(+:m2) reduction(+:m3) reduction(max:se)) {
    double w = rhoGv_S[]*porosity[]*dv();
    m1 += w*Y1[], m2 += w*Y2[], m3 += w*Y3[];
    double sY = 0.;
    for (scalar Y in YGList_S)
      sY += Y[];
    se = max (se, fabs (sY - f[]));
  }
  fprintf (fm, "%g %.12g %.12g %.12g %.3e\n", t, m1, m2, m3, se);
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
  fprintf (fp, " TS u\n");
  foreach()
    if (y < L0/(1 << LEVEL)) { // the first row: the solution is 1-D
    fprintf (fp, "%.10g", x);
    for (int jj = 0; jj < NGS; jj++) {
      scalar YS = YGList_S[jj];
      fprintf (fp, " %.12g", YS[]);
    }
    fprintf (fp, " %.12g %.12g\n", TS[], u.x[]);
    }
  fclose (fp);
  fprintf (stderr, "eps %g capacity %d ", eps0, CORRECTIVE_CAPACITY);
  fprintf (stderr, "level %d coupled %d: steps %d, t %g, dt min %g,"
           " corrective limit min %g\n", LEVEL, COUPLED_DIFFUSION, i, t,
           dtmin, cdtmin);
  return 1;
}
