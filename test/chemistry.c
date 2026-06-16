/**
# Gas-phase chemistry integration test

Minimal test of the **gas-phase** branch of `chemistry.h` with no flow, no
transport and no velocity. The domain is a grid of independent 0-D batch
reactors (one per cell). Every cell is pure gas (`f = 0`), so the solid-gas
branch is skipped and only `gas_batch_nonisothermal_constantpressure` runs.

The initial composition is a CH4 / O2 / N2 mixture that varies along Y:

* CH4 mass fraction: 0.8 at the bottom (y = 0)  ->  0.2 at the top (y = L0)
* O2  mass fraction: 0.1 at the bottom          ->  0.4 at the top
* N2  fills the remainder (0.1 at the bottom    ->  0.4 at the top)

so each row of cells ignites a different fuel/oxidiser ratio. The gas-phase
kinetics come from the `biomass/dummy-solid-gas` scheme, whose gas reactions
include CH4 oxidation, the water-gas shift, etc.
*/

#define F_ERR 1e-10
#define SOLVE_TEMPERATURE 1      // gas batch always integrates TG[]
#define GAS_PHASE_REACTIONS 1    // enable the gas-phase reaction branch

scalar f[];           // solid volume fraction (= 0 everywhere here -> all gas)

#include "run.h"

/**
The species bookkeeping in `memoryallocation-varprop.h` uses the VOF
`tracers`/`inverse` field attributes (normally declared by `vof.h`, which we
do not include because it would drag in the flow solver). */

attribute {
  scalar * tracers, c;
  bool inverse;
}

scalar p[];           // pressure: kept at 0, so the reactor sees P = Pref
scalar porosity[];    // unused in the gas branch, but referenced by chemistry.h
scalar zeta[];        // idem

double rhoG = 0.3;    // gas density   [kg/m3] (constant, no VARPROP)
double rhoS = 1500.;  // solid density [kg/m3] (unused, f = 0)

#include "memoryallocation-varprop.h"

/**
`run.h` provides no CFL constraint, so we drive the timestep ourselves at a
fixed `DT`. Declared before `chemistry.h` so it runs first each iteration. */

event timestep (i++) {
  dtnext (DT);
}

#include "chemistry.h"

/**
Species indices, resolved once the kinetic scheme is loaded. */
int iCH4, iO2, iN2;

/**
Output a vertical profile of the composition and temperature, sampled at the
centre column (x = L0/2). */

void write_profile (const char * fname) {
  FILE * fp = fopen (fname, "w");
  fprintf (fp, "#y CH4 O2 N2 CO H2O CO2 H2 TG\n");
  int n = 1 << 4;
  for (int j = 0; j < n; j++) {
    double yy = Y0 + L0*(j + 0.5)/n, xx = L0/2.;
    fprintf (fp, "%g", yy);
    const char * sp[] = {"CH4","O2","N2","CO","H2O","CO2","H2"};
    for (int k = 0; k < 7; k++) {
      scalar YG = YGList_G[OpenSMOKE_IndexOfSpecies (sp[k])];
      fprintf (fp, " %g", interpolate (YG, xx, yy));
    }
    fprintf (fp, " %g\n", interpolate (TG, xx, yy));
  }
  fclose (fp);
}

int main() {
  TG0 = 900.; TS0 = 300.;          // gas temperature [K] (lowered to slow the kinetics)
  cpS = 1500.; cpG = 1100.;        // constant specific heats [J/kg/K]
  kinfolder = "biomass/Red-gas-2507";
  DT = 1e-4;                       // macro step (chemistry sub-steps internally)
  init_grid (1 << 4);              // 16 x 16 independent batch reactors
  run();
}

/**
We initialise every field by hand and set `restarted = true` so that the
uniform initialisation in `memoryallocation-varprop.h` is skipped: this makes
the spatial gradient the true initial condition regardless of event order. */

event init (i = 0) {
  restarted = true;

  iCH4 = OpenSMOKE_IndexOfSpecies ("CH4");
  iO2  = OpenSMOKE_IndexOfSpecies ("O2");
  iN2  = OpenSMOKE_IndexOfSpecies ("N2");

  foreach() {
    f[] = 0.;                      // pure gas cell
    porosity[] = 0.;
    zeta[] = 0.;
    p[] = 0.;

    double s = (y - Y0)/L0;        // 0 at the bottom, 1 at the top
    double yCH4 = 0.8  - 0.6*s;    // 0.80 -> 0.20
    double yO2  = 0.15 + 0.25*s;   // 0.15 -> 0.40
    double yN2  = 1. - yCH4 - yO2; // remainder (0.05 -> 0.40)

    // Gas-phase mass fractions (external gas list); all others start at 0.
    for (int jj = 0; jj < NGS; jj++) {
      scalar YG = YGList_G[jj];
      YG[] = (jj == iCH4) ? yCH4 :
             (jj == iO2)  ? yO2  :
             (jj == iN2)  ? yN2  : 0.;
    }

    // Fields the gas branch / framework still expect to be defined.
    for (int jj = 0; jj < NGS; jj++) {
      scalar YGs   = YGList_S[jj];   YGs[]  = 0.;
      scalar YGi   = YGList_Int[jj]; YGi[]  = 0.;
      scalar sSexp = sSexpList[jj];  sSexp[] = 0.;
      scalar sGexp = sGexpList[jj];  sGexp[] = 0.;
    }
    for (int jj = 0; jj < NSS; jj++) {
      scalar YS = YSList[jj]; YS[] = 0.;
    }

    TG[] = TG0;                    // gas temperature
    TS[] = 0.;                     // no solid
    T[]  = TG0;
    TInt[] = 0.;
  }

  write_profile ("profile-initial");   // true initial condition (pre-reaction)
}

event profile1 (t = end) { write_profile ("profile-final"); }

/**
Track the bottom (fuel-rich) and top (oxidiser-rich) cells over time. */

event monitor (i++) {
  scalar YCH4 = YGList_G[OpenSMOKE_IndexOfSpecies ("CH4")];
  fprintf (stderr, "%g %g %g %g %g\n", t,
           interpolate (YCH4, L0/2., Y0 + 0.5*L0/16.),   // CH4 bottom
           interpolate (TG,   L0/2., Y0 + 0.5*L0/16.),   // T   bottom
           interpolate (YCH4, L0/2., Y0 + 15.5*L0/16.),  // CH4 top
           interpolate (TG,   L0/2., Y0 + 15.5*L0/16.)); // T   top
}

event stop (i = 600);
