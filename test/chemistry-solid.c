/**
# Solid-phase chemistry integration test

Analogue of `chemistry.c` (the gas-phase test) for the **solid-gas** branch of
`chemistry.h`. There is no flow, no transport and no velocity: the domain is a
grid of independent 0-D batch reactors (one per cell). Every cell is fully
solid (`f = 1`), so the gas-phase branch is skipped and only
`solid_batch_nonisothermal_constantpressure` runs.

The kinetics are the `biomass/dummy-solid-gas` scheme, whose solid reactions are

    BIOMASS => 0.8 TAR + H2O + 1.2 CHAR     (pyrolysis, endothermic)
    MOIST   => H2O                          (drying)
    H2O     => MOIST                         (re-condensation)

The initial *solid* composition varies along Y, so each row of cells pyrolyses a
different biomass loading:

* BIOMASS mass fraction: 0.85 at the bottom (y = 0)  ->  0.45 at the top (y = L0)
* ASH     mass fraction: 0.10 at the bottom          ->  0.50 at the top (inert filler)
* MOIST   mass fraction: 0.05 everywhere (mild drying)

(Moisture is kept low and fixed: above ~15% the dummy scheme's MOIST<->H2O heat
release drives a thermal runaway, which would swamp the pyrolysis signal.)

The pores are initially filled with inert N2 gas (`YGList_S`), into which the
released TAR / H2O accumulate.

In addition to the composition and temperature, we track the **source term**
predicted by the ODE right-hand side: `chemistry.h` stores
`data.sources[NGS+NSS]` (the net solid reaction rate, i.e. the gas mass
generation rate per unit solid volume) into the cell field `omega[]`.
*/

#define F_ERR 1e-10
#define SOLVE_TEMPERATURE 1      // integrate the solid temperature TS[]
#define STORE_SOURCES 1          // copy the full data.sources vector into sourcesList
// no GAS_PHASE_REACTIONS -> pure solid-gas pyrolysis

scalar f[];           // solid volume fraction (= 1 everywhere here -> all solid)

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
scalar porosity[];    // solid-matrix porosity (evolves with pyrolysis)
scalar zeta[];        // shrinkage / gas-release partitioning (0 -> pure gas release)

double rhoG = 0.3;    // gas density   [kg/m3] (constant, no VARPROP)
double rhoS = 1500.;  // solid density [kg/m3]

#include "memoryallocation-varprop.h"

/**
Per-cell copy of the reactor `data.sources` vector (length NEQ = NGS+NSS+2 in
the solid branch). Allocated once the kinetic scheme is loaded. */

scalar * sourcesList = NULL;

/**
`run.h` provides no CFL constraint, so we drive the timestep ourselves at a
fixed `DT`. Declared before `chemistry.h` so it runs first each iteration. */

event timestep (i++) {
  dtnext (DT);
}

#include "chemistry.h"

/**
Species indices, resolved once the kinetic scheme is loaded. */
int iBIO, iMOIST, iASH, iCHAR, iN2;

/**
The reactor's `data.sources` vector (now mirrored in `sourcesList`) is laid out
as, in the solid branch:

* `[0 .. NGS-1]`        gas-species mass sources  (gas-phase reactions only -> 0 here)
* `[NGS .. NGS+NSS-1]`  solid-species mass sources (never set by the reactor -> 0)
* `[NGS+NSS]`           porosity-equation source = -(net solid reaction rate) == omega[]
* `[NGS+NSS+1]`         temperature-equation source (heat release rate)

`source_labels()` writes a matching header so the columns can be identified. */

void source_labels (FILE * fp) {
  for (int jj = 0; jj < NGS; jj++)
    fprintf (fp, " s_%s", OpenSMOKE_NamesOfSpecies (jj));
  for (int jj = 0; jj < NSS; jj++)
    fprintf (fp, " s_%s", OpenSMOKE_NamesOfSolidSpecies (jj));
  fprintf (fp, " s_poros s_temp");
}

/**
Output a vertical profile of the solid composition, temperature, porosity and
the **full** `data.sources` vector, sampled at the centre column (x = L0/2).
Since `f = 1`, the solid mass fractions are simply `YS[]`. */

void write_profile (const char * fname) {
  FILE * fp = fopen (fname, "w");
  fprintf (fp, "#y BIOMASS CHAR MOIST ASH porosity TS");
  source_labels (fp);
  fprintf (fp, "\n");
  int n = 1 << 4;
  for (int j = 0; j < n; j++) {
    double yy = Y0 + L0*(j + 0.5)/n, xx = L0/2.;
    fprintf (fp, "%g", yy);
    const char * sp[] = {"BIOMASS","CHAR","MOIST","ASH"};
    for (int k = 0; k < 4; k++) {
      scalar YS = YSList[OpenSMOKE_IndexOfSolidSpecies (sp[k])];
      fprintf (fp, " %g", interpolate (YS, xx, yy));
    }
    fprintf (fp, " %g %g", interpolate (porosity, xx, yy), interpolate (TS, xx, yy));
    for (scalar src in sourcesList)
      fprintf (fp, " %g", interpolate (src, xx, yy));
    fprintf (fp, "\n");
  }
  fclose (fp);
}

int main() {
  TS0 = 800.; TG0 = 800.;          // solid/gas temperature [K] (hot enough to pyrolyse)
  cpS = 1500.; cpG = 1100.;        // constant specific heats [J/kg/K]
  kinfolder = "biomass/dummy-solid-gas";
  DT = 1e-3;                       // macro step (chemistry sub-steps internally)
  init_grid (1 << 4);              // 16 x 16 independent batch reactors
  run();
}

/**
We initialise every field by hand and set `restarted = true` so that the
uniform initialisation in `memoryallocation-varprop.h` is skipped: this makes
the spatial gradient the true initial condition regardless of event order. */

event init (i = 0) {
  restarted = true;

  /**
  Allocate one field per `data.sources` entry. NGS/NSS are set by the header's
  `defaults` event (i = 0), which runs before this `init`. */

  int NEQ = NGS + NSS + 2;
  for (int jj = 0; jj < NEQ; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[64];
    snprintf (name, sizeof(name), "src%d", jj);
    a.name = strdup (name);
    sourcesList = list_append (sourcesList, a);
  }
  reset (sourcesList, 0.);

  iBIO   = OpenSMOKE_IndexOfSolidSpecies ("BIOMASS");
  iCHAR  = OpenSMOKE_IndexOfSolidSpecies ("CHAR");
  iMOIST = OpenSMOKE_IndexOfSolidSpecies ("MOIST");
  iASH   = OpenSMOKE_IndexOfSolidSpecies ("ASH");
  iN2    = OpenSMOKE_IndexOfSpecies ("N2");

  foreach() {
    f[] = 1.;                      // fully solid cell
    porosity[] = 0.4;              // initial pore fraction
    zeta[] = 0.;                   // gas release only, no shrinkage
    p[] = 0.;

    double s = (y - Y0)/L0;        // 0 at the bottom, 1 at the top
    double yBIO   = 0.85 - 0.4*s;   // 0.85 -> 0.45
    double yMOIST = 0.05;           // mild, fixed moisture
    double yASH   = 1. - yBIO - yMOIST; // inert filler (0.10 -> 0.50)

    // Solid mass fractions (x f, with f = 1); all others start at 0.
    for (int jj = 0; jj < NSS; jj++) {
      scalar YS = YSList[jj];
      YS[] = (jj == iBIO)   ? yBIO   :
             (jj == iMOIST) ? yMOIST :
             (jj == iASH)   ? yASH   : 0.;
    }

    // Pore gas: pure inert N2; released TAR / H2O accumulate here.
    for (int jj = 0; jj < NGS; jj++) {
      scalar YGs = YGList_S[jj];
      YGs[] = (jj == iN2) ? 1. : 0.;
    }

    // Fields the framework still expects to be defined.
    for (int jj = 0; jj < NGS; jj++) {
      scalar YGg   = YGList_G[jj];   YGg[]  = 0.;
      scalar YGi   = YGList_Int[jj]; YGi[]  = 0.;
      scalar sSexp = sSexpList[jj];  sSexp[] = 0.;
      scalar sGexp = sGexpList[jj];  sGexp[] = 0.;
    }

    TS[] = TS0;                     // solid temperature (f = 1)
    TG[] = 0.;                      // no free gas
    T[]  = TS0;
    TInt[] = 0.;
  }

  write_profile ("profile-initial");   // true initial condition (pre-reaction)
}

event profile1 (t = end) { write_profile ("profile-final"); }

/**
Track the bottom (biomass-rich) cell over time: BIOMASS fraction, solid
temperature, and the **full** `data.sources` vector. The column layout is
`t BIOMASS TS` followed by `source_labels()` (the i=0 header is printed once). */

event monitor (i++) {
  double xb = L0/2., yb = Y0 + 0.5*L0/16.;
  scalar YBIO = YSList[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")];
  if (i == 0) {
    fprintf (stderr, "#t BIOMASS TS");
    source_labels (stderr);
    fprintf (stderr, "\n");
  }
  fprintf (stderr, "%g %g %g", t,
           interpolate (YBIO, xb, yb), interpolate (TS, xb, yb));
  for (scalar src in sourcesList)
    fprintf (stderr, " %g", interpolate (src, xb, yb));
  fprintf (stderr, "\n");
}

event stop (i = 300);
