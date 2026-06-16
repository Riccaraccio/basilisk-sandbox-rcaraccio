const double tend = 1., P = 101325.;
double TG0 = 600;
int NGS;

#include "common.h"
#include "run.h"
#include "opensmoke.h"
#include "binning.h"
#include "view.h"

#define R_GAS 8.31446261815324

typedef struct {
  double rhos;
  double rhog;
  double cps;
  double cpg;
  double P;
  double T;
  double zeta;
  double* sources;
} UserDataODE;

void gas_batch_nonisothermal_constantpressure (const double * y, const double dt, double * dy, void * args) {

  UserDataODE data = *(UserDataODE *)args;
  // Clamp to the thermo/kinetics validity range: the stiff Gear corrector can
  // probe unphysical temperature iterates during a sharp ignition, and an
  // out-of-range T overflows the Arrhenius/equilibrium evaluations to Inf/NaN
  // (raised as SIGFPE inside OpenSMOKE). Bounding the RHS evaluation keeps the
  // Newton iteration finite so it can converge to the physical solution.
  double Temperature = clamp (y[NGS], 280., 3500.);

  OpenSMOKE_GasProp_SetTemperature (Temperature);
  OpenSMOKE_GasProp_SetPressure (data.P);

  // Unpack mass fractions
  double gasmassfracs[NGS], gasmolefracs[NGS];
  for (int jj=0; jj<NGS; jj++)
    gasmassfracs[jj] = y[jj] < 0. ? 0. : y[jj];

  // Calculate mole fractions and mixture molecular weight
  double MWMix;
  OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, &MWMix, gasmassfracs);

  // Calculate concentrations and reaction rates
  double ctot = data.P/(R_GAS*1000*Temperature); // kmol/m3
  double cgas[NGS], rgas[NGS];
  for (int jj=0; jj<NGS; jj++) {
    cgas[jj] = ctot*gasmolefracs[jj];
    rgas[jj] = 0.;
  }
  data.rhog = ctot*MWMix;
  data.cpg = OpenSMOKE_GasProp_HeatCapacity (gasmolefracs);
  OpenSMOKE_GasProp_ReactionRates (cgas);
  OpenSMOKE_GasProp_FormationRates (rgas); //[kmol/m3_gas/s]

  double QRgas = OpenSMOKE_GasProp_HeatRelease (rgas);

  for (int jj=0; jj<NGS; jj++) {
    dy[jj] = OpenSMOKE_MW(jj)*rgas[jj]/data.rhog;
    if (data.sources)
      data.sources[jj] = dy[jj]*data.rhog; // source term for gas expansion
  }

  //Temperature equation
  dy[NGS] = QRgas/(data.rhog*data.cpg);
  if (data.sources)
    data.sources[NGS] = dy[NGS]*data.rhog*data.cpg;
}

int main() {
  init_grid (1 << 5);
  DT = 1e-2;
  run();
}

scalar T[], *YList = NULL;          // binning solution
scalar Tref[], *YListRef = NULL;    // reference (per-cell) solution
double eps = 1e-2;
scalar binid[];
event defaults (i = 0) {
  char kinfolder_root[128];
  sprintf (kinfolder_root, "%s/kinetics/skeletal/methanol/kinetics",
      getenv ("OPENSMOKE_INTERFACE"));

  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder_root);
  NGS = OpenSMOKE_NumberOfSpecies();

  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[64];
    snprintf (name, sizeof(name), "Y_%s",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    YList = list_append (YList, a);
  }
  reset (YList, 0.);

  // Independent species list for the reference (per-cell) solution.
  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[64];
    snprintf (name, sizeof(name), "Yref_%s", OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    YListRef = list_append (YListRef, a);
  }
  reset (YListRef, 0.);

  OpenSMOKE_InitODESolver();
}

#define radial(r,r1,r2,T1,T2) (-r1*r2/(r*(r1-r2))*(T1 - T2) + (r1*T1 - r2*T2)/(r1-r2))
#define gaussian(x, y, sigma) ( exp(-((sq(x - 0.) + sq(y - 0.)) / (2 * sq(sigma)))) )
scalar rho[], cp[], mask[], rhoRef[], cpRef[], * fields = NULL, * targets = NULL;
void initial_conditions () {
  foreach() {
    double r = sqrt (sq(x) + sq(y));
    mask[] = (r <= 0.2) ? 0. : 1.;
    T[] = radial (r, 0.2, 0.8, 750, TG0);
    T[] = (r <= 0.2) ? 750 : (r >= 0.8) ? TG0 : T[];

    scalar fuel = YList[OpenSMOKE_IndexOfSpecies("CH3OH")];
    fuel[] = gaussian (x, y, 0.2)*0.1;

    // O2 is depleted near the origin and approaches an air-like value far away
    scalar oxy = YList[OpenSMOKE_IndexOfSpecies("O2")];
    oxy[] = (1. - gaussian (x, y, 0.2))*0.23;

    // N2 makes the mass fractions sum to 1 everywhere
    scalar inert = YList[OpenSMOKE_IndexOfSpecies("N2")];
    inert[] = 1. - fuel[] - oxy[];

    T[] *= mask[];
    for (int jj = 0; jj < NGS; jj++) {
      scalar Y = YList[jj];
      Y[] *= mask[];
    }

    rho[] = 0.;
    cp[] = 0.;
    if (mask[] > 0) {
      double yg[NGS];
      for (int jj = 0; jj < NGS; jj++) {
        scalar Y = YList[jj];
        yg[jj] = Y[];
      }

      double xg[NGS], MWMix;
      OpenSMOKE_MoleFractions_From_MassFractions(xg, &MWMix, yg);
      rho[] = OpenSMOKE_GasProp_Density_IdealGas(T[], P, MWMix);

      OpenSMOKE_GasProp_SetTemperature (T[]);
      OpenSMOKE_GasProp_SetPressure (P);
      cp[] = OpenSMOKE_GasProp_HeatCapacity(xg);
    }
  }
}

event init (i = 0) {
  fields = list_concat (YList, {T});
  initial_conditions();

  /**
  The binning *targets* are the thermochemical state variables used to group
  cells: the temperature and the fuel mass fraction. Cells whose (T, Y_CH3OH)
  fall in the same eps-cell of this normalized 2D space are agglomerated into
  one bin and integrated together. We deliberately keep this set small: the bin
  id is built as a mixed-radix number over the targets, so too many target
  dimensions would overflow. */

  targets = list_append (targets, T);
  targets = list_append (targets, YList[OpenSMOKE_IndexOfSpecies("CH3OH")]);
  //targets = list_append (targets, YList[OpenSMOKE_IndexOfSpecies("O2")]);

  /**
  Copy the initial state into the reference solution so the binning and
  non-binning approaches start from identical conditions. */

  foreach() {
    Tref[] = T[];
    rhoRef[] = rho[];
    cpRef[] = cp[];
    for (int jj = 0; jj < NGS; jj++) {
      scalar Y = YList[jj], Yr = YListRef[jj];
      Yr[] = Y[];
    }
  }
}

event timestep (i++) {
  dtnext (DT);
}

size_t nactive_bins = 0;

/**
Reference ("non-binning") solution: integrate the chemistry ODE directly in
every active cell — one stiff solve per cell. */

void chemistry_percell (scalar * YL, scalar Tf, scalar rhof, scalar cpf) {
  foreach() {
    if (mask[] > 0) {
      double y0ode[NGS + 1];
      for (int jj = 0; jj < NGS; jj++) {
        scalar Y = YL[jj];
        y0ode[jj] = Y[];
      }
      y0ode[NGS] = Tf[];

      UserDataODE data;
      data.P = P;
      data.sources = NULL;
      data.rhog = rhof[];
      data.cpg = cpf[];

      OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
                           NGS + 1, dt, y0ode, &data);

      for (int jj = 0; jj < NGS; jj++) {
        scalar Y = YL[jj];
        Y[] = fmax (0., y0ode[jj]);
      }
      rhof[] = data.rhog;
      cpf[] = data.cpg;
      Tf[] = y0ode[NGS];
    }
  }
}

event chemistry (i++) {

  /**
  Reference solution: one direct stiff solve per active cell. */

  chemistry_percell (YListRef, Tref, rhoRef, cpRef);

  /**
  Agglomerate the masked cells into bins of similar thermochemical state. The
  bin carries the mass-averaged field values (`fields` = species + T, in that
  order) together with the mass-averaged density and heat capacity. */

  BinTable * table = binning (fields, targets, (double[]){eps, eps, eps},
                              rho, cp, mask);

  /**
  Integrate the stiff chemistry ODE **once per bin** instead of once per cell.
  `bin->phi[j]` holds the averaged value of `fields[j]`, so the ODE state vector
  `y0ode` maps directly onto it: entries `[0..NGS-1]` are the gas species and
  entry `[NGS]` is the temperature. */

  foreach_bin (table) {
    double y0ode[NGS + 1];
    for (size_t j = 0; j < bin->nfields; j++)
      y0ode[j] = bin->phi[j];

    UserDataODE data;
    data.P = P;
    data.sources = NULL;
    data.rhog = bin->rho;
    data.cpg = bin->cp;

    OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
                         NGS + 1, dt, y0ode, &data);

    for (size_t j = 0; j < bin->nfields; j++)
      bin->phi[j] = (j < (size_t)NGS) ? fmax (0., y0ode[j]) : y0ode[j];

    bin->rho = data.rhog;
    bin->cp = data.cpg;
  }

  /**
  Store the partition for visualization, then conservatively map the bin-level
  increments back onto the individual cells and free the table. */

  binning_ids (table, binid);
  nactive_bins = binning_stats (table).nactive;
  binning_remap (table, fields, rho, cp);
  binning_cleanup (table);
}

event log (i++) {
  int iFuel = OpenSMOKE_IndexOfSpecies("CH3OH");
  int iH2O  = OpenSMOKE_IndexOfSpecies("H2O");

  // Global mass of the key species, binning vs. reference.
  double F_bin   = statsf(YList[iFuel]).sum,  F_ref   = statsf(YListRef[iFuel]).sum;
  double H2O_bin = statsf(YList[iH2O]).sum,   H2O_ref = statsf(YListRef[iH2O]).sum;

  // Pointwise discrepancy (max over active cells) between the two solutions.
  scalar fb = YList[iFuel], fr = YListRef[iFuel];
  double errT = 0., errFuel = 0.;
  foreach (reduction(max:errT) reduction(max:errFuel))
    if (mask[] > 0) {
      errT    = max (errT,    fabs (T[]   - Tref[]));
      errFuel = max (errFuel, fabs (fb[]  - fr[]));
    }

  // columns: t  CH3OH_bin CH3OH_ref  H2O_bin H2O_ref  maxErrT[K] maxErrFuel  nbins
  fprintf (stderr, "%g %g %g %g %g %g %g %zu\n",
           t, F_bin, F_ref, H2O_bin, H2O_ref, errT, errFuel, nactive_bins);

  clear();
  view (tx=-0.5, ty=-0.5, width=1080, height=1080);
  squares ("binid");
  labels ("binid");
  save ("binid.mp4");
}


event cleanup (t = tend) {
  free (targets), targets = NULL;
  delete (YList), free (YList), YList = NULL;
  delete (YListRef), free (YListRef), YListRef = NULL;
  OpenSMOKE_CleanODESolver();
}
