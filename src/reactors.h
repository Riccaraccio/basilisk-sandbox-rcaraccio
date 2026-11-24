/*
# Ideal reactors
This header file implements different functions for the calculation
of the differential equation system for chemistry integration.
Each function calculates the differential dy at the current time t
given the current variables status.
The system of ODEs is solved using the *batch* method.
*/

#include "radiation.h"

/**
## The gas-phase reactions

`TURN_OFF_GAS_REACTIONS` at 1 removes every gas-phase reaction. Two sites
read it. The solid reactor below skips the reactions of the pore gas, and
the `chemistry` event of `chemistry.h` skips the integration of the gas
cells. The solid reactions and the gas that they release stay.

The default is 0. The test is `#if`, so `-DTURN_OFF_GAS_REACTIONS=0` means
the gas reactions are on.

The default lives here and not in `chemistry.h`, because `chemistry.h`
includes this file before its `TURN_OFF_REACTIONS` block. A case can then
read the flag in every build, also in a build without chemistry.

`GAS_PHASE_REACTIONS` does not replace it. No source file reads that name,
so the gas kinetics always run unless this flag is 1.

The flag covers the per-cell path and the `BINNING` path of `chemistry.h`. */

#ifndef TURN_OFF_GAS_REACTIONS
# define TURN_OFF_GAS_REACTIONS 0
#endif

/**
## The reduction of the gas mechanism (DRG)

`GAS_DRG` at 1 reduces the gas mechanism for each stiff solve of the gas
reactor, with the Directed Relation Graph. `chemistry.h` then includes
`drg.h`, which holds the settings and the solve. The default is 0, and the
test is `#if`.

The default and `gas_drg_active` live here, because the gas right-hand side
below reads them, and `chemistry.h` includes this file before `drg.h`.
`gas_drg_active` is true only during a reduced solve of `drg.h`. The
frozen-cell gate and the solid reactor call the right-hand side outside of
that solve, so they always use the full mechanism. */

#ifndef GAS_DRG
# define GAS_DRG 0
#endif

#if GAS_DRG
static bool gas_drg_active = false;
#endif

#define R_GAS 8.31446261815324

/**
The temperature handed to the RHS is an *unconstrained* iterate of the ODE
integrator: the Gear predictor and the finite-difference Jacobian of
`FindCorrection` can evaluate the system at a negative (or NaN) temperature,
even when the initial and final states are perfectly physical. OpenSMOKE's
`Set{Gas,Sol}Prop_SetTemperature` takes `log(T)` (Arrhenius/PLOG), so such an
iterate raises FE_INVALID and SIGFPEs on machines with traps armed. The
kinetics are meaningful only inside a physical window anyway, so clamp before
every property call. The comparisons are written so that a NaN falls back to
`TMIN_ODE` rather than propagating. */

#define TMIN_ODE 200.
#define TMAX_ODE 5000.

static inline double clamp_temperature (double T) {
  T = (T > TMIN_ODE) ? T : TMIN_ODE;
  return (T < TMAX_ODE) ? T : TMAX_ODE;
}

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

/**
## Isothermal function
Deprecated, to be removed
*/

void solid_batch_isothermal_constantpressure (const double* y, const double dt, double* dy, void* args) {

  /**
  Unpack data for the ODE system.
  */

  UserDataODE data = *(UserDataODE *)args;
  double rhos = data.rhos;
  double Pressure = clamp_pressure (data.P);
  double z = data.zeta;
  double* sources = data.sources;

  double epsilon = y[NGS+NSS];
  double Temperature = clamp_temperature (data.T);

  epsilon = clamp(epsilon, 0., 1.);

  OpenSMOKE_SolProp_SetTemperature (Temperature);
  OpenSMOKE_SolProp_SetPressure (Pressure);

  double gasmass[NGS]; double totgasmass = 0.;
  for (int jj=0; jj<NGS; jj++) {
    gasmass[jj] = y[jj];
    gasmass[jj] = (gasmass[jj] < 0.) ? 0. : gasmass[jj];
    totgasmass += gasmass[jj];
  }

  double gasmassfracs[NGS], gasmolefracs[NGS];
  for (int jj=0; jj<NGS; jj++) {
    gasmassfracs[jj] =  gasmass[jj]/totgasmass;
  }
  double MWmix;
  OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, &MWmix, gasmassfracs);

  double ctot = Pressure/R_GAS/Temperature; //ideal gases, mol/m3
  double cgas[NGS], rgas[NGS];
  for (int jj=0; jj<NGS; jj++) {
    cgas[jj] = ctot*gasmolefracs[jj];
  }

  double solidmass[NSS]; double totsolidmass = 0.;
  for (int jj=0; jj<NSS; jj++) {
    solidmass[jj] = y[jj+NGS];
    solidmass[jj] = (solidmass[jj] < F_ERR) ? 0. : solidmass[jj];
    totsolidmass += solidmass[jj];
  }

  double solmassfracs[NSS];
  double csolid[NSS], rsolid[NSS];
  for (int jj=0; jj<NSS; jj++) {
    solmassfracs[jj] = (totsolidmass < F_ERR) ? 0. : solidmass[jj]/totsolidmass;
    csolid[jj] = rhos*solmassfracs[jj]/sol_MWs[jj];
  }

  OpenSMOKE_SolProp_KineticConstants ();
  OpenSMOKE_SolProp_ReactionRates (cgas,csolid);
  OpenSMOKE_SolProp_FormationRates (rgas, rsolid); //rates are given per m3 of solid

  // update the gas species rates
  for (int jj=0; jj<NGS; jj++) {
    dy[jj] = gas_MWs[jj]*rgas[jj]*(1-epsilon);
  }

  // update the solid species rates
  for (int jj=0; jj<NSS; jj++) {
    dy[jj+NGS] = sol_MWs[jj]*rsolid[jj]*(1-epsilon);
  }

  double totsolidreaction = 0.;
  for (int jj=0; jj<NSS; jj++) {
    totsolidreaction += (sol_MWs[jj]*rsolid[jj]);
  }

  //epsilon equation
  dy[NGS+NSS] = -totsolidreaction*(1-epsilon)*(1-z)/rhos;
  if (sources)
    sources[NGS+NSS] = -totsolidreaction;
}

/**
## *solid_batch_nonisothermal_constantpressure()*
solve chemical species reactions in the solid phase

* *y*: vector (length = NS+1) containing the initial values of the system
* *dt*: simulation time step
* *dy*: right-hand-side of the batch system of equations
* *args*: structure with additional arguments to be passed to the system
*/

void solid_batch_nonisothermal_constantpressure (const double * y, const double dt, double * dy, void * args) {

  UserDataODE data = *(UserDataODE *)args;
  double epsilon = y[NGS+NSS];
  double Temperature = clamp_temperature (y[NGS+NSS+1]);
  double Pressure = clamp_pressure (data.P);

  epsilon = clamp(epsilon, 0., 1.);

  OpenSMOKE_SolProp_SetTemperature (Temperature);
  OpenSMOKE_SolProp_SetPressure (Pressure);
  OpenSMOKE_GasProp_SetTemperature (Temperature);
  OpenSMOKE_GasProp_SetPressure (Pressure);

  double totgasmass = 0.;
  double gasmassfracs[NGS], gasmolefracs[NGS];
  for (int jj=0; jj<NGS; jj++) {
    double mass = fmax (0., y[jj]);
    totgasmass += mass;
    gasmassfracs[jj] = mass;
  }

  double inv_totgasmass = (totgasmass > 0.) ? 1./totgasmass : 0.;
  for (int jj=0; jj<NGS; jj++)
    gasmassfracs[jj] *= inv_totgasmass;

  double MWgasMix;
  OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, &MWgasMix, gasmassfracs);

  double ctot = Pressure/(R_GAS*1000*Temperature); // kmol/m3
  double cgas[NGS], rgas[NGS];
  for (int jj=0; jj<NGS; jj++) {
    cgas[jj] = ctot*gasmolefracs[jj];
    rgas[jj] = 0.;
  }
#ifdef VARPROP
  data.rhog = ctot*MWgasMix;
  data.cpg = OpenSMOKE_GasProp_HeatCapacity (gasmolefracs);
#endif

  double totsolmass = 0.;
  double solmassfracs[NSS];
  for (int jj=0; jj<NSS; jj++) {
    double mass = fmax (0., y[jj+NGS]);
    totsolmass += mass;
    solmassfracs[jj] = mass;
  }

  double csolid[NSS], rsolid[NSS];
  double inv_totsolmass = (totsolmass > 0.) ? 1./totsolmass : 0.;
  for (int jj=0; jj<NSS; jj++) {
    solmassfracs[jj] *= inv_totsolmass;
    csolid[jj] = data.rhos*solmassfracs[jj]/sol_MWs[jj];
    rsolid[jj] = 0.;
  }

  #ifdef VARPROP
  double solmolefracs[NSS], MWsolidMix;
  OpenSMOKE_SolidMoleFractions_From_SolidMassFractions(solmolefracs, &MWsolidMix, solmassfracs);
  data.cps = OpenSMOKE_SolProp_HeatCapacity (solmolefracs);
  #endif

  OpenSMOKE_SolProp_ReactionRates (cgas,csolid); 
  OpenSMOKE_SolProp_FormationRates (rgas, rsolid); //[kmol/m3_solid/s]

  double rgas_pure[NGS], QRgas = 0.;
  for (int jj=0; jj<NGS; jj++)
    rgas_pure[jj] = 0.;

  /**
  The reactions of the pore gas. `TURN_OFF_GAS_REACTIONS` skips them, so
  `rgas_pure` and `QRgas` stay at 0. The solid reactions and the gas that they
  release stay. */

#if !TURN_OFF_GAS_REACTIONS
  OpenSMOKE_GasProp_ReactionRates (cgas);
  OpenSMOKE_GasProp_FormationRates (rgas_pure); //[kmol/m3_gas/s]
  QRgas = OpenSMOKE_GasProp_HeatRelease (rgas_pure);
#endif

  for (int jj=0; jj<NGS; jj++) {
    dy[jj] = gas_MWs[jj]* (rgas[jj]*(1-epsilon) + rgas_pure[jj]*epsilon);
    if (data.sources)
      data.sources[jj] = gas_MWs[jj]*rgas_pure[jj]*epsilon; //to be checked
  }

  for (int jj=0; jj<NSS; jj++) {
    dy[jj+NGS] = sol_MWs[jj]*rsolid[jj]*(1-epsilon);
    // data.sources[jj+NGS] = dy[jj+NGS]; //unused
  }

  double totsolreaction = 0.;
  for (int jj=0; jj<NSS; jj++)
    totsolreaction += sol_MWs[jj]*rsolid[jj];

  //epsilon equation
  dy[NGS+NSS] = -totsolreaction*(1-epsilon)*(1-data.zeta)/data.rhos;
  if (data.sources)
    data.sources[NGS+NSS] = -totsolreaction;

  //Temperature equation
  double QRsol = OpenSMOKE_SolProp_HeatRelease (rgas, rsolid);

  dy[NGS+NSS+1] = (QRsol*(1-epsilon) + QRgas*epsilon)/((data.rhog*data.cpg*epsilon) + data.rhos*data.cps*(1-epsilon));
  if (data.sources)
    data.sources[NGS+NSS+1] = QRsol*(1-epsilon) + QRgas*epsilon;
  #ifdef TURN_OFF_HEAT_OF_REACTION
  dy[NGS+NSS+1] *= 0.;
  if (data.sources)
    data.sources[NGS+NSS+1] *= 0.;
  #endif
}

/**
## *gas_batch_nonisothermal_constantpressure()*
solve chemical species reactions in the gas phase

* *y*: vector (length = NS+1) containing the initial values of the system
* *dt*: simulation time step
* *dy*: right-hand-side of the batch system of equations
* *args*: structure with additional arguments to be passed to the system

### The clip of the negative mass fractions

`GAS_RHS_CLIP` (value-tested, default 0) controls the clip of the negative
mass fractions in this function. At 0, the function uses the mass fractions
of the solver as they are. At 1, it sets each negative mass fraction to 0.

Keep the flag at 0 for a scheme with integer reaction orders. The clip
makes the right-hand side flat below 0. If the solver computes the Jacobian
at a state with a negative mass fraction, the Jacobian does not contain the
destruction of that species. For a fast radical with a small steady value,
the step size then stays at the explicit limit, and the solver stops at its
maximum number of steps. In the isobutanol flame of `run/laminarflame.c`,
OCHO (k = 1.1e11 1/s) caused all these failures. One failure cost 12-17 s.
Without the clip, the same state needs 107 calls.

Set the flag to 1 for a scheme with a reaction order that is not an
integer (`FORD`, or a fractional coefficient of a reactant), for example
`biomass/dummy-solid-gas`. OpenSMOKE++ computes such an order with
`pow (C, order)` and no guard. A negative C gives FE_INVALID, and the run
stops with SIGFPE.

The solid reactor above always clips its gas and solid masses. The char
oxidation of the biomass schemes has the order 0.78 for O2. */

#ifndef GAS_RHS_CLIP
# define GAS_RHS_CLIP 0
#endif

void gas_batch_nonisothermal_constantpressure (const double * y, const double dt, double * dy, void * args) {

  UserDataODE data = *(UserDataODE *)args;
  double Temperature = clamp_temperature (y[NGS]);
  double Pressure = clamp_pressure (data.P);

  OpenSMOKE_GasProp_SetTemperature (Temperature);
  OpenSMOKE_GasProp_SetPressure (Pressure);

  // Unpack mass fractions
  double gasmassfracs[NGS], gasmolefracs[NGS];
  for (int jj=0; jj<NGS; jj++)
#if GAS_RHS_CLIP
    gasmassfracs[jj] = y[jj] < 0. ? 0. : y[jj];
#else
    gasmassfracs[jj] = y[jj];
#endif

  // Calculate mole fractions and mixture molecular weight
  double MWMix;
  OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, &MWMix, gasmassfracs);

  // Calculate concentrations and reaction rates
  double ctot = Pressure/(R_GAS*1000*Temperature); // kmol/m3
  double cgas[NGS], rgas[NGS];
  for (int jj=0; jj<NGS; jj++) {
    cgas[jj] = ctot*gasmolefracs[jj];
    rgas[jj] = 0.;
  }
#ifdef VARPROP
  data.rhog = ctot*MWMix;
  data.cpg = OpenSMOKE_GasProp_HeatCapacity (gasmolefracs);
#endif
#if GAS_DRG
  if (gas_drg_active)
    OpenSMOKE_DRG_ReducedFormationRates (cgas, rgas); //[kmol/m3_gas/s]
  else
#endif
  {
    OpenSMOKE_GasProp_ReactionRates (cgas);
    OpenSMOKE_GasProp_FormationRates (rgas); //[kmol/m3_gas/s]
  }

  OpenSMOKE_OpticallyThinProperties otp;
  otp.T = Temperature;
  otp.P = Pressure;
  otp.x = gasmolefracs;

  double QRgas = OpenSMOKE_GasProp_HeatRelease (rgas);

  for (int jj=0; jj<NGS; jj++) {
    dy[jj] = gas_MWs[jj]*rgas[jj]/data.rhog;
    if (data.sources)
      data.sources[jj] = dy[jj]*data.rhog; // source term for gas expansion
  }

  //Temperature equation
  dy[NGS] = (QRgas +  divq_rad (&otp))/(data.rhog*data.cpg);
  if (data.sources)
    data.sources[NGS] = dy[NGS]*data.rhog*data.cpg;
#ifdef TURN_OFF_HEAT_OF_REACTION
  dy[NGS] *= 0.;
#endif
}
