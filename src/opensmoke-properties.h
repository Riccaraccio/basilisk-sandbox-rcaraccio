/**
# OpenSMOKE++ Properties

We compute the material properties of a mixture using the
OpenSMOKE++ library.
*/

#include "opensmoke.h"
#include "var-prop.h"

/**
## Properties Functions

Functions for the update of the density, given the thermodynamic
state.
*/

/**
### *opensmoke_gasprop_density()*: gas phase density according to the ideal gas low
*/

double opensmoke_gasprop_density (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  double MWmix = OpenSMOKE_MolecularWeight_From_MoleFractions (ts->x);
  return OpenSMOKE_GasProp_Density_IdealGas (ts->T, ts->P, MWmix);
}

/**
### *opensmoke_gasprop_viscosity()*: gas phase dynamic viscosity
*/

double opensmoke_gasprop_viscosity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_GasProp_DynamicViscosity (ts->x);
}

/**
### *opensmoke_gasprop_thermalconductivity()*: gas phase thermal conductivity
*/

double opensmoke_gasprop_thermalconductivity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_GasProp_ThermalConductivity (ts->x);
}

/**
### *opensmoke_gasprop_heatcapacity()*: gas phase specific heat capacity
*/

double opensmoke_gasprop_heatcapacity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_GasProp_HeatCapacity (ts->x);
}

/**
### *opensmoke_gasprop_heatcapacity_species()*: gas phase species heat capacity
*/

double opensmoke_gasprop_heatcapacity_species (void * p, int i) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  double Cps[OpenSMOKE_NumberOfSpecies()];
  OpenSMOKE_GasProp_HeatCapacity_PureSpecies (Cps);
  return Cps[i];
}

/**
### *opensmoke_gasprop_diff()*: diffusion coefficient of a species in gas phase
*/

double opensmoke_gasprop_diff (void * p, int i) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_GasProp_SetTemperature (ts->T);
  OpenSMOKE_GasProp_SetPressure (ts->P);
  return OpenSMOKE_GasProp_Dmix (ts->x, i);
}

/**
### *opensmoke_liqprop_heatcapacity()*: liquid phase mixture specific heat capacity
*/

double opensmoke_solprop_heatcapacity (void * p) {
  ThermoState * ts = (ThermoState *)p;
  OpenSMOKE_SolProp_SetTemperature (ts->T);
  OpenSMOKE_SolProp_SetPressure (ts->P);
  return OpenSMOKE_SolPropHeatCapacity (ts->x);
}

/**
### *opensmoke_solprop_density*: dummy function, rhos is constant and is not read from OS
*/

double opensmoke_solprop_density (void * p, double rho) {
  return rho;
}

/**
### *opensmoke_solprop_thermalconductivity()*: dummy function, lambdas is constant and is not read from OS
*/

double opensmoke_solprop_thermalconductivity (void * p, double lambda) {
  return lambda;
}

/**
## Thermodynamic Properties

We create the instance of two structures with the
thermodynamic properties, *tp1* for the liquid phase
and *tp2* for the gas phase. The same nomenclature is used
for the thermodynamic states.
*/

ThermoProps tp1, tp2;
ThermoState ts1, ts2;

/**
## Initialization

We set the thermodynamic properties function pointers
to the specific opensmoke functions declared above.
*/

event defaults (i = 0) {

  /**
  We set the thermodynamic properties functions to the
  correct opensmoke functions that compute material
  properties. */

  tp1.cpv     = opensmoke_solprop_heatcapacity;
  tp1.rhov    = opensmoke_solprop_density;
  tp1.lambdav = opensmoke_solprop_thermalconductivity;

  tp2.rhov    = opensmoke_gasprop_density;
  tp2.muv     = opensmoke_gasprop_viscosity;
  tp2.lambdav = opensmoke_gasprop_thermalconductivity;
  tp2.cpv     = opensmoke_gasprop_heatcapacity;
  tp2.diff    = opensmoke_gasprop_diff;
  tp2.cps     = opensmoke_gasprop_heatcapacity_species;
}

