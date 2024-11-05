/*
** TO BE REVISED **
Solves the system of ODEs for the chemical species reactions

The system of ODEs is solved using the *batch* method.
*/

//#include "radiation.h" not included for now

#define R_GAS 8.31446261815324

typedef struct {
  double rhos;
  double cps;
  double cpg;
  double P;
  double T;
  double f;
  double zeta;
  double* sources;
} UserDataODE;


void batch_isothermal_constantpressure (const double* y, const double dt, double* dy, void* args) {

  /**
  Unpack data for the ODE system. */

  UserDataODE data = *(UserDataODE *)args;
  double rhos = data.rhos;
  double Pressure = data.P;
  double z = data.zeta;
  double* sources = data.sources;
  double f = data.f;

  double epsilon = y[NGS+NSS];
  double Temperature = data.T;

  epsilon = clamp(epsilon, 0., 1.);

  OpenSMOKE_SolProp_SetTemperature (Temperature);
  OpenSMOKE_SolProp_SetPressure (Pressure);

  //////////////////////////////////////////////////////////////////////////

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

  OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, gas_MWs, gasmassfracs);
  double MWMix = OpenSMOKE_MolecularWeight_From_MassFractions (gasmassfracs);

  double ctot = Pressure/(R_GAS*1000.)/Temperature; //ideal gases
  double cgas[NGS], rgas[NGS];
  for (int jj=0; jj<NGS; jj++) {
    cgas[jj] = ctot*gasmolefracs[jj];
  }
  //double rhog = ctot*MWMix;
  double rhog = rhoG;

  //////////////////////////////////////////////////////////////////////////

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

  //////////////////////////////////////////////////////////////////////////

  OpenSMOKE_SolProp_KineticConstants ();
  OpenSMOKE_SolProp_ReactionRates (cgas,csolid);
  OpenSMOKE_SolProp_FormationRates (rgas, rsolid); //rates are given per m3 of solid

  for (int jj=0; jj<NGS; jj++) {
    dy[jj] = gas_MWs[jj]*rgas[jj]*(1-epsilon)*f;
  }

  for (int jj=0; jj<NSS; jj++) {
    dy[jj+NGS] = sol_MWs[jj]*rsolid[jj]*(1-epsilon)*f;
  }

  double totsolidreaction = 0.;
  for (int jj=0; jj<NSS; jj++) {
    totsolidreaction += (sol_MWs[jj]*rsolid[jj]);
  }

  //epsilon equation
  dy[NGS+NSS] = -totsolidreaction*(1-epsilon)*(1-z)/rhos;
  sources[NGS+NSS] = -totsolidreaction;
}

/**
## *batch_nonisothermal_constantpressure()*: solve chemical species reactions

* *y*: vector (length = NS+1) containing the initial values of the system
* *dt*: simulation time step
* *dy*: right-hand-side of the batch system of equations
* *args*: structure with additional arguments to be passed to the system
*/

 void batch_nonisothermal_constantpressure (const double * y, const double dt, double * dy, void * args) {

  /**
  Unpack data for the ODE system. */

  UserDataODE data = *(UserDataODE *)args;
  double rhos = data.rhos;
  double cps = data.cps;
  double cpg = data.cpg;
  double Pressure = data.P;
  double z = data.zeta;
  double f = data.f;
  double* sources = data.sources;

  double epsilon = y[NGS+NSS];
  double Temperature = y[NGS+NSS+1];

  epsilon = clamp(epsilon, 0., 1.);

  /**
  Set temperature and pressure of the system,
  for the calculation of the kinetic constants. */

  OpenSMOKE_SolProp_SetTemperature (Temperature);
  OpenSMOKE_SolProp_SetPressure (Pressure);

  //////////////////////////////////////////////////////////////////////////

  double gasmass[NGS]; double totgasmass = 0.;
  for (int jj=0; jj<NGS; jj++) {
    gasmass[jj] = y[jj];
    gasmass[jj] = (gasmass[jj] < 0.) ? 0. : gasmass[jj];
    totgasmass += gasmass[jj];
  }

  double gasmassfracs[NGS], gasmolefracs[NGS];
  for (int jj=0; jj<NGS; jj++) {
    gasmassfracs[jj] = gasmass[jj]/totgasmass;
  }

  OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, gas_MWs, gasmassfracs);
  double MWMix = OpenSMOKE_MolecularWeight_From_MassFractions (gasmassfracs);

  double ctot = Pressure/(R_GAS*1000.)/Temperature; //ideal gases
  double cgas[NGS], rgas[NGS];
  for (int jj=0; jj<NGS; jj++) {
    cgas[jj] = ctot*gasmolefracs[jj];
    rgas[jj] = 0.;
  }
  //double rhog = ctot*MWMix;
  double rhog = rhoG;

  //////////////////////////////////////////////////////////////////////////

  double solidmass[NSS]; double totsolidmass = 0.;
  for (int jj=0; jj<NSS; jj++) {
    solidmass[jj] = y[jj+NGS];
    solidmass[jj] = (solidmass[jj] < 0.) ? 0. : solidmass[jj];
    totsolidmass += solidmass[jj];
  }

  double solmassfracs[NSS];
  double csolid[NSS], rsolid[NSS];
  for (int jj=0; jj<NSS; jj++) {
    solmassfracs[jj] = (totsolidmass < F_ERR) ? 0. : solidmass[jj]/totsolidmass;
    csolid[jj] = rhos*solmassfracs[jj]/sol_MWs[jj];
  }
  //////////////////////////////////////////////////////////////////////////

  OpenSMOKE_SolProp_KineticConstants ();
  OpenSMOKE_SolProp_ReactionRates (cgas,csolid);
  OpenSMOKE_SolProp_FormationRates (rgas, rsolid); //rates are given per m3 of solid

  for (int jj=0; jj<NGS; jj++) {
    dy[jj] = gas_MWs[jj]*rgas[jj]*(1-epsilon)*f;
  }

  for (int jj=0; jj<NSS; jj++) {
    dy[jj+NGS] = sol_MWs[jj]*rsolid[jj]*(1-epsilon)*f;
  }

  double totsolidreaction = 0.;
  for (int jj=0; jj<NSS; jj++) {
    totsolidreaction += (sol_MWs[jj]*rsolid[jj]);
  }

  //epsilon equation
  dy[NGS+NSS] = -totsolidreaction*(1-epsilon)*(1-z)/rhos;
  sources[0] = -totsolidreaction;

  /**
  Get the heat of reaction and compute the equation for the
  temperature. We add non-linear contributions such as the
  heat dissipated for radiation */

  double QRgas = OpenSMOKE_GasProp_HeatRelease (rgas); //should be 0 for now
  double QRsol = OpenSMOKE_SolProp_HeatRelease (rgas, rsolid);

  //dy[OpenSMOKE_NumberOfSpecies()] = (QR + divq_rad (&otp))/rho/cp; //NO RADIATION FOR NOW
  dy[NGS+NSS+1] = (QRsol*(1-epsilon)*f + QRgas*(1-f +epsilon*f))/(rhog*cpg*epsilon + rhos*cps*(1-epsilon));
}
