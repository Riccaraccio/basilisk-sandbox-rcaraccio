/*
** TO BE REVISED **
Solves the system of ODEs for the chemical species reactions

The system of ODEs is solved using the *batch* method.
*/

//#include "radiation.h" not included for now

#define R_GAS 8.31446261815324

typedef struct {
  double rho;
  double cp;
  double P;
  double T;
  int isinterface;  //1 if interface, 0 if not
  double epsi;      //porosity, used when isinterface = 1
} UserDataODE;


void batch_isothermal_constantpressure (const double * y, const double dt, double * dy, void * args) {
  
  /**
  Unpack data for the ODE system. */

  UserDataODE data = *(UserDataODE *)args;
  double rho = data.rho;
  double cp = data.cp;
  double Pressure = data.P;
  double Temperature = data.T;
  double porosity = y[NGS+NSS];;

  //porosity = (porosity < 0.) ? 0. : porosity;
 
  OpenSMOKE_SolProp_SetTemperature (Temperature);
  OpenSMOKE_SolProp_SetPressure (Pressure);

  //////////////////////////////////////////////////////////////////////////

  double gasmass[NGS];
  for (int jj=0; jj<NGS; jj++) {
    gasmass[jj] = y[jj];
    gasmass[jj] = (gasmass[jj] < 0.) ? 0. : gasmass[jj];
  }

  double totgasmass = 0.;
  for (int jj=0; jj<NGS; jj++) {
    totgasmass += gasmass[jj];
  }

  double gasmassfracs[NGS], gasmolefracs[NGS];
  for (int jj=0; jj<NGS; jj++) {
    gasmassfracs[jj] = gasmass[jj]/totgasmass;
  }
  
 
  OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, gas_MWs, gasmassfracs);
  double MWMix = OpenSMOKE_MolecularWeight_From_MassFractions (gasmassfracs);

  double ctot = Pressure/(R_GAS*1000.)/Temperature;
  double cgas[NGS], rgas[NGS];
  for (int jj=0; jj<NGS; jj++) {
    cgas[jj] = ctot*gasmolefracs[jj];
    rgas[jj] = 0.;
  }
  rho = ctot*MWMix;

  //////////////////////////////////////////////////////////////////////////

  double solidmass[NSS];
  for (int jj=0; jj<NSS; jj++) {
    solidmass[jj] = y[jj+NGS];
    solidmass[jj] = (solidmass[jj] < 0.) ? 0. : solidmass[jj];
  }

  double totsolidmass = 0.;
  for (int jj=0; jj<NSS; jj++) {
    totsolidmass += solidmass[jj];
  }

  double solmassfracs[NSS];
  double csolid[NSS], rsolid[NSS];
  for (int jj=0; jj<NSS; jj++){
    solmassfracs[jj] = (totsolidmass < T_ERR) ? 0. : solidmass[jj]/totsolidmass;
    csolid[jj] = rho1*solmassfracs[jj]/sol_MWs[jj];
  }

  //////////////////////////////////////////////////////////////////////////
 
  OpenSMOKE_SolProp_KineticConstants ();
  OpenSMOKE_SolProp_ReactionRates (cgas,csolid);
  OpenSMOKE_SolProp_FormationRates (rgas, rsolid);
 
  //fprintf(stderr, "CELL = % g\n", solmassfracs[0]);

  for (int jj=0; jj<NGS; jj++) {
    dy[jj] = gas_MWs[jj]*rgas[jj]*porosity;
  }

  for (int jj=0; jj<NSS; jj++) {
    dy[jj+NGS] = sol_MWs[jj]*rsolid[jj]*porosity;
  }

  double totsolidreaction = 0.;
  for (int jj=0; jj<NSS; jj++) {
    totsolidreaction += (sol_MWs[jj]*rsolid[jj]); //epsilon equation
  }

  //totsolidreaction = 1e-2;
  if (data.isinterface == 1) {
    //shrinking
    dy[NGS+NSS] = 0;
    omega = totsolidreaction*porosity;
  } else {
    // porosity loss  
    dy[NGS+NSS] = totsolidreaction*porosity/rho1;
    omega = 0;
  }

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
  double rho = data.rho;
  double cp = data.cp;
  double Pressure = data.P;
  double Temperature = data.T;

  double porosity = y[NGS+NSS];
  
  /**
  Set temperature and pressure of the system,
  for the calculation of the kinetic constants. */
  
  OpenSMOKE_SolProp_SetTemperature (Temperature);
  OpenSMOKE_SolProp_SetPressure (Pressure);

  //////////////////////////////////////////////////////////////////////////

  double gasmass[NGS];
  for (int jj=0; jj<NGS; jj++) {
    gasmass[jj] = y[jj];
    gasmass[jj] = (gasmass[jj] < 0.) ? 0. : gasmass[jj];
  }

  double totgasmass = 0.;
  for (int jj=0; jj<NGS; jj++) {
    totgasmass += gasmass[jj];
  }

  double gasmassfracs[NGS], gasmolefracs[NGS];
  for (int jj=0; jj<NGS; jj++) {
    gasmassfracs[jj] = gasmass[jj]/totgasmass;
  }
  
  OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, gas_MWs, gasmassfracs);
  double MWMix = OpenSMOKE_MolecularWeight_From_MassFractions (gasmassfracs);

  double ctot = Pressure/(R_GAS*1000.)/Temperature;
  double cgas[NGS], rgas[NGS];
  for (int jj=0; jj<NGS; jj++) {
    cgas[jj] = ctot*gasmolefracs[jj];
    rgas[jj] = 0.;
  }
  rho = ctot*MWMix;

  //////////////////////////////////////////////////////////////////////////

  double solidmass[NSS];
  for (int jj=0; jj<NSS; jj++) {
    solidmass[jj] = y[jj+NGS];
    solidmass[jj] = (solidmass[jj] < 0.) ? 0. : solidmass[jj];
  }

  // double totsolidmass = 0.;
  // for (int jj=0; jj<NSS; jj++) {
  //   totsolidmass += solidmass[jj];
  // }

  double solmassfracs[NSS];
  double csolid[NSS], rsolid[NSS];
  for (int jj=0; jj<NSS; jj++){
  // solmassfracs[jj] = solidmass[jj]/totsolidmass;
    csolid[jj] = rho1*solmassfracs[jj]/sol_MWs[jj];
  }

  //////////////////////////////////////////////////////////////////////////

  OpenSMOKE_SolProp_KineticConstants ();
  OpenSMOKE_SolProp_ReactionRates (cgas,csolid);
  OpenSMOKE_SolProp_FormationRates (rgas, rsolid);

  for (int jj=0; jj<NGS; jj++) {
    dy[jj] = gas_MWs[jj]*rgas[jj]*porosity;
  }

  for (int jj=0; jj<NSS; jj++) {
    dy[jj+NGS] = sol_MWs[jj]*rsolid[jj]*porosity;
  }

  for (int jj=0; jj<NGS; jj++) {
    dy[NGS+NSS] += (sol_MWs[jj]*rsolid[jj]*porosity)/rho1; //epsilon equation
  }

  /**
  Get the heat of reaction and compute the equation for the
  temperature. We add non-linear contributions such as the
  heat dissipated for radiation */

  double QRgas = OpenSMOKE_GasProp_HeatRelease (rgas);
  double QRsol = OpenSMOKE_SolProp_HeatRelease (rgas, rsolid);

  //dy[OpenSMOKE_NumberOfSpecies()] = (QR + divq_rad (&otp))/rho/cp; //NO RADIATION FOR NOW
  //dy[NGS+NSS+1] = (QRsol*solvolume + QRgas*(1-solvolume))/(rho1*cp1*solvolume + rho*cp2*(1-solvolume));
 
}
