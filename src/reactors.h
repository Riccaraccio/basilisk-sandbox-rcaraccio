/*
Solves the system of ODEs for the chemical species reactions

The system of ODEs is solved using the *batch* method.
*/

//#include "radiation.h" not included for now

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


void solid_batch_isothermal_constantpressure (const double* y, const double dt, double* dy, void* args) {

  /**
  Unpack data for the ODE system. */

  UserDataODE data = *(UserDataODE *)args;
  double rhos = data.rhos;
  double Pressure = data.P;
  double z = data.zeta;
  double* sources = data.sources;

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
  double MWmix;
  OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, &MWmix, gasmassfracs);

  double ctot = Pressure/R_GAS/Temperature; //ideal gases, mol/m3
  double cgas[NGS], rgas[NGS];
  for (int jj=0; jj<NGS; jj++) {
    cgas[jj] = ctot*gasmolefracs[jj];
  }

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
    dy[jj] = gas_MWs[jj]*rgas[jj]*(1-epsilon);
  }

  for (int jj=0; jj<NSS; jj++) {
    dy[jj+NGS] = sol_MWs[jj]*rsolid[jj]*(1-epsilon);
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

void solid_batch_nonisothermal_constantpressure (const double * y, const double dt, double * dy, void * args) {

  UserDataODE data = *(UserDataODE *)args;
  double epsilon = y[NGS+NSS];
  double Temperature = y[NGS+NSS+1];

  epsilon = clamp(epsilon, 0., 1.);

  OpenSMOKE_SolProp_SetTemperature (Temperature);
  OpenSMOKE_SolProp_SetPressure (data.P);
  OpenSMOKE_GasProp_SetTemperature (Temperature);
  OpenSMOKE_GasProp_SetPressure (data.P);

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

  double ctot = data.P/(R_GAS*1000*Temperature); // kmol/m3
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

#ifdef GAS_PHASE_REACTIONS
  OpenSMOKE_GasProp_ReactionRates (cgas);
  OpenSMOKE_GasProp_FormationRates (rgas_pure); //[kmol/m3_gas/s]
  QRgas = OpenSMOKE_GasProp_HeatRelease (rgas_pure);
#endif

  for (int jj=0; jj<NGS; jj++) {
    dy[jj] = gas_MWs[jj]* (rgas[jj]*(1-epsilon) + rgas_pure[jj]*epsilon);
    // data.sources[jj] = dy[jj]*epsilon; //to be checked
  }

  for (int jj=0; jj<NSS; jj++) {
    dy[jj+NGS] = sol_MWs[jj]*rsolid[jj]*(1-epsilon);
    // data.sources[jj+NGS] = dy[jj+NGS]*(1-epsilon); //to be checked
  }

  double totsolreaction = 0.;
  for (int jj=0; jj<NSS; jj++)
    totsolreaction += sol_MWs[jj]*rsolid[jj];

  //epsilon equation
  dy[NGS+NSS] = -totsolreaction*(1-epsilon)*(1-data.zeta)/data.rhos;
  data.sources[NGS+NSS] = -totsolreaction;

  //Temperature equation
  double QRsol = OpenSMOKE_SolProp_HeatRelease (rgas, rsolid);

  dy[NGS+NSS+1] = (QRsol*(1-epsilon) + QRgas*epsilon)/((data.rhog*data.cpg*epsilon) + data.rhos*data.cps*(1-epsilon));
  data.sources[NGS+NSS+1] = QRsol*(1-epsilon) + QRgas*epsilon;
  #ifdef TURN_OFF_HEAT_OF_REACTION
  dy[NGS+NSS+1] *= 0.;
  data.sources[NGS+NSS+1] *= 0.;
  #endif
}

void gas_batch_nonisothermal_constantpressure (const double * y, const double dt, double * dy, void * args) {

  /**
  Unpack data for the ODE system. */

  UserDataODE data = *(UserDataODE *)args;
  double rhog = data.rhog;
  double cpg = data.cpg;
  double Pressure = data.P;

  double Temperature = y[NGS];

  OpenSMOKE_GasProp_SetTemperature (Temperature);
  OpenSMOKE_GasProp_SetPressure (Pressure);

  double gasmassfracs[NGS], gasmolefracs[NGS];
  for (int jj=0; jj<NGS; jj++)
    gasmassfracs[jj] = y[jj] < 0. ? 0. : y[jj];

  double MWMix;
  OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, &MWMix, gasmassfracs);


  double ctot = Pressure/(R_GAS*1000)/Temperature; // kmol/m3
  double cgas[NGS], rgas[NGS];
  for (int jj=0; jj<NGS; jj++) {
    cgas[jj] = ctot*gasmolefracs[jj];
    rgas[jj] = 0.;
  }
#ifdef VARPROP
  rhog = ctot*MWMix;
  cpg = OpenSMOKE_GasProp_HeatCapacity (gasmolefracs);
#endif
  OpenSMOKE_GasProp_KineticConstants ();
  OpenSMOKE_GasProp_ReactionRates (cgas);
  OpenSMOKE_GasProp_FormationRates (rgas); //[kmol/m3_gas/s]

  double QRgas = OpenSMOKE_GasProp_HeatRelease (rgas);

  for (int jj=0; jj<NGS; jj++)
    dy[jj] = gas_MWs[jj]*rgas[jj];

  //Temperature equation
  dy[NGS] = QRgas/(rhog*cpg);
  double* sources = data.sources;
  sources[NGS] = QRgas;
#ifdef TURN_OFF_HEAT_OF_REACTION
  dy[NGS] *= 0.;
#endif
}
