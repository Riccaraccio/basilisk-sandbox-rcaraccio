/**
## Update Properties

Update the thermodynamic properties for the multicomponent phase
change model, and compute the lagrangian derivative of the density,
which is used as a sorce term for the velocity divergence, to
describe low Mach compressibility effects. */

#ifdef VARPROP

enum solid_thermal_conductivity_model {
  L_CONST,
  L_CORBETTA,
  L_HUANG,
  L_ANCACOUCE,
  L_KK
};

enum solid_thermal_conductivity_model lambdaSmodel;

scalar rhoGv0[];
extern scalar porosity;
scalar DTDtS[], DTDtG[];
scalar * DYDtG = NULL;

void update_properties_initial (void) {
  foreach() {
    ThermoState tsGh, tsSh;
    double Diff_coeff[NGS];
    // Update the properties of the gas phase
    double xG[NGS];
    double MWmixGh;
    OpenSMOKE_MoleFractions_From_MassFractions(xG, &MWmixGh, gas_start);
    MWmixG[] = MWmixGh;

    tsGh.T = TG0 * (1. - f[]) + TS0 * f[];
    tsGh.P = Pref;
    tsGh.x = xG;

    rhoGv[] = tpG.rhov(&tsGh);
    cpGv[] = tpG.cpv(&tsGh);
    muGv[] = tpG.muv(&tsGh);
    if (f[] > F_ERR)
      muGv[] = muGv[] / (porosity[] / f[]);

    lambdaGv[] = tpG.lambdav(&tsGh);
    tpG.diff(&tsGh, Diff_coeff);
#ifdef MASS_DIFFUSION_ENTHALPY
    double cpG[NGS];
    tpG.cpvs(&tsGh, cpG);
    for (int jj = 0; jj < NGS; jj++) {
      scalar cpGv = cpGList[jj];
      cpGv[] = cpG[jj];
    }
#endif // MASS_DIFFUSION_ENTHALPY

    for (int jj = 0; jj < NGS; jj++) {
      scalar DmixGv = DmixGList[jj];
#ifdef CONST_DIFF
      DmixGv[] = CONST_DIFF;
#else
      DmixGv[] = Diff_coeff[jj];
      if (f[] > F_ERR)
        DmixGv[] *= pow(porosity[] / f[], 4. / 3.); // effect of solid
#endif
    }

    if (f[] > F_ERR) {
      //Update the properties of the solid phase
      double xS[NSS];
      double MWmixS;
      OpenSMOKE_SolidMoleFractions_From_SolidMassFractions (xS, &MWmixS, sol_start);

      tsSh.T = TS0;
      tsSh.P = Pref;
      tsSh.x = xS;

      rhoSv[] = tpS.rhov (&tsSh);
      cpSv[] = tpS.cpv (&tsSh);

      switch (lambdaSmodel) {
        case L_CONST:
          foreach_dimension()
            lambda1v.x[] = (1. - porosity[] / f[]) * lambdaS + porosity[] / f[] * lambdaGv[];
          break;

        case L_CORBETTA: {
          double char_cond = 0.1405;
          double bio_cond = 0.1937;
          scalar char_field = YSList[OpenSMOKE_IndexOfSolidSpecies("CHAR")];
          double char_fraction = char_field[] / f[];
          foreach_dimension()
              lambda1v.x[] = char_cond * char_fraction + bio_cond * (1. - char_fraction) + porosity[]/f[] * lambdaGv[];
          break;
        }
        
        case L_HUANG: {
          double char_cond = 0.071;
          double bio_cond = 0.21;
          scalar char_field = YSList[OpenSMOKE_IndexOfSolidSpecies("CHAR")];
          double char_fraction = char_field[] / f[];
          foreach_dimension()
              lambda1v.x[] = (char_cond * char_fraction + bio_cond * (1. - char_fraction))*(1. - porosity[]/f[])
                              + 13.5 * 5.67e-8 * pow(TS[]/f[], 3) * 80e-06 / 0.95 + porosity[]/f[] * lambdaGv[];
          break;
        }

        case L_ANCACOUCE: {
          double char_cond= 0.125;
          double bio_cond = 0.056 + 2.6e-4 * tsSh.T;
          scalar char_field = YSList[OpenSMOKE_IndexOfSolidSpecies("CHAR")];
          double char_fraction = char_field[] / f[];

          foreach_dimension()
              lambda1v.x[] = char_cond * char_fraction + bio_cond * (1. - char_fraction);
          break;
        }

        case L_KK: {
          double lS_per = 0.430;
          double lS_par = 0.766;
          double leff_per = 1 / ((1. - porosity[] / f[]) / lS_per + porosity[] / f[] / lambdaGv[]);
          double leff_par = (1. - porosity[] / f[]) * lS_par + porosity[] / f[] * lambdaGv[];
          // longitudinal direction theta = 1.0
          lambda1v.x[] = leff_par;
          // transversal direction theta = 0.58
          lambda1v.y[] = 0.58 * leff_par + (1. - 0.58) * leff_per;
          break;
        }

        default:
          fprintf(stderr, "WARNING: Unknown solid thermal conductivity model\n");
          break;
      }
    }
  }
}

trace
void update_properties (void) {

  foreach()
    rhoGv0[] = rhoGv[];

  // Reset all the properties fields
  reset ({rhoGv, rhoSv, muGv,
          lambdaGv, lambdaSv, MWmixG,
          cpGv, cpSv}, 0.);
  reset (DmixGList, 0.);

  // Update the gas properties
  foreach() {
    ThermoState tsGh, tsSh;
    double Diff_coeff[NGS];
    double xG[NGS], yG[NGS];
    double MWmixGh;
    // Update internal gas properties
    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList[jj];
      yG[jj] = YG[];
    }
    OpenSMOKE_MoleFractions_From_MassFractions (xG, &MWmixGh, yG);
    MWmixG[] = MWmixGh;

    tsGh.T = T[];
    tsGh.P = Pref+p[];
    tsGh.x = xG;

    rhoGv[] = tpG.rhov (&tsGh);
    cpGv[] = tpG.cpv (&tsGh);
    lambdaGv[] = tpG.lambdav (&tsGh);
    muGv[] = tpG.muv (&tsGh);
    if (f[] > F_ERR)
      muGv[] = muGv[] / (porosity[] / f[]);
    tpG.diff (&tsGh, Diff_coeff);
#ifdef MASS_DIFFUSION_ENTHALPY
    double cpG[NGS];
    tpG.cpvs (&tsGh, cpG);
    for(int jj=0; jj<NGS; jj++) {
      scalar cpGv = cpGList[jj];
      cpGv[] = cpG[jj];
    }
#endif // MASS_DIFFUSION_ENTHALPY

    for(int jj=0; jj<NGS; jj++) {
      scalar DmixGv = DmixGList[jj];
      #ifdef CONST_DIFF
      DmixGv[] = CONST_DIFF;
      #else
      DmixGv[] = Diff_coeff[jj];
      #endif
      if (f[] > F_ERR)
        DmixGv[] *= pow(porosity[]/f[], 4./3.); // effect of solid
    }
 
    if (f[] > F_ERR) {
      // Update internal solid properties
      double xS[NSS], yS[NSS];
      double MWmixS;
      for (int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        yS[jj] = YS[]/f[];
      }
      OpenSMOKE_SolidMoleFractions_From_SolidMassFractions (xS, &MWmixS, yS);

      tsSh.T = TS[]/f[];
      tsSh.P = Pref+p[];
      tsSh.x = xS;

      rhoSv[] = tpS.rhov (&tsSh);
      cpSv[] = tpS.cpv (&tsSh);
    
      switch (lambdaSmodel) {
        case L_CONST:
          foreach_dimension()
              lambda1v.x[] = (1. - porosity[] / f[]) * lambdaS + porosity[] / f[] * lambdaGv[];
          break;

        case L_CORBETTA: {
          double char_cond = 0.1405;
          double bio_cond = 0.1937;
          scalar char_field = YSList[OpenSMOKE_IndexOfSolidSpecies("CHAR")];
          double char_fraction = char_field[] / f[];
          foreach_dimension()
              lambda1v.x[] = char_cond * char_fraction + bio_cond * (1. - char_fraction) + porosity[]/f[] * lambdaGv[];
          break;
        }
        
        case L_HUANG: {
          double char_cond = 0.071;
          double bio_cond = 0.21;
          scalar char_field = YSList[OpenSMOKE_IndexOfSolidSpecies("CHAR")];
          double char_fraction = char_field[] / f[];
          foreach_dimension()
              lambda1v.x[] = (char_cond * char_fraction + bio_cond * (1. - char_fraction))*(1. - porosity[]/f[])
                              + 13.5 * 5.67e-8 * pow(TS[]/f[], 3) * 80e-06 / 0.95 + porosity[]/f[] * lambdaGv[];
          break;
        }

        case L_ANCACOUCE: {
          double char_cond = 0.125;
          double bio_cond = 0.056 + 2.6e-4 * tsSh.T;
          scalar char_field = YSList[OpenSMOKE_IndexOfSolidSpecies("CHAR")];
          double char_fraction = char_field[] / f[];
          foreach_dimension()
              lambda1v.x[] = char_cond* char_fraction + bio_cond * (1. - char_fraction);
          break;
        }

        case L_KK: {
          double lS_per = 0.430;
          double lS_par = 0.766;
          double leff_per = 1 / ((1. - porosity[] / f[]) / lS_per + porosity[] / f[] / lambdaGv[]);
          double leff_par = (1. - porosity[] / f[]) * lS_par + porosity[] / f[] * lambdaGv[];
          // longitudinal direction theta = 1.0
          lambda1v.x[] = leff_par;
          // transversal direction theta = 0.58
          lambda1v.y[] = 0.58 * leff_par + (1. - 0.58) * leff_per;
          break;
        }

        default:
          fprintf(stderr, "WARNING: Unknown solid thermal conductivity model\n");
          break;
      }
    }
  }
}

event init (i = 0) //Should be done in the default event but is executed before OS++ initialization otherwise
{
  update_properties_initial();
  
  DYDtG = NULL;

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "DYDtG_%s", OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    a.nodump = true;
    DYDtG = list_append (DYDtG, a);
  }
  reset (DYDtG, 0.);

  MWmixG.dirty = true;

#if TREE
  for (scalar s in {drhodt}) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
    s.dirty = true; // boundary conditions need to be updated
  }
#endif
}

event cleanup (t = end)
{
  delete (DYDtG), free (DYDtG), DYDtG = NULL;
}

event reset_sources (i++) {
  foreach() {
    DTDtG[] = 0.;
    DTDtS[] = 0.;
  }

  reset (DYDtG, 0.);
}

// event properties (i++) {
//   update_properties();
// }

void update_divergence (void) {

//   // ENSURE THAT THE TRACER FORM IS LOST
//   /**
//   We define the variables used to compute the lagrangian derivative
//   on each level. */

  restriction ({T,TS,TG});
  restriction (YSList);
  restriction (YGList);
#ifdef MOLAR_DIFFUSION
  restriction (XGList);
#endif

  /**
  We calculate the Lagrangian derivative of the temperature fields. */

  face vector lambdagradTS[], lambdagradTG[];
  foreach_face() {
    lambdagradTS.x[] = face_value(lambda1v.x, 0)*face_gradient_x (TS, 0)*fm.x[]*fsS.x[];
    lambdagradTG.x[] = face_value(lambda2v.x, 0)*face_gradient_x (TG, 0)*fm.x[]*fsG.x[];
  }

  foreach() {
    foreach_dimension()
      DTDtS[] += (lambdagradTS.x[1] - lambdagradTS.x[])/Delta;
    DTDtS[] += sST[];

    foreach_dimension()
      DTDtG[] += (lambdagradTG.x[1] - lambdagradTG.x[])/Delta;
    DTDtG[] += sGT[];
  }

  /**
  We calculate the Lagrangian derivative for the chemical species mass
  fractions. */ 

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList[jj];
    scalar DmixGv = DmixGList[jj];
    scalar DYDtGjj = DYDtG[jj];

    face vector rhoDmixYGjj[];
    foreach_face() {
      double rhoGf = face_value (rhoGv, 0);
      double DmixGf = face_value (DmixGv, 0);
      rhoDmixYGjj.x[] = rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
    }

    scalar sgexp = sGexpList[jj];

    foreach() {
      foreach_dimension()
        DYDtGjj[] += (rhoDmixYGjj.x[1] - rhoDmixYGjj.x[])/Delta;
      DYDtGjj[] += sgexp[];
    }
  }

  /**
  We add diffusion correction contributions to the chemical species
  mass fraction derivatives. */

  face vector phicGtot[];
  foreach_face() {
    phicGtot.x[] = 0.;
#ifdef FICK_CORRECTED
    for (int jj=0; jj<NGS; jj++) {
      scalar DmixGv = DmixGList[jj];

      double rhoGf = face_value (rhoGv, 0);
      double DmixGf = face_value (DmixGv, 0);
# ifdef MOLAR_DIFFUSION
      double MWmixGf = face_value (MWmixG, 0);

      scalar XG = XGList[jj];
      phicGtot.x[] += (MWmixGf > 0.) ?
        rhoGf*DmixGf*gas_MWs[jj]/MWmixGf*face_gradient_x (XG, 0)*fm.x[]*fsG.x[] : 0.;
# else
      scalar YG = YGList[jj];
      phicGtot.x[] += rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
# endif // MOLAR_DIFFUSION
    }
#endif  // FICK_CORRECTED
  }

  for (int jj=0; jj<NGS; jj++) {
    face vector phicGjj[];
    foreach_face() {
      phicGjj.x[] = phicGtot.x[];
#ifdef MOLAR_DIFFUSION
      scalar DmixGv = DmixGList[jj];

      double rhoGf = face_value (rhoGv, 0);
      double DmixGf = face_value (DmixGv, 0);
      double MWmixGf = face_value (MWmixG, 0);

      phicGjj.x[] -= (MWmixGf > 0.) ?
        rhoGf*DmixGf/MWmixGf*face_gradient_x (MWmixG, 0)*fm.x[]*fsG.x[] : 0.;
#endif

      scalar YG = YGList[jj];
      phicGjj.x[] *= face_value (YG, 0);
    }

    scalar DYDtGjj = DYDtG[jj];

    foreach()
      foreach_dimension()
        DYDtGjj[] -= (phicGjj.x[1] - phicGjj.x[])/Delta;
  }

  // We calculate the one-field divergence by volume-averaging the liquid and the
  // gas-phase contributions.

  foreach() {
    double divu1 = 0., divu2 = 0., divu = 0.;

    // Add internal gas temperature contribution
    divu1 += (TS[]*rhoGv[]*cpGv[] > 0.) ?
      1./(TS[]*(rhoGv[]*cpGv[]*porosity[]/f[] + rhoSv[]*cpSv[]*(1-porosity[]/f[])))*DTDtS[] : 0.;

    // Add external gas temperature contribution
    divu2 += (TG[]*rhoGv[]*cpGv[] > 0.) ?
      1./(TG[]*rhoGv[]*cpGv[])*DTDtG[] : 0.;

    // Add gas chemical species contribution
    double divu1species = 0.;
    for (int jj=0; jj<NGS; jj++) {
      scalar DYDtGjj = DYDtG[jj];
      divu1species += 1./gas_MWs[jj]*DYDtGjj[];
    }
    divu += (rhoGv[] > 0.) ? MWmixG[]/rhoGv[]*divu1species : 0.;

    // Volume averaged contributions
    drhodt[] = divu1*f[] + divu2*(1. - f[]) + divu;

    // Adjust sign for internal convention
    drhodt[] *= -1.;
  }
}

#if 0
void update_divergence_density_u (void) {

  scalar rhot[];
  foreach()
    rhot[] = rhoGv_S[]*f[] + rhoGv_G[]*(1. - f[]);
  
  scalar eps[];
  foreach()
    eps[] = porosity[] + (1. - f[]);
  
  scalar DrhoDt[];
  foreach()
    DrhoDt[] = (rhot[] - rhoGv_S0[])*eps[]/dt;
  
  vector grho[];
  gradients ({rhot}, {grho});

  foreach()
    foreach_dimension()
      DrhoDt[] += u.x[]*grho.x[];


  foreach(){
    DrhoDt[] = DrhoDt[]*cm[];

    double one_over_rho = (rhot[] > 0.) ? 1./rhot[] : 0.;

    if (iter > 1) {
      drhodt[] = DrhoDt[]*one_over_rho;
    }
  }
}

void update_divergence_density_uf (void) {
  
  scalar rhot[];
  foreach()
    rhot[] = rhoGv_S[]*f[] + rhoGv_G[]*(1. - f[]);
  
  scalar eps[];
  foreach()
    eps[] = porosity[] + (1. - f[]);
  scalar DrhoDt[];

  foreach()
    DrhoDt[] = (rhot[] - rhoGv_S0[])*eps[]/dt;
  
  face vector rhoGflux[];
  tracer_fluxes (rhot, uf, rhoGflux, dt, zeroc);

  foreach() 
    foreach_dimension()
      DrhoDt[] += (rhoGflux.x[1] - rhoGflux.x[] - rhot[]*(uf.x[1] - uf.x[]))/Delta;

  foreach(){
    DrhoDt[] = DrhoDt[]*cm[];

    double one_over_rho = (rhot[] > 0.) ? 1./rhot[] : 0.;

    if (iter > 1) {
      drhodt[] = DrhoDt[]*one_over_rho;
    }
  }
}

void update_divergence_density (void) {
  // Choose one of the two methods to compute drhodt
  update_divergence_density_u();
  // update_divergence_density_uf();
}
#endif
#endif // VARPROP
