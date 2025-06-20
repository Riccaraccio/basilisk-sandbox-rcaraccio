/**
## Update Properties

Update the thermodynamic properties for the multicomponent phase
change model, and compute the lagrangian derivative of the density,
which is used as a sorce term for the velocity divergence, to
describe low Mach compressibility effects. */

#ifndef T_PROP
// # define T_PROP 1e-5
# define T_PROP F_ERR
// # define T_PROP 0.1
#endif

#ifdef VARPROP

enum solid_thermal_conductivity_model {
  L_CONST,
  L_CORBETTA,
  L_ANCACOUCE,
  L_KK
};

enum solid_thermal_conductivity_model lambdaSmodel;

scalar rhoGv_G0[], rhoGv_S0[];
extern scalar porosity;
// scalar DTDtS[], DTDtG[];
// scalar * DYDtG_G = NULL;    // [NSS]
// scalar * DYDtG_S = NULL;    // [NGS]

// scalar Hcheck[];
//scalar rhoSvInt[], rhoGvInt_G[], rhoGvInt_S[];

/*
void update_properties_constant (void) {
  foreach() {
    rhoSv[] = rhoS;
    rhoSv0[] = rhoS;
    rho1vInt[] = rhoS;
    cpSv[] = cpS;
    lambdaSv[] = lambdaS;

    rhoGv_G[] = rhoG;
    rhoGv0_G[] = rhoG;
    muGv_G[] = muG;
    cpGv_G[] = cpG;
    lambdaGv_G[] = lambdaG;

    foreach_elem (Dmix2List_G, jj) {
      scalar DmixGv = DmixGList_G[jj];
      DmixGv[] = DmixG;
    }
    
    rhoGv_S[] = rhoG;
    rhoGv0_S[] = rhoG;
    muGv_S[] = muG;
    cpGv_S[] = cpG;
    lambdaGv_S[] = lambdaG;

    foreach_elem (Dmix2List_S, jj) {
      scalar DmixGv = DmixGList_S[jj];
      DmixGv[] = DmixG;
    }
  }

  boundary ({rhoSv,rhoSv0,rhoSvInt,muSv,cpSv,lambdaSv,
      rhoGv_G, rhoGv0_G, muGv_G, cpGv_G, lambdaGv_G,
      rhoGv_S, rhoGv0_S, muGv_S, cpGv_S, lambdaGv_S});
  boundary (DmixGList_G);
  boundary (DmixGList_S);
}
*/

void update_properties_initial (void) {
  foreach() {
    ThermoState tsGh, tsSh;
    double Diff_coeff[NGS];
    if (f[] > T_ERR) {
      // Update the properties of the internal gas phase
      double xG[NGS];
      double MWmixG;
      OpenSMOKE_MoleFractions_From_MassFractions (xG, &MWmixG, gas_start);
      MWmixG_S[] = MWmixG;

      tsGh.T = TS0;
      tsGh.P = Pref;
      tsGh.x = xG;

      rhoGv_S[] = tpG.rhov (&tsGh);
      cpGv_S[] = tpG.cpv (&tsGh);
      muGv_S[] = tpG.muv (&tsGh);
      lambdaGv_S[] = tpG.lambdav (&tsGh);
      tpG.diff (&tsGh, Diff_coeff);
#ifdef MASS_DIFFUSION_ENTHALPY
      double cpG[NGS];
      tpG.cpvs (&tsGh, cpG);
      for(int jj=0; jj<NGS; jj++) {
        scalar cpGv = cpGList_S[jj];
        cpGv[] = cpG[jj];
      }
#endif // MASS_DIFFUSION_ENTHALPY

      for(int jj=0; jj<NGS; jj++) {
        scalar DmixGv = DmixGList_S[jj];
        #ifdef CONST_DIFF
        DmixGv[] = CONST_DIFF;
        #else
        DmixGv[] = Diff_coeff[jj];
        DmixGv[] *= pow(porosity[]/f[], 4./3.); //effect of solid, to be revised
        #endif
      }
      
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
            lambda1v.x[] = (1. - porosity[] / f[]) * lambdaS + porosity[] / f[] * lambdaGv_S[];
          break;

        case L_CORBETTA: {
          double char_cond = 0.1405;
          double bio_cond = 0.1937;
          scalar char_field = YSList[OpenSMOKE_IndexOfSolidSpecies("CHAR")];
          double char_fraction = char_field[] / f[];
          foreach_dimension()
              lambda1v.x[] = char_cond * char_fraction + bio_cond * (1. - char_fraction);
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
          double leff_per = 1 / ((1. - porosity[] / f[]) / lS_per + porosity[] / f[] / lambdaGv_S[]);
          double leff_par = (1. - porosity[] / f[]) * lS_par + porosity[] / f[] * lambdaGv_S[];
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
    
    if (f[] < 1. - T_PROP) {
      // Update the properties of the external gas phase
      double xG[NGS];
      double MWmixG;
      OpenSMOKE_MoleFractions_From_MassFractions (xG, &MWmixG, gas_start);
      MWmixG_G[] = MWmixG;

      tsGh.T = TG0;
      tsGh.P = Pref;
      tsGh.x = xG;

      rhoGv_G[] = tpG.rhov (&tsGh);
      muGv_G[] = tpG.muv (&tsGh);
      cpGv_G[] = tpG.cpv (&tsGh);
      lambdaGv_G[] = tpG.lambdav (&tsGh);
      tpG.diff (&tsGh, Diff_coeff);
#ifdef MASS_DIFFUSION_ENTHALPY
      double cpG[NGS];
      tpG.cpvs (&tsGh, cpG);
      for(int jj=0; jj<NGS; jj++) {
        scalar cpGv = cpGList_G[jj];
        cpGv[] = cpG[jj];
        // fprintf(stderr, "cpGv[%d] = %g\n", jj, cpGv[]);
      }
#endif // MASS_DIFFUSION_ENTHALPY

      for(int jj=0; jj<NGS; jj++) {
        scalar DmixGv = DmixGList_G[jj];
        #ifdef CONST_DIFF
        DmixGv[] = CONST_DIFF;
        #else
        DmixGv[] = Diff_coeff[jj];
        #endif
      }
    }
  }
}

trace
void update_properties (void) {

  foreach() {
    // rhoGv_S0[] = rhoGv_S[];
    // rhoGv_G0[] = rhoGv_G[];
    rhoGv_S0[] = rhoGv_S[]*f[] + (1.-f[])*rhoGv_G[]; //field looks nicer done in one field
  }

  // Reset all the properties fields
  // reset ({rhoGv_S, rhoGv_G, rhoSv,
  //         muGv_S, muGv_G,
  //         lambdaGv_S, lambdaGv_G, lambdaSv,
  //         cpGv_S, cpGv_G, cpSv}, 0.);
  // reset (DmixGList_S, 0.);
  // reset (DmixGList_G, 0.);
  // reset ({MWmixG_S, MWmixG_G}, 0.);
  
  const double Pref_const = 1.01325e5; // Pa
  foreach() {

    double f_val = f[];
    double one_minus_f = 1. - f_val;
    ThermoState tsGh, tsSh;
    double xG[NGS], yG[NGS];
    double MWmix;
    double P_local = Pref_const + p[];
    double Diff_coeff[NGS];

    if (f_val > T_PROP) {
      double inv_f = 1./f_val;
      double porosity_val = porosity[]*inv_f;
      double TS_val = TS[]*inv_f;
      double por_f_pow = pow(porosity_val, 4./3.);

      // Update internal gas properties
      foreach_elem (YGList_S, jj) {
        scalar YG = YGList_S[jj];
        yG[jj] = YG[]*inv_f;
      }
      OpenSMOKE_MoleFractions_From_MassFractions (xG, &MWmix, yG);
      MWmixG_S[] = MWmix;

      tsGh.T = TS_val;
      tsGh.P = P_local;
      tsGh.x = xG;

      rhoGv_S[] = tpG.rhov (&tsGh);
      cpGv_S[] = tpG.cpv (&tsGh);
      lambdaGv_S[] = tpG.lambdav (&tsGh);
      muGv_S[] = tpG.muv (&tsGh);
      tpG.diff (&tsGh, Diff_coeff);

#ifdef MASS_DIFFUSION_ENTHALPY
      double cpG[NGS];
      tpG.cpvs (&tsGh, cpG);
      for(int jj=0; jj<NGS; jj++) {
        scalar cpGv = cpGList_S[jj];
        cpGv[] = cpG[jj];
      }
#endif // MASS_DIFFUSION_ENTHALPY

      for(int jj=0; jj<NGS; jj++) {
        scalar DmixGv = DmixGList_S[jj];
        #ifdef CONST_DIFF
        DmixGv[] = CONST_DIFF*por_f_pow;
        #else
        DmixGv[] = Diff_coeff[jj]*por_f_pow;
        #endif
      }
      
      // Update internal solid properties
      double xS[NSS], yS[NSS];
      foreach_elem (YSList, jj) {
        scalar YS = YSList[jj];
        yS[jj] = YS[]*inv_f;
      }
      OpenSMOKE_SolidMoleFractions_From_SolidMassFractions (xS, &MWmix, yS);

      tsSh.T = TS_val;
      tsSh.P = Pref+p[];
      tsSh.x = xS;

      rhoSv[] = tpS.rhov (&tsSh);
      cpSv[] = tpS.cpv (&tsSh);
    
      switch (lambdaSmodel) {
        case L_CONST:
          foreach_dimension()
              lambda1v.x[] = (1. - porosity[] / f[]) * lambdaS + porosity[] / f[] * lambdaGv_S[];
          break;

        case L_CORBETTA: {
          double char_cond = 0.1405;
          double bio_cond = 0.1937;
          scalar char_field = YSList[OpenSMOKE_IndexOfSolidSpecies("CHAR")];
          double char_fraction = char_field[] / f[];
          foreach_dimension()
              lambda1v.x[] = char_cond * char_fraction + bio_cond * (1. - char_fraction);
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
          double leff_per = 1 / ((1. - porosity[] / f[]) / lS_per + porosity[] / f[] / lambdaGv_S[]);
          double leff_par = (1. - porosity[] / f[]) * lS_par + porosity[] / f[] * lambdaGv_S[];
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

    if (one_minus_f > T_PROP) {
      double inv_one_minus_f = 1./one_minus_f;
      double TG_val = TG[]*inv_one_minus_f;

      // Update external gas properties
      foreach_elem (YGList_G, jj) {
        scalar YG = YGList_G[jj];
        yG[jj] = YG[]*inv_one_minus_f;
      }
      OpenSMOKE_MoleFractions_From_MassFractions (xG, &MWmix, yG);
      MWmixG_G[] = MWmix;

      tsGh.T = TG_val;
      tsGh.P = P_local;
      tsGh.x = xG;

      rhoGv_G[] = tpG.rhov (&tsGh);
      muGv_G[] = tpG.muv (&tsGh);
      cpGv_G[] = tpG.cpv (&tsGh);
      lambdaGv_G[] = tpG.lambdav (&tsGh);
      tpG.diff (&tsGh, Diff_coeff);

#ifdef MASS_DIFFUSION_ENTHALPY
      double cpG[NGS];
      tpG.cpvs (&tsGh, cpG);
      for(int jj=0; jj<NGS; jj++) {
        scalar cpGv = cpGList_G[jj];
        cpGv[] = cpG[jj];
      }
#endif // MASS_DIFFUSION_ENTHALPY

      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix2v = DmixGList_G[jj];
# ifdef CONST_DIFF
        Dmix2v[] = CONST_DIFF;
# else
        Dmix2v[] = Diff_coeff[jj];
# endif
      }

      double lambda_g = lambdaGv_G[];
      foreach_dimension()
        lambda2v.x[] = lambda_g;
    }
  }
}

event init (i = 0) //Should be done in the default event but is executed before OS++ initialization otherwise
{
  update_properties_initial();
  
  // DYDtG_G = NULL;
  // DYDtG_S = NULL;

  // for (int jj=0; jj<NGS; jj++) {
  //   scalar a = new scalar;
  //   free (a.name);
  //   char name[20];
  //   sprintf (name, "DYDtG_%s_G", OpenSMOKE_NamesOfSpecies(jj));
  //   a.name = strdup (name);
  //   a.nodump = true;
  //   DYDtG_G = list_append (DYDtG_G, a);
  // }
  // reset (DYDtG_G, 0.);
  
  // for (int jj=0; jj<NGS; jj++) {
  //   scalar a = new scalar;
  //   free (a.name);
  //   char name[20];
  //   sprintf (name, "DYDtG_%s_S", OpenSMOKE_NamesOfSpecies(jj));
  //   a.name = strdup (name);
  //   a.nodump = true;
  //   DYDtG_S = list_append (DYDtG_S, a);
  // }
  // reset (DYDtG_S, 0.);

  MWmixG_G.dirty = true;
  MWmixG_S.dirty = true;
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

// event cleanup (t = end)
// {
//   delete (DYDtG_G), free (DYDtG_G), DYDtG_G = NULL;
//   delete (DYDtG_S), free (DYDtG_S), DYDtG_S = NULL;
// }

// event reset_sources (i++)
// {
//   foreach() {
//     DTDtG[] = 0.;
//     DTDtS[] = 0.;


//     for (int jj=0; jj<NGS; jj++) {
//       scalar DYDtGjj = DYDtG_G[jj];
//       DYDtGjj[] = 0.;
//     }
//     for (int jj=0; jj<NGS; jj++) {
//       scalar DYDtGjj = DYDtG_S[jj];
//       DYDtGjj[] = 0.;
//     }
//   }
// }

// event properties (i++) {
//   update_properties();
// }

// void update_divergence (void) { // Unused, but kept for reference

//   // ENSURE THAT THE TRACER FORM IS LOST
//   /**
//   We define the variables used to compute the lagrangian derivative
//   on each level. */

//   restriction ({T,TS,TG});
//   restriction (YSList);
//   restriction (YGList_G);
//   restriction (YGList_S);

//   /**
//   We calculate the Lagrangian derivative of the temperature fields. */

//   face vector lambdagradTS[], lambdagradTG[];
//   foreach_face() {
//     lambdagradTS.x[] = face_value(lambdaGv_S, 0)*face_gradient_x (TS, 0)*fm.x[]*fsS.x[];
//     lambdagradTG.x[] = face_value(lambdaGv_G, 0)*face_gradient_x (TG, 0)*fm.x[]*fsG.x[];
//   }

//   foreach() {
//     foreach_dimension()
//       DTDtS[] += (lambdagradTS.x[1] - lambdagradTS.x[])/Delta;
//     // DTDtS[] += sST[];

//     foreach_dimension()
//       DTDtG[] += (lambdagradTG.x[1] - lambdagradTG.x[])/Delta;
//     // DTDtG[] += sGT[];
//   }

//   // EXTERNAL GAS PHASE
//   /**
//   We calculate the Lagrangian derivative for the chemical species mass
//   fractions. */ 

//   for (int jj=0; jj<NGS; jj++) {
//     scalar YG = YGList_G[jj];
//     scalar DmixGv = DmixGList_G[jj];
//     scalar DYDtGjj = DYDtG_G[jj];

//     face vector rhoDmixYGjj[];
//     foreach_face() {
//       double rhoGf = face_value (rhoGv_G, 0);
//       double DmixGf = face_value (DmixGv, 0);
//       rhoDmixYGjj.x[] = rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
//     }

//     scalar sgexp = sGexpList[jj];
//     // scalar sgimp = sgimpList[jj]; //not used in my case
//     // scalar YGInt = YGList_Int[jj];

//     foreach() {
//       foreach_dimension()
//         DYDtGjj[] += (rhoDmixYGjj.x[1] - rhoDmixYGjj.x[])/Delta;
//       // DYDtGjj[] += (sgexp[] + sgimp[]*YGInt[]);
//       DYDtGjj[] += sgexp[];
//     }
//   }

//   /**
//   We add diffusion correction contributions to the chemical species
//   mass fraction derivatives. */ //TODO, FICK_CORRECTED is not defined

//   face vector phicGtot[];
//   foreach_face() {
//     phicGtot.x[] = 0.;
// #ifdef FICK_CORRECTED
//     for (int jj=0; jj<NGS; jj++) {
//       scalar Dmix2v = Dmix2List[jj];

//       double rho2f = face_value (rho2v, 0);
//       double Dmix2f = face_value (Dmix2v, 0);
// # ifdef MOLAR_DIFFUSION
//       double MW2mixf = face_value (MW2mix, 0);

//       scalar XG = XGList[jj];
//       phicGtot.x[] += (MW2mixf > 0.) ?
//         rho2f*Dmix2f*inMW[jj]/MW2mixf*face_gradient_x (XG, 0)*fm.x[]*fsG.x[] : 0.;
// # else
//       scalar YG = YGList[jj];
//       phicGtot.x[] += rho2f*Dmix2f*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
// # endif // MOLAR_DIFFUSION
//     }
// #endif  // FICK_CORRECTED
//   }

//   for (int jj=0; jj<NGS; jj++) {
//     face vector phicGjj[];
//     foreach_face() {
//       phicGjj.x[] = phicGtot.x[];
// #ifdef MOLAR_DIFFUSION
//       scalar Dmix2v = Dmix2List[jj];

//       double rho2f = face_value (rho2v, 0);
//       double Dmix2f = face_value (Dmix2v, 0);
//       double MW2mixf = face_value (MW2mix, 0);

//       phicGjj.x[] -= (MW2mixf > 0.) ?
//         rho2f*Dmix2f/MW2mixf*face_gradient_x (MW2mix, 0)*fm.x[]*fsG.x[] : 0.;
// #endif

//       scalar YG = YGList_G[jj];
//       phicGjj.x[] *= face_value (YG, 0);
//     }

//     scalar DYDtGjj = DYDtG_G[jj];

//     foreach()
//       foreach_dimension()
//         DYDtGjj[] -= (phicGjj.x[1] - phicGjj.x[])/Delta;
//   }

//   // INTERNAL GAS PHASE
//   /**
//   We calculate the Lagrangian derivative for the chemical species mass
//   fractions. */ 

//   for (int jj=0; jj<NGS; jj++) {
//     scalar YG = YGList_S[jj];
//     scalar DmixGv = DmixGList_S[jj];
//     scalar DYDtGjj = DYDtG_S[jj];

//     face vector rhoDmixYGjj[];
//     foreach_face() {
//       double rhoGf = face_value (rhoGv_S, 0);
//       double DmixGf = face_value (DmixGv, 0);
//       rhoDmixYGjj.x[] = rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsS.x[];
//     }

//     scalar ssexp = sSexpList[jj];
//     // scalar sgimp = sgimpList[jj]; //not used in my case
//     // scalar YGInt = YGList_Int[jj];

//     foreach() {
//       foreach_dimension()
//         DYDtGjj[] += (rhoDmixYGjj.x[1] - rhoDmixYGjj.x[])/Delta;
//       // DYDtGjj[] += (sgexp[] + sgimp[]*YGInt[]);
//       DYDtGjj[] += ssexp[];
//     }
//   }

//   /**
//   We add diffusion correction contributions to the chemical species
//   mass fraction derivatives. */ //TODO, FICK_CORRECTED is not defined

//   foreach_face() {
//     phicGtot.x[] = 0.;
// #ifdef FICK_CORRECTED
//     for (int jj=0; jj<NGS; jj++) {
//       scalar Dmix2v = DmixGList_S[jj];

//       double rho2f = face_value (rho2v, 0);
//       double Dmix2f = face_value (Dmix2v, 0);
// # ifdef MOLAR_DIFFUSION
//       double MW2mixf = face_value (MW2mix, 0);

//       scalar XG = XGList[jj];
//       phicGtot.x[] += (MW2mixf > 0.) ?
//         rho2f*Dmix2f*inMW[jj]/MW2mixf*face_gradient_x (XG, 0)*fm.x[]*fsG.x[] : 0.;
// # else
//       scalar YG = YGList[jj];
//       phicGtot.x[] += rho2f*Dmix2f*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
// # endif // MOLAR_DIFFUSION
//     }
// #endif  // FICK_CORRECTED
//   }

//   for (int jj=0; jj<NGS; jj++) {
//     face vector phicGjj[];
//     foreach_face() {
//       phicGjj.x[] = phicGtot.x[];
// #ifdef MOLAR_DIFFUSION
//       scalar Dmix2v = DmixGList_S[jj];

//       double rho2f = face_value (rho2v, 0);
//       double Dmix2f = face_value (Dmix2v, 0);
//       double MW2mixf = face_value (MW2mix, 0);

//       phicGjj.x[] -= (MW2mixf > 0.) ?
//         rho2f*Dmix2f/MW2mixf*face_gradient_x (MW2mix, 0)*fm.x[]*fsG.x[] : 0.;
// #endif

//       scalar YG = YGList_S[jj];
//       phicGjj.x[] *= face_value (YG, 0);
//     }

//     scalar DYDtGjj = DYDtG_S[jj];

//     foreach()
//       foreach_dimension()
//         DYDtGjj[] -= (phicGjj.x[1] - phicGjj.x[])/Delta;
//   }
  
//   /**
//   We calculate the one-field divergence by volume-averaging the liquid and the
//   gas-phase contributions. */

//   foreach() {
//     double divu1 = 0., divu2 = 0.;

//     // Add internal gas temperature contribution
//     divu1 += (TS[]*rhoGv_S[]*cpGv_S[] > 0.) ?
//       1./(TS[]*rhoGv_S[]*cpGv_S[])*DTDtS[] : 0.;

//     // Add external gas temperature contribution
//     divu2 += (TG[]*rhoGv_G[]*cpGv_G[] > 0.) ?
//       1./(TG[]*rhoGv_G[]*cpGv_G[])*DTDtG[] : 0.;

//     // Add internal gas chemical species contribution
//     double divu1species = 0.;
//     for (int jj=0; jj<NGS; jj++) {
//       scalar DYDtGjj = DYDtG_S[jj];
//       divu1species += 1./gas_MWs[jj]*DYDtGjj[];
//     }
//     // divu1 += (rhoGv_S[] > 0.) ? MWmixG_S[]/rhoGv_S[]*divu1species : 0.; //TODO

//     // Add external gas chemical species contribution
//     double divu2species = 0.;
//     for (int jj=0; jj<NGS; jj++) {
//       scalar DYDtGjj = DYDtG_G[jj];
//       divu2species += 1./gas_MWs[jj]*DYDtGjj[];
//     }
//     // divu2 += (rhoGv_G[] > 0.) ? MWmixG_G[]/rhoGv_G[]*divu2species : 0.; //TODO
    
//     // Volume averaged contributions
//     drhodt[] = divu1*f[] + divu2*(1. - f[]);

//     // Adjust sign for internal convention
//     drhodt[] *= -1.;
//   }
// }

void update_divergence_density (void) {
  
  // vector grhoG[], grhoS[];
  // gradients ({rhoGv_S, rhoGv_G}, {grhoS, grhoG});

  // scalar DrhoDtS[], DrhoDtG[];
  // foreach() {
  //   DrhoDtS[] = (rhoGv_S[] - rhoGv_S0[])/dt;
  //   DrhoDtG[] = (rhoGv_G[] - rhoGv_G0[])/dt;
  // } 

  // // face vector rhoG_Sflux[], rhoG_Gflux[];
  // // tracer_fluxes (rhoGv_S, uf, rhoG_Sflux, dt, zeroc);
  // // tracer_fluxes (rhoGv_G, uf, rhoG_Gflux, dt, zeroc);

  // // foreach() {
  // //   foreach_dimension() {
  // //     DrhoDtS[] += (rhoG_Sflux.x[] - rhoG_Sflux.x[1] - rhoGv_S[]*(uf.x[] - uf.x[1]))/Delta;
  // //     DrhoDtG[] += (rhoG_Gflux.x[] - rhoG_Gflux.x[1] - rhoGv_G[]*(uf.x[] - uf.x[1]))/Delta;
  // //   }
  // // }

  // foreach() {
  //   foreach_dimension()
  //     DrhoDtS[] += u.x[]*grhoS.x[];
  //     DrhoDtG[] += u.x[]*grhoG.x[];
  // }

  // foreach() {
  //   DrhoDtS[] = DrhoDtS[]*cm[];
  //   DrhoDtG[] = DrhoDtG[]*cm[];

  //   double one_over_rhoS = (rhoGv_S[] > 0.) ? 1./rhoGv_S[] : 0.;
  //   double one_over_rhoG = (rhoGv_G[] > 0.) ? 1./rhoGv_G[] : 0.;

  //   if (iter > 1) {
  //     drhodt[] = (one_over_rhoS*DrhoDtS[]*f[] + one_over_rhoG*DrhoDtG[]*(1. - f[]));
  //   }
  // }
  
  scalar rhot[];
  foreach()
    rhot[] = rhoGv_S[]*f[] + rhoGv_G[]*(1. - f[]);
  
  vector grho[];
  gradients ({rhot}, {grho});
  
  scalar DrhoDt[];
  foreach() {
    DrhoDt[] = (rhot[] - rhoGv_S0[])/dt;

  // //This does not work, idk why
  // face vector rhoGflux[];
  // tracer_fluxes (rhot, uf, rhoGflux, dt, zeroc);

  // foreach()
  //   foreach_dimension()
  //     DrhoDt[] += (rhoGflux.x[] - rhoGflux.x[1] - rhot[]*(uf.x[] - uf.x[1]))/Delta;

    foreach_dimension()
      DrhoDt[] += u.x[]*grho.x[];

    DrhoDt[] = DrhoDt[]*cm[];

    double one_over_rho = (rhot[] > 0.) ? 1./rhot[] : 0.;

    if (iter > 1) {
      drhodt[] = DrhoDt[]*one_over_rho;
    }
  }
}
#endif
