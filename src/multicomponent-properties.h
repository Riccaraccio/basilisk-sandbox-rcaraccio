/**
## Update Properties

Update the thermodynamic properties for the multicomponent phase
change model, and compute the lagrangian derivative of the density,
which is used as a sorce term for the velocity divergence, to
describe low Mach compressibility effects. */

#ifndef T_PROP
# define T_PROP 1e-5
// # define T_PROP 0.1
#endif

#ifdef VARPROP

extern scalar porosity;
scalar DTDt1[], DTDt2[];
scalar * DYDt1 = NULL;    // [NSS]
scalar * DYDt2 = NULL;    // [NGS]

scalar Hcheck[];
scalar betaexpG_G[], betaexpG_S[];
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
    if (f[] < 1. - T_PROP) {

      // Update the properties of the external gas phase
      tsGh.T = TG0;
      tsGh.P = Pref;
      tsGh.x = gas_start;

      rhoGv_G[] = tpG.rhov (&tsGh);
      //rhoGvInt_G[] = rhoGv_G[];
      rhoGv0_G[] = rhoGv_G[];
      muGv_G[] = tpG.muv (&tsGh);
      cpGv_G[] = tpG.cpv (&tsGh);
      lambdaGv_G[] = tpG.lambdav (&tsGh);

      for(int jj=0; jj<NGS; jj++) {
        scalar DmixGv = DmixGList_G[jj];
        #ifdef CONST_DIFF
        DmixGv[] = CONST_DIFF;
        #else
        DmixGv[] = tpG.diff (&tsGh, jj);
        #endif
      }
    }

    if (f[] > T_ERR) {
      // Update the properties of the internal gas phase
      tsGh.T = TS0;
      tsGh.P = Pref;
      tsGh.x = gas_start;

      rhoGv_S[] = tpG.rhov (&tsGh);
      //rhoGvInt_S[] = rhoGv_S[];
      rhoGv0_S[] = rhoGv_S[];
      muGv_S[] = tpG.muv (&tsGh);
      cpGv_S[] = tpG.cpv (&tsGh);
      lambdaGv_S[] = tpG.lambdav (&tsGh);

      for(int jj=0; jj<NGS; jj++) {
        scalar DmixGv = DmixGList_S[jj];
        #ifdef CONST_DIFF
        DmixGv[] = CONST_DIFF;
        #else
        DmixGv[] = tpG.diff (&tsGh, jj);
        DmixGv[] *= pow(porosity[]/f[], 3./2.); //effect of solid, to be revised
        #endif
      }

      //Update the properties of the solid phase
      tsSh.T = TS0;
      tsSh.P = Pref;
      tsSh.x = sol_start;

      rhoSv[] = tpS.rhov (&tsSh);
      //rhoSvInt[] = rhoSv[];
      rhoSv0[] = rhoSv[];
      cpSv[] = tpS.cpv (&tsSh);
      lambdaSv[] = tpS.lambdav (&tsSh);
    }
  }

  boundary ({rhoSv, rhoSv0, cpSv, lambdaSv,
            rhoGv_G, rhoGv0_G, muGv_G, cpGv_G, lambdaGv_G,
            rhoGv_S, rhoGv0_S, muGv_S, cpGv_S, lambdaGv_S});
  boundary (DmixGList_G);
  boundary (DmixGList_S);
}

event defaults (i = 0)
{
  // for (int jj=0; jj<NSS; jj++) { //TODO: check if this is needed for the solid phase
  //   scalar a = new scalar;
  //   free (a.name);
  //   char name[20];
  //   sprintf (name, "DYDt1_%s", liq_species[jj]);
  //   a.name = strdup (name);
  //   a.nodump = true;
  //   DYDt1 = list_append (DYDt1, a);
  // }

  // for (int jj=0; jj<NGS; jj++) {
  //   scalar a = new scalar;
  //   free (a.name);
  //   char name[20];
  //   sprintf (name, "DYDtG_%s", OpenSMOKE_NamesOfSpecies(jj));
  //   a.name = strdup (name);
  //   a.nodump = true;
  //   DYDt2 = list_append (DYDt2, a);
  // }
}

event init (i = 0)
{
  // update_properties_initial(); //TODO: it seems that this is not needed, update_properties works fine

#if TREE
  //for (scalar s in {drhodt, drhodtext}) { //TODO: check if this is correct
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
  // delete (DYDt1), free (DYDt1), DYDt1 = NULL;
  // delete (DYDt2), free (DYDt2), DYDt2 = NULL;
}

event reset_sources (i++)
{
  foreach() {
    // DTDt1[] = 0.;
    DTDt2[] = 0.;

    // for (int jj=0; jj<NLS; jj++) {
    //   scalar DYDt1jj = DYDt1[jj];
    //   DYDt1jj[] = 0.;
    // }
    // for (int jj=0; jj<NGS; jj++) {
    //   fprintf(stderr, "DYDt2[%d] = %g\n", jj, DYDt2[jj]);
    //   scalar DYDt2jj = DYDt2[jj];
    //   DYDt2jj[] = 0.;
    // }
  }
}

void update_properties (void) {
  foreach() {
    rhoSv0[] = rhoSv[];
    rhoGv0_G[] = rhoGv_G[];
    rhoGv0_S[] = rhoGv_S[];
  }

  foreach() {
    if (f[] > T_PROP) {
      // Update internal solid properties
      double xS[NSS], yS[NSS];
      foreach_elem (YSList, jj) {
        scalar YS = YSList[jj];
        yS[jj] = (NSS == 1) ?  1. : YS[]/f[];
      }
      double MWmixS;
      OpenSMOKE_SolidMoleFractions_From_SolidMassFractions (xS, &MWmixS, yS);

      ThermoState tsSh;
      tsSh.T = TS[]/f[];
      tsSh.P = Pref;
      tsSh.x = xS;

      rhoSv[] = tpS.rhov (&tsSh);
      cpSv[] = tpS.cpv (&tsSh);
      lambdaSv[] = tpS.lambdav (&tsSh);

      // Update internal gas properties
      double xG[NGS], yG[NGS];
      foreach_elem (YGList_S, jj) {
        scalar YG = YGList_S[jj];
        yG[jj] = YG[]/f[];
      }
      double MWmixG;
      OpenSMOKE_MoleFractions_From_MassFractions (xG, &MWmixG, yG);

      ThermoState tsGh;
      tsGh.T = TS[]/f[];
      tsGh.P = Pref;
      tsGh.x = xG;

      rhoGv_S[] = tpG.rhov (&tsGh);
      cpGv_S[] = tpG.cpv (&tsGh);
      lambdaGv_S[] = tpG.lambdav (&tsGh);
      muGv_S[] = tpG.muv (&tsGh);
      betaexpG_S[] = gasprop_thermal_expansion (&tsGh);

      for(int jj=0; jj<NGS; jj++) {
        scalar DmixGv = DmixGList_S[jj];
        DmixGv[] = tpG.diff (&tsGh, jj);
        DmixGv[] *= pow(porosity[]/f[], 3./2.);
      }
    }

    if ((1. - f[]) > T_PROP) {
      // Update external gas properties
      double xG[NGS], yG[NGS];
      foreach_elem (YGList_G, jj) {
        scalar YG = YGList_G[jj];
        yG[jj] = YG[]/(1. - f[]);
      }
      double MWmixG;
      OpenSMOKE_MoleFractions_From_MassFractions (xG, &MWmixG, yG);

      ThermoState tsGh;
      tsGh.T = TG[]/(1. - f[]);
      tsGh.P = Pref;
      tsGh.x = xG;

      rhoGv_G[] = tpG.rhov (&tsGh);
      muGv_G[] = tpG.muv (&tsGh);
      cpGv_G[] = tpG.cpv (&tsGh);
      lambdaGv_G[] = tpG.lambdav (&tsGh);
      betaexpG_G[] = gasprop_thermal_expansion (&tsGh);

      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix2v = DmixGList_G[jj];
        Dmix2v[] = tpG.diff (&tsGh, jj);
      }
    }
  }

  boundary ({rhoSv, rhoSv0, cpSv, lambdaSv,
            rhoGv_G, rhoGv0_G, muGv_G, cpGv_G, lambdaGv_G,
            rhoGv_S, rhoGv0_S, muGv_S, cpGv_S, lambdaGv_S});
  boundary (DmixGList_G);
  boundary (DmixGList_S);
}

event properties (i++) {
  update_properties();
}

// void update_divergence (void) {

//   /**
//   We define the variables used to compute the lagrangian derivative
//   on each level. */

//   restriction ({T,TS,TG});
//   restriction (YSList);
//   restriction (YGList_G);
//   restriction (YGList_S);

//   scalar dYdt[];
//   foreach()
//     dYdt[] = 0.;

//   face vector phicGtot[];
//   foreach_face() {
//     phicGtot.x[] = 0.;
//     for (int jj=0; jj<NGS; jj++) {
//       scalar Dmix2v = Dmix2List[jj];
//       double rho2f = 0.5*(rho2v[] + rho2v[-1]);
//       double Dmix2f = 0.5*(Dmix2v[] + Dmix2v[-1]);
// #ifdef FICK_CORRECTED
// # ifdef MOLAR_DIFFUSION
//       scalar XG = XGList[jj];
//       double MW2mixf = 0.5*(MW2mix[] + MW2mix[-1]);
//       phicGtot.x[] += (MW2mixf > 0.) ?
//         rho2f*Dmix2f*inMW[jj]/MW2mixf*face_gradient_x (XG, 0)*fsG.x[]*fm.x[] : 0.;
// # else
//       scalar YG = YGList[jj];
//       phicGtot.x[] += rho2f*Dmix2f*face_gradient_x (YG, 0)*fsG.x[]*fm.x[];
// #endif  // MOLAR_DIFFUSION
// #else
//       phicGtot.x[] = 0.;
// #endif  // FICK_CORRECTED
//     }
//   }

//   for (int jj=0; jj<NGS; jj++) {
//     face vector phicGjj[];
//     scalar YG = YGList[jj];
//     //scalar DYDt2jj = DYDt2[jj];
//     scalar Dmix2v = Dmix2List[jj];
//     foreach_face() {
//       double rho2f = 0.5*(rho2v[] + rho2v[-1]);
//       double Dmix2f = 0.5*(Dmix2v[] + Dmix2v[-1]);
// #ifdef MOLAR_DIFFUSION
//       scalar XG = XGList[jj];
//       double MW2mixf = 0.5*(MW2mix[] + MW2mix[-1]);
//       phicGjj.x[] = (MW2mixf > 0.) ?
//         rho2f*Dmix2f*inMW[jj]/MW2mixf*face_gradient_x (XG, 0)*fsG.x[]*fm.x[] : 0.;
// #else
//       phicGjj.x[] = rho2f*Dmix2f*face_gradient_x (YG, 0)*fsG.x[]*fm.x[];
// #endif  // MOLAR_DIFFUSION

//       double YGf = 0.5*(YG[] + YG[-1]);
//       phicGjj.x[] -= YGf*phicGtot.x[];
//     }

//     scalar sgexp = sgexpList[jj];
//     scalar sgimp = sgimpList[jj];

//     foreach() {
//       foreach_dimension()
//         dYdt[] += (phicGjj.x[1] - phicGjj.x[])/Delta;
//       dYdt[] += (sgexp[] + sgimp[]*YG[]);
//       dYdt[] *= 1./inMW[jj];
//     }
//   }

//   scalar DrhoDt1[], DrhoDt2[];

//   foreach() {

//     double DYDt2sum = dYdt[];
//     for (int jj=0; jj<NGS; jj++) {
//       scalar DYDt2jj = DYDt2[jj];
//       DYDt2sum += 1./inMW[jj]*DYDt2jj[];
//     }
//     DYDt2sum *= (rho2v[] > 0.) ? MW2mix[]/rho2v[] : 0.;

//     // Compute temperature contribution
//     double laplT1 = 0.;
//     foreach_dimension() {
//       double lambdafr = 0.5*(lambda1v[1] + lambda1v[]);
//       double lambdafl = 0.5*(lambda1v[] + lambda1v[-1]);
//       laplT1 += (fm.x[1]*fsL.x[1]*lambdafr*face_gradient_x (TL, 1) -
//           fm.x[]*fsL.x[]*lambdafl*face_gradient_x (TL, 0));
//     }
//     laplT1 /= Delta;

//     double laplT2 = 0.;
//     foreach_dimension() {
//       double lambdafr = 0.5*(lambda2v[1] + lambda2v[]);
//       double lambdafl = 0.5*(lambda2v[] + lambda2v[-1]);
//       laplT2 += (fm.x[1]*fsG.x[1]*lambdafr*face_gradient_x (TG, 1) -
//           fm.x[]*fsG.x[]*lambdafl*face_gradient_x (TG, 0));
//     }
//     laplT2 /= Delta;

//     DTDt1[] += (laplT1 + slT[]);
//     DTDt2[] += (laplT2 + sgT[]);

//     //double DrhoDt1 = 0.;
//     //double DrhoDt2 = 0.;
//     DrhoDt1[] = 0.;
//     DrhoDt2[] = 0.;

//     // Add liquid compressibility due to temperature
//     DrhoDt1[] += (rho1v[]*cp1v[] > 0.) ?
//       -betaexp1[]/(rho1v[]*cp1v[])*DTDt1[] : 0.;

//     DrhoDt2[] += (TG[]*rho2v[]*cp2v[] > 0.) ?
//       -1./(TG[]*rho2v[]*cp2v[])*DTDt2[] : 0.;

//     // Add gas compressibility due to composition
//     //DrhoDt2 += ((1. - f[]) > F_ERR) ? -DYDt2sum : 0.;
//     DrhoDt2[] += (f[] == 0.) ? -DYDt2sum : 0.;

//     //drhodt[] = DrhoDt1*f[] + DrhoDt2*(1. - f[])*(f[] < F_ERR);
//     //drhodtext[] = DrhoDt1;
//   }
//   //shift_field (DrhoDt1, f, 1);
//   //shift_field (DrhoDt2, f, 0);

//   foreach() {
//     drhodt[] = DrhoDt1[]*f[] + DrhoDt2[]*(1. - f[]);
//     drhodtext[] = DrhoDt1[];
//   }
//   boundary ({drhodt, drhodtext});
// }

// void update_divergence_density (void) {
//   vector grho1[], grho2[];
//   gradients ({rho1v, rho2v}, {grho1, grho2});

//   scalar DrhoDt1[], DrhoDt2[];
//   foreach() {
//     DrhoDt1[] = (rho1v[] - rho1v0[])/dt;
//     DrhoDt2[] = (rho2v[] - rho2v0[])/dt;

//     foreach_dimension() {
// #ifdef VELOCITY_JUMP
//       DrhoDt1[] += u1.x[]*grho1.x[];
//       DrhoDt2[] += u2.x[]*grho2.x[];
// #else
//       DrhoDt1[] += uext.x[]*grho1.x[];
//       DrhoDt2[] += u.x[]*grho2.x[];
// #endif
//     }

//     DrhoDt1[] = DrhoDt1[]*cm[];
//     DrhoDt2[] = DrhoDt2[]*cm[];

//     double one_over_rho1 = (rho1v[] > 0.) ? 1./rho1v[] : 0.;
//     double one_over_rho2 = (rho2v[] > 0.) ? 1./rho2v[] : 0.;

//     if (iter > 1) {
//       drhodt[] = (one_over_rho1*DrhoDt1[]*f[] + one_over_rho2*DrhoDt2[]*(1. - f[]));
//       drhodtext[] = one_over_rho1*DrhoDt1[]*f[];
//     }
//   }
// }

#endif