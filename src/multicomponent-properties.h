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

extern scalar porosity;
scalar DTDtS[], DTDtG[];
scalar * DYDtG_G = NULL;    // [NSS]
scalar * DYDtG_S = NULL;    // [NGS]

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

event init (i = 0) //Should be done in the default event but is executed before OS++ initialization
{
  DYDtG_G = NULL;
  DYDtG_S = NULL;

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "DYDtG_%s_G", OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    a.nodump = true;
    DYDtG_G = list_append (DYDtG_G, a);
  }
  reset (DYDtG_G, 0.);
  
  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "DYDtG_%s_S", OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    a.nodump = true;
    DYDtG_S = list_append (DYDtG_S, a);
  }
  reset (DYDtG_S, 0.);
}

event init (i = 0)
{
  // update_properties_initial(); //TODO: it seems that this is not needed, update_properties works fine

  MWmixG_G.dirty = true;
  MWmixG_S.dirty = true;
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
  delete (DYDtG_G), free (DYDtG_G), DYDtG_G = NULL;
  delete (DYDtG_S), free (DYDtG_S), DYDtG_S = NULL;
}

event reset_sources (i++)
{
  foreach() {
    DTDtG[] = 0.;
    DTDtS[] = 0.;


    for (int jj=0; jj<NGS; jj++) {
      scalar DYDtGjj = DYDtG_G[jj];
      DYDtGjj[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar DYDtGjj = DYDtG_S[jj];
      DYDtGjj[] = 0.;
    }
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
      MWmixG_S[] = MWmixG;

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
        #ifdef CONST_DIFF
        DmixGv[] = CONST_DIFF;
        #endif
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
      MWmixG_G[] = MWmixG;

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
# ifdef CONST_DIFF
        Dmix2v[] = CONST_DIFF;
# endif
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

void update_divergence (void) {

  //TODO ENSURE THAT THE TRACER FORM IS LOST
  /**
  We define the variables used to compute the lagrangian derivative
  on each level. */

  restriction ({T,TS,TG});
  restriction (YSList);
  restriction (YGList_G);
  restriction (YGList_S);

  /**
  We calculate the Lagrangian derivative of the temperature fields. */

  face vector lambdagradTS[], lambdagradTG[];
  foreach_face() {
    double ef = face_value(porosity, 0);
    double lambda1vh = ef*face_value(lambdaGv_S, 0) + (1. - ef)*face_value(lambdaSv, 0);
    double lambda2vh = face_value(lambdaGv_G, 0);
    lambdagradTS.x[] = lambda1vh*face_gradient_x (TS, 0)*fm.x[]*fsS.x[];
    lambdagradTG.x[] = lambda2vh*face_gradient_x (TG, 0)*fm.x[]*fsG.x[];
  }

  foreach() {
    foreach_dimension()
      DTDtS[] += (lambdagradTS.x[1] - lambdagradTS.x[])/Delta;
    DTDtS[] += sST[];

    foreach_dimension()
      DTDtG[] += (lambdagradTG.x[1] - lambdagradTG.x[])/Delta;
    DTDtG[] += sGT[];
  }

  //EXTERNAL GAS PHASE
  /**
  We calculate the Lagrangian derivative for the chemical species mass
  fractions. */ 

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    scalar DmixGv = DmixGList_G[jj];
    scalar DYDtGjj = DYDtG_G[jj];

    face vector rhoDmixYGjj[];
    foreach_face() {
      double rhoGf = face_value (rhoGv_G, 0);
      double DmixGf = face_value (DmixGv, 0);
      rhoDmixYGjj.x[] = rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
    }

    scalar sgexp = sGexpList[jj];
    // scalar sgimp = sgimpList[jj]; //not used in my case
    // scalar YGInt = YGList_Int[jj];

    foreach() {
      foreach_dimension()
        DYDtGjj[] += (rhoDmixYGjj.x[1] - rhoDmixYGjj.x[])/Delta;
      // DYDtGjj[] += (sgexp[] + sgimp[]*YGInt[]);
      DYDtGjj[] += sgexp[];
    }
  }

  /**
  We add diffusion correction contributions to the chemical species
  mass fraction derivatives. */ //TODO, FICK_CORRECTED is not defined

  face vector phicGtot[];
  foreach_face() {
    phicGtot.x[] = 0.;
#ifdef FICK_CORRECTED
    for (int jj=0; jj<NGS; jj++) {
      scalar Dmix2v = Dmix2List[jj];

      double rho2f = face_value (rho2v, 0);
      double Dmix2f = face_value (Dmix2v, 0);
# ifdef MOLAR_DIFFUSION
      double MW2mixf = face_value (MW2mix, 0);

      scalar XG = XGList[jj];
      phicGtot.x[] += (MW2mixf > 0.) ?
        rho2f*Dmix2f*inMW[jj]/MW2mixf*face_gradient_x (XG, 0)*fm.x[]*fsG.x[] : 0.;
# else
      scalar YG = YGList[jj];
      phicGtot.x[] += rho2f*Dmix2f*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
# endif // MOLAR_DIFFUSION
    }
#endif  // FICK_CORRECTED
  }

  for (int jj=0; jj<NGS; jj++) {
    face vector phicGjj[];
    foreach_face() {
      phicGjj.x[] = phicGtot.x[];
#ifdef MOLAR_DIFFUSION
      scalar Dmix2v = Dmix2List[jj];

      double rho2f = face_value (rho2v, 0);
      double Dmix2f = face_value (Dmix2v, 0);
      double MW2mixf = face_value (MW2mix, 0);

      phicGjj.x[] -= (MW2mixf > 0.) ?
        rho2f*Dmix2f/MW2mixf*face_gradient_x (MW2mix, 0)*fm.x[]*fsG.x[] : 0.;
#endif

      scalar YG = YGList_G[jj];
      phicGjj.x[] *= face_value (YG, 0);
    }

    scalar DYDtGjj = DYDtG_G[jj];

    foreach()
      foreach_dimension()
        DYDtGjj[] -= (phicGjj.x[1] - phicGjj.x[])/Delta;
  }

  //INTERNAL GAS PHASE
  /**
  We calculate the Lagrangian derivative for the chemical species mass
  fractions. */ 

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_S[jj];
    scalar DmixGv = DmixGList_S[jj];
    scalar DYDtGjj = DYDtG_S[jj];

    face vector rhoDmixYGjj[];
    foreach_face() {
      double rhoGf = face_value (rhoGv_S, 0);
      double DmixGf = face_value (DmixGv, 0);
      rhoDmixYGjj.x[] = rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsS.x[];
    }

    scalar ssexp = sSexpList[jj];
    // scalar sgimp = sgimpList[jj]; //not used in my case
    //scalar YGInt = YGList_Int[jj];

    foreach() {
      foreach_dimension()
        DYDtGjj[] += (rhoDmixYGjj.x[1] - rhoDmixYGjj.x[])/Delta;
      // DYDtGjj[] += (sgexp[] + sgimp[]*YGInt[]);
      DYDtGjj[] += ssexp[];
    }
  }

  /**
  We add diffusion correction contributions to the chemical species
  mass fraction derivatives. */ //TODO, FICK_CORRECTED is not defined

  foreach_face() {
    phicGtot.x[] = 0.;
#ifdef FICK_CORRECTED
    for (int jj=0; jj<NGS; jj++) {
      scalar Dmix2v = DmixGList_S[jj];

      double rho2f = face_value (rho2v, 0);
      double Dmix2f = face_value (Dmix2v, 0);
# ifdef MOLAR_DIFFUSION
      double MW2mixf = face_value (MW2mix, 0);

      scalar XG = XGList[jj];
      phicGtot.x[] += (MW2mixf > 0.) ?
        rho2f*Dmix2f*inMW[jj]/MW2mixf*face_gradient_x (XG, 0)*fm.x[]*fsG.x[] : 0.;
# else
      scalar YG = YGList[jj];
      phicGtot.x[] += rho2f*Dmix2f*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
# endif // MOLAR_DIFFUSION
    }
#endif  // FICK_CORRECTED
  }

  for (int jj=0; jj<NGS; jj++) {
    face vector phicGjj[];
    foreach_face() {
      phicGjj.x[] = phicGtot.x[];
#ifdef MOLAR_DIFFUSION
      scalar Dmix2v = DmixGList_S[jj];

      double rho2f = face_value (rho2v, 0);
      double Dmix2f = face_value (Dmix2v, 0);
      double MW2mixf = face_value (MW2mix, 0);

      phicGjj.x[] -= (MW2mixf > 0.) ?
        rho2f*Dmix2f/MW2mixf*face_gradient_x (MW2mix, 0)*fm.x[]*fsG.x[] : 0.;
#endif

      scalar YG = YGList_S[jj];
      phicGjj.x[] *= face_value (YG, 0);
    }

    scalar DYDtGjj = DYDtG_S[jj];

    foreach()
      foreach_dimension()
        DYDtGjj[] -= (phicGjj.x[1] - phicGjj.x[])/Delta;
  }
  
  /**
  We calculate the one-field divergence by volume-averaging the liquid and the
  gas-phase contributions. */

  foreach() {
    double divu1 = 0., divu2 = 0.;

    // Add internal gas temperature contribution
    divu1 += (TS[]*rhoGv_S[]*cpGv_S[] > 0.) ?
      1./(TS[]*rhoGv_S[]*cpGv_S[])*DTDtS[] : 0.;

    // Add external gas temperature contribution
    divu2 += (TG[]*rhoGv_G[]*cpGv_G[] > 0.) ?
      1./(TG[]*rhoGv_G[]*cpGv_G[])*DTDtG[] : 0.;

    // Add internal gas chemical species contribution
    double divu1species = 0.;
    for (int jj=0; jj<NGS; jj++) {
      scalar DYDtGjj = DYDtG_S[jj];
      divu1species += 1./gas_MWs[jj]*DYDtGjj[];
    }
    divu1 += (rhoGv_S[] > 0.) ? MWmixG_S[]/rhoGv_S[]*divu1species : 0.;

    // Add external gas chemical species contribution
    double divu2species = 0.;
    for (int jj=0; jj<NGS; jj++) {
      scalar DYDtGjj = DYDtG_G[jj];
      divu2species += 1./gas_MWs[jj]*DYDtGjj[];
    }
    divu2 += (rhoGv_G[] > 0.) ? MWmixG_G[]/rhoGv_G[]*divu2species : 0.;
    
    // Volume averaged contributions
    drhodt[] = divu1*f[] + divu2*(1. - f[]);

    // Adjust sign for internal convention
    drhodt[] *= -1.;
  }
}

// void update_divergence_density (void) { //TODO
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