#ifndef MULTICOMPONENT
  #define MULTICOMPONENT 1
#endif

#include "intgrad.h"
#include "fracface.h"

#ifdef EXPLICIT_DIFFUSION
  #include "diffusion-explicit.h"
#else
  #include "diffusion.h"
#endif

#include "common-evaporation.h"
#include "memoryallocation.h"
#include "chemistry.h"
#include "int-temperature.h"
#include "int-concentration.h"

event reset_sources (i++) {
#ifdef SOLVE_TEMPERATURE
  foreach() {
    sST[] = 0.;
    sGT[] = 0.;
  }
#endif

  foreach() {
    for (int jj=0; jj<NGS; jj++) {
      scalar sSexp = sSexpList[jj];
      scalar sGexp = sGexpList[jj];
      sSexp[] = 0.;
      sGexp[] = 0.;
    }
  }
}

extern face vector ufsave;
face vector u_prime[];
#ifndef STOP_TRACER_ADVECTION
event tracer_advection (i++) {
  // lose tracer form and extrapolate fields
  foreach() {
    porosity[] = (f[] > F_ERR) ? porosity[]/f[] : 0.;
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]/f[] : 0.;
    TG[] = ((1.-f[]) > F_ERR) ? TG[]/(1.-f[]) : 0.;

    TS[] = (f[] > F_ERR) ? TS[] : TG[];
    TG[] = (f[] < 1.-F_ERR) ? TG[] : TS[];
#endif

    for (int jj=0; jj<NGS; jj++) { 
      scalar YG_S = YGList_S[jj];
      scalar YG_G = YGList_G[jj];

      YG_S[] = (f[] > F_ERR) ? YG_S[]/f[] : 0.;
      YG_G[] = (f[] < 1.-F_ERR) ? YG_G[]/(1.-f[]) : 0.;

      YG_S[] = (f[] > F_ERR) ? YG_S[] : YG_G[];
      YG_G[] = (f[] < 1.-F_ERR) ? YG_G[] : YG_S[];
    }
  }

  advection_div(YGList_S, ufsave, dt);
  advection_div(YGList_G, ufsave, dt);

#ifdef SOLVE_TEMPERATURE
  foreach_face() {
    double ef = face_value(porosity, 0);
    double ff = face_value(f, 0);
    u_prime.x[] = (ff > F_ERR) ? ufsave.x[]*(rhoG*cpG*ef)/(rhoG*cpG*ef + rhoS*cpS*(1.-ef)) : ufsave.x[];
  }

  advection_div({TS}, u_prime, dt);
# ifndef TEMPERATURE_PROFILE
  advection_div({TG}, u_prime, dt);
# endif
#endif

  // Reset the velocity field, just to be sure
  foreach_face()
    uf.x[] = ufsave.x[];
  
  // recover tracer form
  foreach() {
    porosity[] = (f[] > F_ERR) ? porosity[]*f[] : 0.;
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]*f[] : 0.;
    TG[] = ((1.-f[]) > F_ERR) ? TG[]*(1.-f[]) : 0.;
#endif

    for (int jj=0; jj<NGS; jj++) { 
      scalar YG_S = YGList_S[jj];
      scalar YG_G = YGList_G[jj];

      YG_S[] = (f[] > F_ERR) ? YG_S[]*f[] : 0.;
      YG_G[] = (f[] < 1.-F_ERR) ? YG_G[]*(1.-f[]) : 0.;
    }
  }
}
#endif

event tracer_diffusion (i++) {
  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    fS[] = f[]; fG[] = 1. - f[];
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]/f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]/(1. - f[]) : 0.;
#endif

    for (scalar YG in YGList_S)
      YG[] = (f[] > F_ERR) ? YG[]/f[] : 0.;
    
    for (scalar YG in YGList_G)
      YG[] = ((1. - f[]) > F_ERR) ? YG[]/(1. - f[]) : 0.;
  }

  //Compute face gradients
  face_fraction (fS, fsS);
  face_fraction (fG, fsG);

#ifdef SOLVE_TEMPERATURE
  //interface temperature first guess
  foreach() {
    TInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      TInt[] = avg_neighbor (point, TS, f);
  }

  #ifdef FIXED_INT_TEMP //Force interface temperature = TG0
  foreach()
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      TInt[] = TG0;

  #elif TEMPERATURE_PROFILE
  double tv = TemperatureProfile_GetT(t);
  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      TInt[] = tv;
    
    if (f[] < 1.-F_ERR)
      TG[] = tv;
  }

  #else //default: solve for interface temperature
  ijc_CoupledTemperature();
  #endif
#endif

  // first guess for species interface concentration
  foreach() {
    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList_S[jj];
      scalar YGInt = YGList_Int[jj];
      YGInt[] = 0.;
      if (f[] > F_ERR && f[] < 1.-F_ERR)
        YGInt[] = avg_neighbor (point, YG, f);
        YGInt[] = clamp (YGInt[], 0., 1.);
    }
  }

  //find the interface concentration foreach species
  intConcentration();

  //Calculate the source therm
  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord n = facet_normal (point, fS, fsS), p;
      double alpha = plane_alpha (fS[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      //Solid side
      for (int jj=0; jj<NGS; jj++) {
        scalar YGInt = YGList_Int[jj];
        scalar YG    = YGList_S[jj];
        scalar sSexp = sSexpList[jj];
        scalar Dmix2 = DmixGList_S[jj];

        double bc = YGInt[];
        double Strgrad = ebmgrad (point, YG, fS, fG, fsS, fsG, false, bc, &success);

        double jS = rhoG*Dmix2[]*Strgrad;
#ifdef AXI
        sSexp[] += jS*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
        sSexp[] += jS*area/Delta*cm[];
#endif
      }

      //Gas side
      for (int jj=0; jj<NGS; jj++) {
        scalar YGInt = YGList_Int[jj];
        scalar YG    = YGList_G[jj];
        scalar sGexp = sGexpList[jj];
        scalar Dmix2 = DmixGList_G[jj];

        double bc = YGInt[];
        double Gtrgrad = ebmgrad (point, YG, fS, fG, fsS, fsG, true, bc, &success);

        double jG = rhoG*Dmix2[]*Gtrgrad;
#ifdef AXI
        sGexp[] += jG*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
        sGexp[] += jG*area/Delta*cm[];
#endif
      }

#ifdef SOLVE_TEMPERATURE
      double bc = TInt[];
      double Strgrad = ebmgrad (point, TS, fS, fG, fsS, fsG, false, bc, &success);
      double Gtrgrad = ebmgrad (point, TG, fS, fG, fsS, fsG, true, bc, &success);

      double Sheatflux = lambda1v[]*Strgrad;
      double Gheatflux = lambda2v[]*Gtrgrad;

  #ifdef AXI
      sST[] += Sheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sGT[] += Gheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
  #else
      sST[] += Sheatflux*area/Delta*cm[];
      sGT[] += Gheatflux*area/Delta*cm[];
  #endif
#endif
    }
  }

  //Recover tracer form 
  foreach() {
#ifdef SOLVE_TEMPERATURE
    TS[] *= f[];
    TG[] *= (1. - f[]);
    T[] = TS[] + TG[];
#endif
    for (scalar YG in YGList_S)
      YG[] *= f[];
    for (scalar YG in YGList_G)
      YG[] *= (1. - f[]);
    for (scalar YS in YSList)
      YS[] *= f[];
  }

  scalar theta1[], theta2[];

#if TREE
  theta1.refine = theta1.prolongation = fraction_refine;
  theta2.refine = theta2.prolongation = fraction_refine;
  theta1.dirty = true;
  theta2.dirty = true;
#endif

  //internal diffusion
  for (int jj=0; jj<NGS; jj++) {
    face vector Dmixf[];
    scalar Dmix2 = DmixGList_S[jj];
    foreach_face()
      Dmixf.x[] = face_value(Dmix2, 0)*rhoG*fsS.x[]*fm.x[];

    foreach()
      theta1[] = cm[]*max(fS[]*rhoG, F_ERR);

    scalar YG = YGList_S[jj];
    scalar sSexp = sSexpList[jj];

#ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (YG, dt, D=Dmixf, theta=theta1);
#else
    diffusion (YG, dt, D=Dmixf, r=sSexp, theta=theta1);
#endif
  }

  //external diffusion
  for (int jj=0; jj<NGS; jj++) {
    face vector Dmixf[];
    scalar Dmix2 = DmixGList_G[jj];
    foreach_face()
      Dmixf.x[] = face_value(Dmix2, 0)*rhoG*fsG.x[]*fm.x[];

    foreach()
      theta2[] = cm[]*max(fG[]*rhoG, F_ERR);

    scalar YG = YGList_G[jj];
    scalar sGexp = sGexpList[jj];
#ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (YG, dt, D=Dmixf, theta=theta2);
#else
    diffusion (YG, dt, D=Dmixf, r=sGexp, theta=theta2);
#endif
  }

#ifdef SOLVE_TEMPERATURE
  foreach_face() {
    lambda1f.x[] = face_value(lambda1v, 0)*fsS.x[]*fm.x[];
    lambda2f.x[] = face_value(lambda2v, 0)*fsG.x[]*fm.x[];
  }

  foreach() {
    theta1[] = cm[]*max(fS[]*rhocp1v[], F_ERR);
    theta2[] = cm[]*max(fG[]*rhocp2v[], F_ERR);
  }
  #ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (TS, dt, D=lambda1f, r=sST, theta=theta1);
    diffusion_explicit (TG, dt, D=lambda2f, r=sGT, theta=theta2);
  #else
    diffusion (TS, dt, D=lambda1f, r=sST, theta=theta1);
    diffusion (TG, dt, D=lambda2f, r=sGT, theta=theta2);
  #endif
#endif

  //check mass and recover tracer form
  foreach() {
    double totmassgas = 0.;
    for (scalar YG in YGList_S)
      totmassgas += YG[];

    for (scalar YG in YGList_S) {
      YG[] = (totmassgas > 0.) ? YG[]/totmassgas : 0.;
      YG[] = clamp(YG[], 0., 1.);
      YG[] = (YG[] > 1e-10) ? YG[]*f[] : 0.;
    }
  }

  foreach() {
    double totmassgas = 0.;
    for (scalar YG in YGList_G)
      totmassgas += YG[];

    for (scalar YG in YGList_G) {
      YG[] = (totmassgas > 0.) ? YG[]/totmassgas : 0.;
      YG[] = clamp(YG[], 0., 1.);
      YG[] = (YG[] > 1e-10) ? YG[]*(1.-f[]) : 0.;
    }
  }

  foreach() {
#ifdef SOLVE_TEMPERATURE
    TS[] *= f[];
    TG[] *= (1. - f[]);
    T[] = TS[] + TG[];
#endif
  }
}

// calculate diff coefficents, to be moved to properties (eventually)
event properties (i++) {
#ifdef CONST_DIFF
  // OpenSMOKE_GasProp_SetTemperature (500.);
  // OpenSMOKE_GasProp_SetPressure (101325.);
  // double fixedcomp[NGS];
  // fixedcomp[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  // double Dmixv = OpenSMOKE_GasProp_Dmix(fixedcomp, OpenSMOKE_IndexOfSpecies ("N2"));
  // fprintf(stderr, "Dmixv = %g\n", Dmixv);

  double Dmixv =  2.05e-5; //Diff of CO in N2 at 500K, 1 atm
  //double Dmixv = 0.;

  foreach() {
    if (f[] < 1.-F_ERR) {
    //set the same for all species
      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix2 = DmixGList_G[jj];
        Dmix2[] = Dmixv;
      }
    } 
    if (f[] > F_ERR) {
      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix2 = DmixGList_S[jj];
        Dmix2[] = Dmixv*pow(porosity[], 3./2.); //effect of solid, to be revised
      }
    }
  }
#else
  foreach() {
    if (f[]>F_ERR) { //for the solid phase
      //set T and P
  #ifdef SOLVE_TEMPERATURE
      OpenSMOKE_GasProp_SetTemperature (TS[]);
  #else
      OpenSMOKE_GasProp_SetTemperature (T[]);
  #endif
      OpenSMOKE_GasProp_SetPressure (p[]+Pref);

      // reconstruct the mass fracs
      double gasmassfracs [NGS];
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList_S[jj];
        gasmassfracs[jj] = YG[];
      }

      //convert to mole fracs
      double gasmolefracs [NGS];
      double MWmixG;
      OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, MWmixG, gasmassfracs);

      //calculate diff coeff
      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix = DmixGList_S[jj];
        Dmix[] = OpenSMOKE_GasProp_Dmix(gasmolefracs, jj);
        Dmix[] *= pow(porosity[], 3./2.); //effect of solid, to be revised
      }
    }

    if (f[] < 1-F_ERR) { //for the gas phase
      //set T and P
  #ifdef SOLVE_TEMPERATURE
      OpenSMOKE_GasProp_SetTemperature (TG[]);
  #else
      OpenSMOKE_GasProp_SetTemperature (T[]);
  #endif
      OpenSMOKE_GasProp_SetPressure (p[]+Pref);

      // reconstruct the mass fracs
      double gasmassfracs [NGS];
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList_G[jj];
        gasmassfracs[jj] = YG[];
      }

      //convert to mole fracs
      double gasmolefracs [NGS];
      double MWmixG;
      OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, MWmixG, gasmassfracs);

      //calculate diff coeff
      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix = DmixGList_G[jj];
        Dmix[] = OpenSMOKE_GasProp_Dmix(gasmolefracs, jj);
      }
    }
  }

  boundary (DmixGList_G);
  boundary (DmixGList_S);
#endif
}