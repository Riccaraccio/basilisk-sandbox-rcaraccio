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
#include "memoryallocation-varprop.h"
#include "multicomponent-properties.h"
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

    double rhoGvh_S, rhoSvh;
    double cpGvh_S, cpSvh;
    
    #ifdef VARPROP
    rhoGvh_S = face_value(rhoGv_S, 0); rhoSvh = face_value(rhoSv, 0);
    cpGvh_S = face_value(cpGv_S, 0); cpSvh = face_value(cpSv, 0);
    #else
    rhoGvh_S = rhoG; rhoSvh = rhoS;
    cpGvh_S = cpG; cpSvh = cpS;
    #endif

    u_prime.x[] = (ff > F_ERR) ? ufsave.x[]*(rhoGvh_S*cpGvh_S)/(rhoGvh_S*cpGvh_S*ef + rhoSvh*cpSvh*(1.-ef)) : ufsave.x[];
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

  update_properties();

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

  //Can be removed for performance
  //Check the mass fractions
  foreach() {
    double sum = 0.;
    for (scalar YG in YGList_S)
      sum += YG[];

    for (scalar YG in YGList_S) {
      YG[] = (sum > 1e-10) ? YG[]/sum : 0.;
      YG[] = clamp(YG[], 0., 1.);
    }
  }

  foreach() {
    double sum = 0.;
    for (scalar YG in YGList_G)
      sum += YG[];

    for (scalar YG in YGList_G) {
      YG[] = (sum > 1e-10) ? YG[]/sum : 0.;
      YG[] = clamp(YG[], 0., 1.);
    }
  }

  foreach() {
    double sum = 0.;
    for (scalar YS in YSList)
      sum += YS[];
    
    for (scalar YS in YSList) {
      YS[] = (sum > 1e-10) ? YS[]/sum : 0.;
      YS[] = clamp(YS[], 0., 1.);
    }
  }

  //Compute face gradients
  face_fraction (fS, fsS);
  face_fraction (fG, fsG);

#ifdef SOLVE_TEMPERATURE
  //interface temperature first guess
  foreach() {
    TInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      TInt[] = f[]*TS[] + (1.-f[])*TG[];
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
      scalar YG_S = YGList_S[jj];
      scalar YG_G = YGList_G[jj];
      scalar YGInt = YGList_Int[jj];
      YGInt[] = 0.;
      if (f[] > F_ERR && f[] < 1.-F_ERR)
        YGInt[] = f[]*YG_S[] + (1.-f[])*YG_G[];
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
        scalar DmixG = DmixGList_S[jj];

        double bc = YGInt[];
        double Strgrad = ebmgrad (point, YG, fS, fG, fsS, fsG, false, bc, &success);

        double rhoGvh_S;
        #ifdef VARPROP
        rhoGvh_S = rhoGv_S[];
        #else
        rhoGvh_S = rhoG;
        #endif

        double jS = rhoGvh_S*DmixG[]*Strgrad;
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
        scalar DmixG = DmixGList_G[jj];

        double bc = YGInt[];
        double Gtrgrad = ebmgrad (point, YG, fS, fG, fsS, fsG, true, bc, &success);

        double rhoGvh_G;
#ifdef VARPROP
        rhoGvh_G = rhoGv_G[];
#else
        rhoGvh_G = rhoG;
#endif

        double jG = rhoGvh_G*DmixG[]*Gtrgrad;
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

      double Sheatflux = (n.x*lambda1v.x[] + n.y*lambda1v.y[])*Strgrad;
      double Gheatflux = (n.x*lambda2v.x[] + n.y*lambda2v.y[])*Gtrgrad;

# ifdef AXI
      sST[] += Sheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sGT[] += Gheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
# else
      sST[] += Sheatflux*area/Delta*cm[];
      sGT[] += Gheatflux*area/Delta*cm[];
# endif
#endif
    }
  }

#if defined VARPROP && !defined NO_EXPANSION
  update_divergence();
#endif

  scalar theta1[], theta2[];

#if TREE
  theta1.refine = theta1.prolongation = fraction_refine;
  theta2.refine = theta2.prolongation = fraction_refine;
  theta1.dirty = true;
  theta2.dirty = true;
#endif

  // Internal gas diffusion
  for (int jj=0; jj<NGS; jj++) {
    face vector DmixGf[];
    scalar DmixG = DmixGList_S[jj];
    foreach_face() {
      double rhoGvh_S;
#ifdef VARPROP
      rhoGvh_S = face_value(rhoGv_S, 0);
#else
      rhoGvh_S = rhoG;
#endif
      DmixGf.x[] = face_value(DmixG, 0)*rhoGvh_S*fsS.x[]*fm.x[];
    }

    foreach() {
#ifdef VARPROP
      theta1[] = cm[]*max(fS[]*rhoGv_S[], F_ERR);
#else
      theta1[] = cm[]*max(fS[]*rhoG, F_ERR);
#endif
    }

    scalar YG = YGList_S[jj];
    scalar sSexp = sSexpList[jj];

#ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (YG, dt, D=DmixGf, theta=theta1);
#else
    diffusion (YG, dt, D=DmixGf, r=sSexp, theta=theta1);
#endif
  }

  //external diffusion
  for (int jj=0; jj<NGS; jj++) {
    face vector DmixGf[];
    scalar DmixG = DmixGList_G[jj];
    foreach_face() {
      double rhoGvh_G;
#ifdef VARPROP
      rhoGvh_G = face_value(rhoGv_G, 0);
#else
      rhoGvh_G = rhoG;
#endif
      DmixGf.x[] = face_value(DmixG, 0)*rhoGvh_G*fsG.x[]*fm.x[];
    }
    foreach() {
#ifdef VARPROP
      theta2[] = cm[]*max(fG[]*rhoGv_G[], F_ERR);
#else
      theta2[] = cm[]*max(fG[]*rhoG, F_ERR);
#endif
    }

    scalar YG = YGList_G[jj];
    scalar sGexp = sGexpList[jj];
#ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (YG, dt, D=DmixGf, theta=theta2);
#else
    diffusion (YG, dt, D=DmixGf, r=sGexp, theta=theta2);
#endif
  }

#ifdef SOLVE_TEMPERATURE
  foreach_face() {
    lambda1f.x[] = face_value(lambda1v.x, 0)*fsS.x[]*fm.x[];
    lambda2f.x[] = face_value(lambda2v.x, 0)*fsG.x[]*fm.x[];
  }

  foreach() {
    double theta1vh, theta2vh;
# ifdef VARPROP
    theta1vh = porosity[]*rhoGv_S[]*cpGv_S[] + (1. - porosity[])*rhoSv[]*cpSv[];
    theta2vh = rhoGv_G[]*cpGv_G[];
# else
    theta1vh = porosity[]*rhoG*cpG + (1. - porosity[])*rhoS*cpS;
    theta2vh = rhoG*cpG;
# endif

    theta1[] = cm[]*max(fS[]*theta1vh, F_ERR);
    theta2[] = cm[]*max(fG[]*theta2vh, F_ERR);
  }
# ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (TS, dt, D=lambda1f, r=sST, theta=theta1);
    diffusion_explicit (TG, dt, D=lambda2f, r=sGT, theta=theta2);
# else
    diffusion (TS, dt, D=lambda1f, r=sST, theta=theta1);
    diffusion (TG, dt, D=lambda2f, r=sGT, theta=theta2);
# endif
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
