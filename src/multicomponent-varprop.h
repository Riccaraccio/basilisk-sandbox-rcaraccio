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

// function to check the variables and dump the simulation if they are out of range
void check_variables() {
  bool status = true;

  //check f
  foreach()
    if (f[] < 0. || f[] > 1.) {
      fprintf(stderr, "f out of range: %g\n", f[]);
      status = false;
    }

  //check solid temperature
  if (status) {
    foreach() {
      if (f[] > F_ERR) {
        if (TS[] < 100. || TS[] > 1000.) {
          fprintf(stderr, "Solid temperature out of range: %g\n", TS[]);
          status = false;
        }
      }
    }
  }
    
  //check gas temperature
  if (status) {
    foreach() {
      if (f[] < 1.- F_ERR) {
        if (TG[] < 100. || TG[] > 1000.) {
          fprintf(stderr, "Gas temperature out of range: %g\n", TG[]);
          status = false;
        }
      }
    }
  }

  //check solid mass fractions
  if (status) {
    foreach() {
      if (f[] > F_ERR) {
        for (int jj=0; jj<NGS; jj++) {
          scalar YG_S = YGList_S[jj];
          if (YG_S[] < 0. || YG_S[] > 1.) {
            fprintf(stderr, "Solid mass fraction out of range: %g\n", YG_S[]);
            status = false;
          }
        }
      }
    }
  }

  //check external gas mass fractions
  if (status) {
    foreach() {
      if (f[] < 1.- F_ERR) {
        for (int jj=0; jj<NGS; jj++) {
          scalar YG_G = YGList_G[jj];
          if (YG_G[] < 0. || YG_G[] > 1.) {
            fprintf(stderr, "Gas mass fraction out of range: %g\n", YG_G[]);
            status = false;
          }
        }
      }
    }
  }

  //check internal gas mass fractions
  if (status) {
    foreach() {
      if (f[] > F_ERR) {
        for (int jj=0; jj<NGS; jj++) {
          scalar YG= YGList_S[jj];
          if (YG[] < 0. || YG[] > 1.) {
            fprintf(stderr, "Gas mass fraction out of range: %g\n", YG[]);
            status = false;
          }
        }
      }
    }
  }
 
  //check porosity
  if (status) {
    foreach() {
      if (f[] > F_ERR) {
        if (porosity[] < 0. || porosity[] > 1.) {
          fprintf(stderr, "Porosity out of range: %g\n", porosity[]);
          status = false;
        }
      }
    }
  }
  
  if (!status) {
    dump();
    exit(1);
  }
}

event reset_sources (i++) {
#ifdef SOLVE_TEMPERATURE
  foreach() {
    sST[] = 0.;
    sGT[] = 0.;
  }
#endif

  reset (sGexpList, 0.);
  reset (sSexpList, 0.);
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
    // u_prime.x[] = (ff > F_ERR) ? ufsave.x[]*(rhoGvh_S*cpGvh_S)/(rhoGvh_S*cpGvh_S*ef + rhoSvh*cpSvh*(1.-ef))*ff + (1-ff)*ufsave.x[] : ufsave.x[];
  }

  advection_div({TS}, u_prime, dt);
# ifndef TEMPERATURE_PROFILE
  advection_div({TG}, ufsave, dt);
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

void check_and_correct_fractions (scalar* YList, int n, bool inverse, char* name) { //YList in tracer form
  bool warning = false;
  foreach() {
    double sum = 0.;

    for (int jj = 0; jj < n; jj++) {
      scalar Y = YList[jj];
      if (!inverse) {
        Y[] = (f[] > F_ERR) ? Y[]/f[] : 0.;
      } else {
        Y[] = ((1.-f[]) > F_ERR) ? Y[]/(1. - f[]) : 0.;
      }
      if (Y[] < 1e-10) Y[] = 0.;
      sum += Y[];
    }

    if ((fabs(sum - 1.) > 1.e-1) && !warning && (sum > 1e-10)) {
      fprintf(stderr, "Warning: sum of mass fractions is not equal to 1 for %s: %g\n", name, sum);

      for (int jj = 0; jj < n; jj++) {
        scalar Y = YList[jj];
        fprintf(stderr, "Y[%d] = %g\n", jj, Y[]);
      }
      fprintf(stderr, "f = %g\n", f[]);
      fprintf(stderr, "Inverse = %d\n", inverse);
      fprintf(stderr, "Point: (%g, %g, %g)\n", x, y, z);

      // warning = true; // set a flag to avoid multiple warnings
    }

    for (int jj = 0; jj < n; jj++) {
      scalar Y = YList[jj];
      Y[] = (sum > 1e-10) ? Y[]/sum : 0.;
      if (!inverse) {
        Y[] = (f[] > F_ERR) ? Y[]*f[] : 0.;
      } else {
        Y[] = ((1. - f[]) > F_ERR) ? Y[]*(1. - f[]) : 0.;
      }
    }
  }
}

event tracer_diffusion (i++) {
  
  //Check the mass fractions Can be removed for performance
  check_and_correct_fractions (YGList_S, NGS, false, "YGList_S-1");
  check_and_correct_fractions (YGList_G, NGS, true,  "YGList_G-1");
  check_and_correct_fractions (YSList,   NSS, false, "YSList-1");

#ifdef VARPROP
  update_properties();
#endif

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
      TInt[] = (TS[] + TG[])/2;
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
      scalar YGInt = YGList_Int[jj];
      YGInt[] = 0.;
      if (f[] > F_ERR && f[] < 1.-F_ERR) {
        scalar YG_S = YGList_S[jj];
        scalar YG_G = YGList_G[jj];
        YGInt[] = (YG_G[] + YG_S[])/2;
        YGInt[] = clamp (YGInt[], 0., 1.);
      }
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

      n.x = fabs(n.x); n.y = fabs(n.y);

      double lambda1vh = n.x/(n.x+n.y)*lambda1v.x[] + n.y/(n.x+n.y)*lambda1v.y[];
      double lambda2vh = n.x/(n.x+n.y)*lambda2v.x[] + n.y/(n.x+n.y)*lambda2v.y[];

      double Sheatflux = lambda1vh*Strgrad;
      double Gheatflux = lambda2vh*Gtrgrad;

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
  update_divergence_density();
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
      theta1[] = cm[]*max(fS[]*rhoGv_S[]*porosity[], F_ERR);
#else
      theta1[] = cm[]*max(fS[]*rhoG*porosity[], F_ERR);
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
#  ifndef TEMPERATURE_PROFILE
    diffusion (TG, dt, D=lambda2f, r=sGT, theta=theta2);
#  endif
# endif
#endif

  //recover tracer form
  foreach() {
    for (scalar YG in YGList_S)
      YG[] = (f[] > F_ERR) ? YG[]*f[] : 0.;

    for (scalar YG in YGList_G)
      YG[] = ((1. - f[]) > F_ERR) ? YG[]*(1. - f[]) : 0.;
    
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]*f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]*(1. - f[]) : 0.;
    T[] = TS[] + TG[];
#endif
  }

  check_and_correct_fractions (YGList_S, NGS, false, "YGList_S-2");
  check_and_correct_fractions (YGList_G, NGS, true,  "YGList_G-2");
  check_and_correct_fractions (YSList,   NSS, false, "YSList-2");
}
