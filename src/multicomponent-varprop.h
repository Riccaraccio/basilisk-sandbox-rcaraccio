#ifndef MULTICOMPONENT
  #define MULTICOMPONENT 1
#endif

#include "intgrad.h"

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

void check_and_correct_fractions(scalar *YList, int n, bool inverse) {
  foreach () {
    double sum = 0.;
    double temp[n];
    for (int jj = 0; jj < n; jj++) {
      scalar Y = YList[jj];
      double val = Y[];
      temp[jj] = (val < F_ERR) ? 0. : val;
      sum += temp[jj];
    }

    if (sum > F_ERR) {
      double scale = 1. / sum;
      for (int jj = 0; jj < n; jj++) {
        scalar Y = YList[jj];
        Y[] = temp[jj] * scale;
      }
    }
    else {
      for (int jj = 0; jj < n; jj++) {
        scalar Y = YList[jj];
        Y[] = 0.;
      }
    }
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

foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    f[] = (f[] < 1.-F_ERR) ? f[] : 1.;
    fS[] = f[]; fG[] = 1. - f[];
  }

  check_and_correct_fractions (YGList, NGS, true);
  check_and_correct_fractions (YSList,   NSS, false);

#ifdef VARPROP
  update_properties();
#else
  update_properties_constant();
#endif

  // lose tracer form and extrapolate fields
  foreach() {
    porosity[] = (f[] > F_ERR) ? porosity[]/f[] : 0.;
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]/f[] : 0.;
    TG[] = ((1.-f[]) > F_ERR) ? TG[]/(1.-f[]) : 0.;

    TS[] = (f[] > F_ERR) ? TS[] : TG[];
    TG[] = (f[] < 1.-F_ERR) ? TG[] : TS[];
#endif
  }

  advection_div(YGList, ufsave, dt);

#ifdef SOLVE_TEMPERATURE
  foreach_face() {
    double ef = face_value(porosity, 0);
    double ff = face_value(f, 0);

    double rhoGvh, rhoSvh;
    double cpGvh, cpSvh;
    
    #ifdef VARPROP
    rhoGvh = face_value(rhoGv, 0); rhoSvh = face_value(rhoSv, 0);
    cpGvh = face_value(cpGv, 0); cpSvh = face_value(cpSv, 0);
    #else
    rhoGvh = rhoG; rhoSvh = rhoS;
    cpGvh = cpG; cpSvh = cpS;
    #endif

    u_prime.x[] = (ff > F_ERR) ? ufsave.x[]*(rhoGvh*cpGvh)/(rhoGvh*cpGvh*ef + rhoSvh*cpSvh*(1.-ef)) : ufsave.x[];
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
  }
}
#endif

void update_mole_fields() {
  #ifdef MOLAR_DIFFUSION
  foreach() {
    double xG[NGS], yG[NGS];
    double MWmix;
    for (int jj = 0; jj < NGS; jj++) {
      scalar YG = YGList[jj];
      yG[jj] = YG[];
    }
    OpenSMOKE_MoleFractions_From_MassFractions (xG, &MWmix, yG);
    for (int jj = 0; jj < NGS; jj++) {
      scalar XG = XGList[jj];
      XG[] = xG[jj];
    }
  }
  #endif
}

event tracer_diffusion (i++) {

  check_and_correct_fractions (YGList, NGS, true);
  check_and_correct_fractions (YSList, NSS, false);

#ifdef SOLVE_TEMPERATURE
  foreach() {
    TS[] = (f[] > F_ERR) ? TS[]/f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]/(1. - f[]) : 0.;
  }
#endif

#ifdef MOLAR_DIFFUSION
  update_mole_fields();
#endif

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

#ifdef SOLVE_TEMPERATURE
  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord n = facet_normal (point, fS, fsS), p;
      double alpha = plane_alpha (fS[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);
      
      double bc = TInt[];
      double Strgrad = ebmgrad (point, TS, fS, fG, fsS, fsG, false, bc, &success);
      double Gtrgrad = ebmgrad (point, TG, fS, fG, fsS, fsG, true , bc, &success);

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

#ifdef MASS_DIFFUSION_ENTHALPY
  foreach() {
    if (f[] > 1. - F_ERR) { //Internal gas phase
      double mdeGS = 0.;
      coord gTS = {0., 0., 0.};
      coord gYGj_S = {0., 0., 0.};
      coord gYGsum_S = {0., 0., 0.};

      foreach_dimension()
        gTS.x = (TS[1] - TS[-1])/(2.*Delta);
      
      foreach_dimension() {
        for (int jj=0; jj<NGS; jj++) {
          scalar Dmixv = DmixGList[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList[jj];
          gYGsum_S.x -= (MWmixG[] > 0.) ?
            rhoGv[]*Dmixv[]*gas_MWs[jj]/MWmixG[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
  # else
          scalar YG = YGList[jj];
          gYGsum_S.x -= rhoGv[]*Dmixv[]*(YG[1] - YG[-1])/(2.*Delta);
  # endif
        }

        for (int jj=0; jj<NGS; jj++) {
          scalar YG = YGList[jj];
          scalar cpGv = cpGList[jj];
          scalar Dmixv = DmixGList[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList[jj];
          gYGj_S.x = (MWmixG[] > 0.) ?
            -rhoGv[]*Dmixv[]*gas_MWs[jj]/MWmixG[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
  # else
          gYGj_S.x = -rhoGv[]*Dmixv[]*(YG[1] - YG[-1])/(2.*Delta);
  # endif
          mdeGS += cpGv[]*(gYGj_S.x - YG[]*gYGsum_S.x)*gTS.x;
        }
      }
    sST[] -= mdeGS*cm[];
    }
   
    if (f[] < F_ERR) { //Internal gas phase
      double mdeGG = 0.;
      coord gTG = {0., 0., 0.};
      coord gYGj_G = {0., 0., 0.};
      coord gYGsum_G = {0., 0., 0.};

      foreach_dimension()
        gTG.x = (TG[1] - TG[-1])/(2.*Delta);

      foreach_dimension() {
        for (int jj=0; jj<NGS; jj++) {
          scalar Dmixv = DmixGList[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList[jj];
          gYGsum_G.x -= (MWmixG[] > 0.) ?
            rhoGv[]*Dmixv[]*gas_MWs[jj]/MWmixG[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
  # else
          scalar YG = YGList[jj];
          gYGsum_G.x -= rhoGv[]*Dmixv[]*(YG[1] - YG[-1])/(2.*Delta);
  # endif
        }

        for (int jj=0; jj<NGS; jj++) {
          scalar YG = YGList[jj];
          scalar cpGv = cpGList[jj];
          scalar Dmixv = DmixGList[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList[jj];
          gYGj_G.x = (MWmixG[] > 0.) ?
            -rhoGv[]*Dmixv[]*gas_MWs[jj]/MWmixG[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
  # else
          gYGj_G.x = -rhoGv[]*Dmixv[]*(YG[1] - YG[-1])/(2.*Delta);
  # endif
          mdeGG += cpGv[]*(gYGj_G.x - YG[]*gYGsum_G.x)*gTG.x;
        }
      }
    sGT[] -= mdeGG*cm[];
    }
  }
#endif //MASS_DIFFUSION_ENTHALPY

#if defined VARPROP && !defined NO_EXPANSION
  update_divergence();
  // update_divergence_density();
#endif

#ifdef FICK_CORRECTED
  face vector phicGtot[];
  foreach_face() {
    phicGtot.x[] = 0.;
    for (int jj=0; jj<NGS; jj++) {
      scalar DmixG = DmixGList[jj];
      double DmixGf = face_value(DmixG, 0);
      double rhoGf;
# ifdef VARPROP
      rhoGf = face_value(rhoGv, 0);
# else
      rhoGf = rhoG;
# endif

# ifdef MOLAR_DIFFUSION
      scalar XG = XGList[jj];
      double MWmixf = face_value(MWmixG, 0);
      phicGtot.x[] += (MWmixf > 0.) ? rhoGf*DmixGf*face_gradient_x (XG, 0)*gas_MWs[jj]/MWmixf : 0.;
# else
      scalar YG = YGList[jj];
      phicGtot.x[] += rhoGf*DmixGf*face_gradient_x (YG, 0);
# endif
    }
    phicGtot.x[] *= fsG.x[]*fm.x[];
  }

  //Apply the Fick's law correction
  for (int jj=0; jj<NGS; jj++) {
    face vector phicjj[];
    foreach_face() {
      phicjj.x[] = phicGtot.x[];
#ifdef MOLAR_DIFFUSION
      scalar DmixG = DmixGList[jj];
      double DmixGf = face_value(DmixG, 0);
      double MWmixf = face_value(MWmixG, 0);

      double rhoGf;
# ifdef VARPROP
      rhoGf = face_value(rhoGv, 0);
# else
      rhoGf = rhoG;
# endif
      phicjj.x[] -= (MWmixf > 0.) ? rhoGf*DmixGf/MWmixf*face_gradient_x (MWmixG, 0)*fsG.x[]*fm.x[] : 0.;
#endif
    }

    scalar YG = YGList[jj];
    double (* gradient_backup)(double, double, double) = YG.gradient; // we need to backup the gradient function
    YG.gradient = NULL; //reset the gradient
    face vector flux[];
    tracer_fluxes (YG, phicjj, flux, dt, zeroc); //calculate the fluxes using the corrective velocity
    YG.gradient = gradient_backup; // restore the gradient function
    
    // apply the corrective fluxes
    foreach()
      foreach_dimension()
        YG[] += (rhoGv[] > 0.) ? dt/(rhoGv[])*(flux.x[] - flux.x[1])/(Delta*cm[]) : 0.; 
  }
  #endif //FICK_CORRECTED

  scalar theta1[], theta2[];

#if TREE
  theta1.refine = theta1.prolongation = fraction_refine;
  theta2.refine = theta2.prolongation = fraction_refine;
  theta1.dirty = true;
  theta2.dirty = true;
#endif

  // gas diffusion
  for (int jj=0; jj<NGS; jj++) {
    face vector DmixGf[];
    scalar DmixG = DmixGList[jj];
    foreach_face() {
      double rhoGvh;
#ifdef VARPROP
      rhoGvh = face_value(rhoGv, 0);
#else
      rhoGvh = rhoG;
#endif
      DmixGf.x[] = face_value(DmixG, 0)*rhoGvh*fsS.x[]*fm.x[];
    }

    scalar eps[];
    foreach()
      eps[] = porosity[] + (1. - f[]);

    foreach() {
#ifdef VARPROP
      theta1[] = cm[]*max(rhoGv[]*eps[], F_ERR); // porosity is already multiplied by fS
#else
      theta1[] = cm[]*max(rhoG*eps[], F_ERR); // porosity is already multiplied by fS
#endif
    }

    scalar YG = YGList[jj];
    scalar sSexp = sSexpList[jj];

#ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (YG, dt, D=DmixGf, theta=theta1);
#else
    diffusion (YG, dt, D=DmixGf, r=sSexp, theta=theta1);
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
    theta1vh = fS[] > F_ERR ? porosity[]/fS[]*rhoGv[]*cpGv[] + (1. - porosity[]/fS[])*rhoSv[]*cpSv[] : 0.;
    theta2vh = rhoGv[]*cpGv[];
# else
    theta1vh = fS[] > F_ERR ? porosity[]/fS[]*rhoG*cpG + (1. - porosity[]/fS[])*rhoS*cpS : 0.;
    theta2vh = rhoG*cpG;
# endif

    theta1[] = cm[]*max(fS[]*theta1vh, F_ERR);
    theta2[] = cm[]*max(fG[]*theta2vh, F_ERR);
  }

#ifdef VARCOEFF
  foreach()
    porosity[] = (f[] > F_ERR) ? porosity[]/f[] : 0;

  foreach_face() {
    double ef = face_value(porosity, 0);
    lambda1f.x[] = (ef > F_ERR) ? lambda1f.x[] / (rhoG*cpG*ef + rhoS*cpS*(1. - ef)) : 0.;
    lambda2f.x[] = lambda2f.x[] / (rhoG*cpG);
  }

  foreach() {
    theta1[] = cm[] * max(fS[], F_ERR);
    theta2[] = cm[] * max(fG[], F_ERR);
  }

  foreach() {
    sST[] = (f[] > F_ERR) ? sST[] / (rhoG*cpG*porosity[] + rhoS*cpS*(1. - porosity[])) : 0.;
    sGT[] = sGT[] / (rhoG*cpG);
  }

  foreach()
    porosity[] *= f[];
#endif

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
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]*f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]*(1. - f[]) : 0.;
    T[] = TS[] + TG[];
#endif
  }

  check_and_correct_fractions (YGList, NGS, true);
}
