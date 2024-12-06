#include "intgrad.h"
#include "fracface.h"

#ifdef EXPLICIT_DIFFUSION
  #include "diffusion-explicit.h"
#else
  #include "diffusion.h"
#endif

#include "common-evaporation.h"
#include "memoryallocation-t.h"
//#include "reactions.h"
#include "int-temperature-v.h"
#include "int-condition.h"

event reset_sources (i++) {
  foreach() {
    sST[] = 0.;
    sGT[] = 0.;
  }
}

// extern face vector ufsave;
// event tracer_advection (i++) { //performed in temperature-vt.h
//   // reconstruct darcy velocity
//   face vector darcyv[];
//   foreach_face() {
//     darcyv.x[] = f[] > F_ERR ? ufsave.x[]*porosity[]/f[] : ufsave.x[];
//   }
//
//   //advection of YGList
//   advection_div (YGList, ufsave, dt);
//
//   // ensure that sum(YG) = 1
//   foreach() {
//     double totmassgas = 0.;
//     for (scalar YG in YGList)
//       totmassgas += YG[];
//
//     for (scalar YG in YGList) {
//       YG[] = (totmassgas > 0.) ? YG[]/totmassgas : 0.;
//       YG[] = clamp(YG[], 0., 1.);
//       YG[] = (YG[] > 1e-10) ? YG[] : 0.;
//     }
//   }
// }

// event tracer_diffusion (i++) {
//
//   foreach() {
//     f[] = clamp (f[], 0., 1.);
//     f[] = (f[] > F_ERR) ? f[] : 0.;
//     fs[] = f[]; fG[] = 1. - f[];
//     if (fu[] > F_ERR)
//       for (int jj = 0; jj < NGS; jj++) {
//         scalar YG = YGList_S[jj];
//         YG[] /= fu[]
//       }
//
//     if (fu[] < 1-F_ERR)
//       for (int jj = 0; jj < NGS; jj++) {
//         scalar YG = YGList_G[jj];
//         YG[] /= (1. - fu[])
//       }
//   }
//
//   face_fraction (fS, fsS);
//   face_fraction (fS, fsS);
//
//   foreach() {
//     if (f[] > F_ERR && f[] < 1.-F_ERR)
//       for (int jj = 0; jj<NGS; jj++) {
//         scalar YG = YGList_S[jj];
//         scalar YGInt = YGList_Int[jj];
//         YGInt[] = avg_neighbor (point, YG, f);
//       }
//   }
//
//   //calculate interface concentration
//   int_Concentration();
//
//   foreach() {
//     if (f[]>F_ERR && f[]<1.-F_ERR){
//       coord n = facet_normal (point, fS, fsS), p;
//       double alpha = plane_alpha (fS[], n);
//       double area = plane_area_center (n, alpha, &p);
//       normalize (&n);
//
//       double
//     }
//   }
//
//
//   scalar theta1[], theta2[];
// #if TREE
//   theta1.refine = theta1.prolongation = fraction_refine;
//   theta2.refine = theta2.prolongation = fraction_refine;
//   theta1.dirty = true;
//   theta2.dirty = true;
// #endif
//
//   forech() {
//     YGInt[] = 0.;
//     for (int jj = 0; jj < NGS; jj++) {
//       scalar YG = YGList_S[jj];
//       scalar YGInt = YGList_Int[jj];
//     }
//   }
//   // calculate diff coeffincet in each cell for each species
//   foreach() {
//     //set T and P
//     OpenSMOKE_GasProp_SetTemperature (T[]);
//     OpenSMOKE_GasProp_SetPressure (p[]+Pref);
//
//     // reconstruct the mass fracs
//     double gasmassfracs [NGS];
//     for (int jj=0; jj<NGS; jj++) {
//       scalar YG = YGList[jj];
//       gasmassfracs[jj] = YG[];
//     }
//
//     //convert to mole fracs
//     double gasmolefracs [NGS];
//     OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, gas_MWs, gasmassfracs);
//
//     //calculate diff coeff
//     for (int jj=0; jj<NGS; jj++) {
//       scalar Dmix = Dmix2List[jj];
//       Dmix[] = OpenSMOKE_GasProp_Dmix(gasmolefracs, jj);
//       Dmix[] *= f[] > F_ERR ? pow(porosity[], 3./2.) : 1.; //effect of solid, to be revised
//     }
//   }
//
//   //diffuse gas species
//   for (int jj=0; jj<NGS; jj++) {
//
//     face vector Dmix2f[];
//     scalar Dmix = Dmix2List[jj];
//     foreach_face()
//       Dmix2f.x[] = Dmix[]*fm.x[];
//
//     foreach()
//       thetaG[] = cm[];
//
//     scalar YG = YGList[jj];
// #ifdef EXPLICIT_DIFFUSION
//     diffusion_explicit (YG, dt, D=Dmix2f, theta=thetaG);
// #else
//     diffusion (YG, dt, D=Dmix2f, theta=thetaG);
// #endif
//   }
//
//   foreach() {
//     double totmassgas = 0.;
//     for (scalar YG in YGList)
//       totmassgas += YG[];
//
//     for (scalar YG in YGList) {
//       YG[] = (totmassgas > 0.) ? YG[]/totmassgas : 0.;
//       YG[] = clamp(YG[], 0., 1.);
//       YG[] = (YG[] > 1e-10) ? YG[] : 0.;
//     }
//   }
// }

event tracer_diffusion (i++) {

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    fS[] = f[]; fG[] = 1. - f[];
#ifdef SOLVE_TEMPERATURE
    TS[] = fu[] > F_ERR ? TS[]/fu[] : 0.;
    TG[] = ((1. - fu[]) > F_ERR) ? TG[]/(1. - fu[]) : 0.;
#endif
    if (fu[] > F_ERR)
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_S[jj];
        YG[] /= fu[];
      }

    if (fu[] < 1-F_ERR)
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        YG[] /= (1. - fu[]);
      }
  }

  //Compute face gradients
  face_fraction (fS, fsS);
  face_fraction (fG, fsG);

#ifdef SOLVE_TEMPERATURE
  //Assign interface temperature first guess
  foreach() {
    TInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      TInt[] = avg_neighbor (point, TS, f);
  }

  #ifdef FIXED_INT_TEMP //Force interface temperature = TG0
  foreach()
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      TInt[] = TG0;
  #else //default: solve interface balance
  //REVIEW THIS PART
  ijc_CoupledTemperature();
  #endif
#endif

  //find the interface concentration foreach species
  intConcentration();





  //Compute source terms for the energy balance 
  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord n = facet_normal (point, fS, fsS), p;
      double alpha = plane_alpha (fS[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      double bc = TInt[];
      double Gtrgrad = ebmgrad (point, TG, fS, fG, fsS, fsG, true, bc, &success);
      double Strgrad = ebmgrad (point, TS, fS, fG, fsS, fsG, false, bc, &success);

      double Sheatflux = lambda1v[]*Strgrad;
      double Gheatflux = lambda2v[]*Gtrgrad;

#ifdef AXI
      sST[] += Sheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sGT[] += Gheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      sST[] += Sheatflux*area/Delta*cm[];
      sGT[] += Gheatflux*area/Delta*cm[];
#endif
    }
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
    scalar Dmix2 = Dmix2List_S[jj];
    foreach_face()
      Dmixf.x[] = face_value(Dmix2, 0)*rhoG*fsS.x[]*fm.x[];

    foreach()
      theta1[] = cm[]*max(fS[]*rhoG, F_ERR);
    
    scalar YG = YGList_S[jj];
    scalar sSexp = zeroc; //REVISE
    scalar sSimp = zeroc;

    diffusion (YG, dt, D=Dmixf, r=sSexp, beta=sSimp, theta=theta1);
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

  diffusion (TS, dt, D=lambda1f, r=sST, theta=theta1);
  diffusion (TG, dt, D=lambda2f, r=sGT, theta=theta2);
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

// calculate diff coefficents, to be moved to properties
event properties (i++) {
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
      OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, gas_MWs, gasmassfracs);

      //calculate diff coeff
      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix = Dmix2List_S[jj];
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
      OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, gas_MWs, gasmassfracs);

      //calculate diff coeff
      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix = Dmix2List_G[jj];
        Dmix[] = OpenSMOKE_GasProp_Dmix(gasmolefracs, jj);
      }
    }
  }
}
