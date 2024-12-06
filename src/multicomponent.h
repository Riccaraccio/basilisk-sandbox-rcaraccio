#include "intgrad.h"
#include "fracface.h"

#ifdef EXPLICIT_DIFFUSION
  #include "diffusion-explicit.h"
#else
  #include "diffusion.h"
#endif

#include "memoryallocation.h"
#include "reactions.h"

#ifndef SOLVE_TEMPERATURE
scalar fsS[], fsG[];
#endif

extern scalar T;

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

event tracer_diffusion (i++) {

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    fs[] = f[]; fG[] = 1. - f[];
    if (fu[] > F_ERR)
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_S[jj];
        YG[] /= fu[]
      }

    if (fu[] < 1-F_ERR)
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        YG[] /= (1. - fu[])
      }
  }

  face_fraction (fS, fsS);
  face_fraction (fS, fsS);

  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      for (int jj = 0; jj<NGS; jj++) {
        scalar YG = YGList_S[jj];
        scalar YGInt = YGList_Int[jj];
        YGInt[] = avg_neighbor (point, YG, f);
      }
  }

  //calculate interface concentration
  int_Concentration();

  foreach() {
    if (f[]>F_ERR && f[]<1.-F_ERR){
      coord n = facet_normal (point, fS, fsS), p;
      double alpha = plane_alpha (fS[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      double
    }
  }


  scalar theta1[], theta2[];
#if TREE
  theta1.refine = theta1.prolongation = fraction_refine;
  theta2.refine = theta2.prolongation = fraction_refine;
  theta1.dirty = true;
  theta2.dirty = true;
#endif

  forech() {
    YGInt[] = 0.;
    for (int jj = 0; jj < NGS; jj++) {
      scalar YG = YGList_S[jj];
      scalar YGInt = YGList_Int[jj];
    }
  }
  // calculate diff coeffincet in each cell for each species
  foreach() {
    //set T and P
    OpenSMOKE_GasProp_SetTemperature (T[]);
    OpenSMOKE_GasProp_SetPressure (p[]+Pref);

    // reconstruct the mass fracs
    double gasmassfracs [NGS];
    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList[jj];
      gasmassfracs[jj] = YG[];
    }

    //convert to mole fracs
    double gasmolefracs [NGS];
    OpenSMOKE_MoleFractions_From_MassFractions(gasmolefracs, gas_MWs, gasmassfracs);

    //calculate diff coeff
    for (int jj=0; jj<NGS; jj++) {
      scalar Dmix = Dmix2List[jj];
      Dmix[] = OpenSMOKE_GasProp_Dmix(gasmolefracs, jj);
      Dmix[] *= f[] > F_ERR ? pow(porosity[], 3./2.) : 1.; //effect of solid, to be revised
    }
  }

  //diffuse gas species
  for (int jj=0; jj<NGS; jj++) {

    face vector Dmix2f[];
    scalar Dmix = Dmix2List[jj];
    foreach_face()
      Dmix2f.x[] = Dmix[]*fm.x[];

    foreach()
      thetaG[] = cm[];

    scalar YG = YGList[jj];
#ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (YG, dt, D=Dmix2f, theta=thetaG);
#else
    diffusion (YG, dt, D=Dmix2f, theta=thetaG);
#endif
  }

  foreach() {
    double totmassgas = 0.;
    for (scalar YG in YGList)
      totmassgas += YG[];

    for (scalar YG in YGList) {
      YG[] = (totmassgas > 0.) ? YG[]/totmassgas : 0.;
      YG[] = clamp(YG[], 0., 1.);
      YG[] = (YG[] > 1e-10) ? YG[] : 0.;
    }
  }
}
