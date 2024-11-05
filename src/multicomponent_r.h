#include "intgrad.h"
#include "fracface.h"
#include "diffusion.h"
#include "memoryallocation.h"
#include "reactions.h"

extern scalar T;

extern face vector ufsave;
event tracer_advection (i++) {
  // reconstruct darcy velocity
  face vector darcyv[];
  foreach_face() {
    double ff = face_value (f,0);
    darcyv.x[] = ff > F_ERR ? ufsave.x[]*face_value(porosity, 0)/ff : ufsave.x[];
  }

  //advection of YGList
  advection (YGList, darcyv, dt);

  // ensure that sum(YG) = 1
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

event tracer_diffusion (i++) {

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
  }

  scalar thetaG[];
#if TREE
  thetaG.refine = thetaG.prolongation = fraction_refine;
  thetaG.dirty = true;
#endif

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
      Dmix[] = OpenSMOKE_GasProp_Dmix(gasmolefracs, jj); // *0 TEMP 
      Dmix[] *= f[] > F_ERR ? pow(porosity[], 3./2.) : 1.; //effect of solid, to be revised
    }
  }

  //diffuse each species
  for (int jj=0; jj<NGS; jj++) {

    face vector Dmix2f[];
    scalar Dmix = Dmix2List[jj];
    foreach_face()
      Dmix2f.x[] = Dmix[]*fm.x[];

    foreach()
      thetaG[] = cm[];

    scalar YG = YGList[jj];

    diffusion (YG, dt, D=Dmix2f, theta=thetaG);
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
