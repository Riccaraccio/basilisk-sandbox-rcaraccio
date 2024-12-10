#include "osreactors.h"

extern scalar zeta;
extern scalar T;
extern scalar porosity;
extern scalar rho1v, cp1v;

//TOBEDONE ADD cpsv, cpgv, rhosv, rhogv

void OpenSMOKE_ODESolverEXP (odefunction ode, unsigned int neq, double dt, double* y, void* args) {

  double dy[neq];
  ode(y, dt, dy, args);

  for (int jj=0; jj<neq; jj++)
    y[jj] += dt*dy[jj];
}

event init (i = 0) {
  OpenSMOKE_InitODESolver ();
}

event cleanup (t = end) {
  OpenSMOKE_CleanODESolver ();
}

event chemistry (i++) {

  //reset omega
  foreach()
    omega[] = 0.;

#ifdef SOLVE_TEMPERATURE
  odefunction batch = &batch_nonisothermal_constantpressure;
  unsigned int NEQ = NGS + NSS + 1 + 1; //NGS + NSS + porosity + T
#else
  odefunction batch = &batch_isothermal_constantpressure;
  unsigned int NEQ = NGS + NSS + 1;
#endif

  //Solve solid phase chemistry
  foreach ()
    if (f[] > F_ERR) {

      for(int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        YS[] /= f[];
      }

      for(int jj=0; jj<NGS; jj++) {
        scalar YG = YGList_S[jj];
        YG[] /= f[];
      }
      porosity[] /= f[];

      double y0ode[NEQ];
      UserDataODE data;
      data.P = Pref + p[];
      data.T = TS[]/f[];
      data.rhos = rhoS;
      data.zeta = zeta[];
#ifdef SOLVE_TEMPERATURE
      data.cps = cpS;
      data.cpg = cpG;
#endif
      double sources[NEQ];
      data.sources = sources;

      double gasmass[NGS];
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList_S[jj];
        gasmass[jj] = YG[]*rhoG*porosity[];
        y0ode[jj] = gasmass[jj];
      }

      double solidmass[NSS];
      for (int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        solidmass[jj] = YS[]*rhoS*(1-porosity[]);
        y0ode[jj+NGS] = solidmass[jj];
      }

      y0ode[NGS+NSS] = porosity[];

#ifdef SOLVE_TEMPERATURE
      y0ode[NGS+NSS+1] = TS[]/f[];
#endif

#ifdef EXPLICIT_REACTIONS
    OpenSMOKE_ODESolverEXP (batch, NEQ, dt, y0ode, &data);
#else //default
    OpenSMOKE_ODESolver (batch, NEQ, dt, y0ode, &data); 
#endif

      double totgasmass = 0;
      for (int jj=0; jj<NGS; jj++) {
        totgasmass += y0ode[jj];
      }

      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList_S[jj];
        YG[] = y0ode[jj]/totgasmass*f[];
      }

      double totsolidmass = 0;
      for (int jj=0; jj<NSS; jj++) {
        totsolidmass += y0ode[jj+NGS];
      }

      for (int jj=0; jj<NSS; jj++) { 
        scalar YS = YSList[jj];
        YS[] = (totsolidmass < 1e-5) ? 0. : y0ode[jj+NGS]/totsolidmass*f[];
      }

      porosity[] = y0ode[NGS+NSS]*f[];

#ifdef SOLVE_TEMPERATURE
      TS[] = y0ode[NGS+NSS+1]*f[];
#endif
      omega[] = sources[NGS+NSS];
    }
}
