#include "osreactors.h"

extern scalar zeta;
extern scalar TS;
extern scalar porosity;
extern scalar rho1v, cp1v;

//TOBEDONE ADD cpsv, cpgv, rhosv, rhogv

event init (i = 0) {
  OpenSMOKE_InitODESolver ();
}

event cleanup (t = end) {
  OpenSMOKE_CleanODESolver ();
}

event chemistry (i++) {

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

      // YS are the only one in tracer form
      for(int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        YS[] /= f[];
      }

      porosity[] /= f[];

      double y0ode[NEQ];
      UserDataODE data;
      data.P = Pref + p[];
      data.T = TS[]/f[];
      data.rhos = rhoS;
      data.zeta = zeta[];
      data.f = f[];
#ifdef SOLVE_TEMPERATURE
      data.cps = cpS;
      data.cpg = cpG;
#endif
      double sources[NEQ];
      data.sources = sources;

      double gasmass[NGS];
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList[jj];
        gasmass[jj] = YG[]*rhoG*(1-f[] +porosity[]*f[]);
        y0ode[jj] = gasmass[jj];
      }

      double solidmass[NSS];
      for (int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        solidmass[jj] = YS[]*rhoS*(1-porosity[])*f[];
        y0ode[jj+NGS] = solidmass[jj];
      }

      y0ode[NGS+NSS] = porosity[];

#ifdef SOLVE_TEMPERATURE
      y0ode[NGS+NSS+1] = TS[];
#endif

      OpenSMOKE_ODESolver (batch, NEQ, dt, y0ode, &data); 

      double totgasmass = 0;
      for (int jj=0; jj<NGS; jj++) {
        totgasmass += y0ode[jj];
      }

      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList[jj];
        YG[] = y0ode[jj]/totgasmass;
      }

      double totsolidmass = 0;
      for (int jj=0; jj<NSS; jj++) {
        totsolidmass += y0ode[jj+NGS];
      }

      for (int jj=0; jj<NSS; jj++) { 
        scalar YS = YSList[jj];
        YS[] = (totsolidmass < 1e-5) ? 0. : y0ode[jj+NGS]/totsolidmass*f[];
      }

      data.zeta = 0.;
      porosity[] = y0ode[NGS+NSS]*f[];

#ifdef SOLVE_TEMPERATURE
      TS[] = y0ode[NGS+NSS+1]*f[];
#endif
      omega[] = sources[NGS+NSS];
    }
}
