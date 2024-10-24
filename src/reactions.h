#include "osreactors.h"

scalar omega[];
extern scalar zeta;
extern scalar porosity;
extern scalar rho1v, cp1v;

//TOBEDONE ADD cps, cpg, rhos, rhog

event init (i = 0) {
  OpenSMOKE_InitODESolver ();
}

event cleanup (t = end) {
  OpenSMOKE_CleanODESolver ();
}

event chemisrty (i++) {

#ifdef SOVE_TEMPERATURE
  odefunction batch = &batch_nonisothermal_constantpressure;
  unsigned int NEQ = NGS + NSS + 1 + 1; //NGS + NSS + porosity + T
#else
  odefunction batch = &batch_isothermal_constantpressure;
  unsigned int NEQ = NGS + NSS + 1;
#endif

  //Solve solid phase chemistry
  foreach ()
    if (f[] > F_ERR) {

      porosity[] /= f[];

      double y0ode[NEQ];
      UserDataODE data;
      data.P = Pref + p[];
      data.T = T[];
      data.rhos = rhoS;
      data.rhog = rhoG;
      data.zeta = zeta[];
      data.f = f[];
#ifdef SOLVE_TEMPERATURE
      data.cps = cpS;
      data.cpg = cpG;
#endif

      for(int jj=0; jj<NSS; jj++){
        scalar YS = YSList[jj];
        YS[] /= f[];
      }

      for(int jj=0; jj<NGS; jj++){
        scalar YG = YGList[jj];
        YG[] /= f[]; //Div by f[]??
      }
      double gasmass[NGS];
      for (int jj=0; jj<NGS; jj++){
        gasmass[jj] = YG[]*rhoG*(1-f[]+porosity[]*f[]);
        y0ode[jj] = gasmass[jj];
      }

      double solidmass[NSS];
      for (int jj=0; jj<NSS; jj++){
        scalar YS = YSList[jj];
        solidmass[jj] = YS[]*rhoS*(1-porosiy[])*f[];
        y0ode[jj+NGS] = solidmass[jj];
      }

      y0ode[NGS+NSS] = porosity[];
#ifdef SOLVE_TEMPERATURE
      y0ode[NGS+NSS+1] = T[];
#endif

      OpenSMOKE_ODESolver (batch, NEQ, dt, y0ode, &data);

      double totgasmass = 0;
      for (int jj=0; jj<NGS; jj++){
        totgasmass += y0ode[jj];
      }
      for (int jj=0; jj<NGS; jj++){
        scalar YG = YGList[jj];
        YG[] = y0ode[jj]/totgasmass;
      }

      double totsolidmass = 0;
      for (int jj=0; jj<NSS; jj++){
        totsolidmass += y0ode[jj+NGS];
      }

      for (int jj=0; jj<NSS; jj++){
        scalar YS = YSList[jj];
        YS[] = (totsolidmass<1e-5) ? 0. : y0ode[jj+NGS]/totsolidmass;
      }

      porosity[] = y0ode[NGS+NSS];

#ifdef SOLVE_TEMPERATURE
      T[] = y0ode[NGS+NSS+1];
#endif
      //recover tracer form
      porosity[] *= f[];
    }
}
