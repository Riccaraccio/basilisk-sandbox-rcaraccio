#include "intgrad.h"
#include "fracface.h"
#include "diffusion.h"
#include "memoryallocation.h"
#include "osreactors.h"

//TOBEDONE ADD cps, cpg, rhos, rhog

event init (i = 0) {
  OpenSMOKE_InitODESolver ();
}

event cleanup (t = end) {
  OpenSMOKE_CleanODESolver ();
}

event chemisrty (i++) {

  odefunction batch = &batch_nonisothermal_constantpressure;
  unsigned int NEQ = NGS + NSS + 1 + 1; //NGS + NSS + porosity + T

  //Solve solid phase chemistry
  foreach () 
    if (f[] > F_ERR) {

      for (int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        YS[] /= f[];
      }

      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList[jj];
        YG[] /= f[];
      }

      epsilon[] /= f[];

      double y0ode[NEQ];
      UserDataODE data;
      data.P = Pref + p[];
      data.T = T[];
      data.rho = rhos*(1-porosity[]) + rhog*porosity[];
      data.cp = cps*(1-porosity[]) + cpg*porosity[];

      double gasmass[NGS];
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList[jj];
        gasmass[jj] = YG[]*rhog*(1 -f[] +porosity[]*f[]);
        y0ode[jj] = gasmass[jj];
      }

      double solidmass[NSS];
      for (int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        solidmass[jj] = YS[]*rhos*(f[] - porosity[]*f[]);
        y0ode[jj+NGS] = solidmass[jj];
      }

      y0ode[NGS+NSS] = porosity[];
      y0ode[NGS+NSS+1] = T[];

      OpenSMOKE_ODESolver (batch, NEQ, dt, y0ode, &data);

      double totgasmass = 0;
      for (int jj=0; jj<NGS; jj++) {
        totgasmass += y0ode[jj];
      }

      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList[jj];
        YG[] = totgasmass < 1e-5 ? 0. : y0ode[jj]/totgasmass;
      }

      double totsolidmass = 0;
      for (int jj=0; jj<NSS; jj++) {
        totsolidmass += y0ode[jj+NGS];
      }

      for (int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        YS[] = totsolidmass < 1e-5 ? 0. : y0ode[jj+NGS]/totsolidmass;
      }

      epsilon[] = y0ode[NGS+NSS];
      T[] = y0ode[NGS+NSS+1];


      //recover tracer form

      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YSList[jj];
        YG[] *= f[];
      }

      for (int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        YS[] *= f[];
      }

      porosity[] *= f[];

    }


}
