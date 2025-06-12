#include "reactors.h"

extern scalar zeta;
extern scalar T;
extern scalar porosity;

#ifndef TURN_OFF_REACTIONS
void ODESolverEXP (odefunction ode, unsigned int neq, double dt, double* y, void* args) {

  double dy[neq];
  ode(y, dt, dy, args);

  for (int jj=0; jj<neq; jj++)
    y[jj] += dt*dy[jj];
}

void RungeKutta45EXP(odefunction ode, unsigned int neq, double dt, double *y, void *args) {

  // Allocate arrays for the k values and temporary y values
  double k1[neq], k2[neq], k3[neq], k4[neq], k5[neq], k6[neq];
  double ytmp[neq];

  // Coefficients for the RK45 method
  const double a2 = 1.0 / 4.0;
  const double a3 = 3.0 / 8.0;
  const double a4 = 12.0 / 13.0;
  const double a6 = 1.0 / 2.0;

  const double b31 = 3.0 / 32.0;
  const double b32 = 9.0 / 32.0;

  const double b41 = 1932.0 / 2197.0;
  const double b42 = -7200.0 / 2197.0;
  const double b43 = 7296.0 / 2197.0;

  const double b51 = 439.0 / 216.0;
  const double b52 = -8.0;
  const double b53 = 3680.0 / 513.0;
  const double b54 = -845.0 / 4104.0;

  const double b61 = -8.0 / 27.0;
  const double b62 = 2.0;
  const double b63 = -3544.0 / 2565.0;
  const double b64 = 1859.0 / 4104.0;
  const double b65 = -11.0 / 40.0;

  // Coefficients for the 5th order solution
  const double c1 = 16.0 / 135.0;
  const double c3 = 6656.0 / 12825.0;
  const double c4 = 28561.0 / 56430.0;
  const double c5 = -9.0 / 50.0;
  const double c6 = 2.0 / 55.0;

  // Step 1: Calculate k1 = f(t, y)
  ode(y, 0, k1, args);

  // Step 2: Calculate k2 = f(t + a2*dt, y + a2*k1*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * a2 * k1[j];
  ode(ytmp, a2 * dt, k2, args);

  // Step 3: Calculate k3 = f(t + a3*dt, y + b31*k1*dt + b32*k2*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * (b31 * k1[j] + b32 * k2[j]);
  ode(ytmp, a3 * dt, k3, args);

  // Step 4: Calculate k4 = f(t + a4*dt, y + b41*k1*dt + b42*k2*dt + b43*k3*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * (b41 * k1[j] + b42 * k2[j] + b43 * k3[j]);
  ode(ytmp, a4 * dt, k4, args);

  // Step 5: Calculate k5 = f(t + a5*dt, y + b51*k1*dt + b52*k2*dt + b53*k3*dt + b54*k4*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * (b51 * k1[j] + b52 * k2[j] + b53 * k3[j] + b54 * k4[j]);
  ode(ytmp, dt, k5, args);

  // Step 6: Calculate k6 = f(t + a6*dt, y + b61*k1*dt + b62*k2*dt + b63*k3*dt + b64*k4*dt + b65*k5*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * (b61 * k1[j] + b62 * k2[j] + b63 * k3[j] + b64 * k4[j] + b65 * k5[j]);
  ode(ytmp, a6 * dt, k6, args);

  // Update y using the 5th order solution
  for (int j = 0; j < neq; j++) {
    y[j] += dt * (c1 * k1[j] + c3 * k3[j] + c4 * k4[j] + c5 * k5[j] + c6 * k6[j]);
    y[j] = y[j] < 0 ? 0 : y[j]; // Ensure non-negativity
  }

  y[NGS+NSS] = clamp (y[NGS+NSS], 0., 1.); // Ensure boundness for porosity
}

event init (i = 0) {
  OpenSMOKE_InitODESolver ();
}

event cleanup (t = end) {
  OpenSMOKE_CleanODESolver ();
}

event reset_sources (i++) {
  foreach()
    omega[] = 0.;
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
#ifdef VARPROP
      data.rhos = rhoSv[];
      data.rhog = rhoGv_S[];
#else
      data.rhos = rhoS;
      data.rhog = rhoG;
#endif
      data.zeta = zeta[];
#ifdef SOLVE_TEMPERATURE
# ifdef VARPROP
      data.cps = cpSv[];
      data.cpg = cpGv_S[];
# else
      data.cps = cpS;
      data.cpg = cpG;
# endif
#endif
      double sources[NEQ];
      data.sources = sources;

      double gasmass[NGS];
      double rhoGvh;
      #ifdef VARPROP
      rhoGvh = rhoGv_S[];
      #else
      rhoGvh = rhoG;
      #endif

      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList_S[jj];
        gasmass[jj] = YG[]*rhoGvh*porosity[];
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
    // ODESolverEXP (batch, NEQ, dt, y0ode, &data);
    RungeKutta45EXP (batch, NEQ, dt, y0ode, &data);
#else //default
    OpenSMOKE_ODESolver (batch, NEQ, dt, y0ode, &data); 
#endif

      double totgasmass = 0;
      for (int jj=0; jj<NGS; jj++) {
        totgasmass += y0ode[jj];
      }

      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList_S[jj];
        YG[] = (totgasmass < 1e-8) ? 0. : y0ode[jj]/totgasmass*f[];
      }

      double totsolidmass = 0;
      for (int jj=0; jj<NSS; jj++) {
        totsolidmass += y0ode[jj+NGS];
      }

      for (int jj=0; jj<NSS; jj++) { 
        scalar YS = YSList[jj];
        YS[] = (totsolidmass < 1e-8) ? 0. : y0ode[jj+NGS]/totsolidmass*f[];
      }

      porosity[] = y0ode[NGS+NSS]*f[];

// #ifdef VARPROP
//       for (int jj=0; jj<NGS; jj++) {
//         scalar DYDtGjj = DYDtG_S[jj];
//         DYDtGjj[] += sources[jj]*cm[];
//       }
// #endif

#ifdef SOLVE_TEMPERATURE
      TS[] = y0ode[NGS+NSS+1]*f[];
// # ifdef VARPROP
//       DTDtS[] += sources[NGS+NSS+1]*cm[];
// # endif
#endif
      omega[] = sources[NGS+NSS];
    }
}
#endif