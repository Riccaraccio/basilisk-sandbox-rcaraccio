#define NO_ADVECTION_DIV    1
#define SOLVE_TEMPERATURE   1
#define VARCOEFF            1

#ifndef GAS_VELOCITY
# define GAS_VELOCITY 0.5
#endif

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#define MULTICOMPONENT
#include "constant-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
#include "view.h"

const double Uin = GAS_VELOCITY; // 50 cm/s, avg inlet velocity from Nobili
const double H0 = 4e-3; // 4 mm, from picuture in paper

u.n[top] = dirichlet (0.)*f[] + (1-f[])*neumann (0.);
u.t[top] = dirichlet (0.)*f[] + (1-f[])*neumann (0.);
p[top]   = neumann (0.)*f[] + (1-f[])*dirichlet (0.);
pf[top]  = neumann (0.)*f[] + (1-f[])*dirichlet (0.);
psi[top] = neumann (0.);

u.n[right]   = dirichlet(-Uin);
u.t[right]   = dirichlet(0.);
p[right]     = neumann(0.);
psi[right]   = dirichlet(0.);

int maxlevel = 6; int minlevel = 3;

int main() {
  lambdaS = 0.5*2; lambdaG = 0.10931; //N2 at 2000 K and 1 atm
  cpS = 2000;     cpG = 1284; //N2 at 2000 K and 1 atm

  TS0 = 600.; TG0 = 2000.;
  rhoS = 1078; rhoG = 0.1745; //rhoG equimolar mix N2 and CH2O at 2000 K and 1 atm
  muG = 6.54e-5; // N2 at 2000 K and 1 atm
  eps0 = 0.1;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  L0 = (H0 + 0.01); //10 mm distance of nozzel from fuel surface
  Da = (coord) {1e-14, 1e-14};

  zeta_policy = ZETA_SHRINK;

  shift_prod = true;
  DT = 1e-1;

  kinfolder = "biomass/PE";
  init_grid(1 << maxlevel);

  run();
}

event init(i=0) {
  fraction (f, H0-x);

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("P-HDPE-P")] = 1.;

  foreach()
    porosity[] = eps0*f[];

  foreach()
    u.y[] = Uin*(1. - f[]);

#ifdef SOLVE_TEMPERATURE
  TG[right]   = dirichlet (TG0);
#endif

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2"))
      YG[right] = dirichlet (1.);
    else
      YG[right] = dirichlet (0.);
    
    YG[top] = neumann (0.);
    YG[bottom] = neumann (0.);
  }
}

double delta = H0, RR;
event end_timestep (i++) {
  double delta0 = delta;
  delta = 2./sq(L0)*statsf(f).sum;
  RR = fabs (delta - delta0)/dt*1e3; //mm/s
}

event output (t += 0.01) {
  fprintf (stderr, "%g\n", t);

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double Tint_avg = avg_interface (TInt, f);
  
  double Tliq_avg = 0., fSum = 0.;
  foreach(reduction(+:Tliq_avg) reduction(+:fSum)) {
    Tliq_avg += TS[];
    fSum += f[];
  }
  Tliq_avg /= fSum;

  double OmegaTot = statsf(omega).sum;

  fprintf (fp, "%g %g %g %g %g %g\n", t, delta, RR, Tint_avg, Tliq_avg, OmegaTot);
  fflush(fp);
}

event movie (t += 0.1) {
  clear();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f", lw = 2.5);
  squares ("u.x", spread = -1);
  vectors ("u", scale = 1e-4, level = 5);
  isoline ("u.x", lw = 1.8);
  save ("velocity.mp4");

  clear();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f", lw = 2.5);
  squares ("T", spread = -1, min = 300, max = TG0);
  save ("temperature.mp4");

  clear();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f", lw = 2.5);
  squares ("omega", spread = -1);
  save ("omega.mp4");
}

event adapt (i++) {
  #if TREE
  adapt_wavelet_leave_interface ({T, u.x, u.y}, {f},
    (double[]){1.e-2, 1.e-2, 1.e-2}, maxlevel, minlevel, 1);
  #endif
}

event stop (t = 500);

/** 
~~~gnuplot
~~~
**/