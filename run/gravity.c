#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define RADIATION_INTERFACE 0.9

#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "gravity.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"

double tend = 600.; //simulation time
double Uin = 0.3; //inlet velocity

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = neumann (0.);
pf[right]     = dirichlet (0.);

psi[top]     = dirichlet (0.);

u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
pf[left]     = neumann (0.);
psi[left]    = neumann (0.);

int maxlevel = 8; int minlevel = 4;
double D0 = 1e-1;
double solid_mass0 = 0., moisture0 = 0.;

int main() {
  
  lambdaS = 0.1987;
  lambdaSmodel = L_HUANG;
  TS0 = 600.; TG0 = 900.;
  rhoS = 920;
  eps0 = 0.4;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_SHRINK;

  G.x = -9.81;

  L0 = 10*D0;
  origin (-L0/2, 0);

  DT = 1e-1;

  kinfolder = "biomass/dummy-solid";
  shift_prod = true;
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

double r0;
event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  foreach()
    porosity[] = eps0*f[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  TG[right] = neumann (0.);
  TG[top] = dirichlet (TG0);
  TG[left] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[right] = dirichlet (1.);
      YG[top] = dirichlet (1.);
    } else {
      YG[right] = dirichlet (0.);
      YG[top] = dirichlet (0.);
    }
  }
}

#if TREE
event adapt (i++) {
  fprintf(stderr, "%g\n", t);
  scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];
  adapt_wavelet_leave_interface ({T, u.x, u.y, inert}, {f},
    (double[]){1.e-2, 1.e-2, 1.e-2, 1e-2}, maxlevel, minlevel, 1);

  unrefine (x < -L0/3);
}
#endif

event stop (t = tend);

/** 
~~~gnuplot
~~~
**/