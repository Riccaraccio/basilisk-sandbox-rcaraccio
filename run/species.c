#define NO_ADVECTION_DIV 1
#define FSOLVE_ABSTOL 1.e-3

//#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "var-prop"
#include "two-phase.h"
#include "temperature-v.h"
#include "multicomponent_r.h"
#include "shrinking.h"
//#include "darcy.h"
#include "view.h"

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
psi[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
psi[right] = dirichlet (0.);

int maxlevel = 6; int minlevel = 2;
double D0 = 1e-3;

scalar omega[];

double lambda1 = 0.124069;
double lambda2 = 0.0295641;

double cp1 = 2244.92;
double cp2 = 1041.52;

double TG0 = 1000.;
double TS0 = 300.;

int main() {
  rho1 = 681.042, rho2 = 9.75415;
  mu1 = 0.00037446, mu2 = 2.02391e-5;
  L0 = 3.5*D0;

  for (maxlevel = 7; maxlevel <=7; maxlevel++) {
    kinfolder = "biomass/Solid-only-2003";
    init_grid(1 << maxlevel);
    run();
  }
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 1.;

  foreach()
    porosity[] = eps0*f[];
}

event bcs (i=0) {
  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
}

event phasechange (i++){
  foreach()
    omega[] = 0.; //fixed interface
}

//event adapt (i++) {
//  adapt_wavelet_leave_interface ({T, u.x, u.y, ubf.x, ubf.y}, {f},
//      (double[]){1.e0,1.e-1,1.e-1,1.e-1,1.e-1}, maxlevel, minlevel, 1);
//}

event stop (i = 1) {
  return 1;
}
