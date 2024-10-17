#define NO_ADVECTION_DIV 1
#define FSOLVE_ABSTOL 1.e-3

#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "var-prop.h"
#include "two-phase.h"
#include "shrinking.h"
//#include "darcy.h"
#include "view.h"

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
psi[top] = dirichlet (0.);
//ubf.t[top] = neumann (0.);
//ubf.n[top] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
psi[right] = dirichlet (0.);
//ubf.t[right] = neumann (0.);
//ubf.n[right] = neumann (0.);

int maxlevel = 6; int minlevel = 2;
double D0 = 1e-3;

scalar omega[];

int main() {
  rho1 = 1.; rho2 = 1.;
  mu1 = 1.; mu2 = 1.;
  L0 = 3.5*D0;
  DT = 5e-3;
  for (maxlevel = 7; maxlevel <=7; maxlevel++) {
    init_grid(1 << maxlevel);
    run();
  }
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  foreach()
     porosity[] = eps0*f[];

  zeta_policy = ZETA_SWELLING;
}

event phasechange (i++){
  foreach()
    omega[] = 1.;
}

//event adapt (i++) {
//  adapt_wavelet_leave_interface ({T, u.x, u.y, ubf.x, ubf.y}, {f},
//      (double[]){1.e0,1.e-1,1.e-1,1.e-1,1.e-1}, maxlevel, minlevel, 1);
//}

event stop (t=1) {
  return 1;
}
