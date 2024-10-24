#define NO_ADVECTION_DIV 1
#define FSOLVE_ABSTOL 1.e-3

#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "var-prop.h"
#include "two-phase.h"
#include "shrinking.h"
//#include "darcy.h"

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
psi[top] = dirichlet (0.);
ubf.t[top] = neumann (0.);
ubf.n[top] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
psi[right] = dirichlet (0.);
ubf.t[right] = neumann (0.);
ubf.n[right] = neumann (0.);

int maxlevel = 5; int minlevel = 2;
double D0 = 1;
scalar omega[];

int main() {

  rhoS = 681.042; rhoG = 9.75415;
  muG = 2.02391e-5;

  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;
  L0 = 1.*D0;
  eps0 = 0.;
  DT = 5e-3;
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  foreach()
     porosity[] = f[]*eps0;

  zeta_policy = ZETA_SHRINK;
}

event phasechange (i++){
  foreach()
    omega[] = f[] > F_ERR ? 100.: 0.;//fixed interface
}

event logprofile (t += 0.01) {
  fprintf(stderr, "t = %g\n", t);
}
event stop (t=20) {
  return 1;
}
