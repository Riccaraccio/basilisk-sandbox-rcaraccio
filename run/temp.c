#define NO_ADVECTION_DIV 1

#ifndef MAXLEVEL
# define MAXLEVEL 9
#endif

#include "navier-stokes/centered-phasechange.h"
#include "constant-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "view.h"
#include "superquadric.h"

u.n[top] = neumann(0.);
u.t[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);
psi[top] = dirichlet(0.);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
psi[right] = dirichlet(0.);

u.n[left] = neumann(0.);
u.t[left] = neumann(0.);
p[left] = dirichlet(0.);
pf[left] = dirichlet(0.);
psi[left] = dirichlet(0.);

int maxlevel = MAXLEVEL, minlevel = 4;
double D0 = 1;
scalar omega[];

scalar fS[];
face vector fsS[];

int main() {
  eps0 = 0.4;
  rho1 = rho2 = 1.;
  mu1 = mu2 = 1.;
  muG = 1.e-3;
  L0 = 10*D0;
  DT = 1e-3;

  rhoS = 100.;
  rhoG = 1.;

  origin (-L0/2, 0.);
  init_grid (1 << maxlevel);
  run();
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))
event init(i = 0) {
  //fraction (f, circle(x, y, 0.5*D0));
  fraction (f, superquadric (x, y, 20, 0.5*D0, 0.5*D0));

  foreach()
    porosity[] = eps0*f[];

  clear();
  view (tx = -L0/2);
  cells();
  box();
  draw_vof("f", "k", lw = 2);
  save ("initial.png");

}

event chemistry(i++) {
  foreach()
    omega[] = rhoS/10.;
  foreach() {
    if (f[] > F_ERR) {
      porosity[] = porosity[]/f[];
      porosity[] += (omega[]*(1.-porosity[])*(1.-zeta[])/rhoS)*dt;
      porosity[] *= f[];
    }
  }
}

event adapt (i++) {
    adapt_wavelet_leave_interface({porosity, u.x, u.y}, {f}, 
      (double[]){1e-2, 1e-2, 1e-2}, maxlevel, minlevel, 2);
}

event log (t += 0.1) {
  fprintf(stderr, "%g\n", t);
}


event stop (t = 10);
