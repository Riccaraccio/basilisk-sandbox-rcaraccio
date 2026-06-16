#define NO_ADVECTION_DIV 1

#include "navier-stokes/centered-phasechange.h"
#include "constant-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "darcy.h"
#include "view.h"

const double u0 = 0.02;

u.n[left]    = dirichlet (u0);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = neumann (0.);

int maxlevel = 9, minlevel = 3;
double D0 = 1e-3;
scalar omega[];

scalar fS[];
face vector fsS[];

int main() {
  eps0 = 0.4;
  rho1 = rho2 = 1.;
  mu1 = mu2 = 1.;
  muG = 1.e-5;
  L0 = 8*D0;
  DT = 1e-3;

  rhoS = 1000.;
  rhoG = 1.;

  zeta_policy = ZETA_SHRINK;
  origin (-L0/2, 0);

  init_grid(1 << maxlevel);
  run();
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))
event init(i = 0) {

  vertex scalar phi[];
  foreach_vertex() {
    double theta = atan2(y, x + D0), r = sqrt (sq(x + D0) + sq(y));
    phi[] = (0.25 + 0.05*cos(6*theta))*D0*0.75/0.3 - r;
    phi[] *= -circle(x - D0, y, 0.5*D0);
  }
  fractions (phi, f);

  foreach() {
    u.x[] = u0*(1.- f[]);
    porosity[] = eps0*f[];
  }
}

event chemistry(i++) {
  foreach()
    omega[] = rhoS/10.;
  foreach() {
    if (f[] > F_ERR) {
      porosity[] = porosity[]/f[];
      porosity[] += (omega[]*(1. - porosity[])*(1. - zeta[])/rhoS)*dt;
      porosity[] *= f[];
    }
  }
}

event adapt (i++) {
    adapt_wavelet_leave_interface({porosity, u.x, u.y}, {f}, 
      (double[]){1e-2, 1e-2, 1e-2}, maxlevel, minlevel);
}

event log (t = 0.05; t += 0.1) {
  fprintf(stderr, "%g\n", t);
;
  clear();
  view (ty=-0.5, height = 1080, width = 1080);
  draw_vof("f", "k", lw=1.5);
  cells();
  //vectors("u", "k", scale=3e-4);
  save("movie.mp4");
}

event stop (t = 5);
