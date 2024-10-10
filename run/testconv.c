#define ufext uf
//#define VARPROP
#define NO_ADVECTION_DIV 1

#include "navier-stokes/centered-evaporation.h"
#include "two-phase.h"
#include "temperature.h"
#include "evaporation.h"
#include "shrinking.h"
//#include "darcy.h"
#include "view.h"

double TS0 = 300.;
double TG0 = 600.;

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
psi[top] = dirichlet (0.);

//u.n[bottom] = neumann (0.);
//u.t[bottom] = neumann (0.);
//p[bottom] = dirichlet (0.);
//psi[bottom] = dirichlet (0.);
//
//u.n[left] = neumann (0.);
//u.t[left] = neumann (0.);
//p[left] = dirichlet (0.);
//psi[left]  = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
psi[right] = dirichlet (0.);

int maxlevel = 7; int minlevel = 2;
double D0 = 0.7;

scalar omega[];
double solid_mass0;

double lambda1 = 2e-1;
double cp1 =1.e3;

double lambda2 = 1e-1;
double cp2 =1.e3;

int main() {
  rho1 = 100., rho2 = 1.;
  mu1 = 1e-5, mu2 = 1e-5;
  L0 = 1e-1;
  TOLERANCE = 1e-4;
  DT = 1e-1;
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0*L0));

  foreach()
    feps[] = eps0*f[]; //TODO: move to shrinking.h

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0)) {
    solid_mass0 += feps[]*rhos*dv();
  }

  zeta_policy = ZETA_SHARP;

}

event bcs (i=0) {
  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
}

event phasechange (i++){
  foreach()
    omega[] = f[] > F_ERR ? TS[]/f[]/TG0 - 0.5: 0.;
}

event logprofile (t+=10) {
  fprintf (stderr, "t = %g\n", t);
}

event stability (i++,last) {
  face vector us[];
  foreach_face()
    us.x[] = ubf.x[];
  dt = dtnext (timestep (us, dtmax));
}

event movie (t += 10) {
  clear();
  view (ty=-0.5, width=1400.);
  draw_vof ("f", lw=2);
  squares ("T", min=TS0, max=TG0, linear=true);
  mirror ({1.,0.}) {
    vectors ("ubf", scale=10);
    draw_vof ("f", lw=2);
    squares ("feps", min=0., max=eps0, linear=true);
  }
  save ("movie.mp4");
}

event adapt (i++) {
  adapt_wavelet_leave_interface ({T, u.x, u.y, ubf.x, ubf.y}, {f},
      (double[]){1.e-1,1.e-1,1.e-1,1.e-1,1.e-1}, maxlevel, minlevel, 1);
}

#if TRACE> 1
event profiling (i += 1) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif 

event stop (t = 1000) {
  return 1;
}
