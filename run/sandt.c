#define ufext uf
//#define VARPROP

#include "axi.h"
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

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
psi[right] = dirichlet (0.);

int maxlevel = 7; int minlevel = 5;
double D0 = 1e-3;

scalar omega[];
double solid_mass0;

double lambda1 = 0.124069;
double lambda2 = 0.0295641;

double cp1 = 2244.92;
double cp2 = 1041.52;

double TG0 = 1000.;
double TS0 = 300.;

int main() {
  rho1 = 100., rho2 = 1.;
  mu1 = 1e-5, mu2 = 1e-5;
  L0 = 1e-1;
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

  foreach(){
    //omega[] = circle (x-0.5, y, 0.5*D0) > 0 ? 2. : 0.;
    //omega[] = 2.5*x +2.5;
    omega[] = 1.;
  }

  zeta_policy = ZETA_LEVELSET;

}

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
}

event phasechange (i++){
  foreach()
    omega[] = 10*T[]/TG0;
    //omega[] = T[]/TG0-0.5;
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

event movie (t += 1) {
  clear();
  view (ty=-0.5, width=1400.);
  draw_vof ("f", lw=1.5);
  squares ("T", linear = true, min=TS0, max=TG0);
  mirror ({1.,0.}) {
    vectors ("ubf");
    draw_vof ("f", lw=1.5);
    squares ("zeta*omega");
  }
  save ("movie.mp4");
}

//event adapt (i++) {
// scalar fa[];
// foreach()
//   fa[] = f[];
// adapt_wavelet ({fa,T,u.x,u.y}, (double[]){1.e-3,1e0,1e-3,1e-3}, maxlevel, 4);
//}

event stop (t = 500) {
  return 1;
} 
