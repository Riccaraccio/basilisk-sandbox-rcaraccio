#define NO_ADVECTION_DIV 1

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "prop.h"
#include "two-phase.h"
#include "shrinking.h"
#include "view.h"

u.n[top]      = neumann (0.);
u.t[top]      = neumann (0.);
p[top]        = dirichlet (0.);
psi[top]      = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = dirichlet (0.);

u.n[left]    = neumann (0.);
u.t[left]    = neumann (0.);
p[left]      = dirichlet (0.);
psi[left]    = dirichlet (0.);

int maxlevel = 7; int minlevel = 2;
scalar omega[];
double D0 = 1e-2;
double H0 = 2e-2;

int main() {
  eps0 = 0.4;
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;
  muG = 1.e-3;
  TOLERANCE = 1e-4;
  DT = 1e-2;

  rhoS = 100.;
  rhoG = 1.;

  L0 = 6*D0;
  DT = 1e-2;
  origin(-L0/2, 0);

  //0: SHRINK, 
  //1: SWELLING, 
  //2: SMOOTH, 
  //3: SHARP, 
  //4: LEVELSET
  zeta_policy = ZETA_SHRINK;
  init_grid(1 << maxlevel);
  run();
}

#define rect(x,y)(fabs(x) < 0.5*H0 && fabs(y) < 0.5*D0)

event init(i=0) {
  fraction (f, rect(x, y));
  mask (y > 0.5*L0 ? top : none);
  
  foreach()
    porosity[] = eps0*f[];
}

event output (t += 0.1) {
  fprintf(stderr, "%g\n", t);
}

//update the porosity field
event chemistry(i++) {  
  foreach()
    omega[] = 10.;
  foreach() {
    if (f[] > F_ERR) {
      porosity[] = porosity[]/f[];
      porosity[] += (omega[]*(1.-porosity[])*(1.-zeta[])/rhoS)*dt;
      porosity[] *= f[];
    }
  }
}

event movie (t += 0.1) {
  clear();
  view (tx = -0.5*L0, ty = 0);
  draw_vof ("f");
  squares("p", spread = -1);
  mirror ({0, 1}) {
    draw_vof("f");
    squares("(u.x^2 + u.y^2)^0.5", spread = -1);
  }
  save ("movie.mp4");
}

event stop (t = 10);