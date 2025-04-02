double D0 = 2e-2;

#define NO_ADVECTION_DIV 1
#define NO_EXPANSION 1
#define SOLVE_TEMPERATURE 1
#define CONST_DIFF 2.05e-5
#define FSOLVE_ABSTOL 1.e-3
#define RADIATION_INTERFACE 1
#define MULTICOMPONENT 1

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "const-prop.h"
// #include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "view.h"
#include "superquadric.h"

u.n[top] = neumann(0.);
u.t[top] = neumann(0.);
p[top] = dirichlet(0.);
psi[top] = dirichlet(0.);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right] = dirichlet(0.);
psi[right] = dirichlet(0.);

int maxlevel = 6, minlevel = 2;
double tend = 500.; //800s

int main() {
  eps0 = 0.4;
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;
  
  lambdaS = 0.6; lambdaG = 0.04;
  cpS = 1800; cpG = 1000.;
  TS0 = 300.; TG0 = 723.;
  muG = 1.e-3;
  L0 = 2.5*D0;
  DT = 1e-2;

  rhoS = 1000.;
  rhoG = 1.;

  zeta_policy = ZETA_SWELLING;

  kinfolder = "biomass/dummy-solid";
  for (maxlevel = 6; maxlevel <= 8; maxlevel++) {
    fprintf(stderr, "maxlevel = %d\n", maxlevel);
    init_grid(1 << maxlevel);
    run();
  }
}

event init(i = 0) {
  fraction (f, superquadric(x, y, 2, 0.5*D0, 0.5*D0));

  foreach()
    porosity[] = eps0*f[];

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj != OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[top] = dirichlet (0.);
      YG[right] = dirichlet (0.);
    }
  }
}

event output (t += 1) {
  fprintf (stderr, "%g\n", t);

  scalar of[];
  foreach() {
    of[] = f[] > F_ERR ? omega[]*f[] : 0.;
  }

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  fprintf (fp, "%g %g %g\n", t, statsf(omega).max, statsf(of).max);
  fflush (fp);
}

event movie(t += 5) {

  scalar of[];
  foreach()
    of[] = f[] > F_ERR ? omega[]*f[] : 0.;

  scalar ze[];
  foreach() {
    ze[] = f[] > F_ERR ? omega[]/statsf(of).max : 0.;
    ze[] = clamp(ze[], 0., 1.);
  }

  clear();
  box();
  view (ty=-0.5, width = 1400.);
  draw_vof("f", lw=2);
  squares ("T", min=300, max=800);
  mirror ({1.,0.}) {
    draw_vof ("f", lw=2);
    squares ("ze", min=0., max=1);
  }

  char name[80];
  sprintf (name, "movie-%d.mp4", maxlevel);
  save (name);
}

event adapt (i++) {
  adapt_wavelet_leave_interface ({T, u.x, u.y, porosity}, {f},
    (double[]){1.e0, 1.e-1, 1.e-1, 1e-1}, maxlevel, minlevel, padding=1);
}

event stop (t = tend);

/*
~~~gnuplot ReactionRate
reset
set xlabel "t [s]"
set ylabel "max(Omega)"
set xrange [0:500]
set key top right box width 1

plot  "OutputData-6" u 1:2 w l lw 2 dt 1 lc "red"       t "LVL 6", \
      "OutputData-7" u 1:2 w l lw 2 dt 1 lc "web-green" t "LVL 7", \
      "OutputData-8" u 1:2 w l lw 2 dt 1 lc "web-blue"  t "LVL 8", \
      "OutputData-6" u 1:3 w l lw 2 dt 2 lc "red"       t "LVL 6", \
      "OutputData-7" u 1:3 w l lw 2 dt 2 lc "web-green" t "LVL 7", \
      "OutputData-8" u 1:3 w l lw 2 dt 2 lc "web-blue"  t "LVL 8"
~~~
*/