#define NO_ADVECTION_DIV 1
#define NO_EXPANSION 1
#define SOLVE_TEMPERATURE 1
#define CONST_DIFF 2.05e-5
#define FSOLVE_ABSTOL 1.e-3
#define RADIATION_INTERFACE 1
#define KK_CONDUCTIVITY 1

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
#include "view.h"
#include "superquadric.h"

u.n[top]      = neumann (0.);
u.t[top]      = neumann (0.);
p[top]        = dirichlet (0.);
psi[top]      = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = dirichlet (0.);

int maxlevel = 8; int minlevel = 2;
double D0 = 2e-2; //2cm
double H0 = 3e-2; //3cm
double tend = 800.; //800s

int main() {
  eps0 = 0.4;
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  lambdaS = 0.2;
  TS0 = 300.; TG0 = 723.;
  rhoS = 1200.;

  L0 = 5*H0;
  DT = 1e-1;
  // origin(-L0/2, 0);

  zeta_policy = ZETA_REACTION;
  kinfolder = "biomass/Solid-only-2407";
  init_grid(1 << maxlevel);
  run();
}

double solid_mass0 = 0.;
double r0, z0;

event init(i=0) {

  fraction (f, superquadric(x, y, 20, 0.5*H0, 0.5*D0));
  // mask (y > 0.5*L0 ? top : none);
  
  foreach()
    porosity[] = eps0*f[];
  
  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 0.4169;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")] = 0.3147;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")] = 0.1039;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")] = 0.0595;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")] = 0.0005;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")] = 0.0616;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]  = 0.0349;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]  = 0.0080;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")]= 0.0000;

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[top] = dirichlet (1.);
      YG[right] = dirichlet (1.);
      // YG[left] = dirichlet (1.);
    } else {
      YG[top] = dirichlet (0.);
      YG[right] = dirichlet (0.);
      // YG[left] = dirichlet (0.);
    }
  }

  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  boundary({TG});
  
  //radial shrinking
  double step = L0/(1 << maxlevel);
  for (double y=step/2; y<L0; y+=step) {
    double f_cell = interpolate(f, 0, y);
    if (f_cell < 1.-F_ERR && f_cell > F_ERR) {
      r0 = y+step*(f_cell - 0.5);
      break;
    }
  }

  //axial shrinking
  for (double x=step/2; x<L0; x+=step) {
    double f_cell = interpolate(f, x, 0);
    if (f_cell < 1.-F_ERR && f_cell > F_ERR) {
      z0 = x+step*(f_cell-0.5);
      break;
    }
  }

  fprintf(stderr, "r0 = %g, z0 = %g\n", r0, z0);
}

event adapt (i++) {
  adapt_wavelet_leave_interface ({T, u.x, u.y}, {f},
    (double[]){1.e0, 1.e-1, 1.e-1}, maxlevel, minlevel, 1);
}

event output (t += 1) {
  fprintf(stderr, "%g\n", t);

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();
  
  //radial shrinking
  double step = L0/(1 << maxlevel);
  double r = 0.;
  for (double y=step/2; y<L0; y+=step) {
    double f_cell = interpolate(f, 0, y);
    if (f_cell < 1.-F_ERR && f_cell > F_ERR) {
      r = y+step*(f_cell - 0.5);
      break;
    }
  }

  //axial shrinking
  double z = 0.;
  for (double x=step/2; x<L0; x+=step) {
    double f_cell = interpolate(f, x, 0);
    if (f_cell < 1.-F_ERR && f_cell > F_ERR) {
      z = x+step*(f_cell-0.5);
      break;
    }
  }

  //print on file
  fprintf(fp, "%g %g %g %g\n", t, solid_mass/solid_mass0, r/r0, z/z0);
  fflush(fp);
}

// event movie (t += 0.1) {
//   clear();
//   view (tx = -0.5*L0, ty = 0);
//   draw_vof ("f");
//   squares("p", spread = -1);
//   mirror ({0, 1}) {
//     draw_vof("f");
//     squares("(u.x^2 + u.y^2)^0.5", spread = -1);
//   }
//   save ("movie.mp4");
// }

event stop (t = tend);

/*
~~~gnuplot Mass profile
reset
set xlabel "t [s]"
set ylabel "M/M_0"
set yrange [0:1]
set key top right box width 1

plot  "OutputData-7" u 1:2 w l lw 2 lc "black" t "Mass profile", \
      "data/mass-exp" u 1:2 w p pt 4 lc "black" t "Exp mass", \
      "data/mass-gentile" u 1:2 w l dt 2 lw 2 lc "black" t "Gentile mass"
~~~

~~~gnuplot Shrinking
reset
set xlabel "t [s]"
set ylabel "Shrinking factor"
set key bottom right box width 1
set yrange [0.5:1]

plot "OutputData-7" u 1:3 w l lw 2 lc "red" t "Radial shrinking", \
     "OutputData-7" u 1:4 w l lw 2 lc "web-green" t "Axial shrinking",\
     "data/radial-exp" u 1:2 w p pt 4 lc "red" t "Radial exp", \
     "data/axial-exp" u 1:2 w p pt 4 lc "web-green" t "Axial exp",\
     "data/radial-gentile" u 1:2 w l dt 2 lw 2 lc "red" t "Radial Gentile", \
     "data/axial-gentile" u 1:2 w l dt 2 lw 2 lc "web-green" t "Axial Gentile"
~~~
*/