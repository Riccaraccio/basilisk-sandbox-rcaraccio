#define NO_ADVECTION_DIV 1
#define NO_EXPANSION 1
#define SOLVE_TEMPERATURE 1
#define CONST_DIFF 2.05e-5
#define FSOLVE_ABSTOL 1.e-3
#define RADIATION_INTERFACE n.y+0.25*n.x
// #define TURN_OFF_HEAT_OF_REACTION 1
#define KK_CONDUCTIVITY 1

double D0 = 2e-2; //2cm
double H0 = 3e-2; //3cm

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

int maxlevel = 7; int minlevel = 2;
double tend = 600.; //800s

int main() {
  eps0 = 0.4;
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  lambdaS = 0.2;
  TS0 = 300.; TG0 = 723.;
  rhoS = 1200.;

  L0 = 3*H0;
  DT = 1.e-1;

  zeta_policy = ZETA_SHRINK;
  kinfolder = "biomass/Solid-only-2407";
  init_grid(1 << maxlevel);
  run();
}

double solid_mass0 = 0.;
double r0, h0;
FILE *fp;

#define rect(x,y)(fabs(x) < 0.5*H0 && fabs(y) < 0.5*D0)

scalar ls[];
event init (i = 0) {

  fraction (f, superquadric(x, y, 20, 0.5*H0, 0.5*D0));

  if (zeta_policy == ZETA_LEVELSET) {
    double d = min(D0, H0)*0.5*0.25;
    foreach()
      ls[] = -superquadric(x, y, 20, 0.5*H0-d, 0.5*D0-d)*f[];

    f.tracers = list_append (f.tracers, ls);
  }
  
  // fraction (f, rect(x, y));
  //
  // scalar f0[];
  // foreach() {
  //   f0[] = f[];
  //   foreach_dimension()
  //     if (f[-1] == 1 && f[1] == 0 && f[] == 1)
  //       f0[] = 0.9999;
  // }
  // foreach()
  //   f[] = f0[];
  
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
    } else {
      YG[top] = dirichlet (0.);
      YG[right] = dirichlet (0.);
    }
  }

  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  //radial shrinking
  r0 = -HUGE;
  foreach_boundary (left, reduction(max:r0)) 
    if (f[] < 1.-F_ERR && f[] > F_ERR)
      r0 = y + Delta*(f[] - 0.5);

  //axial shrinking
  h0 = -HUGE.;
  foreach_boundary (bottom, reduction(max:h0)) 
    if (f[] < 1.-F_ERR && f[] > F_ERR)
      h0 = x + Delta*(f[] - 0.5);

  if (pid() == 0) {
    char name[80];
    sprintf (name, "OutputData-%d", maxlevel);
    fp = fopen(name, "w");
  }
}

event adapt (i++) {
  adapt_wavelet_leave_interface ({T, u.x, u.y, porosity}, {f},
    (double[]){1.e-1, 1.e-2, 1.e-2, 1e-2}, maxlevel, minlevel, 1);
}

event output (t += 1) {
  if (pid() == 0)
    fprintf (stderr, "%g\n", t);

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();
  
  //radial shrinking
  double r = -HUGE;
  foreach_boundary (left, reduction(max:r)) 
    if (f[] < 1.-F_ERR && f[] > F_ERR)
      r = y + Delta*(f[] - 0.5);

  //axial shrinking
  double h = -HUGE;
  foreach_boundary (bottom, reduction(max:h))
    if (f[] < 1.-F_ERR && f[] > F_ERR)
      h = x + Delta*(f[] - 0.5);

  if (pid() == 0) {
    fprintf(fp, "%g %g %g %g\n", t, solid_mass/solid_mass0, r/r0, h/h0);
    fflush(fp);
  }
}

// event movie (t += 1) {
//   clear();
//   box();
//   view (ty = -0.5, width=1400.);
//   draw_vof ("f");
//   squares("T", min=300, max=750, linear=true, spread=-1);
//   mirror ({1, 0}) {
//     draw_vof("f");
//     squares("zeta", min=0, max=1, linear=true, spread=-1);
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
set yrange [0.5:1.0]

plot "OutputData-6" u 1:3 w l lw 2 lc "red" t "Radial shrinking", \
     "OutputData-6" u 1:4 w l lw 2 lc "web-green" t "Axial shrinking",\
     "data/radial-exp" u 1:2 w p pt 4 lc "red" t "Radial exp", \
     "data/axial-exp" u 1:2 w p pt 4 lc "web-green" t "Axial exp",\
     "data/radial-gentile" u 1:2 w l dt 2 lw 2 lc "red" t "Radial Gentile", \
     "data/axial-gentile" u 1:2 w l dt 2 lw 2 lc "web-green" t "Axial Gentile"
~~~
*/
