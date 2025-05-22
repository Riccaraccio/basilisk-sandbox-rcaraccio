#define NO_ADVECTION_DIV 1
#define NO_EXPANSION 1
#define SOLVE_TEMPERATURE 1
#define CONST_DIFF 2.05e-5
#define FSOLVE_ABSTOL 1.e-3
#define RADIATION_INTERFACE 0.9

double D0 = 8e-3; //2cm
double H0 = 19e-3; //3cm

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
#include "view.h"
#include "superquadric.h"
// #include "balances.h"

u.n[top]     = dirichlet (0.);
u.t[top]     = dirichlet (0.);
p[top]       = neumann (0.);
pf[top]      = neumann (0.);
psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);
psi[right]    = dirichlet (0.);

int maxlevel = 7; int minlevel = 2;
double tend = 300.; //300s

int main() {
  eps0 = 0.2;
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  lambdaS = 0.21;
  lambdaSmodel = L_ANCACOUCE;

  TS0 = 300.; TG0 = 823.;
  rhoS = 1500.;

  L0 = 3.*H0;
  DT = 1e-2; 
  Da = 1e-10;

  zeta_policy = ZETA_SWELLING;
  // kinfolder = "biomass/Solid-only-2407";
  kinfolder = "biomass/dummy-solid";
  init_grid(1 << maxlevel);
  run();
}

double solid_mass0 = 0.;
double r0, h0;
FILE *fp;

#define rect(x,y)(fabs(x) < 0.5*H0 && fabs(y) < 0.5*D0)

event init (i = 0) {

  fraction (f, superquadric(x, y, 20, 0.5*H0, 0.5*D0));

  foreach()
    porosity[] = eps0*f[];
  
  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;

  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 0.3990;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")] = 0.2001;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")] = 0.1077;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")] = 0.0917;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")] = 0.0236;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")] = 0.0112;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]  = 0.0436;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]  = 0.0800;
  
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")]= 1.0000;

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[right] = dirichlet (1.);
      YG[top] = dirichlet (1.);
    } else {
      YG[right] = dirichlet (0.);
      YG[top] = dirichlet (0.);
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
  h0 = -HUGE;
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
    (double[]){1.e-1, 1.e0, 1.e0, 1e0}, maxlevel, minlevel, padding=1);
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

  //core temperature
  double Tcore  = interpolate (TS, 0., 0.);

  //surf T: computed as average of the surface temperature
  double Tsurf  = 0.;
  int count = 0;
  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      Tsurf += TInt[];
      count++;
    }
  } 
  Tsurf /= count;

  if (pid() == 0) {
    fprintf(fp, "%g %g %g %g %g %g\n", t, solid_mass/solid_mass0, r/r0, h/h0, Tcore, Tsurf);
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
set terminal svg size 490,415
set output "mass.svg"
set xlabel "t [s]"
set ylabel "M/M_0"
set yrange [0:1]
set key top right box width 1
set xrange [0:800]
set yrange [0:1.0]
set xtics (0, 200, 400, 600, 800)
set grid
err_mass_time_r = 20
err_mass_time_l = 0.

plot  "OutputData-7-2407-Const" u 1:2 w l lw 2 lc "black" t "Mass", \
      "data/mass-exp" u 1:2:($1-err_mass_time_l):($1+err_mass_time_r) w xerrorbars pt 4 lc "black" t "Mass exp", \
      "data/mass-gentile" u 1:2 w l dt 2 lw 2 lc "black" t "Mass Gentile"
~~~

~~~gnuplot Shrinking
reset
set terminal svg size 490,415
set output "shrinking.svg"
set xlabel "t [s]"
set ylabel "Shrinking factor"
set key bottom left box width 1
set yrange [0.5:1.0]
set xrange [0:800]
set xtics (0, 200, 400, 600, 800)
set grid

err_shrink_time_r = 25
err_shrink_time_l = 0.
err_shrink_value_u = 0.
err_shrink_value_d = 0.05
err_mass_time = 20

plot "OutputData-7-2407-Const" u 1:3 w l lw 2 lc "dark-green" t "Radial", \
     "OutputData-7-2407-Const" u 1:4 w l lw 2 lc "black" t "Axial", \
     "data/radial-exp" u 1:2:($1-err_shrink_time_l):($1+err_shrink_time_r):($2-err_shrink_value_d):($2+err_shrink_value_u) w xyerrorbars pt 4 lw 1.5 lc "dark-green" t "Radial exp", \
     "data/axial-exp"  u 1:2:($1-err_shrink_time_l):($1+err_shrink_time_r):($2-err_shrink_value_d):($2+err_shrink_value_u) w xyerrorbars pt 4 lw 1.5 lc "black" t "Axial exp", \
     "data/radial-gentile" u 1:2 w l dt 2 lw 2 lc "dark-green" t "Radial Gentile", \
     "data/axial-gentile" u 1:2 w l dt 2 lw 2 lc "black" t "Axial Gentile"
~~~
*/
