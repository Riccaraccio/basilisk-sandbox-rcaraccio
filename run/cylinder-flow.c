#define NO_ADVECTION_DIV 1
#define NO_EXPANSION 1
#define SOLVE_TEMPERATURE 1
#define CONST_DIFF 2.05e-5
#define FSOLVE_ABSTOL 1.e-3
#define RADIATION_INTERFACE 1
#define KK_CONDUCTIVITY 1
// #define EXPLICIT_REACTIONS 1

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

double Uin = 0.11; //0.11 m/s or 2Nl/min over 8e-4 m^2 area

u.n[left]     = dirichlet (Uin);
u.t[left]     = dirichlet (0.);
p[left]       = neumann (0.);
psi[left]     = dirichlet (0.);

u.n[top]     = dirichlet (0.);
u.t[top]     = dirichlet (0.);
p[top]       = neumann (0.);
psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = dirichlet (0.);

int maxlevel = 7; int minlevel = 2;
double tend = 800.; //800s

int main() {
  eps0 = 0.4;
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  TS0 = 300; TG0 = 723.;
  rhoS = 1200.;

  L0 = 3.5*H0;
  DT = 1e-1;
  origin(-L0/3, 0.);
  Da = 1e-10;

  zeta_policy = ZETA_SHRINK;
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
  mask (y > 6e-3+D0/2 ? top : none);

  fraction (f, superquadric(x, y, 20, 0.5*H0, 0.5*D0));

  foreach()
    porosity[] = eps0*f[];
  
  foreach()
    u.x[] = f[] > F_ERR ? 0. : Uin;

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;

  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 0.4169;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")] = 0.3147;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")] = 0.1039;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")] = 0.0595;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")] = 0.0005;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")] = 0.0616;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]  = 0.0349;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]  = 0.0080;
  
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")]= 1.0000;

  TG[top] = dirichlet (TG0);
  TG[left] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[left] = dirichlet (1.);
    } else {
      YG[left] = dirichlet (0.);
    }
  }

  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  //radial shrinking
  r0 = 0;
  coord pr;
  coord regionr[2] = {{0,0},{0,1.6e-2}};
  coord samplingr = {1, 1.6e-2/(L0/(1<<maxlevel))};
  foreach_region (pr, regionr, samplingr, reduction(+:r0))
    r0 += f[];

  //axial shrinking
  h0 = 0;
  coord ph;
  coord regionh[2] = {{-L0/3,0},{L0/3*2,0}};
  coord samplingh = {(1<<maxlevel), 1};
  foreach_region (ph, regionh, samplingh, reduction(+:h0))
    h0 += f[];

  if (pid() == 0) {
    char name[80];
    sprintf (name, "OutputData-%d", maxlevel);
    fp = fopen(name, "w");
  }
}

// event adapt (i++) {
//   scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];
//   adapt_wavelet_leave_interface ({T, u.x, u.y, porosity, inert}, {f},
//     (double[]){1.e-2, 1.e-2, 1.e-2, 1e0, 1e-1}, maxlevel, minlevel, padding=1);
// }

event output (t += 1) {
  if (pid() == 0)
    fprintf (stderr, "%g\n", t);

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();
  
  //radial shrinking
  double r = 0.;
  coord p;
  coord regionr[2] = {{0,0},{0,1.6e-2}};
  coord samplingr = {1, 1.6e-2/(L0/(1<<maxlevel))};
  foreach_region (p, regionr, samplingr, reduction(+:r))
    r += f[];
  
  //axial shrinking
  double h = 0.;
  coord regionh[2] = {{-L0/3,0},{L0/3*2,0}};
  coord samplingh = {(1<<maxlevel), 1};
  foreach_region (p, regionh, samplingh, reduction(+:h))
    h += f[];

  if (pid() == 0) {
    fprintf(fp, "%g %g %g %g\n", t, solid_mass/solid_mass0, r/r0, h/h0);
    fflush(fp);
  }
}

event movie (t += 10) {
  clear();
  // box();
  view (width=1000., height=900.);
  draw_vof ("f");
  squares("T", min=300, max=750, linear=true, spread=-1);
  mirror ({0, 1}) {
    draw_vof("f");
    // squares("zeta", min=0, max=1, linear=true, spread=-1);
    vectors ("u", scale=5e-4);
  }
  save ("movie.mp4");
}

event stop (t = tend);

/*
~~~gnuplot Mass profile
reset
set terminal svg size 400,350
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

plot  "OutputData-7-red-o" u 1:2 w l lw 2 lc "black" t "Mass", \
      "data/mass-exp" u 1:2:($1-err_mass_time_l):($1+err_mass_time_r) w xerrorbars pt 4 lc "black" t "Mass exp"
      #"data/mass-gentile" u 1:2 w l dt 2 lw 2 lc "black" t "Mass Gentile"
~~~

~~~gnuplot Shrinking
reset
set terminal svg size 400,350
set output "shrinking.svg"
set xlabel "t [s]"
set ylabel "Shrinking factor"
set key bottom left box width 1
set yrange [0:1.0]
set xrange [0:800]
set xtics (0, 200, 400, 600, 800)
set grid

err_shrink_time_r = 25
err_shrink_time_l = 0.
err_shrink_value_u = 0.
err_shrink_value_d = 0.05
err_mass_time = 20

plot "OutputData-7-fail" u 1:3 w l lw 2 lc "dark-green" t "Radial", \
     "OutputData-7-fail" u 1:4 w l lw 2 lc "black" t "Axial", \
     "data/radial-exp" u 1:2:($1-err_shrink_time_l):($1+err_shrink_time_r):($2-err_shrink_value_d):($2+err_shrink_value_u) w xyerrorbars pt 4 lw 1.5 lc "dark-green" t "Radial exp", \
     "data/axial-exp"  u 1:2:($1-err_shrink_time_l):($1+err_shrink_time_r):($2-err_shrink_value_d):($2+err_shrink_value_u) w xyerrorbars pt 4 lw 1.5 lc "black" t "Axial exp"
     #"data/radial-gentile" u 1:2 w l dt 2 lw 2 lc "dark-green" t "Radial Gentile", \
     #"data/axial-gentile" u 1:2 w l dt 2 lw 2 lc "black" t "Axial Gentile"
~~~
*/
