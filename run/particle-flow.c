#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define RADIATION_INTERFACE 0.8
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1

#ifndef T_ENV
# define T_ENV 673
#endif

#include "grid/multigrid.h"
#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"

double Uin = 0.0221; //inlet velocity, 300 NL/min over 8.04e-4 m^2 area
double tend = 500.; //simulation time

u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
pf[left]     = neumann (0.);
psi[left]    = dirichlet (0.);

u.n[top]      = dirichlet (0.);
u.t[top]      = dirichlet (0.);
p[top]        = neumann (0.);
pf[top]       = neumann (0.);
psi[top]      = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);
psi[right]    = dirichlet (0.);

int maxlevel = 7; int minlevel = 2;
double D0 = 2e-2;
double solid_mass0 = 0.;

int main() {
  
  lambdaS = 0.1987;
  lambdaSmodel = L_CORBETTA;
  TS0 = 300.; TG0 = T_ENV;
  rhoS = 920;
  eps0 = 0.39;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_REACTION;

#if TREE
  L0 = 4e-2;
#else
  int n_proc = 3;
  size((4e-2)*n_proc);
  dimensions(nx=n_proc, ny=1);
#endif
  origin (-L0/2, 0);

  DT = 1e-1;

  kinfolder = "biomass/dummy-solid";
  shift_prod = true;
  // kinfolder = "biomass/Solid-only-2407";
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  // No ultimate analysis was found, average from smilar chinese hardwood
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.4229;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")]  = 0.1830;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.1759;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0364;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0081;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0362;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0245;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0130;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.1000;

  foreach()
    porosity[] = eps0*f[];

  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  TG[left] = dirichlet (TG0);
  TG[right] = neumann (0);
  TG[top] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[left] = dirichlet (1.);
    } else {
      YG[left] = dirichlet (0.);
    }
  }

  foreach()
    u.x[] = f[] > F_ERR ? 0. : Uin; 
}

event output (t+=1) {
  fprintf (stderr, "%g\n", t);

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();

  //calculate radius
  double radius = cbrt (3./2.*statsf(f).sum);

  //save temperature profile
  double Tcore  = interpolate (T, 0., 0.);
  double Tr2    = interpolate (T, radius/2., 0.);
  double Tsurf  = interpolate (T, radius, 0.);

  fprintf (fp, "%g %g %g %g %g %g\n", t, solid_mass/solid_mass0, Tcore, Tr2, Tsurf, radius/(D0/2.));
  fflush(fp);
}

#if TREE
event adapt (i++) {
  scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];
  adapt_wavelet_leave_interface ({T, u.x, u.y, inert}, {f},
    (double[]){1.e-1, 1.e-1, 1.e-1, 1e-1}, maxlevel, minlevel, 1);
}
#endif


event stop (t = tend);

/** 
~~~gnuplot temperature profiles
reset
set xlabel "t [s]"
set ylabel "temperature [K]"
set key bottom right box width 1
set xrange [0:1000]
set yrange [300:800]
set grid

plot  "OutputData-6" u 1:3 w l lw 2 lc "red" t "Core", \
      "OutputData-6" u 1:4 w l lw 2 lc "web-green" t "R/2", \
      "OutputData-6" u 1:5 w l lw 2 lc "web-blue" t "Surface", \
      "data/corbetta-core.txt"  u 1:2 w p pt 7 lc "red" t "Corbetta core", \
      "data/corbetta-r2.txt"    u 1:2 w p pt 7 lc "web-green" t "Corbetta R/2", \
      "data/corbetta-surf.txt"  u 1:2 w p pt 7 lc "web-blue" t "Corbetta surface"
~~~

~~~gnuplot temperature profiles
reset
set xlabel "t [s]"
set ylabel "temperature [K]"
set key bottom right box width 1
set xrange [0:1000]
set yrange [300:800]
set grid

plot  "OutputData-7" u 1:3 w l lw 2 lc "red" t "Core", \
      "OutputData-7" u 1:4 w l lw 2 lc "web-green" t "R/2", \
      "OutputData-7" u 1:5 w l lw 2 lc "web-blue" t "Surface", \
      "data/biosmoke-core.txt"  u 1:2 w l lw 2 dt 2 lc "red" t "BioSMOKE core", \
      "data/biosmoke-r2.txt"    u 1:2 w l lw 2 dt 2 lc "web-green" t "BioSMOKE R/2", \
      "data/biosmoke-surf.txt"  u 1:2 w l lw 2 dt 2 lc "web-blue" t "BioSMOKE surface"

~~~
~~~gnuplot temperature profiles
reset
set xlabel "t [s]"
set ylabel "temperature [K]"
set key bottom right box width 1
set xrange [0:1000]
set yrange [300:800]
set grid

plot  "OutputData-7" u 1:3 w l lw 2 lc "red" t "Core", \
      "OutputData-7" u 1:4 w l lw 2 lc "web-green" t "R/2", \
      "OutputData-7" u 1:5 w l lw 2 lc "web-blue" t "Surface", \
      "data/corbetta-core.txt"  u 1:2 w p pt 7 lc "red" t "Corbetta core", \
      "data/corbetta-r2.txt"    u 1:2 w p pt 7 lc "web-green" t "Corbetta R/2", \
      "data/corbetta-surf.txt"  u 1:2 w p pt 7 lc "web-blue" t "Corbetta surface", \
      "data/biosmoke-core.txt"  u 1:2 w l lw 2 dt 2 lc "red" t "BioSMOKE core", \
      "data/biosmoke-r2.txt"    u 1:2 w l lw 2 dt 2 lc "web-green" t "BioSMOKE R/2", \
      "data/biosmoke-surf.txt"  u 1:2 w l lw 2 dt 2 lc "web-blue" t "BioSMOKE surface"
~~~
~~~gnuplot temperature profiles
reset
set xlabel "t [s]"
set ylabel "temperature [K]"
set key bottom right box width 1
set xrange [0:1000]
set yrange [300:800]
set grid

plot  "OutputData-7" u 1:3 w l lw 2 lc "red" t "Core", \
      "OutputData-7" u 1:4 w l lw 2 lc "web-green" t "R/2", \
      "OutputData-7" u 1:5 w l lw 2 lc "web-blue" t "Surface", \
      "data/corbetta-core.txt"  u 1:2 w p pt 7 lc "red" t "Corbetta core", \
      "data/corbetta-r2.txt"    u 1:2 w p pt 7 lc "web-green" t "Corbetta R/2", \
      "data/corbetta-surf.txt"  u 1:2 w p pt 7 lc "web-blue" t "Corbetta surface", \
      "data/temperature.ris"    u 1:2 w l lw 2 dt 2 lc "red" t "BioSMOKE core", \
      "data/temperature.ris"    u 1:3 w l lw 2 dt 2 lc "web-green" t "BioSMOKE R/2", \
      "data/temperature.ris"  u 1:11 w l lw 2 dt 2 lc "web-blue" t "BioSMOKE surface"
~~~
**/