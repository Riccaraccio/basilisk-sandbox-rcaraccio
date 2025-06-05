#define NO_ADVECTION_DIV    1
#define SOLVE_TEMPERATURE   1
#define CONST_DIFF 2.05e-5

#ifndef TSTEEL
  #define TSTEEL 500 //C
#endif

#ifndef HBLOCK
  #define HBLOCK 400 // W/m^2/K
#endif

#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"

u.n[top]      = neumann (0.);
u.t[top]      = neumann (0.);
p[top]        = dirichlet (0.);
pf[top]       = dirichlet (0.);
psi[top]      = dirichlet (0.);

int maxlevel = 6; int minlevel = 2;
double H0 = 1e-3;

int main() {
  lambdaS = 0.1987; lambdaG = 0.076;
  lambdaSmodel = L_KK;
  cpS = 1600; cpG = 1167;
 
  TS0 = 300.; TG0 = 300.;
  rhoS = 850; rhoG = 0.674;
  muG = 3.53e-5;
  eps0 = 0.3;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  L0 = 2.1415*H0;

  zeta_policy = ZETA_REACTION;

  DT = 5e-4;

  // kinfolder = "biomass/dummy-solid";
  kinfolder = "biomass/Solid-only-2407";
  init_grid(1 << maxlevel);

  // double Tarray[4] = {500, 550, 600, 700};
  // double harray[4] = {400, 400, 400, 400}; // Paulsen
  // double harray[4] = {550, 650, 850, 950}; // Gentile
  // for (int i=0; i<4; i++) {
  //   TSteel = Tarray[i] + 273.15;
  //   hblock = harray[i];
  //   run();
  // }

  fprintf(stderr, "Running for TSteel = %d C, hblock = %d W/m^2/K\n", TSTEEL, HBLOCK);
  run();
}

double solid_mass0, h0;
event init(i=0) {
  fraction (f, H0-y);

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 0.4176;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")] = 0.2102;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")] = 0.0782;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")] = 0.1264;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")] = 0.1053;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")] = 0.0184;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]  = 0.0415;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]  = 0.0024;

  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  foreach()
    porosity[] = eps0*f[];

#ifdef SOLVE_TEMPERATURE
  TG[top] = neumann (0.);
  TG[bottom] = dirichlet (0.);

  TS[top] = dirichlet (0.);
  TS[bottom] = lambda1v.y[] > 0. ? neumann (HBLOCK*(TSTEEL + 273.15 - TS[])/lambda1v.y[]) : neumann (0.);
#endif

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    YG[top] = neumann (0.);
    YG[bottom] = dirichlet (0.);
  }

  solid_mass0 = 0.;
  foreach(reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv();

  h0 = 0.;
  coord p;
  coord regionh[2] = {{L0/2, 0},{L0/2, L0}};
  coord samplingh = {1, (1<<maxlevel)};
  foreach_region (p, regionh, samplingh, reduction(+:h0))
    h0 += f[]; 

  fprintf (stderr, "h0 = %g\n", h0*(L0/(1<<maxlevel)));
}

event output (t+=0.01) {

  fprintf (stderr, "%g\n", t);

  char name[80];
  sprintf(name, "OutputData-%d-%d", maxlevel, TSTEEL);
  static FILE * fp = fopen (name, "w");

  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv(); 

  //also try calculating the height
  double h = 0.;
  coord p;
  coord regionh[2] = {{L0/2, 0},{L0/2, L0}};
  coord samplingh = {1, (1<<maxlevel)};
  foreach_region (p, regionh, samplingh, reduction(+:h))
    h += f[];

  fprintf (fp, "%g %g %g\n", t, solid_mass/solid_mass0, h/h0);
  fflush(fp);
}

event adapt (i++) {
  // scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];

  adapt_wavelet_leave_interface ({T, u.x, u.y}, {f},
    (double[]){1.e-2, 1.e-2, 1.e-2, 1e-2}, maxlevel, minlevel, 1);
}

event stop (t = 0.25);

/** 
~~~gnuplot temperature profiles
reset
set terminal svg size 450,400
set xlabel "Time [s]"
set ylabel "Normalized Volume"
set key top right box width 1
set xrange [0:25]
set yrange [0.4:1.0]
set grid

plot  "OutputData-6-700" u 1:2 w l lw 2 lc "dark-green" t "Sym 700 C", \
      "OutputData-6-600" u 1:2 w l lw 2 lc "blue" t "Sym 600 C", \
      "OutputData-6-550" u 1:2 w l lw 2 lc "orange" t "Sym 550 C", \
      "OutputData-6-500" u 1:2 w l lw 2 lc "black" t "Sym 500 C", \
      "data/exp700" u 1:2 w p pt 4 lc "dark-green" t "Exp. 700 C", \
      "data/exp600" u 1:2 w p pt 4 lc "blue" t "Exp. 600 C", \
      "data/exp550" u 1:2 w p pt 4 lc "orange" t "Exp. 550 C", \
      "data/exp500" u 1:2 w p pt 4 lc "black" t "Exp. 500 C"


~~~
**/