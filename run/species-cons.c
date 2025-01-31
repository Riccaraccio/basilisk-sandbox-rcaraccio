#define NO_ADVECTION_DIV    1
#define FSOLVE_ABSTOL       1.e-3       //tolerance of fsolve, used in the interface condition
#define SOLVE_TEMPERATURE   1           //wheter to solve the temperature field
#define TURN_OFF_HEAT_OF_REACTION 1     //turn off the heat of reaction
//#define STOP_TRACER_ADVECTION 1         //stop the advection of tracers
//#define NO_2D_COMPRESSION   1
#define CONST_DIFF          1         //constant diffusion coefficient
//#define EXPLICIT_REACTIONS  1         //explicit reactions
//#define EXPLICIT_DIFFUSION  1         //explicit diffusion
#define FIXED_INT_TEMP    1           //fixed interface temperature

//#include "temperature-profile.h"
//#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "prop.h"
#include "two-phase.h"
#include "multicomponent-t.h"
#include "shrinking.h"
//#include "darcy.h"
#include "balances.h"

u.n[top]      = neumann (0.);
u.t[top]      = neumann (0.);
p[top]        = dirichlet (0.);
psi[top]      = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = dirichlet (0.);

int maxlevel = 6; int minlevel = 2;
double D0 = 2*1.27e-2;
int zetamodel = 0;
int amr = 0;

int main() {
  lambdaS = 0.2; lambdaG = 0.08;
  cpS = 1600; cpG = 1200;
#ifdef TEMPERATURE_PROFILE
  TS0 = 300.; TG0 = 300.;
#else
  TS0 = 600; TG0 = 600.;
#endif
  rhoS = 1000; rhoG = 1;
  muG = 3.e-5;
  eps0 = 0.4;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  L0 = 3.5*D0;

#ifdef EXPLICIT_DIFFUSION
  fprintf(stderr, "Using EXPLICIT_DIFFUSION\n");
  DT = 1e-1;
#else   
  fprintf(stderr, "Using IMPLICIT_DIFFUSION\n");
  DT = 1e-1;
#endif

  kinfolder = "biomass/dummy-solid";

// amr = 0;
// for (maxlevel=5; maxlevel<=7; maxlevel++) {
//   for (zetamodel=0; zetamodel <= 1; zetamodel++) {
//         fprintf(stderr, "Running maxlevel %d zetamodel %d, amr %d\n",
//           maxlevel, zetamodel, amr);
//         init_grid(1 << maxlevel);
//         run();
//     }
//   }

  maxlevel = 7;
  zetamodel = 1;
  amr = 0;
  
  //0: SHRINK, 
  //1: SWELLING, 
  //2: SMOOTH, 
  //3: SHARP, 
  //4: LEVELSET
  zeta_policy = zetamodel;

  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i = 0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  foreach()
    porosity[] = eps0*f[];

#ifdef SOLVE_TEMPERATURE
#ifndef TEMPERATURE_PROFILE
  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
#endif
#endif

  scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];
  scalar tar = YGList_G[OpenSMOKE_IndexOfSpecies ("TAR")];
  scalar water = YGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];

  inert[top] = dirichlet (1.);
  tar[top]  = dirichlet (0.);
  water[top] = dirichlet (0.);
  
  inert[right] = dirichlet (1.);
  tar[right]  = dirichlet (0.);
  water[right] = dirichlet (0.);

#ifdef BALANCES
  //specify the output file name
  sprintf(mb.name, "test");
#endif
}

event output (t+=1) {
  fprintf (stderr, "%g\n", t);
}

event adapt (i++) {
  if (amr > 0) {
  scalar tar = YGList_G[OpenSMOKE_IndexOfSpecies ("TAR")];
  adapt_wavelet_leave_interface ({tar, T, u.x, u.y}, {f},
     (double[]){1.e-2, 1.e-1, 1.e-2, 1.e-2}, maxlevel, minlevel, 1);
  }
}

#if TRACE > 1
  event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif

event stop (t = 500);

/** 
~~~gnuplot Total Mass Balance
reset
set terminal pdfcairo enhanced color font ',10'
set title "Total Mass Balance"
set output 'plot0.pdf'
set xlabel 't'
set ylabel "M_{s}/M_{s0}"
set grid
plot "balances-7" u 1:2 w l t "Total solid mass", \
      "balances-7" u 1:3 w l t "Total gas mass", \
      "balances-7" u 1:4 w l t "Total mass"
~~~

~~~gnuplot Solid Mass Balance
reset
set terminal pdfcairo enhanced color font ',10'
set title "Solid Mass Balance"
set output 'plot1.pdf'
set xlabel 't'
set ylabel "M_{s}/M_{s0}"
set grid
plot "balances-7" u 1:2 w l t "Total solid mass", \
      "balances-7" u 1:8 w l t "Biomass", \
      "balances-7" u 1:9 w l t "Char", \
      "balances-7" u 1:($8+$9) w l t "Bio+Char"
~~~

~~~gnuplot Gas Mass Balance
reset
set terminal pdfcairo enhanced color font ',10'
set title "Gas Mass Balance"
set output 'plot2.pdf'
set xlabel 't'
set ylabel "M_{s}/M_{s0}"
set grid
set key top left
plot "balances-7" u 1:3 w l t "Total gas exited", \
      "balances-7" u 1:6 w l t "TAR", \
      "balances-7" u 1:7 w l t "H2O", \
      "balances-7" u 1:($6+$7+$5) w l t "TAR+H2O+N2"
~~~
**/