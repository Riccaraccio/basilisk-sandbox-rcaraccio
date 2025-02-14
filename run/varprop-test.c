#define NO_ADVECTION_DIV    1
#define FSOLVE_ABSTOL       1.e-3       //tolerance of fsolve, used in the interface condition
#define SOLVE_TEMPERATURE   1           //wheter to solve the temperature field
#define TURN_OFF_HEAT_OF_REACTION 1     //turn off the heat of reaction
//#define STOP_TRACER_ADVECTION 1         //stop the advection of tracers
//#define EXPLICIT_REACTIONS  1         //explicit reactions
//#define EXPLICIT_DIFFUSION  1         //explicit diffusion
#define FIXED_INT_TEMP    1           //fixed interface temperature

//#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking-varprop.h"

#include "multicomponent-varprop.h"

u.n[top]      = neumann (0.);
u.t[top]      = neumann (0.);
p[top]        = dirichlet (0.);
psi[top]      = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = dirichlet (0.);

int maxlevel = 6; int minlevel = 2;
double D0 = 2.e-2;

int main() {
  lambdaS = 0.2; lambdaG = 0.08;
  cpS = 1600; cpG = 1200;
  rhoS = 1000; rhoG = 1;
  muG = 3.e-5;

  TS0 = 600; TG0 = 600.;
  eps0 = 0.4;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  L0 = 3.5*D0;
  DT = 1e-1;

  kinfolder = "biomass/dummy-solid";
  // kinfolder = "biomass/Solid-only-2003";
  
  //0: SHRINK, 
  //1: SWELLING, 
  //2: SMOOTH, 
  //3: SHARP, 
  //4: LEVELSET
  zeta_policy = ZETA_SHRINK;

  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i = 0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 1.;

  foreach()
    porosity[] = eps0*f[];

#ifdef SOLVE_TEMPERATURE
  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
#endif

  scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];
  scalar tar = YGList_G[OpenSMOKE_IndexOfSpecies ("TAR")];
  // scalar tar = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  scalar water = YGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];

  inert[top] = dirichlet (1.);
  tar[top]  = dirichlet (0.);
  water[top] = dirichlet (0.);
  
  inert[right] = dirichlet (1.);
  tar[right]  = dirichlet (0.);
  water[right] = dirichlet (0.);
}

event output (t+=1) {
  fprintf (stderr, "%g\n", t);
}

event stop (t = 1000);

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