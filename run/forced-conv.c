#define NO_ADVECTION_DIV    1
#define FSOLVE_ABSTOL       1.e-3       //tolerance of fsolve, used in the interface condition
#define SOLVE_TEMPERATURE   1           //wheter to solve the temperature field
#define TURN_OFF_HEAT_OF_REACTION 1     //turn off the heat of reaction
//#define TURN_OFF_REACTIONS 1     //turn off the reactions
//#define STOP_TRACER_ADVECTION 1         //stop the advection of tracers
//#define NO_1D_COMPRESSION   1
//#define CONST_DIFF          1         //constant diffusion coefficient
//#define EXPLICIT_REACTIONS  1         //explicit reactions
//#define EXPLICIT_DIFFUSION  1         //explicit diffusion
//#define FIXED_INT_TEMP    1           //fixed interface temperature

//#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "prop.h"
#include "two-phase.h"
#include "multicomponent-t.h"
#include "shrinking.h"
#include "darcy.h"
//#include "balances.h"

//conditions: inflow from the left, outflow to the right
//simmetry at the bottom and at the top
double U0 = 0.2;
u.n[left] = dirichlet(U0);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);
psi[left] = neumann(0.);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right] = dirichlet(0.);
psi[right] = dirichlet(0.);

int maxlevel = 7; int minlevel = 2;
double D0 = 1e-2;

int main() {
  lambdaS = 0.2; lambdaG = 0.08;
  cpS = 1600; cpG = 1200;
  TS0 = 300; TG0 = 600.;
  rhoS = 1000; rhoG = 1;
  muG = 1.e-4;
  eps0 = 0.4;
  Da = 5e-3;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  L0 = 6*D0;
  DT = 1e-2;

  kinfolder = "biomass/dummy-solid";
  origin(-L0/2.5, 0);
  //0: SHRINK, 
  //1: SWELLING, 
  //2: SMOOTH, 
  //3: SHARP, 
  //4: LEVELSET
  zeta_policy = ZETA_SWELLING;
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));
  mask (y > 0.5*L0 ? top : none);

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  foreach()
    porosity[] = eps0*f[];

#ifdef SOLVE_TEMPERATURE
  TG[left] = dirichlet (TG0);
  TG[right] = neumann (0.);
#endif

  scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];
  scalar tar = YGList_G[OpenSMOKE_IndexOfSpecies ("TAR")];
  scalar water = YGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];

  inert[left] = dirichlet (1.);
  tar[left]  = dirichlet (0.);
  water[left] = dirichlet (0.);
  
  inert[right] = neumann (0.);
  tar[right]  = neumann (0.);
  water[right] = neumann (0.);
  
  foreach()
    u.x[] = f[] > F_ERR ? 0. :U0;

#ifdef BALANCES
  mb.print_iter = 100;
#endif
}

event output (t += 0.1) {
  fprintf(stderr, "%g\n", t);
}

event adapt (i++) {
  adapt_wavelet_leave_interface ({T, u.x, u.y}, {f},
     (double[]){1.e-1, 1.e-2, 1.e-2}, maxlevel, minlevel, 1);
}

event stop (t = 100);