#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1
#define GAS_PHASE_REACTIONS 1
#define RADIATION_TEMP 300

#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "gravity.h"
#include "superquadric.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
#include "view.h"
#include "flame.h"

const double Uin = 0.28; //inlet velocity
u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = neumann (0.);

const double tend = 60 ; //simulation time 60 s
int maxlevel = 9; int minlevel = 2;

double D0 = 8e-3, H0 = 4e-3;
double solid_mass0 = 0.;

int main() {

  lambdaSmodel = L_LU;
  TS0 = 300.; TG0 = 1408.;
  rhoS = 1200; // not specified in the paper, from Swedish softwood density 
  eps0 = 0.4; // guessed

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_CONST;

  DT = 1;

  G.x = -9.81;

  // kinfolder = "biomass/dummy-solid-gas";
  kinfolder = "biomass/Red-gas-2507";
  shift_prod = true;

  L0 = 20*D0;
  origin (-L0/2, 0);
  init_grid(1 << maxlevel);

  emissivity = emissivity_lu;

  run();
}

double r0;

event init (i= 0) {
  scalar f0[];
  fraction (f0, superquadric (x, y, 20, 0.5*H0, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.0272;
  gas_start[OpenSMOKE_IndexOfSpecies ("CO2")] = 0.8371;
  gas_start[OpenSMOKE_IndexOfSpecies ("H2O")] = 0.1357;

  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 0.935; // 93.5% biomass
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.061; // 6.1% moisture
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.004; // 0.4% ash

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.5690;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("GMSW")]  = 0.2493;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.0722;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0022;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0000;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0166;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0025;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0017;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.0865;

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += f0[]*(1. - eps0)*rhoS*dv();

  TG[left] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);
  TG[right] = neumann (0.);
  TG[bottom] = neumann (0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("CO2")) {
      YG[left] = dirichlet (0.8371);
      YG[top] = dirichlet (0.8371);
    } else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      YG[left] = dirichlet (0.0272);
      YG[top] = dirichlet (0.0272);
    } else if (jj == OpenSMOKE_IndexOfSpecies ("H2O")) {
      YG[left] = dirichlet (0.1357);
      YG[top] = dirichlet (0.1357);
    } else {
      YG[left] = dirichlet (0.);
      YG[top] = dirichlet (0.);
    }
  }

  if (restore (file = "last-snapshot", list = all)) {
    fprintf (stderr, "Restart file found!\n");
    restarted = true;
  } else {
    fprintf (stderr, "No restart file found, starting from scratch!\n");

    foreach() {
      f[] = f0[];
      porosity[] = eps0*f[];
    }
  }
}

event output (t += 0.1) {
  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");
  fprintf(fp, "%g %g %g\n", t, solid_mass/solid_mass0, statsf(T).max);
  fflush(fp);
}

#if TREE
event adapt (i++) {
  scalar fuel = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  scalar oxidiser = YGList_G[OpenSMOKE_IndexOfSpecies ("O2")];

  adapt_wavelet_leave_interface ({T, u.x, u.y, fuel, oxidiser, porosity}, {f},
    (double[]){1.e-1, 1.e-1, 1.e-1, 1e-1, 1e-1}, maxlevel, minlevel, 2);

  // Unrefine for outflow condition
  unrefine (x > L0*0.4);
}
#endif

event movie (t += 1) {
  scalar XCO_G = XGList_G[OpenSMOKE_IndexOfSpecies ("CO")];
  scalar XCO_S = XGList_S[OpenSMOKE_IndexOfSpecies ("CO")];
  scalar XCO[];
  foreach()
    XCO[] = XCO_S[]*f[] + XCO_G[]*(1. - f[]);

  clear();
  view (theta=0, phi=0, psi=-pi/2., width = 1080, height = 1080);
  squares ("T", min = 300, max = 2000, spread = -1, linear = true);
  isoline ("T", val = statsf(T).max);
  isoline ("zmix - zsto", lw = 1.5, lc = {1., 1., 1.});
  draw_vof ("f", lw = 1.5);
  mirror ({0, 1}) {
    squares ("XCO", min = 0, max = 0.5, spread = -1, linear = true);
    isoline ("zmix - zsto", lw = 1.5, lc = {1., 1., 1.});
    draw_vof ("f", lw = 1.5);
  }
  save ("movie.mp4");
}

event dump (t += 1) {
  dump ("last-snapshot");
}

event stop (t = tend);

/** 
~~~gnuplot Mass
set terminal svg size 450, 450
set output "mass.svg"
set xlabel "Time [s]"
set ylabel "Normalized solid mass [-]"
set xrange [0:60]
set yrange [0:1.05]

plot "cluster/F3/rho1100/OutputData-9"  u 1:2 w l lw 2 lc "black" t "rho 1100", \
     "cluster/F3/rho1500/OutputData-9"  u 1:2 w l lw 2 lc "red" t "rho  1500", \
     "cluster/gasification/OutputData-9"  u 1:2 w l lw 2 lc "blue" t "rho 1250", \
     "../../data/gasification/F3-mass" u 1:2 w p pt 4 lc "black" notitle
~~~
**/
