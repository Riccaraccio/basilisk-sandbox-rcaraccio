#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1
#define GAS_PHASE_REACTIONS 1
#define RADIATION_TEMP 1273

#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "gravity.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
#include "view.h"

const double Uin = 0.6; //inlet velocity
u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = neumann (0.);

const double tend = 100; //simulation time 300 s
int maxlevel = 9; int minlevel = 3;

const double D0 = 9.5e-3;
double solid_mass0 = 0.;

int main() {

  lambdaSmodel = L_HUANG;
  TS0 = 300.; TG0 = 1050.;
  rhoS = 950.; // biomass density
  eps0 = 0.4; // low, compressed pellet

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_REACTION;

  DT = 1;

  G.x = -9.81;

  kinfolder = "biomass/dummy-solid-gas";
  // kinfolder = "biomass/Red-gas-2507";
  shift_prod = true;

  L0 = 15*D0;
  origin (-L0/2, 0);
  init_grid(1 << maxlevel);

  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

double r0;

event init (i= 0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 0.594; // 93.5% biomass
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")]   = 0.400; // 6.1% moisture
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]     = 0.006; // 0.6% ash

  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.3104;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("GMSW")]  = 0.1175;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.1245;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0004;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0000;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0407;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0004;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0060;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.4000;

  foreach() {
    porosity[] = eps0*f[];
    u.x[] = Uin*(1. - f[]);
  }

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f
  
  TG[left] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);
  TG[right] = neumann (0.);
  TG[bottom] = neumann (0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[left] = dirichlet (0.765);
      YG[left] = dirichlet (0.765);
    } else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      YG[left] = dirichlet (0.235);
      YG[left] = dirichlet (0.235);
    }     
    else {
      YG[left] = dirichlet (0.);
      YG[left] = dirichlet (0.);
    }
  }

  divq_rad = opensmoke_optically_thin;
}

event output (t += 1) {
  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();

  double T_center = interpolate (T, 0., 0.);
  
  // calculate surface temperature as average over interface cells
  double T_surface = avg_interface (TInt, f, F_ERR);
  double TS_avg = avg_interface (TS, f, F_ERR);

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");
  fprintf(fp, "%g %g %g %g %g %g\n", t, solid_mass / solid_mass0, statsf(T).max, T_center, T_surface, TS_avg);
  fflush(fp);
}

#if TREE
event adapt (i++) {
  // scalar fuel = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  scalar fuel = YGList_G[OpenSMOKE_IndexOfSpecies ("TAR")];

  adapt_wavelet_leave_interface ({T, u.x, u.y, fuel, porosity}, {f},
    (double[]){1.e0, 1.e-1, 1.e-1, 1e-1}, maxlevel, minlevel, 2);

  // Unrefine for outflow condition
  unrefine (x > L0*0.4);
}
#endif

event movie (t += 1) {
  scalar XH2O_G = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar XH2O_S = XGList_S[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar XH2O[];
  foreach()
    XH2O[] = XH2O_S[]*f[] + XH2O_G[]*(1. - f[]);

  clear();
  view (theta=0, phi=0, psi=-pi/2., width = 1080, height = 1080);
  squares ("T", min = 300, max = 2000, spread = -1, linear = true);
  isoline ("T", val = statsf(T).max);
  draw_vof ("f", lw = 1.5);
  mirror ({0, 1}) {
    squares ("XH2O", min = 0, max = 1., spread = -1, linear = true);
    draw_vof ("f", lw = 1.5);
  }
  save ("movie.mp4");
}

event stop (t = tend);
/** 
~~~gnuplot Mass
reset
set terminal svg size 400, 400
set output "lu-mass.svg"
set xlabel "Time [s]"
set ylabel "Normalized Mass [-]"
set yrange [0.:1.1]
set xrange [0:100]

plot  "../../data/lu/mass" u 1:(1 - ($2 + $3)/2):(abs(($2 + $3)/2 -$2)) w yerrorbars pt 4 ps 0.8 lc "black" notitle, \
      "OutputData-9" u 1:2 w l lw 2 lc "black" notitle

~~~

~~~gnuplot Temperature evolution
reset
set terminal svg size 400, 400
set output "T-lu.svg"
set xlabel "Time [s]"
set ylabel "Temperature [K]"
set yrange [200:1800]
set xrange [0:100]
set key bottom right

plot  "../../data/lu/T-particle-core" u 1:(($2 + $3)/2):(abs(($2 + $3)/2 -$2)) w yerrorbars pt 4 ps 0.8 lc "black" notitle, \
      "../../data/lu/T-particle-surf" u 1:(($2 + $3)/2):(abs(($2 + $3)/2 -$2)) w yerrorbars pt 4 ps 0.8 lc "red" notitle, \
      "OutputData-9" u 1:4 w l lw 2 lc "black" title "Core", \
      "OutputData-9" u 1:5 w l lw 2 lc "red" title "Surface"


~~~

*/