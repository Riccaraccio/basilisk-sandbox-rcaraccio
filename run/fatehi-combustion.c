#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define RADIATION_INTERFACE 0.9
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1
#define GAS_PHASE_REACTIONS 1

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

const double Uin = 0.035; //inlet velocity
u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
psi[left]    = neumann (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = dirichlet (0.);

const double tend = 300; //simulation time 300 s
int maxlevel = 9; int minlevel = 3;

const double D0 = 8e-3, H0 = 8e-3;
double solid_mass0 = 0.;

int main() {
  
  lambdaSmodel = L_HUANG;
  TS0 = 300.; TG0 = 1123.;
  rhoS = 1550;
  eps0 = 0.2; // low, compressed pellet

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_CONST;

  DT = 1;

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
  fraction (f, superquadric (x, y, 20, 0.5*H0, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 0.919;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.061; // 6.1% moisture
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.02; // 2% ash

  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.4358;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("GMSW")]  = 0.2115;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.1352;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0788;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0179;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0168;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0420;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0010;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.0610;

  foreach()
    porosity[] = eps0*f[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f
  
  fprintf(stderr, "Initial solid mass: %g kg\n", solid_mass0);

  TG[left] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);
  TG[right] = neumann (0.);

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
  // Interpolate temperature
  static FILE * fpT = fopen ("TemperatureProfile.dat", "w");
  double T_2mm = interpolate(T, H0/2 + 2e-3, 0.); // 2mm from surface
  double T_4mm = interpolate(T, H0/2 + 4e-3, 0.); // 4mm from surface
  double T_8mm = interpolate(T, H0/2 + 8e-3, 0.); // 8mm from surface
  double T_11mm = interpolate(T, H0/2 + 11e-3, 0.); // 11mm from surface
  double T_15mm = interpolate(T, H0/2 + 15e-3, 0.); // 15mm from surface
  fprintf (fpT, "%g %g %g %g %g %g\n", t, T_2mm, T_4mm, T_8mm, T_11mm, T_15mm);
  fflush(fpT);

  // Interpolate Water vapor mass fraction
  static FILE * fpYH2O = fopen ("YH2OProfile.dat", "w");
  scalar YH2O = YGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
  double YH2O_2mm = interpolate(YH2O, H0/2 + 2e-3, 0.); // 2mm from surface
  double YH2O_4mm = interpolate(YH2O, H0/2 + 4e-3, 0.); // 4mm from surface
  double YH2O_8mm = interpolate(YH2O, H0/2 + 8e-3, 0.); // 8mm from surface
  double YH2O_11mm = interpolate(YH2O, H0/2 + 11e-3, 0.); // 11mm from surface
  fprintf (fpYH2O, "%g %g %g %g %g\n", t, YH2O_2mm, YH2O_4mm, YH2O_8mm, YH2O_11mm);
  fflush(fpYH2O);

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();

  fprintf (fp, "%g %g %g\n", t, solid_mass/solid_mass0, statsf(T).max);
  fflush(fp);
}

#if TREE
event adapt (i++) {
  // scalar fuel = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  scalar fuel = YGList_G[OpenSMOKE_IndexOfSpecies ("TAR")];

  adapt_wavelet_leave_interface ({T, u.x, u.y, fuel, porosity}, {f},
    (double[]){1.e0, 1.e-1, 1.e-1, 1e-1}, maxlevel, minlevel, 2);
}
#endif

event movie (t += 1) {
  scalar YCO2_G = YGList_G[OpenSMOKE_IndexOfSpecies ("CO2")];
  scalar YCO2_S = YGList_S[OpenSMOKE_IndexOfSpecies ("CO2")];
  scalar YCO2[];
  foreach()
    YCO2[] = YCO2_S[] + YCO2_G[];

  clear();
  view (theta=0, phi=0, psi=-pi/2., width = 1080, height = 1080);
  squares ("T", min=300, max=2500, spread=-1, linear=true);
  isoline ("T", val=statsf(T).max);
  draw_vof("f", lw=1.5);
  mirror ({0, 1}) {
    squares ("YCO2", min=0, max=1., spread=-1, linear=true);
    draw_vof("f", lw=1.5);
  }
  save ("movie.mp4");
}

event stop (t = tend);

/** 
~~~gnuplot
~~~
**/
