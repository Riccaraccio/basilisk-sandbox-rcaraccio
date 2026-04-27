#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
// #define MOLAR_DIFFUSION 1
// #define FICK_CORRECTED 1
// #define MASS_DIFFUSION_ENTHALPY 1
#define GAS_PHASE_REACTIONS 1
#define NO_DARCY_CORRECTION 1

#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "gravity.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
#include "view.h"

const double Uin = 0.02; //inlet velocity
u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
psi[left]    = dirichlet (0.);
TG[left]     = dirichlet(TG0);

psi[top]     = dirichlet (0.);
TG[top]      = dirichlet(TG0);

u.n[right]   = neumann (0.);
u.t[right]   = neumann (0.);
p[right]     = dirichlet (0.);
psi[right]   = neumann (0.);

double tend = 10;
int maxlevel = 7, minlevel = 3;
double solid_mass0 = 0.;
double D0 = 3e-3;

int main() {
  
  lambdaSmodel = L_CONST;
  lambdaS = 0.21;
  TS0 = 300.; TG0 = 1473.;
  rhoS = 1000;
  eps0 = 0.4;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_CONST;

  DT = 1e-2;

  G.x = -9.81;

  kinfolder = "biomass/dummy-solid-gas";
  //kinfolder = "biomass/Red-gas-2507";
  shift_prod = true;

  L0 = 10*D0;
  origin (-L0/2, 0);
  emissivity = emissivity_lu;
  init_grid(1 << maxlevel);

  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

double r0;
event init (i= 0) {

  scalar f0[];
  fraction (f0, circle (x, y, 0.5*D0));
  
  gas_start[OpenSMOKE_IndexOfSpecies("N2")] = 1.;
  // gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;

  sol_start[OpenSMOKE_IndexOfSolidSpecies("BIOMASS")] = 1.;

  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.4752;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")]  = 0.2039;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.1547;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0202;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0000;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0332;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0168;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0030;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.0930;
    
  solid_mass0 = 0.;
  foreach (reduction(+ : solid_mass0))
    solid_mass0 += f0[]*(1. - eps0)*rhoS*dv(); // Note: (1-e) = (1-ef)!= (1-e)f

  for (int jj = 0; jj < NGS; jj++) {
      scalar YG = YGList_G[jj];
      if (jj == OpenSMOKE_IndexOfSpecies("N2")) {
        YG[left] = dirichlet(1.);
        YG[top] = dirichlet(1.);
      }
      // else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      //   YG[left] = dirichlet (0.235);
      //   YG[top] = dirichlet (0.235);
      // }
      else {
        YG[left] = dirichlet(0.);
        YG[top] = dirichlet(0.);
      }
    }

  if (restore (file = "snapshot", list = all)) {
    fprintf (stderr, "Restart file found!\n");
    restarted = true;
  } else {
    fprintf (stderr, "No restart file found, starting from scratch!\n");

    foreach () {
      f[] = f0[];
      porosity[] = eps0*f[];
    }
  }
}

event output (t += 0.01) {

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[] - porosity[])*rhoS*dv();

  //calculate radius, only meaningful for spherical particles
  double radius = cbrt (3.*statsf(f).sum);

  fprintf (fp, "%g %g %g %g\n", t, solid_mass/solid_mass0, statsf(T).max, radius/(D0/2.));

  fflush(fp);
}

#if TREE
event adapt (i++) {
  // scalar fuel = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  scalar fuel = YGList_G[OpenSMOKE_IndexOfSpecies ("TAR")];
  scalar oxidiser = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];

  adapt_wavelet_leave_interface ({T, u.x, u.y, fuel, oxidiser, porosity}, {f},
    (double[]){1.e-1, 1.e-0, 1.e-0, 1e-1, 1e-1}, maxlevel, minlevel, 2);

  // Unrefine for outflow condition
  unrefine (x > L0*0.4);
}
#endif

event dump (i = 100) {
  dump("snapshot");

  clear();
  squares ("T", linear = true);
  draw_vof ("f", lw = 1.5);
  save ("restart.ppm");
}

event stop (i = 200);

/** 
~~~gnuplot
~~~
**/
