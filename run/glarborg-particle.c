#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1
#define GAS_PHASE_REACTIONS 1
#define NO_DARCY_CORRECTION 1

#ifndef CASE_NUMBER
# define CASE_NUMBER 0
#endif

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

const double Uin = 0.02; //inlet velocity
u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = neumann (0.);

double tend = 10;
int maxlevel = 9, minlevel = 3;
double solid_mass0 = 0.;
double D0, H0;

double D0_arr[] = {3e-3, 2.08e-3, 1.65e-3, 1.44e-3, 1.31e-3};
double H0_arr[] = {0., 4.16e-3, 6.60e-3, 8.65e-3, 10.48e-3}; 

int main() {
  
  lambdaSmodel = L_LU;
  TS0 = 300.; TG0 = 1473.;
  rhoS = 1000;
  eps0 = 0.4;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_CONST;

  DT = 1e-2;

  G.x = -9.81;

  kinfolder = "biomass/Red-gas-2507";
  shift_prod = true;
  
  if (CASE_NUMBER < 0 || CASE_NUMBER >= 5) {
    fprintf(stderr, "Invalid CASE_NUMBER %d. Must be between 0 and 4.\n", CASE_NUMBER);
    return 1;
  } else {
    D0 = D0_arr[CASE_NUMBER];
    H0 = H0_arr[CASE_NUMBER];
  }

  L0 = 20*D0;
  origin (-L0/2, 0);
  emissivity = emissivity_lu;
  init_grid(1 << maxlevel);

  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

double r0;
event init (i= 0) {
  if (CASE_NUMBER == 0) {
    fraction (f, circle(x, y, 0.5 * D0));
  } else {
    fraction (f, superquadric (x, y, 20, 0.5*H0, 0.5*D0));
  }

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.4752;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")]  = 0.2039;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.1547;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0202;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0000;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0332;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0168;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0030;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.0930;

  foreach()
    porosity[] = eps0*f[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[] - porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  TG[left] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[left] = dirichlet (0.765);
      YG[top] = dirichlet (0.765);
    } else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      YG[left] = dirichlet (0.235);
      YG[top] = dirichlet (0.235);
    }     
    else {
      YG[left] = dirichlet (0.);
      YG[top] = dirichlet (0.);
    }
  }

  divq_rad = opensmoke_optically_thin;
}

event output (t += 0.01) {

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();

  //calculate radius, only meaningful for spherical particles
  double radius = cbrt (3.*statsf(f).sum);

  fprintf (fp, "%g %g %g %g\n", t, solid_mass/solid_mass0, statsf(T).max, radius/(D0/2.));

  fflush(fp);
}

#if TREE
event adapt (i++) {
  scalar fuel = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  scalar oxidiser = YGList_G[OpenSMOKE_IndexOfSpecies ("O2")];

  adapt_wavelet_leave_interface ({T, u.x, u.y, fuel, oxidiser, porosity}, {f},
    (double[]){1.e-1, 1.e-0, 1.e-0, 1e-1, 1e-1}, maxlevel, minlevel, 2);

  // Unrefine for outflow condition
  unrefine (x > L0*0.4);
}
#endif

event movie (t += 0.1) {
  clear();
  view (theta=0, phi=0, psi=-pi/2., width = 1080, height = 1080);
  squares ("T", min = 300, max = 2000, spread = -1, linear = true);
  isoline ("T", val = statsf(T).max);
  isoline ("zmix - zsto", lw = 1.5, lc = {1., 1., 1.});
  draw_vof ("f", lw = 1.5);
  mirror ({0, 1}) {
    squares ("O2_G + O2_S", min = 0., max = 0.235, spread = -1, linear = true);
    isoline ("zmix - zsto", lw = 1.5, lc = {1., 0., 0.});
    draw_vof ("f", lw = 1.5);
  }
  save ("movie.mp4");
}

event stop (t = tend);

/** 
~~~gnuplot
~~~
**/
