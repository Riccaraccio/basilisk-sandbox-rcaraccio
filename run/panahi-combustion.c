#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define RADIATION_INTERFACE 0.9
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1
#define GAS_PHASE_REACTIONS 1
#define RADIATION_TEMP 1400.

#ifndef ASPECT_RATIO
# define ASPECT_RATIO 4. // aspect ratio of the particle
#endif

#ifndef IS_SPHERE
# define IS_SPHERE false // true: sphere, false: cylinder
#endif

#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "superquadric.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
#include "view.h"

// double Uin = 0.01; //free fall velocity = 2*g*r^2/9/nu + gas velocity
// 2D shows no difference in results with Uin = 0. So we set Uin = 0.
u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = dirichlet (0.);

u.n[top]    = neumann (0.);
u.t[top]    = neumann (0.);
p[top]      = dirichlet (0.);
psi[top]    = dirichlet (0.);

double tend = 600e-3; //simulation time 300 ms
int maxlevel = 9; int minlevel = 3;
double aspect_ratio = ASPECT_RATIO; // aspect ratio of the particle
bool is_sphere = IS_SPHERE; // true: sphere, false: cylinder
double average_mass = 5.3e-8; //average biomass particle mass 53 mg

double D0, H0;
double solid_mass0 = 0.;

int main() {
  
  lambdaS = 0.1987; lambdaG = 0.076;
  cpS = 2200; cpG = 1167;
  #ifdef VARPROP
  lambdaSmodel = L_HUANG;
  #endif
  TS0 = 300.; TG0 = 1350.;
  rhoS = 1500; rhoG = 0.674;
  muG = 3.53e-5;
  eps0 = 0.7;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_CONST;

  DT = 2e-5;

  kinfolder = "biomass/dummy-solid-gas";
  // kinfolder = "biomass/Red-gas-2507";
  shift_prod = true;

  // Calculate D0 and H0 based on average mass and aspect ratio
  double particle_volume = average_mass / ((1.-eps0)*rhoS);
  if (is_sphere && aspect_ratio != 1.) {
    fprintf (stderr, "Error: aspect ratio must be 1 for a sphere.\n");
    return 1;
  }

  if (is_sphere) {
    D0 = cbrt (particle_volume*6/M_PI);
    fprintf (stderr, "Calculated sphere diameter D0: %e m\n", D0*1e3);
  } else {
    D0 = cbrt (4*particle_volume/(M_PI*aspect_ratio));
    H0 = aspect_ratio * D0;
    fprintf (stderr, "Calculated cylinder diameter D0: %e mm, height H0: %e mm\n", D0*1e3, H0*1e3);
  }

  L0 = 9*max(D0, H0);
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

double r0;
event init (i= 0) {
  if (is_sphere)
    fraction (f, circle (x, y, 0.5*D0));
  else
    fraction (f, superquadric (x, y, 20, 0.5*H0, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 0.93;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.05; // 5% moisture
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.02; // 2% ash

  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.4832;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")]  = 0.1200;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.1652;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0003;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0000;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.1668;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0005;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0190;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.0450;

  foreach()
    porosity[] = eps0*f[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f
  
  fprintf(stderr, "Initial solid mass: %g kg\n", solid_mass0);

  TG[right] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[right] = dirichlet (0.765);
      YG[top] = dirichlet (0.765);
    } else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      YG[right] = dirichlet (0.235);
      YG[top] = dirichlet (0.235);
    }     
    else {
      YG[right] = dirichlet (0.);
      YG[top] = dirichlet (0.);
    }
  }

  divq_rad = opensmoke_optically_thin;
}

event output (t += 1e-3) {

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();

  double BIOMASS_mass = 0.;
  scalar YBIOMASS = YSList[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")];
  foreach (reduction(+:BIOMASS_mass))
    if (f[] > F_ERR)
      BIOMASS_mass += YBIOMASS[]*(1. - porosity[]/f[])*rhoS*dv();

  double CHAR_mass = 0.;
  scalar YCHAR = YSList[OpenSMOKE_IndexOfSolidSpecies ("CHAR")];
  foreach (reduction(+:CHAR_mass))
    if (f[] > F_ERR)
      CHAR_mass += YCHAR[]*(1. - porosity[]/f[])*rhoS*dv();
      
  //calculate radius, only meaningful for spherical particles
  double radius = cbrt (3.*statsf(f).sum);

  fprintf (fp, "%g %g %g %g %g %g\n", t, solid_mass/solid_mass0, radius/(D0/2.), statsf(T).max, BIOMASS_mass/solid_mass0, CHAR_mass/solid_mass0);

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

event movie (t += 5e-3) {
  clear();
  view (ty = -0.5, tx = -0.5, width = 1080, height = 1080);
  squares ("T", min=300, max=2500, spread=-1);
  isoline ("T", val=statsf(T).max);
  draw_vof("f");
  save ("movie.mp4");
}

event stop (t = tend);

/** 
~~~gnuplot
~~~
**/
