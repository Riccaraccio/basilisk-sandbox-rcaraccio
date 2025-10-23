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
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"

double tend = 60.; //simulation time

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = dirichlet (0.);

u.n[top]    = neumann (0.);
u.t[top]    = neumann (0.);
p[top]      = dirichlet (0.);
psi[top]    = dirichlet (0.);

int maxlevel = 8; int minlevel = 3;
double D0 = 9.5e-3;
double solid_mass0 = 0.;

int main() {
  
  lambdaS = 0.1987; lambdaG = 0.076;
  cpS = 2200; cpG = 1167;
  #ifdef VARPROP
  lambdaSmodel = L_HUANG;
  #endif
  TS0 = 300.; TG0 = 1050.;
  rhoS = 970; rhoG = 0.674;
  muG = 3.53e-5;
  eps0 = 0.4;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_SWELLING;

  L0 = 10*D0;
  origin (0, 0);

  DT = 2.5e-2;

  kinfolder = "biomass/dummy-solid-gas";
  shift_prod = true;
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

double r0;
event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  foreach()
    porosity[] = eps0*f[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  TG[right] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[right] = dirichlet (0.79);
      YG[top] = dirichlet (0.79);
    } else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      YG[right] = dirichlet (0.21);
      YG[top] = dirichlet (0.21);
    }     
    else {
      YG[right] = dirichlet (0.);
      YG[top] = dirichlet (0.);
    }
  }

  divq_rad = opensmoke_optically_thin;
}

event end_timestep (i++) {
  fprintf (stderr, "%g\n", t);
}

event output (t += 0.1) {

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();

  //calculate radius
  double radius = cbrt (3./2.*statsf(f).sum);

  fprintf (fp, "%g %g %g\n", t, solid_mass/solid_mass0, radius/(D0/2.));

  fflush(fp);
}

#if TREE
event adapt (i++) {
  scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];
  scalar product = YGList_G[OpenSMOKE_IndexOfSpecies ("TAR")];

  scalar umag[];
  foreach()
    umag[] = norm (u);

  adapt_wavelet_leave_interface ({T, u.x, u.y, inert, product, umag}, {f},
    (double[]){1.e-1, 1.e-1, 1.e-1, 1e-1, 1e-1, 1e-1}, maxlevel, minlevel, 2);
}
#endif

event stop (t = tend);

/** 
~~~gnuplot
~~~
**/
