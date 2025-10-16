#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define RADIATION_INTERFACE 0.9
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1

#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "gravity.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"

double tend = 600.; //simulation time

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = dirichlet (0.);

u.n[left]    = neumann (0.);
u.t[left]    = neumann (0.);
p[left]      = dirichlet (0.);
psi[left]    = dirichlet (0.);

u.n[top]    = neumann (0.);
u.t[top]    = neumann (0.);
p[top]      = dirichlet (0.);
psi[top]    = dirichlet (0.);

u.n[bottom] = dirichlet (0.);
u.t[bottom] = dirichlet (0.);
p[bottom]   = neumann (0.);


int maxlevel = 8; int minlevel = 4;
double D0 = 2e-2;
double solid_mass0 = 0., moisture0 = 0.;

int main() {
  
  lambdaS = 0.1987;
  lambdaSmodel = L_HUANG;
  TS0 = 650.; TG0 = 700.;
  rhoS = 920;
  eps0 = 0.4;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_SWELLING;

  G.x = -9.81;

  L0 = 10*D0;
  origin (-L0/2, 0);

  DT = 1e-1;

  kinfolder = "biomass/dummy-solid";
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
      YG[right] = dirichlet (1.);
      YG[top] = dirichlet (1.);
    } else {
      YG[right] = dirichlet (0.);
      YG[top] = dirichlet (0.);
    }
  }
}

event output (t += 1) {
  fprintf (stderr, "%g\n", t);

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
  adapt_wavelet_leave_interface ({T, u.x, u.y, inert}, {f},
    (double[]){1.e-2, 1.e-2, 1.e-2, 1e-2}, maxlevel, minlevel, 1);

  unrefine (x < -L0/3);
}
#endif


event stop (t = tend);

/** 
~~~gnuplot
~~~
**/