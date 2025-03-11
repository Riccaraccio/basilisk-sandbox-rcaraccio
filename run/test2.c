#define NO_ADVECTION_DIV    1
#define FSOLVE_ABSTOL       1.e-3
#define SOLVE_TEMPERATURE   1
#define MULTICOMPONENT 1

#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
// #include "const-prop.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"

u.n[top]      = neumann (0.);
u.t[top]      = neumann (0.);
p[top]        = dirichlet (0.);
psi[top]      = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = dirichlet (0.);

int maxlevel = 6; int minlevel = 2;
double D0 = 2*1.27e-2;
double solid_mass0 = 0.;

int main() {
  lambdaS = 0.1987; lambdaG = 0.076;
  cpS = 1600; cpG = 1167;
  TS0 = 300.; TG0 = 750.;
  rhoS = 850; rhoG = 0.674;
  muG = 3.53e-5;
  eps0 = 0.4;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  // Da = 1e-4;
  
  L0 = 4*D0;
  zeta_policy = ZETA_SHRINK;

  fprintf(stderr, "Using IMPLICIT_DIFFUSION\n");
  DT = 1e-1;

  kinfolder = "biomass/dummy-solid";
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  foreach()
    porosity[] = eps0*f[];

  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f


#ifdef SOLVE_TEMPERATURE
  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
#endif

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[top] = dirichlet (1.);
      YG[right] = dirichlet (1.);
    } else {
      YG[top] = dirichlet (0.);
      YG[right] = dirichlet (0.);
    }
  }
}

event output (t+=1) {
  fprintf (stderr, "%g\n", t);

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();

  //calculate radius
  double radius = pow (3.*statsf(f).sum, 1./3.);

  //save temperature profile
  double Tcore  = interpolate (T, 0., 0.);
  double Tr2    = interpolate (T, radius/2., 0.);
  double Tsurf  = interpolate (T, radius, 0.);

  #ifdef TEMPERATURE_PROFILE
    Tsurf = TemperatureProfile_GetT(t);
  #endif

  fprintf (fp, "%g %g %g %g %g\n", t, solid_mass/solid_mass0, Tcore, Tr2, Tsurf);
  fflush(fp);
}

event adapt (i++) {
  adapt_wavelet_leave_interface ({T, u.x, u.y}, {f},
    (double[]){1.e0, 1.e-1, 1.e-1}, maxlevel, minlevel, 1);
}

event stop (t = 10);