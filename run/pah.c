#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1
#define TURN_OFF_HEAT_OF_REACTION 1
#define GAS_PHASE_REACTIONS 1

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "view.h"

u.n[top]      = neumann (0.);
u.t[top]      = neumann (0.);
p[top]        = dirichlet (0.);
pf[top]       = dirichlet (0.);
psi[top]      = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);
psi[right]    = dirichlet (0.);

const double tend = 60 ; //simulation time 60 s
int maxlevel = 7; int minlevel = 2;
double D0 = 2e-2;
double solid_mass0 = 0.;

int main () {
    lambdaSmodel = L_LU;
    emissivity = emissivity_lu;

    TS0 = 650.; TG0 = 1173.;
    rhoS = 1500;
    eps0 = 0.4;

    //dummy properties
    rho1 = 1., rho2 = 1.;
    mu1 = 1., mu2 = 1.;

    zeta_policy = ZETA_CONST;
    shift_prod = true;

    DT = 1.e-2;
    // kinfolder = "biomass/full-solid-gas";
    kinfolder = "biomass/dummy-solid-gas";

    L0 = 8*D0;
    init_grid(1 << maxlevel);

    run();
}

// Benzene, naphthalene, acenanthrene, pyrene profiles
FILE* fyC6H6, * fyC10H8, * fyC12H8, * fyC16H10, * fyBIN;

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init (i = 0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;

  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 0.4082;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")] = 0.1991;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")] = 0.0189;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")] = 0.1536;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")] = 0.1633;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")] = 0.0184;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")] = 0.0385;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  foreach()
    porosity[] = eps0*f[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[] - porosity[])*rhoS*dv();

  TG[right] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[right] = dirichlet (1.);
      YG[top] = dirichlet (1.);
    } 
  }

  divq_rad = opensmoke_optically_thin;

  fyC6H6 = fopen("yC6H6.dat", "w");
  fyC10H8 = fopen("yC10H8.dat", "w");
  fyC12H8 = fopen("yC12H8.dat", "w");
  fyC16H10 = fopen("yC16H10.dat", "w");
  fyBIN = fopen("yBIN.dat", "w");
}

// Prints radial profile of a scalar at y=0 location
void print_radial_profile (scalar f, FILE* fp, double time, int n_samples = 100, const double length = L0) {
  double step = length/n_samples;
  for (double xx = 0.; xx <= length; xx += step) {
    double val = interpolate (f, xx, 0.);
    fprintf (fp, "%g %g %g\n", time, xx, val);
  }
}

event output (t += 0.01) {
  // scalar yC6H6 = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H6")];
  scalar yC6H6 = YGList_G[OpenSMOKE_IndexOfSpecies ("CH4")];
  // scalar yC10H8 = YGList_G[OpenSMOKE_IndexOfSpecies ("C10H8")];
  // scalar yC12H8 = YGList_G[OpenSMOKE_IndexOfSpecies ("C12H8")];
  // scalar yC16H10 = YGList_G[OpenSMOKE_IndexOfSpecies ("C16H10")];
  // scalar yBINA = YGList_G[OpenSMOKE_IndexOfSpecies ("BIN1A")];
  // scalar yBINB = YGList_G[OpenSMOKE_IndexOfSpecies ("BIN1B")];

  // scalar yBIN[];
  // foreach() {
  //   yBIN[] = yBINA[] + yBINB[];
  // }

  print_radial_profile (yC6H6, fyC6H6, t);
  // print_radial_profile (yC10H8, fyC10H8, t);
  // print_radial_profile (yC12H8, fyC12H8, t);
  // print_radial_profile (yC16H10, fyC16H10, t);
  // print_radial_profile (yBIN, fyBIN, t);

  fflush(fyC6H6);
  // fflush(fyC10H8);
  // fflush(fyC12H8);
  // fflush(fyC16H10);
  // fflush(fyBIN);

  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[] - porosity[])*rhoS*dv();

  char name[20];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE* fp = fopen (name, "w");
  fprintf (fp, "%g %g\n", t, solid_mass/solid_mass0);
  fflush(fp);
}

event adapt (i++) {
  // scalar yLVG = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  // scalar yC6H6 = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H6")];
  scalar yLVG = YGList_G[OpenSMOKE_IndexOfSpecies ("TAR")];
  scalar yC6H6 = YGList_G[OpenSMOKE_IndexOfSpecies ("CH4")];

  adapt_wavelet_leave_interface ({T, u.x, u.y, yLVG, yC6H6, porosity}, {f},
    (double[]){1.e-1, 1.e-0, 1.e-0, 1e-1, 1e-1, 1e-0}, maxlevel, minlevel, 2);
}

event stop (t = tend) {
  fclose(fyC6H6);
  // fclose(fyC10H8);
  // fclose(fyC12H8);
  // fclose(fyC16H10);
  // fclose(fyBIN);
}

event movie (t += 0.01) {
  view (tx=-0.5, ty=-0.5, width=1080, height=1080);
  clear();
  draw_vof ("f", lw = 2);
  cells();
  save ("movie.mp4");
}

#if TRACE > 1
event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif