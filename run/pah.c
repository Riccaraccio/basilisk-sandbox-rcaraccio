#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1
#define TURN_OFF_HEAT_OF_REACTION 1

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
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

const double tend = 10 ; //simulation time 60 s
int maxlevel = 8; int minlevel = 2;
double D0 = 2e-2;
double solid_mass0 = 0.;

int main () {
    lambdaSmodel = L_LU;
    emissivity = emissivity_lu;

    TS0 = 550.; TG0 = 1173.;
    rhoS = 1500;
    eps0 = 0.4;

    //dummy properties
    rho1 = 1., rho2 = 1.;
    mu1 = 1., mu2 = 1.;

    zeta_policy = ZETA_CONST;
    shift_prod = true;

    DT = 1.e-2;
    kinfolder = "biomass/full-solid-gas";
    L0 = 10*D0;
    init_grid(1 << maxlevel);

    run();
}

// Benzene, naphthalene, acenanthrene, pyrene profiles
FILE* fyC6H6, * fyC10H8, * fyC12H8, * fyC16H10, * fyBIN;

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init (i = 0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  foreach()
    porosity[] = eps0*f[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //

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
void print_radial_profile (scalar f, FILE* fp, int time, int n_samples = 100, const double length = L0) {
  double step = length/n_samples;
  for (double xx = 0.; xx <= length; xx += step) {
    double val = interpolate (f, xx, 0.);
    fprintf (fp, "%d %g %g\n", time, xx, val);
  }
}

event output (t += 0.1) {
  scalar yC6H6 = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H6")];
  scalar yC10H8 = YGList_G[OpenSMOKE_IndexOfSpecies ("C10H8")];
  scalar yC12H8 = YGList_G[OpenSMOKE_IndexOfSpecies ("C12H8")];
  scalar yC16H10 = YGList_G[OpenSMOKE_IndexOfSpecies ("C16H10")];
  scalar yBINA = YGList_G[OpenSMOKE_IndexOfSpecies ("BIN1A")];
  scalar yBINB = YGList_G[OpenSMOKE_IndexOfSpecies ("BIN1B")];

  scalar yBIN[];
  foreach() {
    yBIN[] = yBINA[] + yBINB[];
  }

  print_radial_profile (yC6H6, fyC6H6, (int) round(t));
  print_radial_profile (yC10H8, fyC10H8, (int) round(t));
  print_radial_profile (yC12H8, fyC12H8, (int) round(t));
  print_radial_profile (yC16H10, fyC16H10, (int) round(t));
  print_radial_profile (yBIN, fyBIN, (int) round(t));

  fflush(fyC6H6);
  fflush(fyC10H8);
  fflush(fyC12H8);
  fflush(fyC16H10);
  fflush(fyBIN);
}

#if TREE
event adapt (i++) {
  scalar yLVG = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  scalar yC6H6 = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H6")];

  adapt_wavelet_leave_interface ({T, u.x, u.y, yLVG, yC6H6, porosity}, {f},
    (double[]){1.e-1, 1.e-0, 1.e-0, 1e-1, 1e-1, 1e-0}, maxlevel, minlevel, 2);
}
#endif

event stop (t = tend) {
  fclose(fTprofile_2mm);
  fclose(fTprofile_11mm);
  fclose(fxH2Oprofile_2mm);
  fclose(fxH2Oprofile_11mm);
  fclose(fxOHprofile_2mm);
  fclose(fxOHprofile_11mm);
}

event movie (t += 0.1) {
  view (width = 800, height = 400, fov = 20);
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