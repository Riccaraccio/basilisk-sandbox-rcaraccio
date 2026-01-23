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

const double Uin = 0.03; //inlet velocity
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

#ifndef MOLAR_DIFFUSION
  fprintf (stderr, "MOLAR_DIFFUSION is compulsory for this case\n");
  return 1;
#endif

  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

double r0;
// Output files for profiles
FILE* fTprofile_2mm, * fTprofile_11mm;
FILE* fpxH2Oprofile_2mm, * fpxH2Oprofile_11mm;

event init (i= 0) {
  fraction (f, superquadric (x, y, 20, 0.5*H0, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 0.935; // 93.5% biomass
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.061; // 6.1% moisture
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.004; // 0.4% ash

  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.4344;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("GMSW")]  = 0.2108;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.1347;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0786;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0178;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0167;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0419;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0040;
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

  fTprofile_2mm = fopen ("T_profile_2mm.dat", "w");
  fTprofile_11mm = fopen ("T_profile_11mm.dat", "w");
  fpxH2Oprofile_2mm = fopen ("xH2O_profile_2mm.dat", "w");
  fpxH2Oprofile_11mm = fopen ("xH2O_profile_11mm.dat", "w");
}

// Prints profile of a scalar at a given x location for all y form 0 to L0/2
void print_profile (scalar f, double x_interp, FILE* fp, int time, int n_samples =  100) {
  coord box[2] = {{x_interp, 0.}, {x_interp, L0/2}}
  coord p, n = {1, n_samples};
  double delta = (L0/2)/n_samples;

  foreach_region (p, box, n) {
    double val = interpolate (f, p.x, p.y);
    fprintf (fp, "%d %g %g\n", time, yy, val);
  }
  fflush (fp);
  // double step = (L0/2)/(1<<maxlevel);
  // for (double yy = 0.; yy <= L0/2; yy += step) {
  //   double val = interpolate (f, x, yy);
  // }
}

// Calculates the H2O-density path-averaged temperarure
double T_H2O_weigthed_average (double x_interp, int n_samples = 100) {
  scalar XH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];

  coord box[2] = {{x_interp, 0.}, {x_interp, L0/2}};
  coord p, n = {1, n_samples};
  double delta = (L0/2)/n_samples, numerator = 0., denominator = 0.;
  foreach_region (p, box, n) {
    numerator += interpolate (XH2O, p.x, p.y) * delta;
    denominator += interpolate (XH2O, p.x, p.y) / interpolate (T, p.x, p.y) * delta;
  }

  if (denominator < SEPS) // avoid division by 0
    return TG0;

  return numerator/(denominator + SEPS);
}

// Calculates the path-averaged xH2O mole fractions
double path_average_XH2O (double x_interp, const int n_points =100) {
  scalar XH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];

  coord box[2] = {{x_interp, 0.}, {x_interp, L0/2}};
  coord p, n = {1, n_points};
  double delta = (L0/2)/n_points, numerator = 0.;
  foreach_region (p, box, n)
    numerator += interpolate (XH2O, p.x, p.y) * delta;

  return numerator/(L0/2.);
}

event print_profile (t = 0; t += 10) {
  update_mole_fields();
  int time = (int) round(t);

  // Temperature profiles
  print_profile (T, H0/2 + 2e-3, fTprofile_2mm, time);
  print_profile (T, H0/2 + 11e-3, fTprofile_11mm, time);

  // Water vapor mole fraction profile
  scalar XH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
  print_profile (XH2O, H0/2 + 2e-3, fpxH2Oprofile_2mm, time);
  print_profile (XH2O, H0/2 + 11e-3, fpxH2Oprofile_11mm, time);
}

event output (t += 1) {
  update_mole_fields();

  // Interpolate Water vapor mass fraction
  static FILE * fpxH2O = fopen ("xH2OProfile.dat", "w");
  double xH2O[5], sample_points[5] = {H0/2 + 2e-3, H0/2 + 4e-3, H0/2 + 8e-3, H0/2 + 11e-3, H0/2 + 15e-3};

  for (int i = 0; i < 5; i++)
    xH2O[i] = path_average_XH2O (sample_points[i]);

  fprintf (fpxH2O, "%g %g %g %g %g %g\n", t, xH2O[0], xH2O[1], xH2O[2], xH2O[3], xH2O[4]);
  fflush(fpxH2O);

  // Interpolate temperature
  static FILE * fpT = fopen ("TemperatureProfile.dat", "w");
  double Tavg[5];
  for (int i = 0; i < 5; i++)
    if (xH2O[i] < 1e-4)
      Tavg[i] = TG0; //ambient temperature
    else
      Tavg[i] = T_H2O_weigthed_average (sample_points[i]);

  fprintf (fpT, "%g %g %g %g %g %g\n", t, Tavg[0], Tavg[1], Tavg[2], Tavg[3], Tavg[4]);
  fflush(fpT);

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
  scalar XH2O_G = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar XH2O_S = XGList_S[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar XH2O[];
  foreach()
    XH2O[] = XH2O_S[] + XH2O_G[];

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

event stop (t = tend) {
  fclose (fTprofile_2mm);
  fclose (fTprofile_11mm);
  fclose (fpxH2Oprofile_2mm);
  fclose (fpxH2Oprofile_11mm);
}

/** 
~~~gnuplot Mass
set terminal svg size 400, 400
set output "mass-fatehi.svg"
set xlabel "Time [s]"
set ylabel "Normalized solid mass [-]"
set grid
set xrange [0:150]
set yrange [0:1.05]

plot "OutputData-9" u 1:2 w l lw 3 lc "black" notitle, \
     "../../data/fatehi/mass" u 1:2 w lp pt 64 ps 1 lw 1 lc "black" notitle
~~~

~~~gnuplot Temperature
reset
set terminal svg size 400, 400
set output "temperature-fatehi.svg"
set xlabel "Time [s]"
set ylabel "Temperature [K]"
set grid
set yrange [400:1900]
set xrange [0:80]

plot "TemperatureProfile.dat" u 1:2 w l lw 3 lc "black" t "T-2mm", \
     "../../data/fatehi/T-2mm" u 1:2 w lp pt 1 ps 1 lw 1 lc "black" notitle, \
      "TemperatureProfile.dat" u 1:4 w l lw 3 lc "red" t "T-8mm", \
      "../../data/fatehi/T-8mm" u 1:2 w lp pt 1 ps 1 lw 1 lc "red" notitle

~~~
**/
