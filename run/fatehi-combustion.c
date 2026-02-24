#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
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
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = neumann (0.);

const double tend = 150 ; //simulation time 150 s
int maxlevel = 10; int minlevel = 2;

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

  zeta_policy = ZETA_SHRINK;

  DT = 1;

  G.x = -9.81;

  // kinfolder = "biomass/dummy-solid-gas";
  kinfolder = "biomass/Red-gas-2507";
  shift_prod = true;

  L0 = 20*D0;
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
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 0.935; // 93.5% biomass
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.061; // 6.1% moisture
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.004; // 0.4% ash

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.4344;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("GMSW")]  = 0.2108;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.1347;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0786;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0178;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0167;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0419;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0041;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.0610;

  foreach()
    porosity[] = eps0*f[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f
  
  fprintf(stderr, "Initial solid mass: %g kg\n", solid_mass0);

  TG[left] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);
  TG[right] = neumann (0.);
  TG[bottom] = neumann (0.);

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

  fTprofile_2mm = fopen("T_profile_2mm.dat", "w");
  fTprofile_11mm = fopen("T_profile_11mm.dat", "w");
  fpxH2Oprofile_2mm = fopen("xH2O_profile_2mm.dat", "w");
  fpxH2Oprofile_11mm = fopen("xH2O_profile_11mm.dat", "w");
}

// Prints profile of a scalar at a given x location for all y form 0 to L0/2
void print_profile (scalar f, double x_interp, FILE* fp, int time, int n_samples = 100, const double length = L0/2) {
  double step = length/n_samples;
  for (double yy = 0.; yy <= length; yy += step) {
    double val = interpolate (f, x_interp, yy);
    fprintf (fp, "%d %g %g\n", time, yy, val);
  }
}

// Calculates the H2O-density path-averaged temperarure
double T_H2O_weigthed_average (double x_interp, int n_samples = 100, const double length = L0/2) {
  scalar XH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];

  double step = length/n_samples;
  double numerator = 0., denominator = 0.;
  for (double yy = 0.; yy <= length; yy += step) {
    double xH2O_local = interpolate (XH2O, x_interp, yy);
    numerator += xH2O_local * step;
    denominator += xH2O_local / interpolate (T, x_interp, yy) * step; 
  }

  if (denominator < SEPS) // avoid division by 0
    return TG0;

  return numerator/denominator;
}

// Calculates the path-averaged xH2O mole fractions
// NOTE: if too high, can be weighted by 1/T
double path_average_XH2O (double x_interp, const int n_samples = 100, const double length = L0/2) {
  scalar XH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];

  double step = length/n_samples, numerator = 0.;
  for (double yy = 0.; yy <= length; yy += step)
    numerator += interpolate (XH2O, x_interp, yy) * step;

  return numerator/(length);
}

event print_profile (t += 1; t <= 100) {
  const int time = (int) round(t);

  // Temperature profiles
  print_profile (T, H0/2 + 2e-3, fTprofile_2mm, time, 200, L0/2);
  print_profile (T, H0/2 + 11e-3, fTprofile_11mm, time, 200, L0/2);

  // Water vapor mole fraction profile
  scalar XH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
  print_profile (XH2O, H0/2 + 2e-3, fpxH2Oprofile_2mm, time, 200, L0/2);
  print_profile (XH2O, H0/2 + 11e-3, fpxH2Oprofile_11mm, time, 200, L0/2);

  fflush (fTprofile_2mm);
  fflush (fTprofile_11mm);
  fflush (fpxH2Oprofile_2mm);
  fflush (fpxH2Oprofile_11mm);
}

event output (t += 1) {
  // Interpolate Water vapor mass fraction
  double xH2O[5], sample_points[5] = {H0/2 + 2e-3, H0/2 + 4e-3, H0/2 + 8e-3, H0/2 + 11e-3, H0/2 + 15e-3};

  const double L_flame_exp = 30e-3/2;
  for (int ii = 0; ii < 5; ii++)
    xH2O[ii] = path_average_XH2O (sample_points[ii], 200, L_flame_exp);

  // Interpolate temperature
  double Tavg[5];
  for (int ii = 0; ii < 5; ii++)
    if (xH2O[ii] < 1e-4)
      Tavg[ii] = TG0; //ambient temperature
    else
      Tavg[ii] = T_H2O_weigthed_average (sample_points[ii], 200, L_flame_exp);

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();

  static FILE *fpxH2O = fopen("xH2OProfile.dat", "w");
  fprintf(fpxH2O, "%g %g %g %g %g %g\n", t, xH2O[0], xH2O[1], xH2O[2], xH2O[3], xH2O[4]);
  fflush(fpxH2O);
  
  static FILE * fpT = fopen ("TemperatureProfile.dat", "w");
  fprintf (fpT, "%g %g %g %g %g %g\n", t, Tavg[0], Tavg[1], Tavg[2], Tavg[3], Tavg[4]);
  fflush(fpT);
  
  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");
  fprintf(fp, "%g %g %g\n", t, solid_mass / solid_mass0, statsf(T).max);
  fflush(fp);
}

#if TREE
event adapt (i++) {
  scalar fuel = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  // scalar fuel = YGList_G[OpenSMOKE_IndexOfSpecies ("TAR")];

  adapt_wavelet_leave_interface ({T, u.x, u.y, fuel, porosity}, {f},
    (double[]){1.e0, 1.e-1, 1.e-1, 1e-1}, maxlevel, minlevel, 2);

  // Unrefine for outflow condition
  unrefine (x > L0*0.4);
}
#endif

event movie (t += 1) {
  scalar XH2O_G = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar XH2O_S = XGList_S[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar XH2O[];
  foreach()
    XH2O[] = XH2O_S[]*f[] + XH2O_G[]*(1. - f[]);

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

event save_fields (t += 10) {
  char name[80];
  sprintf (name, "fields-%g", t);

  // H2O
  scalar XH2O_G = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar XH2O_S = XGList_S[OpenSMOKE_IndexOfSpecies ("H2O")];
  scalar XH2O[];
  foreach()
    XH2O[] = XH2O_S[]*f[] + XH2O_G[]*(1. - f[]);

  // CO2
  scalar XCO2_G = XGList_G[OpenSMOKE_IndexOfSpecies ("CO2")];
  scalar XCO2_S = XGList_S[OpenSMOKE_IndexOfSpecies ("CO2")];
  scalar XCO2[];
  foreach()
    XCO2[] = XCO2_S[]*f[] + XCO2_G[]*(1. - f[]);
  
  // OH
  scalar XOH_G = XGList_G[OpenSMOKE_IndexOfSpecies ("OH")];
  scalar XOH_S = XGList_S[OpenSMOKE_IndexOfSpecies ("OH")];
  scalar XOH[];
  foreach()
    XOH[] = XOH_S[]*f[] + XOH_G[]*(1. - f[]);

  // CO
  scalar XCO_G = XGList_G[OpenSMOKE_IndexOfSpecies ("CO")];
  scalar XCO_S = YGList_S[OpenSMOKE_IndexOfSpecies ("CO")];
  scalar XCO[];
  foreach()    
    XCO[] = XCO_S[]*f[] + XCO_G[]*(1. - f[]);

  // LVG
  scalar LVG_G = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  scalar LVG_S = YGList_S[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  scalar LVG[];
  foreach()
    LVG[] = LVG_S[]*f[] + LVG_G[]*(1. - f[]);

  output_field ({T, f, u.x, u.y, XH2O, XCO2, XOH, XCO, LVG, omega}, name);
}

event stop (t = tend) {
  fclose(fTprofile_2mm);
  fclose(fTprofile_11mm);
  fclose(fpxH2Oprofile_2mm);
  fclose(fpxH2Oprofile_11mm);
}

/** 
~~~gnuplot Mass
set terminal svg size 550, 550
set output "mass-fatehi.svg"
set xlabel "Time [s]"
set ylabel "Normalized solid mass [-]"
set xrange [0:150]
set yrange [0:1.05]

plot "cluster/full-9/fatehi-combustion/OutputData-9" u 1:2 w l lw 4 lc "black" notitle, \
     "../../data/fatehi/mass" u 1:2 w p pt 64 ps 1 lw 3 lc "black" notitle
~~~

~~~gnuplot H2O mole fraction
reset
set terminal svg size 550, 350
set output "xH2O-fatehi.svg"
set xlabel "Time [s]"
set ylabel "H2O Mole Fraction [-]"
set yrange [0.:0.5]
set xrange [0:70]

folder = "cluster/full-9/fatehi-combustion/"

error_margin = 0.03
shift = 5

plot  "../../data/fatehi/yH2O-11mm" u 1:2:(error_margin) w yerrorbars pt 64 ps 0.8 lc "gray" notitle  ,\
      "../../data/fatehi/yH2O-2mm" u 1:2:(error_margin) w yerrorbars pt 64 ps 0.8 lc "light-coral" notitle  ,\
      sprintf("%s%s", folder, "xH2OProfile.dat") u ($1-shift):4 w l lw 2 lc "gray" title "11 mm",\
      sprintf("%s%s", folder, "xH2OProfile.dat") u ($1-shift):2 w l lw 2 lc "light-coral" title "2 mm", \
      sprintf("%s%s", folder, "results_2mm/effective_values.dat") u ($1-shift):($3/100) w l lw 2 lc "red" title "2 mm effective", \
      sprintf("%s%s", folder, "results_11mm/effective_values.dat") u ($1-shift):($3/100) w l lw 2 lc "black" title "11 mm effective" 
~~~

~~~gnuplot Temperature
reset
set terminal svg size 550, 350
set output "temperature-fatehi.svg"
set xlabel "Time [s]"
set ylabel "Temperature [K]"
set yrange [800:1900]
set xrange [0:60]

folder = "cluster/full-9/fatehi-combustion/"

error_margin = 50
shift = 5

plot  "../../data/fatehi/T-11mm" u 1:2:(error_margin) w yerrorbars pt 4 ps 0.8 lw 1 lc "black" notitle, \
      "../../data/fatehi/T-2mm" u 1:2:(error_margin) w yerrorbars pt 4 ps 0.8 lw 1 lc "light-coral" notitle, \
      sprintf("%s%s", folder, "TemperatureProfile.dat") u ($1-shift):5 w l lw 2 lc "gray" title "11 mm", \
      sprintf("%s%s", folder, "TemperatureProfile.dat") u ($1-shift):2 w l lw 2 lc "light-coral" title "2 mm", \
      sprintf("%s%s", folder, "results_2mm/effective_values.dat") u ($1-shift):2 w l lw 2 lc "red" title "2 mm effective", \
      sprintf("%s%s", folder, "results_11mm/effective_values.dat") u ($1-shift):2 w l lw 2 lc "black" title "11 mm effective"
~~~

~~~gnuplot Temperature evolution
reset
set terminal svg size 550, 550 
set output "temperature-evolution-fatehi-2mm.svg"
load "/root/gnuplot-palettes/inferno.pal"

set xlabel "y position"
set ylabel "Temperature"
set key right top
set ytics 900,200,2300
set yrange [800:2300]
set xrange [-0.06:0.06]
set title "Temperature 2mm"

end_time = 60 
step_time = 10

set cbrange [0:end_time]
unset colorbox

folder = "cluster/full-9/fatehi-combustion/"

plot for [t=0:end_time:step_time] sprintf("%s%s", folder, "T_profile_2mm.dat") u ($1==t ? $2 : 1/0):3:(t) w l lw 3 lc palette title sprintf("%d s", t), \
     for [t=0:end_time:step_time] sprintf("%s%s", folder, "T_profile_2mm.dat") u ($1==t ? -$2 : 1/0):3:(t) w l lw 3 lc palette notitle
~~~

~~~gnuplot Temperature evolution
reset
set terminal svg size 550, 550
set output "temperature-evolution-fatehi-11mm.svg"
load "/root/gnuplot-palettes/inferno.pal"

set xlabel "Axial distance [m]"
set ylabel "Temperature [K]"
set key right top
set ytics 900,200,2300
set yrange [950:2500]
set xrange [-0.06:0.06]
set title "Temperature 11mm"

folder = "cluster/full-9/fatehi-combustion/"

end_time = 60
step_time = 10

set cbrange [0:end_time]
unset colorbox

plot for [t=0:end_time:step_time] sprintf("%s%s", folder, "T_profile_11mm.dat") u ($1==t ? $2 : 1/0):3:(t) w l lw 3 lc palette title sprintf("%d s", t), \
     for [t=0:end_time:step_time] sprintf("%s%s", folder, "T_profile_11mm.dat") u ($1==t ? -$2 : 1/0):3:(t) w l lw 3 lc palette notitle
~~~

~~~gnuplot xH2O evolution
reset
set terminal svg size 550, 550
set output "xH2O-evolution-fatehi-2mm.svg"
load "/root/gnuplot-palettes/inferno.pal"
set xlabel "y position"
set ylabel "H2O Mole Fraction"
set key right top
set yrange [0:0.5]
set xrange [-0.06:0.06]
set title "xH2O 2mm"

end_time = 60
step_time = 10

set cbrange [0:end_time]
unset colorbox 

plot for [t=0:end_time:step_time] "cluster/full-9/fatehi-combustion/xH2O_profile_2mm.dat" u ($1==t ? $2 : 1/0):3:(t) w l lw 2 lc palette title sprintf("%d s", t), \
     for [t=0:end_time:step_time] "cluster/full-9/fatehi-combustion/xH2O_profile_2mm.dat" u ($1==t ? -$2 : 1/0):3:(t) w l lw 2 lc palette notitle

~~~

~~~gnuplot xH2O evolution
reset
set terminal svg size 550, 550
set output "xH2O-evolution-fatehi-11mm.svg"
load "/root/gnuplot-palettes/inferno.pal"
set xlabel "y position"
set ylabel "H2O Mole Fraction"
set key right top
set yrange [0:0.5]
set xrange [-0.06:0.06]
set title "xH2O 11mm"

end_time = 60
step_time = 10

set cbrange [0:end_time]
unset colorbox 

folder = "cluster/full-9/fatehi-combustion/"

plot for [t=0:end_time:step_time] sprintf("%s%s", folder, "xH2O_profile_11mm.dat") u ($1==t ? $2 : 1/0):3:(t) w l lw 2 lc palette title sprintf("%d s", t), \
     for [t=0:end_time:step_time] sprintf("%s%s", folder, "xH2O_profile_11mm.dat") u ($1==t ? -$2 : 1/0):3:(t) w l lw 2 lc palette notitle

~~~

~~~gnuplot Test plot
reset
set terminal svg size 550, 350
set output "test-plot-fatehi.svg"
set xlabel "Time [s]"
set ylabel "H2O Mole  Fraction[-]"
set yrange [0.:0.5]
set xrange [0:60]

endtime = 70
step = 1
n = endtime/step + 1
error_margin = 0.03
shift = 5

folder_path = "cluster/full-9/fatehi-combustion/"

# Define the upper space limit
space_max = 0.015

# ---- 2mm case ----
system(sprintf("awk '$2 <= %f { sum[$1]+=$3; count[$1]++; if(!($1 in max) || $3>max[$1]) max[$1]=$3 } END { for(t in sum) print t, sum[t]/count[t], max[t] }' %sxH2O_profile_2mm.dat | sort -n > %savg_h2o_2mm.dat", space_max, folder_path, folder_path))
system(sprintf("awk '$2 <= %f { sum[$1]+=$3; count[$1]++; if(!($1 in max) || $3>max[$1]) max[$1]=$3 } END { for(t in sum) print t, sum[t]/count[t], max[t] }' %sxH2O_profile_11mm.dat | sort -n > %savg_h2o_11mm.dat", space_max, folder_path, folder_path))


plot  "../../data/fatehi/yH2O-2mm" u 1:2:(error_margin) w yerrorbars pt 64 ps 0.8 lc "red" notitle, \
      "../../data/fatehi/yH2O-11mm" u 1:2:(error_margin) w yerrorbars pt 64 ps 0.8 lc "black" notitle, \
      sprintf("%savg_h2o_11mm.dat", folder_path) u ($1-shift):3 w l lw 2 lc "black" t "max 11mm", \
      sprintf("%savg_h2o_2mm.dat", folder_path) u ($1-shift):3 w l lw 2 lc "red" t "max 2mm", \
      sprintf("%savg_h2o_11mm.dat", folder_path) u ($1-shift):2 w l lw 2 lc "gray" t "11mm", \
      sprintf("%savg_h2o_2mm.dat", folder_path) u ($1-shift):2 w l lw 2 lc "light-coral" t "2mm"




~~~

**/
