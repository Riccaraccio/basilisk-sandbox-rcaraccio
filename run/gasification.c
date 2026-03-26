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
#include "flame.h"

const double Uin = 0.32; //inlet velocity
u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = neumann (0.);

const double tend = 60 ; //simulation time 60 s
int maxlevel = 9; int minlevel = 2;

double D0 = 8e-3, H0 = 4e-3;
double solid_mass0 = 0.;

int main() {

  lambdaSmodel = L_GALGANO;
  TS0 = 300.; TG0 = 1587.;
  rhoS = 800; // not specified in the paper, from Swedish softwood density 
  eps0 = 0.4; // guessed

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_CONST;

  DT = 1;

  G.x = -9.81;

  // kinfolder = "biomass/dummy-solid-gas";
  kinfolder = "biomass/Red-gas-2507";
  shift_prod = true;

  L0 = 15*D0;
  origin (-L0/2, 0);
  init_grid(1 << maxlevel);

#ifndef MOLAR_DIFFUSION
  fprintf (stderr, "MOLAR_DIFFUSION is compulsory for this case\n");
  return 1;
#endif

  emissivity = emissivity_diblasi;

  run();
}

double r0;

event init (i= 0) {
  fraction (f, superquadric (x, y, 20, 0.5*H0, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.0272;
  gas_start[OpenSMOKE_IndexOfSpecies ("CO2")] = 0.8371;
  gas_start[OpenSMOKE_IndexOfSpecies ("H2O")] = 0.1357;

  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 0.935; // 93.5% biomass
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.061; // 6.1% moisture
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.004; // 0.4% ash

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.5690;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("GMSW")]  = 0.2493;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.0722;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0022;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0000;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0166;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0025;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0017;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.0865;

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
    if (jj == OpenSMOKE_IndexOfSpecies ("CO2")) {
      YG[left] = dirichlet (0.8371);
      YG[left] = dirichlet (0.8371);
    } else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      YG[left] = dirichlet (0.0272);
      YG[left] = dirichlet (0.0272);
    } else if (jj == OpenSMOKE_IndexOfSpecies ("H2O")) {
      YG[left] = dirichlet (0.1357);
      YG[left] = dirichlet (0.1357);
    } else {
      YG[left] = dirichlet (0.);
      YG[left] = dirichlet (0.);
    }
  }

  divq_rad = opensmoke_optically_thin;
}

event output (t += 0.1) {
  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");
  fprintf(fp, "%g %g %g\n", t, solid_mass / solid_mass0, statsf(T).max);
  fflush(fp);
}

#if TREE
event adapt (i++) {
  scalar fuel = YGList_G[OpenSMOKE_IndexOfSpecies ("C6H10O5")];
  scalar oxidiser = YGList_G[OpenSMOKE_IndexOfSpecies ("O2")];
  // scalar fuel = YGList_G[OpenSMOKE_IndexOfSpecies ("TAR")];

  adapt_wavelet_leave_interface ({T, u.x, u.y, fuel, oxidiser, porosity}, {f},
    (double[]){1.e-1, 1.e-1, 1.e-1, 1e-1, 1e-1}, maxlevel, minlevel, 2);

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
  isoline ("zmix - zsto", lw = 1.5, lc = {1., 1., 1.});
  draw_vof ("f", lw = 1.5);
  mirror ({0, 1}) {
    squares ("XH2O", min = 0, max = 1., spread = -1, linear = true);
    isoline ("zmix - zsto", lw = 1.5, lc = {1., 1., 1.});
    draw_vof ("f", lw = 1.5);
  }
  save ("movie.mp4");
}


event stop (t = tend);

/** 
~~~gnuplot Mass
set terminal svg size 450, 450
set output "mass-fatehi.svg"
set xlabel "Time [s]"
set ylabel "Normalized solid mass [-]"
set xrange [0:150]
set yrange [0:1.05]

plot "cluster/10/diblasi/fatehi-combustion/OutputData-10"     u 1:2 w l lw 2 lc "black" t "Diblasi", \
     "cluster/10/const/fatehi-combustion/OutputData-10"       u 1:2 w l lw 2 lc "red" t "Constant", \
     "cluster/10/diblasi-swl/fatehi-combustion/OutputData-10" u 1:2 w l lw 2 lc "blue" t "Diblasi-SWL", \
     "cluster/10/lu/fatehi-combustion/OutputData-10"          u 1:2 w l lw 2 lc "orange" t "Lu", \
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

folder = "cluster/10/diblasi/fatehi-combustion/"

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
set ylabel "H2O Mole Fraction[-]"
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
