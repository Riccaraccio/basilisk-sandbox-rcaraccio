#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1

#include "temperature-profile.h"
#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
#include "view.h"
//#include "balances-interface-rop.h"

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

int maxlevel = 8; int minlevel = 2;
double D0 = 2*1.27e-2;
double solid_mass0 = 0.;

int main() {
  lambdaS = 0.1987; lambdaG = 0.076;
  lambdaSmodel = L_CORBETTA;
  cpS = 1600; cpG = 1167;
#ifdef TEMPERATURE_PROFILE
  TS0 = 300.; TG0 = 300.;
#else
  TS0 = 650.; TG0 = 650.;
#endif
  rhoS = 850; rhoG = 0.674;
  muG = 3.53e-5;
  eps0 = 0.4;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  // Da = (coord){1e-10, 1e-10};
  
  L0 = 2.5*D0;

  shift_prod = true;
  zeta_policy = ZETA_CONST;

  DT = 1.e-2;

  kinfolder = "biomass/Solid-only-2407";
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 0.4807;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")] = 0.2611;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")] = 0.1325;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")] = 0.0957;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")] = 0.0214;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]  = 0.0086;

  foreach()
    porosity[] = eps0*f[];

  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

#ifdef SOLVE_TEMPERATURE
#ifndef TEMPERATURE_PROFILE
  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
#endif
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

#ifdef TEMPERATURE_PROFILE
  double timeprofile[] = {0, 1, 2, 3, 4, 5, 13.8996139, 41.6988417, 88.03088803, 166.7953668, 
    254.8262548, 342.8571429, 454.0540541, 574.5173745, 694.980695, 806.1776062, 
    917.3745174, 1037.837838, 1200};
  double temperatureprofile[] = {300, 309.492891, 366.8388626, 424.1848341, 486.7440758, 559.7298578, 
    611.8625592, 656.1753555, 697.8815166, 723.9478673, 736.9810427, 752.6208531, 
    750.014218, 752.6208531, 752.6208531, 750.014218, 750.014218, 750.014218, 750.014218};
  
  TemperatureProfile_Set(timeprofile, temperatureprofile, sizeof(timeprofile)/sizeof(double));
#endif
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

  fprintf (fp, "%g %g %g %g %g %g\n", t, solid_mass/solid_mass0, Tcore, Tr2, Tsurf, radius/D0/2);
  fflush(fp);
}

event adapt (i++) {
  scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];

  adapt_wavelet_leave_interface ({T, u.x, u.y, inert}, {f},
    (double[]){1.e-2, 1.e-2, 1.e-2, 1e-2}, maxlevel, minlevel, 1);
}

// event movie(t+=5) {
//  clear();
//  box();
//  view (ty=-0.5, width = 1400.);
//  draw_vof("f", lw=2);
//  squares ("T", min=300, max=800, linear=true);
//  mirror ({1.,0.}) {
//    draw_vof ("f", lw=2);
//    squares ("C6H10O5_G+C6H10O5_S", min=0., max=0.3, linear=true);
// }
// save ("movie.mp4");
// }

#if DUMP
int count = 0;
event snapshots (t += 1) {
  // we keep overwriting the last two snapshots
  if (count == 0) {
    dump ("snapshot-0");
    count++;
  } else {
    dump ("snapshot-1");
    count = 0;
  }
}
#endif

// event outputfields (t = {100, 200, 305, 600}) {
//  char name[80];
//  sprintf (name, "OutputFields-%d", (int)(t));
//  static FILE * fs = fopen (name, "w");

//  output_field ({T, f, u.x, u.y, omega, zeta, porosity}, fs);
//  fflush (fs);
// }


event stop (i = 1);

/** 
~~~gnuplot temperature profiles
reset
#set terminal svg size 450,400
#set title "Zeta reaction rate"
set terminal epslatex size 3.6, 3.6 color colortext
set output "corbetta-temperature-profile.tex"
datafile = "OutputData-8-rate"

set xlabel "Time [s]"
set size square
set ylabel "Temperature [K]"
set key bottom right box width 1
set xrange [0:1000]
set yrange [300:850]
set grid
plot  "data/corbetta-core.txt"  u 1:2 w p pt 64 ps 2 lw 3 lc "dark-green" notitle, \
      "data/corbetta-r2.txt"    u 1:2 w p pt 64 ps 2 lw 3 lc "blue" notitle, \
      "data/corbetta-surf.txt"  u 1:2 w p pt 64 ps 2 lw 3 lc "black" notitle, \
      datafile u 1:3 w l lw 6 lc "dark-green" t "core", \
      datafile u 1:4 w l lw 6 lc "blue" t "half-radius", \
      datafile u 1:5 w l lw 6 lc "black" t "surface"
~~~

~~~gnuplot species profiles
reset
set terminal svg size 1350,400
set multiplot layout 1,3

# Define a function for centered derivative with multiplier
# centered_diff(file, column, multiplier) = sprintf("< awk -v col=%d -v mult=%g 'NR==2{a=$1;b=$col;next} NR==3{c=$1;d=$col;next} NR>=4{dt=$1-a; dy=$col-b; if(dt!=0) print c,(dy/dt)*mult; a=c; b=d; c=$1; d=$col}' %s", column, multiplier, file)

# Function with sampling parameter (most readable)
centered_diff(file, column, multiplier, interval) = sprintf("< awk -v col=%d -v mult=%g -v skip=%d 'NR%%skip==1{count++; if(count==2){a=$1;b=$col;next} if(count==3){c=$1;d=$col;next} if(count>=4){dt=$1-a; dy=$col-b; if(dt!=0) print c,(dy/dt)*mult; a=c; b=d; c=$1; d=$col}}' %s", column, multiplier, interval, file)

datafile = "balances-7-const"

#set title "Zeta Rate"
set xlabel "Time [s]"
set ylabel "Production rate [mg/s/g_{init}]"
set key top right box width 2.2
set grid
set xrange [0:600]
set yrange [0:0.5]
plot  centered_diff(datafile, 32, 1000, 25) u 1:2 w l lw 2 lc "dark-green" t "sym CO_{2}", \
      "data/expCO2" u 1:2 w p pt 4 lc "dark-green" t "exp CO_{2}", \
      centered_diff(datafile, 31, 1000, 25) u 1:2 w l lw 2 lc "blue" t "sym CO", \
      "data/expCO" u 1:2 w p pt 4 lc "blue" t "exp CO", \
      centered_diff(datafile, 19, 1000, 25) u 1:2 w l lw 2 lc "black" t "sym HCHO", \
      "data/expHCHO" u 1:2 w p pt 4 lc "black" t "exp HCHO"

set xlabel "Time [s]"
#set ylabel "Production rate [mg/s/g_{init}]"
unset ylabel
set key top right box width 2.2
set grid
set xrange [0:600]
set yrange [0:0.4]
plot  centered_diff(datafile, 7, 1000, 25) u 1:2 w l lw 2 lc "dark-green" t "sym CH_{3}COOH", \
      "data/expCH3COOH" u 1:2 w p pt 4 lc "dark-green" t "exp CH_{3}COOH", \
      centered_diff(datafile, 18, 1000, 25) u 1:2 w l lw 2 lc "blue" t "sym CH_{3}OH", \
      "data/expCH3OH" u 1:2 w p pt 4 lc "blue" t "exp CH_{3}OH"

set xlabel "Time [s]"
#set ylabel "Production rate [mg/s/g_{init}]"

set key top right box width 2.2
set grid
set xrange [0:600]
set yrange [0:0.05]
plot  centered_diff(datafile, 33, 1000*10, 25) u 1:2 w l lw 2 lc "dark-green" t "sym H_{2}*10", \
      "data/expH2" u 1:($2*10) w p pt 4 lc "dark-green" t "exp H_{2}*10", \
      centered_diff(datafile, 34, 1000, 25) u 1:2 w l lw 2 lc "blue" t "sym CH_{4}", \
      "data/expCH4" u 1:2 w p pt 4 lc "blue" t "exp CH_{4}", \
      centered_diff(datafile, 20, 1000, 25) u 1:2 w l lw 2 lc "black" t "sym HCOOH", \
      "data/expHCOOH" u 1:2 w p pt 4 lc "black" t "exp HCOOH"

unset multiplot
~~~

~~~gnuplot species profiles
reset
set terminal epslatex size 10, 3.6 color colortex
set output "corbetta-species-profile.tex"
set multiplot layout 1,3

# Function with sampling parameter (most readable)
centered_diff(file, column, multiplier, interval) = sprintf("< awk -v col=%d -v mult=%g -v skip=%d 'NR%%skip==1{count++; if(count==2){a=$1;b=$col;next} if(count==3){c=$1;d=$col;next} if(count>=4){dt=$1-a; dy=$col-b; if(dt!=0) print c,(dy/dt)*mult; a=c; b=d; c=$1; d=$col}}' %s", column, multiplier, interval, file)

datafile = "balances-7-rate"

#set title "Zeta Rate"
set xlabel "Time [s]"
set ylabel "Production rate [mg/s/g_{init}]"
set key top right box width 2.2
set size square 
set grid
set xrange [0:600]
set xtics 200
set ytics 0.1
set yrange [0:0.5]

plot  "data/expCO2"   u 1:2 w p pt 64 ps 2 lw 3 lc "dark-green" notitle, \
      "data/expCO"    u 1:2 w p pt 64 ps 2 lw 3 lc "blue" notitle, \
      "data/expHCHO"  u 1:2 w p pt 64 ps 2 lw 3 lc "black" notitle, \
      centered_diff(datafile, 32, 1000, 25) u 1:2 w l lw 6 lc "dark-green" t "CO_{2}", \
      centered_diff(datafile, 31, 1000, 25) u 1:2 w l lw 6 lc "blue" t "CO", \
      centered_diff(datafile, 19, 1000, 25) u 1:2 w l lw 6 lc "black" t "HCHO"

set xlabel "Time [s]"
#set ylabel "Production rate [mg/s/g_{init}]"
unset ylabel
set key top right box width 2.2
set size square
set grid
set xrange [0:600]
set xtics 200
set ytics 0.1
set yrange [0:0.4]
      
plot  "data/expCH3COOH" u 1:2 w p pt 64 ps 2 lw 3 lc "dark-green" notitle, \
      "data/expCH3OH"   u 1:2 w p pt 64 ps 2 lw 3 lc "blue" notitle, \
      centered_diff(datafile, 7, 1000, 25)  u 1:2 w l lw 6 lc "dark-green" t "CH_{3}COOH", \
      centered_diff(datafile, 18, 1000, 25) u 1:2 w l lw 6 lc "blue" t "CH_{3}OH"

set xlabel "Time [s]"
#set ylabel "Production rate [$mg/s/g_{init}$]"

set key top right box width 2.2
set grid
set xrange [0:600]
set size square
set xtics 200
set ytics 0.01
set yrange [0:0.05]

plot  "data/expH2"  u 1:($2*10) w p pt 64 ps 2 lw 3 lc "dark-green" notitle, \
      "data/expCH4" u 1:2       w p pt 64 ps 2 lw 3 lc "blue" notitle, \
      "data/expHCOOH" u 1:2     w p pt 64 ps 2 lw 3 lc "black" notitle, \
      centered_diff(datafile, 33, 1000*10, 25) u 1:2 w l lw 6 lc "dark-green" t "H_{2}*10", \
      centered_diff(datafile, 34, 1000, 25)    u 1:2 w l lw 6 lc "blue" t "CH_{4}", \
      centered_diff(datafile, 20, 1000, 25)    u 1:2 w l lw 6 lc "black" t "HCOOH"

unset multiplot
~~~

~~~gnuplot species profiles using rop files
reset
set terminal svg size 1350,400
set multiplot layout 1,3

datafile = "balances-8-const-rop"

set title "Zeta Rate"
set xlabel "Time [s]"
set ylabel "Production rate [mg/s/g_{init}]"
set key top right box width 2.2
set grid
set xrange [0:600]
set yrange [0:0.5]
plot  datafile u 1:($32*1000) every 200 w l lw 2 lc "dark-green" t "sym CO_{2}", \
      "data/expCO2" u 1:2 w p pt 4 lc "dark-green" t "exp CO_{2}", \
      datafile u 1:($31*1000) every 200 w l lw 2 lc "blue" t "sym CO", \
      "data/expCO" u 1:2 w p pt 4 lc "blue" t "exp CO", \
      datafile u 1:($19*1000) every 200 w l lw 2 lc "black" t "sym HCHO", \
      "data/expHCHO" u 1:2 w p pt 4 lc "black" t "exp HCHO"

set xlabel "Time [s]"
set ylabel "Production rate [mg/s/g_{init}]"
set key top right box width 2.2
set grid
set xrange [0:600]
set yrange [0:0.4]
plot  datafile u 1:($7*1000) every 200 w l lw 2 lc "dark-green" t "sym CH_{3}COOH", \
      "data/expCH3COOH" u 1:2 w p pt 4 lc "dark-green" t "exp CH_{3}COOH", \
      datafile u 1:($18*1000) every 200 w l lw 2 lc "blue" t "sym CH_{3}OH", \
      "data/expCH3OH" u 1:2 w p pt 4 lc "blue" t "exp CH_{3}OH"

set xlabel "Time [s]"
set ylabel "Production rate [mg/s/g_{init}]"
set key top right box width 2.2
set grid
set xrange [0:600]
set yrange [0:0.05]
plot  datafile u 1:($33*10*1000) every 200 w l lw 2 lc "dark-green" t "sym H_{2}*10", \
      "data/expH2" u 1:($2*10) w p pt 4 lc "dark-green" t "exp H_{2}*10", \
      datafile u 1:($34*1000) every 200 w l lw 2 lc "blue" t "sym CH_{4}", \
      "data/expCH4" u 1:2 w p pt 4 lc "blue" t "exp CH_{4}", \
      datafile u 1:($20*1000) every 200 w l lw 2 lc "black" t "sym HCOOH", \
      "data/expHCOOH" u 1:2 w p pt 4 lc "black" t "exp HCOOH"

unset multiplot
~~~

~~~gnuplot testing plot
reset
set terminal svg size 450,400

set xlabel "Time [s]"
set grid
# Function with sampling parameter
centered_diff(file, column, multiplier, interval) = sprintf("< awk -v col=%d -v mult=%g -v skip=%d 'NR%%skip==1{count++; if(count==2){a=$1;b=$col;next} if(count==3){c=$1;d=$col;next} if(count>=4){dt=$1-a; dy=$col-b; if(dt!=0) print c,(dy/dt)*mult; a=c; b=d; c=$1; d=$col}}' %s", column, multiplier, interval, file)
 
datafile = "balances-7-old"

plot centered_diff(datafile, 6, 1, 2) u 1:2 w l lw 2 lc "black" t "Tar", \
     centered_diff(datafile, 7, 1, 2) u 1:2 w l lw 2 lc "dark-green" t "H2O", \
     "balances-7" u 1:6 w l dt 2 lw 2 lc "black" t "Tar", \
     "balances-7" u 1:7 w l dt 2 lw 2 lc "dark-green" t "H2O"
~~~
**/