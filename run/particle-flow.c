#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define RADIATION_INTERFACE 0.8
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1

#ifndef T_ENV
# define T_ENV 673
#endif

#include "grid/multigrid.h"
#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"

double Uin = 0.0221; //inlet velocity, 400 NL/min over 5.03e-3 m^2 area
double tend = 600.; //simulation time

u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
pf[left]     = neumann (0.);
psi[left]    = dirichlet (0.);

u.n[top]      = dirichlet (0.);
u.t[top]      = dirichlet (0.);
p[top]        = neumann (0.);
pf[top]       = neumann (0.);
psi[top]      = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);
psi[right]    = dirichlet (0.);

int maxlevel = 7; int minlevel = 2;
double D0 = 2e-2;
double solid_mass0 = 0., moisture0 = 0.;

int main() {
  
  lambdaS = 0.1987;
  lambdaSmodel = L_HUANG;
  TS0 = 300.; TG0 = T_ENV;
  rhoS = 920;
  eps0 = 0.39;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_REACTION;

#if TREE
  L0 = 4e-2;
#else
  int n_proc = 3;
  size((4e-2)*n_proc);
  dimensions(nx=n_proc, ny=1);
#endif
  origin (-L0/2, 0);

  DT = 1e-1;

  // kinfolder = "biomass/dummy-solid";
  shift_prod = true;
  kinfolder = "biomass/Solid-only-2507";
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  // No ultimate analysis was found, average from smilar chinese hardwood
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.4229;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")]  = 0.1830;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.1759;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0364;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0081;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0362;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.0245;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0130;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.1000;

  foreach()
    porosity[] = eps0*f[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  moisture0 = solid_mass0*sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")];

  TG[left] = dirichlet (TG0);
  TG[right] = neumann (0);
  TG[top] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[left] = dirichlet (1.);
    } else {
      YG[left] = dirichlet (0.);
    }
  }

  foreach()
    u.x[] = f[] > F_ERR ? 0. : Uin; 
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
  double radius = cbrt (3./2.*statsf(f).sum);

  //log moisture profile
  double moisture = 0.;
  scalar Ymoist = YSList[OpenSMOKE_IndexOfSolidSpecies ("MOIST")];
  foreach (reduction(+:moisture)) {
    if (f[] > F_ERR)
      moisture += (f[]-porosity[])*rhoS*Ymoist[]/f[]*dv(); //Note: (1-e) = (1-ef)!= (1-e)f
  }

  //save temperature profile
  double Tcore  = interpolate (T, 0., 0.);
  double Tr2    = interpolate (T, radius/2., 0.);
  double Tsurf  = interpolate (T, radius, 0.);

  fprintf (fp, "%g %g %g %g %g %g %g\n", t, solid_mass/solid_mass0, Tcore, 
                                      Tr2, Tsurf, radius/(D0/2.), moisture/moisture0);
  fflush(fp);
}

#if TREE
event adapt (i++) {
  scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];
  adapt_wavelet_leave_interface ({T, u.x, u.y, inert}, {f},
    (double[]){1.e-1, 1.e-1, 1.e-1, 1e-1}, maxlevel, minlevel, 1);
}
#endif


event stop (t = tend);

/** 
~~~gnuplot
reset
set terminal svg size 450, 450
set output "huang-mass.svg"

set xlabel "Time [s]"
set ylabel "Normalized mass [-]"
set grid
set xrange [0:650]
set xtics 100
set yrange [0:1.1]
set ytics 0.2
set key top right box width 2.2
set size square

plot  "OutputData-20-773" u 1:2 w l lw 3 lc "black" t "2 cm", \
      "../../data/huang-particleflow/20/773-mass" u 1:2 w p pt 4 ps 1.2 lc "black" notitle, \
      "OutputData-30-773" u 1:2 w l lw 3 lc "dark-green" t "3 cm", \
      "../../data/huang-particleflow/30/773-mass" u 1:2 w p pt 4 ps 1.2 lc "dark-green" notitle
~~~
~~~gnuplot
reset
set terminal svg size 450, 450
set output "huang-shrink-20.svg"
set xlabel "Time [s]"
set ylabel "Shrinking factor [-]"
set grid
set xrange [0:650]
set xtics 100
set yrange [0.5:1.05]
set ytics 0.1
set key top right box width 2.2
set size square 

plot  "OutputData-20-673" u 1:6 w l lw 1 lc "black" t "673K", \
      "OutputData-20-773" u 1:6 w l lw 1 lc "dark-green" t "773K", \
      "OutputData-20-873" u 1:6 w l lw 1 lc "blue" t "873", \
      "OutputData-20-973" u 1:6 w l lw 1 lc "orange" t "973", \
      "../../data/huang-particleflow/20/673-shrinking" u 1:2 w p pt 4 ps 0.8 lc "black" notitle, \
      "../../data/huang-particleflow/20/773-shrinking" u 1:2 w p pt 4 ps 0.8 lc "dark-green" notitle, \
      "../../data/huang-particleflow/20/873-shrinking" u 1:2 w p pt 4 ps 0.8 lc "blue" notitle, \
      "../../data/huang-particleflow/20/973-shrinking" u 1:2 w p pt 4 ps 0.8 lc "orange" notitle
      #"OutputData-30-673" u 1:6 w l lw 1 lc "dark-green" t "3 cm", \
      #"OutputData-30-773" u 1:6 w l lw 1 lc "dark-green" notitle, \
      #"OutputData-30-873" u 1:6 w l lw 1 lc "dark-green" notitle, \
      #"OutputData-30-973" u 1:6 w l lw 1 lc "dark-green" notitle, \
      #"../../data/huang-particleflow/30/673-shrinking" u 1:2 w p pt 1 ps 0.8 lc "dark-green" notitle, \
      #"../../data/huang-particleflow/30/773-shrinking" u 1:2 w p pt 4 ps 0.8 lc "dark-green" notitle, \
      #"../../data/huang-particleflow/30/873-shrinking" u 1:2 w p pt 3 ps 0.8 lc "dark-green" notitle, \
      #"../../data/huang-particleflow/30/973-shrinking" u 1:2 w p pt 2 ps 0.8 lc "dark-green" notitle
~~~
~~~gnuplot
reset
set terminal svg size 450, 450
set output "huang-shrink-30.svg"
set xlabel "Time [s]"
set ylabel "Shrinking factor [-]"
set grid
set xrange [0:650]
set xtics 100
set yrange [0.5:1.05]
set ytics 0.1
set key top right box width 2.2
set size square 

plot  "OutputData-30-673" u 1:6 w l lw 1 lc "black" t "673K", \
      "OutputData-30-773" u 1:6 w l lw 1 lc "dark-green" t "773K", \
      "OutputData-30-873" u 1:6 w l lw 1 lc "blue" t "873K", \
      "OutputData-30-973" u 1:6 w l lw 1 lc "orange" t "973K", \
      "../../data/huang-particleflow/30/673-shrinking" u 1:2 w p pt 4 ps 0.8 lc "black" notitle, \
      "../../data/huang-particleflow/30/773-shrinking" u 1:2 w p pt 4 ps 0.8 lc "dark-green" notitle, \
      "../../data/huang-particleflow/30/873-shrinking" u 1:2 w p pt 4 ps 0.8 lc "blue" notitle, \
      "../../data/huang-particleflow/30/973-shrinking" u 1:2 w p pt 4 ps 0.8 lc "orange" notitle

~~~
**/