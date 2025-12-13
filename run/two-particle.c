#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define RADIATION_INTERFACE 0.9
//#define MOLAR_DIFFUSION 1
//#define FICK_CORRECTED 1
//#define MASS_DIFFUSION_ENTHALPY 1

#ifndef TWO_PARTICLES
# define TWO_PARTICLES 0
#endif

#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"

double Uin = 0.1;
double tend = 600.; //simulation time

u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
pf[left]     = neumann (0.);
psi[left]    = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);
psi[right]    = dirichlet (0.);

int maxlevel = 8; int minlevel = 2;
double D0 = 2e-2;

int main() {
  
  lambdaS = 0.1987;
#ifdef VARPROP
  lambdaSmodel = L_HUANG;
#endif
  TS0 = 300.; TG0 = 773;
  rhoS = 920;
  eps0 = 0.39;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_CONST;

  L0 = D0*8;
  origin (-L0/2, 0);

  DT = 1e-1;

  shift_prod = true;
  kinfolder = "biomass/dummy-solid";
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

double solid_mass0;
event init(i = 0) {
  
  if (TWO_PARTICLES) {
    fraction (f, circle (x - D0, y, 0.5*D0) *
                 -circle (x + D0, y, 0.75*D0));
  } else {
    fraction (f, circle (x, y, 0.5*D0));
  }

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  foreach()
    porosity[] = eps0*f[];
  
  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0)) {
    if (TWO_PARTICLES) {
      if (x > 0)
        solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f
    } else {
      solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f
    }
  }

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

event output (t += 1) {
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass)) {
    if (TWO_PARTICLES) {
      if (x > 0)
        solid_mass += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f
    } else {
      solid_mass += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f
    }
  }
  fprintf (stderr, "%g %g\n", t, solid_mass/solid_mass0);
}

#if TREE
event adapt (i++) {
  scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];
  adapt_wavelet_leave_interface ({T, u.x, u.y, inert}, {f},
    (double[]){1.e-1, 1.e-1, 1.e-1, 1e-1}, maxlevel, minlevel, 2);
}
#endif

event outputfields (t = {100, 200, 300, 600}) {
  char name[80];
  sprintf (name, "OutputFields-%d", (int)(t));
  static FILE * fs = fopen (name, "w");

  output_field ({T, f, u.x, u.y}, fs);
  fflush (fs);
}

event stop (t = tend);

/** 
~~~gnuplot
reset
#set terminal svg size 450, 450
#set output "huang-mass.svg"

set terminal epslatex size 3.6, 3.6 color colortext
set output "huang-mass.tex"

set xlabel "Time [s]"
set ylabel "Normalized mass [-]"
set grid
set xrange [0:650]
set xtics 100
set yrange [0:1.1]
set ytics 0.2
set key top right box width 2.2
set size square

plot  "OutputData-20-773" u 1:2 w l lw 6 lc "black" t "2 cm", \
      "../../data/huang-particleflow/20/773-mass" u 1:2 w p pt 64 ps 2 lw 3 lc "black" notitle, \
      "OutputData-30-773" u 1:2 w l lw 6 lc "dark-green" t "3 cm", \
      "../../data/huang-particleflow/30/773-mass" u 1:2 w p pt 64 ps 2 lw 3 lc "dark-green" notitle
~~~
~~~gnuplot
reset
set terminal svg size 450, 450
set output "huang-shrink-20.svg"

#set terminal epslatex size 3.6, 3.6 color colortext
#set output "huang-shrink-20.tex"

set xlabel "Time [s]"
set ylabel "Shrinking factor [-]"
set grid
set xrange [0:650]
set xtics 100
set yrange [0.5:1.05]
set ytics 0.1
set key bottom right box width 2.2
set size square 

 plot "../../data/huang-particleflow/20/673-shrinking" u 1:2 w p pt 64 ps 1.2 lw 3 lc "black" notitle, \
      "../../data/huang-particleflow/20/773-shrinking" u 1:2 w p pt 64 ps 1.2 lw 3 lc "dark-green" notitle, \
      "../../data/huang-particleflow/20/873-shrinking" u 1:2 w p pt 64 ps 1.2 lw 3 lc "blue" notitle, \
      "../../data/huang-particleflow/20/973-shrinking" u 1:2 w p pt 64 ps 1.2 lw 3 lc "orange" notitle, \
      "OutputData-20-673" u 1:6 w l lw 3 lc "black" t "673K", \
      "OutputData-20-773" u 1:6 w l lw 3 lc "dark-green" t "773K", \
      "OutputData-20-873" u 1:6 w l lw 3 lc "blue" t "873K", \
      "OutputData-20-973" u 1:6 w l lw 3 lc "orange" t "973K"
~~~

~~~gnuplot
reset
set terminal svg size 450, 450
set output "huang-shrink-30.svg"

#set terminal epslatex size 3.6, 3.6 color colortext
#set output "huang-shrink-30.tex"

set xlabel "Time [s]"
set ylabel "Shrinking factor [-]"
set grid
set xrange [0:650]
set xtics 100
set yrange [0.5:1.05]
set ytics 0.1
set key bottom right box width 2.2
set size square 

plot  "../../data/huang-particleflow/30/673-shrinking" u 1:2 w p pt 64 ps 1.2 lw 3 lc "black" notitle, \
      "../../data/huang-particleflow/30/773-shrinking" u 1:2 w p pt 64 ps 1.2 lw 3 lc "dark-green" notitle, \
      "../../data/huang-particleflow/30/873-shrinking" u 1:2 w p pt 64 ps 1.2 lw 3 lc "blue" notitle, \
      "../../data/huang-particleflow/30/973-shrinking" u 1:2 w p pt 64 ps 1.2 lw 3 lc "orange" notitle, \
      "OutputData-30-673" u 1:6 w l lw 3 lc "black" t "673K", \
      "OutputData-30-773" u 1:6 w l lw 3 lc "dark-green" t "773K", \
      "OutputData-30-873" u 1:6 w l lw 3 lc "blue" t "873K", \
      "OutputData-30-973" u 1:6 w l lw 3 lc "orange" t "973K"
~~~

~~~gnuplot
reset
#set terminal svg size 450, 450
#set output "huang-t-core.svg"

set terminal epslatex size 3.6, 3.6 color colortext
set output "huang-t-core.tex"

set xlabel "Time [s]"
set ylabel "Temperature [K]"
set grid
set xrange [0:650]
set xtics 100
set yrange [200:850]
set ytics 100
set key bottom right box width 2.1
set size square

plot "OutputData-20-773" u 1:3 w l lw 6 lc "black" t "2 cm", \
     "../../data/huang-particleflow/20/773-T-core" u 1:2 w p pt 64 ps 2 lw 5 lc "black" notitle, \
     "OutputData-30-773" u 1:3 w l lw 6 lc "dark-green" t "3 cm", \
     "../../data/huang-particleflow/30/773-T-core" u 1:2 w p pt 64 ps 2 lw 5 lc "dark-green" notitle
~~~

~~~gnuplot
reset
set terminal svg size 450, 450
set output "huang-t-r2.svg"
set xlabel "Time [s]"
set ylabel "Temperature [K]"
set grid
set xrange [0:650]
set xtics 100
set yrange [200:850]
set ytics 100
set key bottom right box width 2.1
set size square

plot "OutputData-20-773" u 1:4 w l lw 3 lc "black" t "2 cm", \
     "../../data/huang-particleflow/20/773-T-r2" u 1:2 w p pt 4 ps 1.2 lc "black" notitle, \
     "OutputData-30-773" u 1:4 w l lw 3 lc "dark-green" t "3 cm", \
     "../../data/huang-particleflow/30/773-T-r2" u 1:2 w p pt 4 ps 1.2 lc "dark-green" notitle
~~~

~~~gnuplot
reset
set terminal svg size 450, 450
set output "huang-t-surf.svg"
set xlabel "Time [s]"
set ylabel "Temperature [K]"
set grid
set xrange [0:650]
set xtics 100
set yrange [200:850]
set ytics 100
set key bottom right box width 2.1
set size square

plot "OutputData-20-773" u 1:5 w l lw 3 lc "black" t "2 cm", \
     "../../data/huang-particleflow/20/773-T-surf" u 1:2 w p pt 4 ps 1.2 lc "black" notitle, \
     "OutputData-30-773" u 1:5 w l lw 3 lc "dark-green" t "3 cm", \
     "../../data/huang-particleflow/30/773-T-surf" u 1:2 w p pt 4 ps 1.2 lc "dark-green" notitle
~~~
**/
