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
#include "view.h"

double Uin = 0.25;
double tend = 800.; //simulation time
//double tend = 5; //simulation time

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
    vertex scalar phi[];
    foreach_vertex() {
      double theta = atan2(y, x + D0), r = sqrt (sq(x + D0) + sq(y));
      phi[] = (0.25 + 0.05*cos(6*theta))*D0*0.75/0.3 - r;
      phi[] *= -circle(x - D0, y, 0.5*D0);
    }
    fractions (phi, f);
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

//event stop (t = tend) {
//  clear();
//  view (ty=-0.5, width = 1400., height=1400.);
//  squares ("T", min=300, max=773, linear=true, spread=-1);
//  draw_vof ("f", lw=2);
//  cells();
//  save ("final-temperature.png");
//}

/** 
~~~gnuplot
reset
set terminal svg size 450, 450
#set terminal epslatex size 3.6, 3.6 color colortext
set output "mass.svg"

set xlabel "Time [s]"
set ylabel "Normalized mass [-]"
set grid
set xrange [0:600]
set xtics 100
set yrange [0:1.1]
set ytics 0.2
set key top right box
set size square

plot  "log-one" u 1:2 w l lw 6 dt 1 lc "black" t "one particle", \
      "log-two" u 1:2 w l lw 6 dt 3 lc "black" t "two particles"
~~~
**/
