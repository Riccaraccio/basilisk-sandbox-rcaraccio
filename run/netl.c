#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define RADIATION_INTERFACE 0
#define MOLAR_DIFFUSION 1
#define FICK_CORRECTED 1
#define MASS_DIFFUSION_ENTHALPY 1

#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
#include "view.h"

double Uin = 0.32; //inlet velocity
double tend = 30; //simulation time

u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
pf[left]     = neumann (0.);
psi[left]    = dirichlet (0.);

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

int maxlevel = 9; int minlevel = 2;
double D0 = 0.7e-3;
double solid_mass0 = 0.;

int main() {
  
  lambdaSmodel = L_CORBETTA;
  TS0 = 300.; TG0 = 773.15;
  rhoS = 1482.5; // 593 / (1-0.6)
  eps0 = 0.6;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_REACTION;

  L0 = 10*D0;
  origin (-L0/2, 0);
  Da = (coord){1e-12, 1e-12};

  DT = 1e-1;

  shift_prod = true;
  kinfolder = "biomass/Solid-only-2507";
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;

  // Residues composition
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")]  = 0.2714;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("GMSW")]  = 0.2061;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")]  = 0.2542;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")]  = 0.0823;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")]  = 0.0054;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")]  = 0.0150;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]   = 0.1018;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]   = 0.0147;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")] = 0.0491;

  foreach()
    porosity[] = eps0*f[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  TG[left] = dirichlet (TG0);
  TG[right] = neumann (0);
  TG[top] = neumann (0);

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

event output (t += 0.1) {
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

  //average temperature of the surface
  double Tsurf_avg = 0.; 
  int count = 0;
  foreach(reduction(+:Tsurf_avg)reduction(+:count)) {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      Tsurf_avg += TInt[];
      count++;
    }
  }
  Tsurf_avg /= count;

  double Tcore  = interpolate (T, 0., 0.);
  double Tr2    = interpolate (T, 0., radius/2.);

  fprintf (fp, "%g %g %g %g %g %g\n", 
            t, solid_mass/solid_mass0, Tcore, Tr2, Tsurf_avg, 
            radius/(D0/2.));

  fflush(fp);
}

#if TREE
event adapt (i++) {
  scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];
  scalar product = YGList_G[OpenSMOKE_IndexOfSpecies ("TAR")];
  adapt_wavelet_leave_interface ({T, u.x, u.y, inert, product}, {f},
    (double[]){1.e-2, 1.e-2, 1.e-2, 1e-2, 1e-2}, maxlevel, minlevel, 2);
  unrefine (x > L0*2/5);
}

#endif

event movie (i += 10) {
  clear();
  box();
  view (ty=-L0*3/2, width=1920, height=640, fov=7);
  draw_vof("f", lw=2);
  squares ("T", min=300, max=800, linear=true);
  vectors ("u", scale=1e-5, level=maxlevel-2);
  mirror ({0.,1.}) {
    draw_vof ("f", lw=2);
    squares ("C6H10O5_G+C6H10O5_S", min=0., max=0.3, linear=true);
  }
  save ("movie.mp4");
}

// event stop (t = tend);
event stop (i = 10);

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
