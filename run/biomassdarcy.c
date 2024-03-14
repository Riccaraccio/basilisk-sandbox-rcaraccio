#define VARPROP
#include "navier-stokes/centered-evaporation.h"
#define ufext uf
#include "two-phase.h"
#include "evaporation.h"
#include "porous-media.h"
#include "darcy.h"
#include "view.h"

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);

u.n[bottom] = neumann (0.);
u.t[bottom] = neumann (0.);
p[bottom] = dirichlet (0.);

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

int maxlevel = 8;
double D0 = 1.;

double rhop1 = 100., rhop2 = 1.;
double eps0 = 0.5;
double omega = 1.;
double solid_mass0 = 0.;

int main (void) {
  rho1 = 100., rho2 = 1.;
  mu1 = 1.e-3, mu2 = 1.e-5;

  DT = 1.e-1;

  L0 = 2.*D0;
  X0 = -0.5*L0;
  Y0 = -0.5*L0;
  init_grid (1 << maxlevel);
  run();

}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init (i = 0) {
  fraction (f, circle (x, y, 0.5*D0));

  // Check back on this
  foreach() {
    feps[] = eps0*f[];
    eps[] = eps0;
  }

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0)) {
    solid_mass0 += feps[]*rhop1*dv();
  }
  fprintf(stderr, "%g\n", solid_mass0);
  zeta_policy = ZETA_SMOOTH;
}

scalar * interfaces_save = interfaces;

event vof (i++) {
  interfaces = NULL;
}

event tracer_advection (i++) {
  interfaces = interfaces_save;
}

// event adapt (i++) {
//   scalar fa[];
//   foreach()
//     if (f[] > F_ERR && f[] < 1.-F_ERR)
//       fa[] = noise();
//   adapt_wavelet ({fa}, (double[]){1.e-3}, maxlevel);
// }

double exact (double t) {
  return exp(-omega/rhop1*t);
}

event output (i++, last) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double solid_mass = 0.;
  foreach (reduction(+:solid_mass)) {
    solid_mass += feps[]*rhop1*dv();
  }
  fprintf (stderr, "%g\n", solid_mass0);
  double relerr = fabs (solid_mass/solid_mass0 - exact (t))/exact (t);
  double radius = sqrt (statsf(f).sum/pi);
  fprintf (fp, "%g %g %g %g\n", t, solid_mass/solid_mass0, relerr, radius/(D0/2));
}

event movie (t += 1) {
  clear();
  draw_vof ("f", lw = 1.5);
  squares ("feps", min = 0., max = eps0);
  save ("movie.mp4");
}

event stop (t = 50);

/**
~~~gnuplot Evolution of the solid mass
reset

set multiplot layout 1,2

set xlabel "time [s]"
set ylabel "mass/mass_0 [-]"
set key top right
set size square
set grid

rhop1 = 1
omega = 1.e-1
f(x) = exp(-omega/rhop1*x)

plot f(x) w p pt 3 t "Analytic", "OutputData-6" u 1:2 w l t "LEVEL 6", "OutputData-7" u 1:2 w l t "LEVEL 7", "OutputData-8" u 1:2 w l t "LEVEL 8", "OutputData-9" u 1:2 w l t "LEVEL 9"

set xlabel "time [s]"
set ylabel "radius/R0 [-]"
set key top right
set size square
set grid

plot "OutputData-6" u 1:4 w l t "LEVEL 6", \
     "OutputData-7" u 1:4 w l t "LEVEL 7", \
     "OutputData-8" u 1:4 w l t "LEVEL 8", \
     "OutputData-9" u 1:4 w l t "LEVEL 9
~~~

~~~gnuplot
reset

stats "OutputData-6" using (last6=$3) nooutput
stats "OutputData-7" using (last7=$3) nooutput
stats "OutputData-8" using (last8=$3) nooutput
stats "OutputData-9" using (last9=$3) nooutput

set print "Errors.csv"

print sprintf ("%d %.12f", 2**6, last6)
print sprintf ("%d %.12f", 2**7, last7)
print sprintf ("%d %.12f", 2**8, last8)
print sprintf ("%d %.12f", 2**9, last9)

unset print

reset
set xlabel "N"
set ylabel "Relative Error"

set logscale x 2
set logscale y

set xr[2**5:2**10]
#set yr[1e-4:10]

set size square
set grid

f(x) = a*x**(-b)
fit f(x) "Errors.csv" via a,b
ftitle(a,b) = sprintf("%.2fx^{-%.3f}", a, b)

plot "Errors.csv" w p pt 8 ps 2 title "Results", \
     f(x) w l lw 1 t ftitle (a, b)
~~~
*/