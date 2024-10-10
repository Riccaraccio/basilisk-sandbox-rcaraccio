#define ufext uf
#define VARPROP

#include "navier-stokes/centered-evaporation.h"
#include "two-phase.h"
#include "evaporation.h"
#include "shrinking.h"
#include "darcy.h"
#include "view.h"

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
psi[top]    = dirichlet (0.);

u.n[bottom] = neumann (0.);
u.t[bottom] = neumann (0.);
p[bottom] = dirichlet (0.);
psi[bottom] = dirichlet (0.);

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);
psi[left]  = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
psi[right] = dirichlet (0.);

int maxlevel = 7; int minlevel = 5;
double D0 = 1.;

double omega = 1.;
double solid_mass0 = 0.;

int main (void) {
  rho1 = 1., rho2 = 1.;
  mu1 = 1.e-5, mu2 = 1.e-5;

  DT = 1.e-1;

  L0 = 2.*D0;
  origin (-0.5*L0, -0.5*L0);

  // for (maxlevel = 6; maxlevel <= 8; maxlevel++) {
  //  init_grid (1 << maxlevel);
  //  run();
  // }

  init_grid (1 << maxlevel);
  run();
  // init_grid (1 << maxlevel);
  // for (GAMMA = 0.4; GAMMA <= 0.9; GAMMA += 0.1) {
  //   run();
  // }
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init (i = 0) {
  fraction (f, circle (x, y, 0.5*D0));

  foreach()
    feps[] = eps0*f[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0)) {
    solid_mass0 += feps[]*rhos*dv();
  }
}

event adapt (i++) {
 scalar fa[];
 foreach()
   fa[] = f[];
 adapt_wavelet ({fa,u.x,u.y}, (double[]){1.e-3,1e-3,1e-3}, maxlevel);
}

double exact (double t) {
  return exp(-omega/rhos*t);
}

event output (t+=1, last) {
  //fprintf (stderr, "t = %g\n", t);
  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double solid_mass = 0.;
  foreach (reduction(+:solid_mass)) {
    solid_mass += feps[]*rhos*dv();
  }
  double relerr = fabs (solid_mass/solid_mass0 - exact (t))/exact (t);
  double radius = sqrt (statsf(f).sum/pi);
  fprintf (fp, "%g %g %g %g\n", t, solid_mass/solid_mass0, exact(t),relerr);
  
  fflush(fp);
}
  
vector ub[];

event movie (t += 1) {

  foreach()
    foreach_dimension()
      ub.x[] = (ubf.x[1,0] + ubf.x[-1,0])/2.;
  
  clear();
  box();
  view(ty = -0.5, width = 1200.);
  draw_vof ("f");
  squares ("feps",min = 0., max = 0.5);
  vectors ("u", scale=0.01);
  mirror ({1.,0.}) {
    draw_vof ("f");
    squares ("feps",min = 0., max = 0.5);
    vectors ("ub", scale = 0.02);
  }
  save ("movie.mp4");
}

event stop (t = 100);

event logfile (t += 1) {
  fprintf (stdout, "%d %g %g\n", i, t, dt);
  fflush (stdout);
}

/**
~~~gnuplot Evolution of the solid mass
reset

set title "Mass"
set xlabel "time [s]"
set ylabel "mass/mass_0 [-]"
set key top right
set size square
set grid

rhop1 = 100.
omega = 1.
f(x) = exp(-omega/rhop1*x)

plot f(x) w p pt 3 t "Analytic", \
      "OutputData-6" u 1:2 w l t "LEVEL 6", \
      "OutputData-7" u 1:2 w l t "LEVEL 7", \
      "OutputData-8" u 1:2 w l t "LEVEL 8"
~~~

~~~gnuplot Convergence rate
reset

stats "OutputData-6" using (last6=$4) nooutput
stats "OutputData-7" using (last7=$4) nooutput
stats "OutputData-8" using (last8=$4) nooutput

set print "Errors.csv"

print sprintf ("%d %.12f", 2**6, last6)
print sprintf ("%d %.12f", 2**7, last7)
print sprintf ("%d %.12f", 2**8, last8)

unset print

reset
set xlabel "N"
set ylabel "Relative Error"

set logscale x 2
set logscale y

set xr[2**5:2**9]

set size square
set grid

f(x) = a*x**(-b)
fit f(x) "Errors.csv" via a,b
ftitle(a,b) = sprintf("%.2fx^{-%.3f}", a, b)

plot "Errors.csv" using 1:2 with points pt 8 ps 2 title "Results", f(x) with lines lw 1 title ftitle(a,b)
~~~
*/
