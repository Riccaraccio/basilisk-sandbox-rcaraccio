#define POROUS_ADVECTION
#define F_ERR 1e-8

#include "navier-stokes/centered-phasechange.h"
#include "fractions.h"
#include "darcy.h"

double uin = 0.8;
u.n[left] = dirichlet (uin);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);
pf[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
pf[right] = dirichlet (0.);

int maxlevel = 10;
double eps0 = 0.6;
scalar porosity[];
double rhoG, muG;

int main () {
  const face vector muv[] = {1e-1, 1e-1};
  mu = muv;
  muG = 1e-1;
  rhoG = 1.;

  DT = 1e-4;
  Da = (coord) {1e-4, 1e-4};
  stokes = true;
  init_grid (1 << maxlevel);
  run();
}

scalar f[];
vector v[];
scalar true_p[];

event init (i = 0) {

  fraction (f, x > 0.4 && x < 0.6);

  foreach()
    porosity[] = f[]*eps0;

  foreach()
    u.x[] = (1. - f[])*uin;
}
event endtimestep (i++) {

  scalar eps[];

  foreach()
    eps[] = (1.-f[]) + porosity[];

  foreach()
    foreach_dimension()
      v.x[] = u.x[]/eps[];

  foreach()
    true_p[] = p[]/eps[];
}

event logprofile (t = end) {
  char name[80];
  sprintf (name, "Output");
  static FILE * fp = fopen (name, "w");

  for (double x = 0.; x < L0; x += L0/(1<<maxlevel))
    fprintf (fp, "%g %g %g %g %g\n", x, interpolate (u.x, x, 0.5), interpolate (v.x, x, 0.5),  interpolate (p, x, 0.5),  interpolate (true_p, x, 0.5));
  fflush (fp);
}

event stop (i = 100);

/** 
~~~gnuplot velocity profile
reset
#set terminal epslatex size 3.6, 3.6 color colortext
#set output "septum-velocity-darcy.tex"

set terminal svg size 350, 350  
set output "septum-velocity-darcy.svg"

set xlabel "x [m]"
set ylabel "Velocity [m/s]"
set yrange [0:2]
set key top right box width 1 opaque
set samples 20
set xtics 0.2

set object 1 rectangle from 0.4,graph 0 to 0.6,graph 1 behind fc "black" fs solid 0.2
set label 1 at 0.5,0.1 "porous region" center front

u_in = 0.8
epsilon = 0.6
analytical_u(x) = u_in
analytical_v(x) = (x < 0.4) ? u_in : (x > 0.6) ? u_in : u_in/epsilon

plot  analytical_u(x) w p pt 64 ps 2 lc "dark-green" notitle, \
      analytical_v(x) w p pt 64 ps 2 lc "blue" notitle, \
      "Output"           u 1:2 w l lw 4 lc "dark-green" t "v", \
      "Output" every ::1 u 1:3 w l lw 4 lc "blue" t "u"
~~~

~~~gnuplot pressure profile
reset
#set terminal epslatex size 3.5, 3.5 color
#set output "septum-pressure-darcy.tex"

set terminal svg size 350, 350
set output "septum-pressure-darcy.svg"

set xlabel "x"
set ylabel "pressure"
#set yrange [-1:3.5]
set samples 20
set key top right box width 1 opaque
set xtics 0.2

set object 1 rectangle from 0.4,graph 0 to 0.6,graph 1 behind fc "black" fs solid 0.2
set label 1 at 0.5, -0.8 "porous region" center front

mu = 1e-1
u_in = 0.8
epsilon = 0.6
K = 1e-4
rho = 1
F = 1.75/sqrt (150*epsilon**3)
B = mu*u_in*epsilon/K + rho*F*epsilon*u_in*u_in/K**0.5
analytical_p(x) = (x < 0.4) ? B*(0.6-0.4) : \
                  (x > 0.6) ? 0 : \
                  rho*u_in*u_in*(1-1/epsilon) + B*(0.4-x + (0.6-0.4))

analytical_p_prime(x) = (x < 0.4) ? analytical_p(x) : analytical_p(x)/epsilon

plot "Output-old" u 1:4 w l lw 4 lc "blue" t "old imp", \
      "Output" u 1:4 w l lw 4 lc "red" t "new imp", \
      analytical_p(x) w p pt 4 ps 1.2 lc "blue" t "analytical p"
      #"Output" u 1:5 w l lw 2 lc "web-green" t "p", \
      #analytical_p_prime(x) w p pt 4 lc "web-green" t "analytical p", \
~~~
*/