#define POROUS_ADVECTION 1
#define F_ERR 1e-10

#include "grid/multigrid.h"
#include "navier-stokes/centered-phasechange.h"
// #include "fractions.h"
#include "two-phase.h"
#include "constant-properties.h"
#include "darcy.h"

double uin = 1.;
double Re = 1.;
int maxlevel = 8;
double epsi0 = 0.7;
double H = 1.;

u.n[left] = dirichlet (y*6/sq(H)*(H-y));
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);
pf[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
pf[right] = dirichlet (0.);

u.n[top] = dirichlet (0.);
u.t[top] = dirichlet (0.);
p[top] = neumann (0.);
pf[top] = neumann (0.);

u.n[bottom] = dirichlet (0.);
u.t[bottom] = dirichlet (0.);
p[bottom] = neumann (0.);
pf[bottom] = neumann (0.);

scalar porosity[];
double rhoG, muG;

int main () {
  muG = 1e-2;
  const face vector muv[] = {muG, muG};
  mu = muv;
  mu1 = mu2 = 1.;
  rho1 = rho2 = 1.;
  rhoG = muG*Re/uin/H;

  size(8*H);
  dimensions(nx=8, ny=1);
  DT = 1e-3;

  Da = (coord) {1e-2/sq(H), 1e-2/sq(H)};
  init_grid (1 << maxlevel);
  run();
}

event init (i = 0) {
  fraction (f, (5+1e-4-x)*(x-3-1e-4));
  fraction (f, (5-x)*(x-3));
  foreach() {
    porosity[] = f[]*epsi0;
    u.x[] = y*6/sq(y)*(H-y);
  }
}

event logprofile (t = end) {
  for (double x = 0.; x < L0; x += L0/(1<<maxlevel))
    fprintf (stderr, "%g %g %g %g %g\n", x/H, 
                                   interpolate (u.x, x, 0.5*H), 
                                   interpolate (p, x, 0.5*H),
                                  //  interpolate (u.x, x, 0.5*H)/interpolate (eps, x, 0.5*H), 
                                  //  interpolate (p, x, 0.5)/interpolate (eps, x, 0.5*H)
                                  );
}

const double max_cfl = 0.1;
event stability (i++) {
  if (CFL > max_cfl)
    CFL = max_cfl;
}

// Avoid the transport of the interface
face vector ufsave[];
event vof (i++) {
  foreach_face() {
    ufsave.x[] = uf.x[];
    uf.x[] = 0.;
  }
}

/**
 * We restore the original velocity field after the interface and tracer 
 * advection has been performed.
 */

event tracer_diffusion (i++) {
  foreach_face()
      uf.x[] = ufsave.x[];
}


event stop (t = 1);

/** 
~~~gnuplot velocity profile
reset

set terminal svg size 750, 350  
set output "plot.svg"

set multiplot layout 1,2

unset key
set size square
set grid
set xtics 2
set xlabel "x [m]"
set ylabel "Velocity [m/s]"
set yrange [1.2:1.6]
set xrange [0:8]
plot "log" u 1:2 w l lc "blue",\
     "../../data/porous-plug/velocity" w p pt 4
#plot "log" u 1:4 w l lc "blue"

set size square
set grid
set xtics 2
set xlabel "x"
set ylabel "Pressure"
set yrange [0:350]
set xrange [0:8]
plot "log" u 1:3 w l lc "blue",\
     "../../data/porous-plug/pressure" w p pt 4
#plot "log" u 1:5 w l lc "blue"

unset multiplot
~~~
*/