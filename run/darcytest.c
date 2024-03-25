#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "darcy.h"
#include "view.h"

// BOUNDARIES: left and right open, top and bottom no-slip

double U0 = 0.1;

u.n[left] = dirichlet (U0);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);

u.n[top] = dirichlet (0.);
u.t[top] = dirichlet (0.);
p[top] = neumann (0.); 

u.n[bottom] = dirichlet (0.);
u.t[bottom] = dirichlet (0.);
p[bottom] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

int maxlevel = 7;

int main() {
  rho1 = 1., rho2 = 1.;
  mu1 = 5e-3, mu2 = 5e-3;
  L0 = 16.;

  for (maxlevel = 7; maxlevel <= 9; maxlevel++){
  init_grid (1 << maxlevel);
  run();
  }

}

event init (i = 0) {
  mask (y >= L0/4 ? top : none); // domain is 16x4
  fraction (f, x>7.5 && x<8.5); // interface centerd in x=8, with 1 unit width

  foreach() {
    u.x[] = U0;
  }
  foreach_face(x){
    uf.x[] = U0;
  }
}

// Avoid tranporting the interface: the solid is fixed
static scalar * interfaces_save = NULL;

event vof (i++) {
  interfaces_save = interfaces; 
  interfaces = NULL;
}

event tracer_advection (i++) {  
  interfaces = interfaces_save;
}

// Movie event
event movie (t += 1) {
  clear();
  view (tx = -0.5);
  draw_vof ("f", lw = 1.5);
  squares ("u.x");
  save ("movie.mp4");
}

// Output data
event logprofile (t = end){
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  FILE * fp = fopen (name, "w");

  double step = L0/(1<<maxlevel);
  for (double x = 0; x<L0; x+=step){
    fprintf(fp, "%g %g %g\n", x, interpolate(u.x, x, L0/8), interpolate(p, x, L0/8));
  }
  fprintf (fp, "\n\n");
}

event end (t = 100) {} 
// PLOTS
/**
~~~gnuplot Velocity on y = 2
reset

set title "Velocity u_x  on y = 2"
set xlabel "Longitudinal position x [m]"
set ylabel "velocity [m/s]"
set key bottom right
set size square
set grid

plot  "OutputData-7" u 1:2 w l t "LEVEL 7", \
      "OutputData-8" u 1:2 w l t "LEVEL 8", \
      "OutputData-9" u 1:2 w l t "LEVEL 9", \
      "Data/OpenFoamValidationData.csv" u 1:8 w l t "pisoFOAM"
~~~

~~~gnuplot Pressure on y = 2
reset
set title "Pressure on y = 2"
set xlabel "longitudianl position [m]"
set ylabel "Pressure "
set key top right
set size square
set grid

plot  "OutputData-7" u 1:3 w l t "LEVEL 7", \
      "OutputData-8" u 1:3 w l t "LEVEL 8", \
      "OutputData-9" u 1:3 w l t "LEVEL 9", \
      "Data/OpenFoamValidationData.csv" u 1:10 w l t "pisoFOAM"

~~~
*/