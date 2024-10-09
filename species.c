#define ufext uf
#define NO_ADVECTION_DIV 1
#define FSOLVE_ABSTOL 1.e-3

//#include "axi.h" 
#include "navier-stokes/centered-evaporation.h"
#include "two-phase.h"
#include "evaporation.h"
#include "temperature.h"
#include "shrinking.h"
//#include "darcy.h"
#include "view.h"

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
psi[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
psi[right] = dirichlet (0.);

int maxlevel = 6; int minlevel = 2;
double D0 = 1e-3;

scalar omega[];
double solid_mass0;

double lambda1 = 0.124069;
double lambda2 = 0.0295641;

double cp1 = 2244.92;
double cp2 = 1041.52;

double TG0 = 1000.;
double TS0 = 300.;

int main() {
  rho1 = 681.042, rho2 = 9.75415;
  mu1 = 0.00037446, mu2 = 2.02391e-5;
  L0 = 3.5*D0;
  for (maxlevel = 7; maxlevel <=7; maxlevel++) {
    init_grid(1 << maxlevel);
    run();
  }
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  foreach()
    feps[] = eps0*f[]; //TODO: move to shrinking.h

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0)) {
    solid_mass0 += feps[]*rhos*dv();
  }
}

event bcs (i=0) {
  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
}

event phasechange (i++){
  foreach()
    omega[] = 0.; //fixed interface
}

event logprofile (t += 0.01) {
  fprintf(stderr, "t = %g\n", t);

  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE* fp = fopen (name, "w");

  coord p = {D0/2, 0}; //interface @ y = 0
  double Tinterface = 0.;
  foreach_point (p.x, p.y)
    Tinterface = TInt[];

  fprintf (fp, "%g %g\n", t, Tinterface);
  fflush (fp);
}

event adapt (i++) {
  adapt_wavelet_leave_interface ({T, u.x, u.y, ubf.x, ubf.y}, {f},
      (double[]){1.e0,1.e-1,1.e-1,1.e-1,1.e-1}, maxlevel, minlevel, 1);
}

//event movie (t += 10) {
//  clear();
//  view (ty=-0.5, width=1400.);
//  draw_vof ("f", lw=2);
//  squares ("T", min=TS0, max=TG0, linear=true);
//  mirror ({1.,0.}) {
//    vectors ("ubf", scale=10);
//    draw_vof ("f", lw=2);
//    squares ("feps", min=0., max=eps0, linear=true);
//  }
//  save ("movie.mp4");
//}

int count = 1;
event profiles (t = {0.056606, 0.541909}) {
  char name[80];
  fprintf (stderr, "%d %g\n", count, t);
  sprintf (name, "T-Profile-%d", count);

  FILE* fpp = fopen (name, "w");
  for (double x = 0.; x < L0; x+= 0.5*L0/(1<<maxlevel)) {
    fprintf (fpp, "%g %g\n", x, interpolate (T, x, 0.));
  }

  fflush (fpp);
  count++;
}


event stop (t=1) {
  return 1;
}

/**
## Results

~~~gnuplot Interface temperature
reset
set terminal pngcairo enhanced size 960, 540
set xlabel "t [s]"
set ylabel "T [k]"
set key bottom right
set xrange [0:1]
set yrange [250:600]
set grid

plot "data/Output-R1000K/Interface.out" u 1:12 w l lw 2 t "microgravity", \
     "OutputData-6" u 1:2 w l lw 2 t "LEVEL 6", \
     "OutputData-7" u 1:2 w l lw 2 t "LEVEL 7", \
     "OutputData-8" u 1:2 w l lw 2 t "LEVEL 8", \
     "../c7pathak/log" u 1:2 w l lw 2 t "Edo, LEVEL 6"
~~~

~~~gnuplot T profiles
reset
set terminal pngcairo enhanced size 960, 960
set multiplot layout 2,1 title "Temperature snapshots"

set xlabel "x [mm]"
set ylabel "T [k]"
set key bottom right
set xrange [0:3]
set yrange [250:1100]
set grid
set title "t = 0.056606"
plot "T-Profile-1" u ($1*1000):2 w l lw 2 t "Basilisk", \
     "data/Output-R1000K/Snapshot-0.056606.txt" every 5 u 3:5 w p ps 2 t "microgravity"

set xlabel "x [mm]"
set ylabel "T [k]"
set key bottom right
set xrange [0:3]
set yrange [250:1100]
set grid
set title "t = 0.541909"
plot "T-Profile-2" u ($1*1000):2 w l lw 2 t "Basilisk", \
     "data/Output-R1000K/Snapshot-0.541909.txt" every 5 u 3:5 w p ps 2 t "microgravity"

unset multiplot

~~~




*/
