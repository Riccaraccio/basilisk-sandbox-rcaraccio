#define NO_ADVECTION_DIV 1
#define FSOLVE_ABSTOL 1.e-3

//#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "var-prop.h"
#include "two-phase.h"
#include "temperature-of.h"
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

T[right] = dirichlet (TG0);
T[top] = dirichlet (TG0);

int maxlevel = 4; int minlevel = 2;
double D0 = 1e-3;

scalar omega[];

int main() {

  lambdaS = 0.124069;
  lambdaG = 0.0295641;

  cpS = 2244.92;
  cpG = 1041.52;

  rhoS = 681.042;
  rhoG = 9.75415;

  muG = 1e-5;

  TG0 = 1000.;
  TS0 = 300.;

  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  L0 = 1*D0;
  eps0 = 0.; //No internal gas
  DT = 5e-3;
  // for (maxlevel = 6; maxlevel <= 6; maxlevel++) {
    init_grid(1 << maxlevel);
    run();
  // }
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

int count = 1; // for log profiles
event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  foreach()
    porosity[] = eps0*f[];
  count = 1;
}

event phasechange (i++){
  foreach()
    omega[] = 100.;
}

event logprofile (t += 0.01) {
  fprintf(stderr, "t = %g\n", t);

  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE* fp = fopen (name, "w");

  coord p = {D0/2, 0}; //interface @ y = 0
  double Tinterface = 0.;
  foreach_point (p.x, p.y)
    Tinterface = T[];

  fprintf (fp, "%g %g\n", t, Tinterface);
  fflush (fp);
}

//event adapt (i++) {
//  adapt_wavelet_leave_interface ({T, u.x, u.y, ubf.x, ubf.y}, {f},
//      (double[]){1.e0,1.e-1,1.e-1,1.e-1,1.e-1}, maxlevel, minlevel, 1);
//}

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

event profiles (t = {0.056606, 0.541909}) {
  char name[80];
  fprintf (stderr, "%d %g\n", count, t);
  sprintf (name, "T-Profile-%d-%d", count, maxlevel);

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

/*, maxlevel*
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

plot "../ebvalidation/data/Output-R1000K/Interface.out" u 1:12 w l lw 2 t "microgravity", \
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
plot "T-Profile-1-6" u ($1*1000):2 w l lw 2 t "LEVEL 6", \
     "T-Profile-1-7" u ($1*1000):2 w l lw 2 t "LEVEL 7", \
     "T-Profile-1-8" u ($1*1000):2 w l lw 2 t "LEVEL 8", \
     "../ebvalidation/data/Output-R1000K/Snapshot-0.056606.txt" every 5 u 3:5 w p pt 6 ps 2 lc rgb "black" t "microgravity"

set xlabel "x [mm]"
set ylabel "T [k]"
set key bottom right
set xrange [0:3]
set yrange [250:1100]
set grid
set title "t = 0.541909"
plot "T-Profile-2-6" u ($1*1000):2 w l lw 2 t "LEVEL 6", \
     "T-Profile-2-7" u ($1*1000):2 w l lw 2 t "LEVEL 7", \
     "T-Profile-2-8" u ($1*1000):2 w l lw 2 t "LEVEL 8", \
     "../ebvalidation/data/Output-R1000K/Snapshot-0.541909.txt" every 5 u 3:5 w p pt 6 ps 2 lc rgb "black" t "microgravity"


unset multiplot

~~~




*/
