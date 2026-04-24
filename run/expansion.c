/**
# Time-reversed expansion of a VOF field

This test case is inspired by the test case 
[Time-reversed VOF advection in a vortex](/src/test/reversed.c).

We want to test the convergence of the error for the transport
of an interface following a non-divergence-free velocity field.
The velocity reverts at half time and the interface should come
back to the original position. */

attribute {
  double sigma;
}

#include "grid/multigrid.h"
#include "advection.h"
#include "vof.h"

scalar f[], cf[];
scalar * interfaces = {f}, * tracers = NULL;
int maxlevel;

int main (void) {
  //DT = 0.001;
  
  f.sigma = 0.03;
  f.tracers = {cf};

  /**
  We run the simulation at different grid refinement levels */
  for (maxlevel = 4; maxlevel <= 9; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}
/**
The initial interface is a circle of radius 0.5 centered on
(0, 0). We use the levelset function `circle()` to define 
this interface.
The period of the stretching cycle is set to 1. */

#define circle(x,y,r) (sq(r) - sq(x) - sq(y))
const double T = 1.;
const double R0 = 0.4;

face vector alpha[];

event init (t = 0) {
  fraction (f, circle(x,y,R0));
  foreach()
    cf[] = f[];

  foreach_face()
    alpha.x[] = 1.;
}

event velocity (i++) {
#if TREE
  adapt_wavelet ({f}, (double[]){1e-3}, maxlevel, list = {f, cf});
#endif

  trash({u});
  const double Umag = sin(2*M_PI/T *t);
  foreach_face() {
    coord o = {x,y,z};
    double mag = sqrt (sq(x) + sq(y));
    u.x[] = o.x*Umag / (mag + 1e-10);
  }
}

double dtmax;

event stability (i++)
{

  /**
  We first compute the minimum and maximum values of $\alpha/f_m =
  1/\rho$, as well as $\Delta_{min}$. */

  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin))
    if (fm.x[] > 0.) {
      if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
      if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
      if (Delta < dmin) dmin = Delta;
    }
  double rhom = (1./amin + 1./amax)/2.;

  /**
  The maximum timestep is set using the sum of surface tension
  coefficients. */

  double sigma = 0.;
  for (scalar c in interfaces)
    sigma += c.sigma;
  if (sigma) {
    double dt = sqrt (rhom*cube(dmin)/(pi*sigma));
    if (dt < dtmax)
      dtmax = dt;
  }
}
event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}

event logfile (t = {0,T}) {
  stats s = statsf (f);

  /**
  We compute the minimum and maximum concentration. They should both
  be equal to one. */
  
  stats sc = statsf (cf);
  double cmin = HUGE, cmax = 0.;
  foreach (reduction(min:cmin) reduction(max:cmax))
    if (f[] > 1e-6) { // round-off errors are a problem
      double c = cf[]/f[];
      if (c < cmin) cmin = c;
      if (c > cmax) cmax = c;
    }
  fprintf (stderr, "# t\t\tf.sum\t\tf.min\t\tf.max\n");
  fprintf (stderr, "# %f %.12f %.f %g\n", t, s.sum, s.min, s.max);
  fprintf (stderr, "# t\t\tcf.sum\t\tc.min - 1\tc.max - 1\n");
  fprintf (stderr, "# %f %.12f %.11f %.11f\n",
	   t, sc.sum, fabs(cmin - 1.), fabs(cmax - 1.));
}

event field (t = T) {
  scalar e[];
  fraction (e, circle(x,y,R0));
  foreach()
    e[] -= f[];
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}

/*
~~~gnuplot Error vs grid size
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)

f(x) = a+b*x
fit f(x) "log" u (log($1)):(log($2)) via a,b

f2(x) = a2+b2*x
fit f2(x) "log" u (log($1)):(log($3)) via a2,b2

f3(x) = a3+b3*x
fit f3(x) "log" u (log($1)):(log($4)) via a3,b3

set xlabel "Maximum resolution"
set ylabel "Error"
set logscale
set xrange [8:1024]
set xtics 8,2,1024
set grid ytics
set key bottom left box
set yrange [0.00001:1]
set format y "10^{%T}"
plot "log" u 1:2 t "avg", exp(f(log(x)))  t ftitle(a,b), \
     "log" u 1:3 t "rms", exp(f2(log(x))) t ftitle(a2,b2), \
     "log" u 1:4 t "max", exp(f3(log(x))) t ftitle(a3,b3)
~~~
*/
