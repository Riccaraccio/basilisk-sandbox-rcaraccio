/**
# Rotation of a circular interface - MPI Load Balancing test
*/

#include "advection.h"
#include "vof.h"

scalar c[];
scalar * interfaces = {c}, * tracers = NULL;
int MAXLEVEL = 8;

int main() {
  origin (-0.5, -0.5);
  init_grid (1 << MAXLEVEL);
  run ();
}

#define circle(x,y) (sq(0.1) - (sq(x-0.25) + sq(y)))

event init (i = 0) {
  fraction (c, circle(x,y));
}

#define end 0.785398

scalar test[];
event velocity (i++) {

 // double w_min;
 // if (weights.i > 0)
 //   w_min = statsf(weights).min;
 // else
 //   w_min = 1.;

 // foreach()
 //   test[] = w_min > 0 ? weights[]/w_min : 0.;

#if TREE
  double cmax = 1e-3;
  adapt_wavelet ({c}, &cmax, MAXLEVEL);
#endif

  double a = -8.;
  trash ({u});
  foreach_face(x) u.x[] = - a*y;
  foreach_face(y) u.y[] =   a*x;
}

event logfile (i++) {

  long ncells[npe()]; double load[npe()];
  balance_score (ncells, load, weights);

  if (pid() == 0 && i > 1) {
    fprintf (stderr, "%g %g ", t, parallel_efficency (load));
    for (int ne = 0; ne < npe(); ne++) fprintf (stderr, "%ld ", ncells[ne]);
    for (int ne = 0; ne < npe(); ne++) fprintf (stderr, "%g ", load[ne]);
    fprintf (stderr, "\n");
  }
}

/**
We simulate a hefty computational cost in the cell with a non-zero
volume fraction
*/

event slowdown (i++) {
  foreach()
    if (c[] > 1e-10) {
      double s = 0.;
      for (int k = 0; k < 10000; k++)
        s += sin (k * 1.000001) * cos (k * 0.999999);
      // force the result to be observed so -O2 can't drop the loop
      if (s == HUGE) fprintf (stderr, "unreachable %g\n", s);
    }
}

#if TRACE > 1
event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif

event stop (t = end);

/**
## Results

~~~gnuplot parallel efficiency
set xlabel "Time"
set ylabel "Parallel efficiency"

set yrange [0:1]

plot "log" u 1:2 w l lw 3 notitle
~~~
~~~gnuplot cells per processor (stacked)

set xlabel "Time"
set ylabel "% cells"
set yrange [0:100]
set key outside
set grid front

# 5 processes (-D_MPI=5): ncells in columns 3..7.
# Draw the cumulative sum 3..i normalised by the total (3..7), largest
# first, filled to the x-axis: each smaller fill is painted on top,
# leaving one stacked band per proc that sums to 100%.
plot for [i=7:3:-1] "log" u 1:(100.*(sum [c=3:i] column(c))/(sum [c=3:7] column(c))) \
     w filledcurves x1 title sprintf("proc %d", i-3)

~~~
~~~gnuplot load per processor (stacked)

set xlabel "Time"
set ylabel "% load"
set yrange [0:100]
set key outside
set grid front

# 5 processes (-D_MPI=5): load values in columns 8..12.
plot for [i=12:8:-1] "log" u 1:(100.*(sum [c=8:i] column(c))/(sum [c=8:12] column(c))) \
     w filledcurves x1 title sprintf("proc %d", i-8)

~~~
*/
