/**
# Rotation of a circular interface - MPI Load Balancing test

This test is a modified version of [rotate.c](/src/test/rotate.c).
The goal is to test the correct split of load between processors, given a user 
assinged scalar field 'weights', a representative of the computational load being
carried out in a specific cell.
This is also a twin to the test [balance-rotate-auto.c](/src/test/balance-rotate-auto.c)
where instead the weights are detected automatically.

We advect a circular interface along a circular path. To simulate an inbalance 
in the load between cells, we add an intensive calculation only in the region
marked by the 'c' field. 

The total load is defined as the sum of the individual weights of each cells.
Ideally, each processor gets total_load/n_proc amout of work.
This is also the configuration that minimze the idle time between processor 
waiting for the other to finish their computations. 
*/

#include "advection.h"
#include "vof.h"

/**
Diagnostic functions needed for intresting statistics.
'balance_score' computes the number of cells and the load assigned to 
each processor.
'parallel_efficency' gives a quantitative measure of the imbalance.
*/
void balance_score (long* counter, double* load, (const) scalar w) {
  for (int ne = 0; ne < npe(); ne++) {
    counter[ne] = 0;
    load[ne] = 0.;
  }

  foreach() {
    counter[pid()] += 1;
    load[pid()] += w[];
  }

  if (pid() == 0) {
    MPI_Reduce(MPI_IN_PLACE, counter, npe(), MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, load, npe(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(counter, NULL, npe(), MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(load, NULL, npe(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
}

double parallel_efficency (double* load) {
  //compute parallel efficency given the load array
  double tot_load = 0., max_load = -1;
  for (int ne = 0; ne < npe(); ne++) {
    tot_load += load[ne];
    if (load[ne] > max_load)
      max_load = load[ne];
  }
  return (tot_load/npe())/max_load;
}

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
  weights = new scalar;
}

#define end 0.785398

event velocity (i++) {
  foreach()
    weights[] = c[]*20. + 1.;

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
volume fraction.
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

plot "log" u 1:2 w l lw 3 lc "red" notitle
~~~
~~~gnuplot cells per processor (stacked)

set title "Realtive cell share in each processor"
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

set title "Realtive load share in each processor"
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
