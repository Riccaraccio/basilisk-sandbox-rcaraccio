/**
# Time-reversed VOF advection in a vortex

This classical test advects and stretches an initially circular
interface in a non-divergent vortical flow. The flow reverts in time
and the interface should come back to its original position. The
difference between the initial and final shapes is a measure of the
errors accumulated during advection.

This also checks that a constant "VOF concentration" remains constant
when it is advected with the VOF tracer.

We will need the advection solver combined with the VOF advection
scheme. 

This is a modified version of reversed.c to test the implementation
of the new load balancing algorithm.
*/

#include "advection.h"
#include "vof.h"
#include "view.h"

scalar f[], cf[];
scalar * interfaces = {f}, * tracers = NULL;
int MAXLEVEL = 7;

int main (int argc, char * argv[])
{
  origin (-0.5, -0.5);
  DT = .1[0,1];

  f.tracers = {cf};

  init_grid (1 << MAXLEVEL);
  run();
}

#define circle(x,y) (sq(0.2) - (sq(x + 0.2) + sq(y + .236338)))
const double T = 15.;

event init (i = 0)
{
  fraction (f, circle(x,y));
  foreach()
    cf[] = f[];

  weights = new scalar;
}

event log (i += 10) {
  long ncells[npe()]; double load[npe()];
  balance_score (ncells, load, weights);

  if (pid() == 0 && i > 1) {
    fprintf (stderr, "%g %g ", t, parallel_efficency (load));
    for (int ne = 0; ne < npe(); ne++) fprintf (stderr, "%ld ", ncells[ne]);
    for (int ne = 0; ne < npe(); ne++) fprintf (stderr, "%g ", load[ne]);
    fprintf (stderr, "\n");
  }
}

event velocity (i++) {

#if TREE
  foreach()
    weights[] = 20*f[] + 1.;
    //weights[] = 1.;

  adapt_wavelet ({f}, {5e-3}, MAXLEVEL, list = {f, cf, weights});
#endif

  vertex scalar psi[];
  double a = 1.5, k = pi;
  foreach_vertex()
    psi[] = - a*sin(2.*pi*t/T)*sin(k*(x + 0.5))*sin(k*(y + 0.5))/pi;
  
  trash ({u});
  coord f = {-1.,1.};
  foreach_face()
    u.x[] = f.x*(psi[0,1] - psi[])/Delta;
}

event slowdown (i++) {
  const struct timespec ts = {0, 1e5}; //sleep for 0.1 ms
  foreach()
    if (f[] > 1e-10)
      nanosleep(&ts, NULL);
}

event movie (i += 10) {

  scalar pid_field[];
  foreach()
    pid_field[] = pid();
 
  clear();
  box();
  draw_vof ("f", lw=2);
  cells();
  squares ("pid_field", min=0, max=npe() - 1);
  save ("mpi.mp4");
}

#if TRACE > 1
event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif

event stop (t = T);

/**
## Results

~~~gnuplot parallel efficiency

set xlabel "Time"
set ylabel "Parallel efficiency"

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

## See also

[Benchmark on GPUs](/src/grid/gpu/Benchmarks.md#time-reversed-vof-advection-in-a-vortex)
*/
