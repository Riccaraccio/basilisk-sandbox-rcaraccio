#define NO_ADVECTION_DIV 1

//#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "prop.h"
#include "two-phase.h"
#include "shrinking.h"
#include "balances.h"

u.n[top] = neumann(0.);
u.t[top] = neumann(0.);
p[top] = dirichlet(0.);
psi[top] = dirichlet(0.);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right] = dirichlet(0.);
psi[right] = dirichlet(0.);

int maxlevel = 7, minlevel = 2;
double D0 = 1;
scalar omega[];

int amr = 0;
int axi = 0;
int zetamodel = 0;

int main() {
  eps0 = 0.4;
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;
  L0 = 1.5*D0;
  DT = 1e-1;

  rhoS = 100.;
  rhoG = 1.;

// for (maxlevel=5; maxlevel<=7; maxlevel++) {
//   for (zetamodel=0; zetamodel <= 2; zetamodel++) {
//     for (amr=0; amr <=1; amr++) {
//         fprintf(stderr, "Running maxlevel %d zetamodel %d, amr %d\n",
//           maxlevel, zetamodel, amr);
//         init_grid(1 << maxlevel);
//         run();
//     }
//   }

amr = 0;
for (maxlevel=5; maxlevel<=8; maxlevel++) {
  for (zetamodel=0; zetamodel <= 1; zetamodel++) {
        fprintf(stderr, "Running maxlevel %d zetamodel %d, amr %d\n",
          maxlevel, zetamodel, amr);
        init_grid(1 << maxlevel);
        run();
    }
  }

  // maxlevel = 7;
  // amr = 0;
  // for (zetamodel=0; zetamodel <= 1; zetamodel++) {
  //   fprintf(stderr, "Running maxlevel %d zetamodel %d, amr %d\n",
  //     maxlevel, zetamodel, amr);
  //   init_grid(1 << maxlevel); 
  //   run();
  // }
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

event init(i = 0) {
  fraction (f, circle(x, y, 0.5*D0));

  foreach()
    porosity[] = eps0*f[];

  //0: SHRINK, 
  //1: SWELLING, 
  //2: SMOOTH, 
  //3: SHARP, 
  //4: LEVELSET
  zeta_policy = zetamodel;

  foreach()
    omega[] = 10.;

#ifdef BALANCES
  //specify the output file name
  sprintf(mb.name, "balances-%d-%d-%d", maxlevel, zetamodel, amr);
#endif
}

// update the porosity field
event chemistry(i++) {
  foreach() {
    if (f[] > F_ERR) {
      porosity[] = porosity[]/f[];
      porosity[] += (omega[]*(1.-porosity[])*(1.-zeta[])/rhoS)*dt;
      porosity[] *= f[];
    }
  }
}

event logfile (t +=1 ) {
  fprintf (stderr, "%g\n", t);
}

event adapt (i++) {
  if (amr > 0) {
  adapt_wavelet_leave_interface({porosity, u.x, u.y}, {f}, 
    (double[]){1e-2, 1e-2, 1e-2}, maxlevel, minlevel);
  }
}

#if TRACE > 1
event profiling (i += 10) {
  static FILE* fp = fopen ("profiling", "w");
  trace_print(fp, 1);
}
#endif

event stop (t = 10);

/** 
~~~gnuplot Mass conservation profiles at level 8
reset
set terminal svg size 1000, 700
set multiplot layout 2,1 title "Mass conservation profiles, LEVEL 8 AXI" \
  margins 0.1, 0.9, 0.1, 0.9 spacing 0.1, 0.1

# Calculate differences in total mass
stats "balances-8-0-0" using 1:4 nooutput
diff_00 = STATS_max_y - STATS_min_y

stats "balances-8-1-0" using 1:4 nooutput
diff_10 = STATS_max_y - STATS_min_y

# SHRINK-NO AMR
set title "SHRINK-NO AMR"
set xlabel "time"
set ylabel "M/M_0"
set key bottom right box width 1
set xrange [0:10]
set yrange [0:1.2]
set grid
plot "balances-8-0-0" u 1:2 w l lw 2 lc "red" t "sol mass", \
     "balances-8-0-0" u 1:3 w l lw 2 lc "web-green" t "gas mass", \
     "balances-8-0-0" u 1:4 w l lw 2 lc "web-blue" t sprintf("total mass (Δ=%.3f)", diff_00)

# SWELLING-NO AMR
set title "SWELLING-NO AMR"
plot "balances-8-1-0" u 1:2 w l lw 2 lc "red" t "sol mass", \
     "balances-8-1-0" u 1:3 w l lw 2 lc "web-green" t "gas mass", \
     "balances-8-1-0" u 1:4 w l lw 2 lc "web-blue" t sprintf("total mass (Δ=%.3f)", diff_10)

unset multiplot
~~~

~~~gnuplot Error vs maxlevel
reset
set terminal svg size 1600, 800
set multiplot layout 2,1 title "Mass errors at different grid refinments AXI"

# Calculate differences in total mass at level 5
stats "balances-5-0-0" u (last_y=$4) nooutput
diff_500 = last_y -1

stats "balances-5-1-0" u (last_y=$4) nooutput
diff_510 = last_y - 1

# Calculate differences in total mass at level 6
stats "balances-6-0-0" u (last_y=$4) nooutput
diff_600 = last_y - 1

stats "balances-6-1-0" u (last_y=$4) nooutput
diff_610 = last_y - 1

# Calculate differences in total mass at level 7
stats "balances-7-0-0" u (last_y=$4) nooutput
diff_700 = last_y - 1

stats "balances-7-1-0" u (last_y=$4) nooutput
diff_710 = last_y - 1

# Calculate differences in total mass at level 8
stats "balances-8-0-0" u (last_y=$4) nooutput
diff_800 = last_y - 1

stats "balances-8-1-0" u (last_y=$4) nooutput
diff_810 = last_y - 1

array x_levels[4] = [2**5, 2**6, 2**7, 2**8]

array y_diff_00[4] = [diff_500, diff_600, diff_700, diff_800]
array y_diff_10[4] = [diff_510, diff_610, diff_710, diff_810]

# create arrays for plot
set title "SHRINK"
set xlabel "maxlevel"
set ylabel "ΔM/M"
set key bottom left box width 1 spacing 1.5
set logscale x 2
set logscale y
set xr [2**4:2**9]
set size square
set grid

plot x_levels u (x_levels[$1]):(y_diff_00[$1]) w p pt 8 ps 2 t "SHRINK",\
      10*x**(-1) lw 2 title "1^{st} order", \
      30*x**(-2) lw 2 title "2^{nd} order"

set title "SWELLING"
plot x_levels u (x_levels[$1]):(y_diff_10[$1]) w p pt 8 ps 2 t "SWELLING",\
      10*x**(-1) lw 2 title "1^{st} order", \
      30*x**(-2) lw 2 title "2^{nd} order"

unset multiplot
~~~

~~~gnuplot Error in solid mass
reset
#set terminal pdfcairo

rhoS = 100.
omega = 10.
analytical = exp(-omega/rhoS*10)

# Calculate differences in solid mass at level 5
stats "balances-5-0-0" u (last_y=$2) nooutput
diff_500 = last_y - analytical

stats "balances-5-1-0" u (last_y=$2) nooutput
diff_510 = last_y - analytical

# Calculate differences in solid mass at level 6
stats "balances-6-0-0" u (last_y=$2) nooutput
diff_600 = last_y - analytical

stats "balances-6-1-0" u (last_y=$2) nooutput
diff_610 = last_y - analytical

# Calculate differences in solid mass at level 7
stats "balances-7-0-0" u (last_y=$2) nooutput
diff_700 = last_y - analytical

stats "balances-7-1-0" u (last_y=$2) nooutput
diff_710 = last_y - analytical

# Calculate differences in solid mass at level 8
stats "balances-8-0-0" u (last_y=$2) nooutput
diff_800 = last_y - analytical

stats "balances-8-1-0" u (last_y=$2) nooutput
diff_810 = last_y - analytical

array x_levels[4] = [2**5, 2**6, 2**7, 2**8]

array y_diff_00[4] = [diff_500, diff_600, diff_700, diff_800]
array y_diff_10[4] = [diff_510, diff_610, diff_710, diff_810]

set title "Solid mass conservation"
set key bottom right box width 1

set xlabel "maxlevel"
set ylabel "error"
set key bottom left box width 1 spacing 1.5
#set logscale x 2
#set logscale y
set xr [2**4:2**9]
set size square
set grid

plot x_levels u (x_levels[$1]):(y_diff_00[$1]) w p pt 8 ps 2 t "SHRINK",\
      x_levels u (x_levels[$1]):(y_diff_10[$1]) w p pt 8 ps 2 t "SWELLING",\
      10*x**(-1) lw 2 title "1^{st} order", \
      30*x**(-2) lw 2 title "2^{nd} order"
~~~

~~~gnuplot Error in gas mass
reset
#set terminal pdfcairo

rhoS = 100.
omega = 10.
analytical = 1-exp(-omega/rhoS*10)

# Calculate differences in gas mass at level 5
stats "balances-5-0-0" u (last_y=$3) nooutput
diff_500 = last_y - analytical

stats "balances-5-1-0" u (last_y=$3) nooutput
diff_510 = last_y - analytical

# Calculate differences in gas mass at level 6
stats "balances-6-0-0" u (last_y=$3) nooutput
diff_600 = last_y - analytical

stats "balances-6-1-0" u (last_y=$3) nooutput
diff_610 = last_y - analytical

# Calculate differences in gas mass at level 7
stats "balances-7-0-0" u (last_y=$3) nooutput
diff_700 = last_y - analytical

stats "balances-7-1-0" u (last_y=$3) nooutput
diff_710 = last_y - analytical

# Calculate differences in gas mass at level 8
stats "balances-8-0-0" u (last_y=$3) nooutput
diff_800 = last_y - analytical

stats "balances-8-1-0" u (last_y=$3) nooutput
diff_810 = last_y - analytical

array x_levels[4] = [2**5, 2**6, 2**7, 2**8]

array y_diff_00[4] = [diff_500, diff_600, diff_700, diff_800]
array y_diff_10[4] = [diff_510, diff_610, diff_710, diff_810]

set title "Gas mass conservation"
set key bottom right box width 1

set xlabel "maxlevel"
set ylabel "Error"
set key bottom left box width 1 spacing 1.5
set logscale x 2
set logscale y
set xr [2**4:2**9]
set size square
set grid


plot x_levels u (x_levels[$1]):(y_diff_00[$1]) w p pt 8 ps 2 t "SHRINK",\
      x_levels u (x_levels[$1]):(y_diff_10[$1]) w p pt 8 ps 2 t "SWELLING",\
      0.5*x**(-1) lw 2 title "1^{st} order", \
      5*x**(-2) lw 2 title "2^{nd} order"
~~~
**/