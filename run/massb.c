#define NO_ADVECTION_DIV 1

#include "axi.h"
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
  DT = 1.;

for (maxlevel=5; maxlevel<=7; maxlevel++) {
  for (zetamodel=0; zetamodel <= 2; zetamodel++) {
    for (amr=0; amr <=1; amr++) {
        fprintf(stderr, "Running maxlevel %d zetamodel %d, amr %d\n",
          maxlevel, zetamodel, amr);
        init_grid(1 << maxlevel);
        run();
    }
  }
}

  // zetamodel = 0;
  // amr = 1;
  // init_grid(1 << maxlevel);
  // run();
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

  //specify the output file name
  sprintf(mb.name, "balances-%d-%d-%d", maxlevel, zetamodel, amr);
}

// update the porosity field
event phasechange(i++) {
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

event stop (t = 10);

/** 
~~~gnuplot Mass conservation profiles at level 7
reset
set terminal svg size 1600, 1200
set multiplot layout 3,2 title "Mass conservation profiles, LEVEL 7 AXI" \
  margins 0.1, 0.9, 0.1, 0.9 spacing 0.1, 0.1

# Calculate differences in total mass
stats "balances-7-0-0" using 1:4 nooutput
diff_00 = STATS_max_y - STATS_min_y

stats "balances-7-0-1" using 1:4 nooutput
diff_01 = STATS_max_y - STATS_min_y

stats "balances-7-1-0" using 1:4 nooutput
diff_10 = STATS_max_y - STATS_min_y

stats "balances-7-1-1" using 1:4 nooutput
diff_11 = STATS_max_y - STATS_min_y

stats "balances-7-2-0" using 1:4 nooutput
diff_20 = STATS_max_y - STATS_min_y

stats "balances-7-2-1" using 1:4 nooutput
diff_21 = STATS_max_y - STATS_min_y

# SHRINK-NO AMR
set title "SHRINK-NO AMR"
set xlabel "time"
set ylabel "M/M_0"
set key bottom right box width 1
set xrange [0:10]
set yrange [0:1.2]
set grid
plot "balances-7-0-0" u 1:2 w l lw 2 lc "red" t "sol mass", \
     "balances-7-0-0" u 1:3 w l lw 2 lc "web-green" t "gas mass", \
     "balances-7-0-0" u 1:4 w l lw 2 lc "web-blue" t sprintf("total mass (Δ=%.3f)", diff_00)

# SHRINK-AMR
set title "SHRINK-AMR"
plot "balances-7-0-1" u 1:2 w l lw 2 lc "red" t "sol mass", \
     "balances-7-0-1" u 1:3 w l lw 2 lc "web-green" t "gas mass", \
     "balances-7-0-1" u 1:4 w l lw 2 lc "web-blue" t sprintf("total mass (Δ=%.3f)", diff_01)

# SWELLING-NO AMR
set title "SWELLING-NO AMR"
plot "balances-7-1-0" u 1:2 w l lw 2 lc "red" t "sol mass", \
     "balances-7-1-0" u 1:3 w l lw 2 lc "web-green" t "gas mass", \
     "balances-7-1-0" u 1:4 w l lw 2 lc "web-blue" t sprintf("total mass (Δ=%.3f)", diff_10)

# SWELLING-AMR
set title "SWELLING-AMR"
plot "balances-7-1-1" u 1:2 w l lw 2 lc "red" t "sol mass", \
     "balances-7-1-1" u 1:3 w l lw 2 lc "web-green" t "gas mass", \
     "balances-7-1-1" u 1:4 w l lw 2 lc "web-blue" t sprintf("total mass (Δ=%.3f)", diff_11)

# SMOOTH-NO AMR
set title "SMOOTH-NO AMR"
plot "balances-7-2-0" u 1:2 w l lw 2 lc "red" t "sol mass", \
     "balances-7-2-0" u 1:3 w l lw 2 lc "web-green" t "gas mass", \
     "balances-7-2-0" u 1:4 w l lw 2 lc "web-blue" t sprintf("total mass (Δ=%.3f)", diff_20)

# SMOOTH-AMR
set title "SMOOTH-AMR"
plot "balances-7-2-1" u 1:2 w l lw 2 lc "red" t "sol mass", \
     "balances-7-2-1" u 1:3 w l lw 2 lc "web-green" t "gas mass", \
     "balances-7-2-1" u 1:4 w l lw 2 lc "web-blue" t sprintf("total mass (Δ=%.3f)", diff_21)

unset multiplot
~~~

~~~gnuplot Error vs maxlevel
reset
set terminal svg size 1600, 800
set multiplot layout 1,3 title "Mass errors at different grid refinments AXI"


# Calculate differences in total mass at level 5
stats "balances-5-0-0" using 1:4 nooutput
diff_500 = STATS_max_y - STATS_min_y

stats "balances-5-0-1" using 1:4 nooutput
diff_501 = STATS_max_y - STATS_min_y

stats "balances-5-1-0" using 1:4 nooutput
diff_510 = STATS_max_y - STATS_min_y

stats "balances-5-1-1" using 1:4 nooutput
diff_511 = STATS_max_y - STATS_min_y

stats "balances-5-2-0" using 1:4 nooutput
diff_520 = STATS_max_y - STATS_min_y

stats "balances-5-2-1" using 1:4 nooutput
diff_521 = STATS_max_y - STATS_min_y


# Calculate differences in total mass at level 6
stats "balances-6-0-0" using 1:4 nooutput
diff_600 = STATS_max_y - STATS_min_y

stats "balances-6-0-1" using 1:4 nooutput
diff_601 = STATS_max_y - STATS_min_y

stats "balances-6-1-0" using 1:4 nooutput
diff_610 = STATS_max_y - STATS_min_y

stats "balances-6-1-1" using 1:4 nooutput
diff_611 = STATS_max_y - STATS_min_y

stats "balances-6-2-0" using 1:4 nooutput
diff_620 = STATS_max_y - STATS_min_y

stats "balances-6-2-1" using 1:4 nooutput
diff_621 = STATS_max_y - STATS_min_y


# Calculate differences in total mass at level 7
stats "balances-7-0-0" using 1:4 nooutput
diff_700 = STATS_max_y - STATS_min_y

stats "balances-7-0-1" using 1:4 nooutput
diff_701 = STATS_max_y - STATS_min_y

stats "balances-7-1-0" using 1:4 nooutput
diff_710 = STATS_max_y - STATS_min_y

stats "balances-7-1-1" using 1:4 nooutput
diff_711 = STATS_max_y - STATS_min_y

stats "balances-7-2-0" using 1:4 nooutput
diff_720 = STATS_max_y - STATS_min_y

stats "balances-7-2-1" using 1:4 nooutput
diff_721 = STATS_max_y - STATS_min_y

array x_levels[3] = [2**5, 2**6, 2**7]

array y_diff_00[3] = [diff_500, diff_600, diff_700]
array y_diff_01[3] = [diff_501, diff_601, diff_701]
array y_diff_10[3] = [diff_510, diff_610, diff_710]
array y_diff_11[3] = [diff_511, diff_611, diff_711]
array y_diff_20[3] = [diff_520, diff_620, diff_720]
array y_diff_21[3] = [diff_521, diff_621, diff_721]

# create arrays for plot
set title "SHRINK"
set xlabel "maxlevel"
set ylabel "ΔM/M"
set key bottom left box width 1 spacing 1.5
set logscale x 2
set logscale y
set xr [2**4:2**8]
set size square
set grid

plot x_levels u (x_levels[$1]):(y_diff_00[$1]) w p pt 8 ps 2 t "SHRINK",\
      10*x**(-1) lw 2 title "1^{st} order", \
      30*x**(-2) lw 2 title "2^{nd} order"

set title "SWELLING"
plot x_levels u (x_levels[$1]):(y_diff_10[$1]) w p pt 8 ps 2 t "SWELLING",\
      10*x**(-1) lw 2 title "1^{st} order", \
      30*x**(-2) lw 2 title "2^{nd} order"

set title "SMOOTH"
plot x_levels u (x_levels[$1]):(y_diff_20[$1]) w p pt 8 ps 2 t "SMOOTH",\
      10*x**(-1) lw 2 title "1^{st} order", \
      30*x**(-2) lw 2 title "2^{nd} order"
unset multiplot
~~~
**/