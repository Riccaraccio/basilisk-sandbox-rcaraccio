#define NO_ADVECTION_DIV 1

// #include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "constant-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "balances.h"
#include "view.h"

u.n[top] = neumann(0.);
u.t[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);
psi[top] = dirichlet(0.);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
psi[right] = dirichlet(0.);

int maxlevel = 7, minlevel = 4;
double D0 = 1;
scalar omega[];

scalar fS[];
face vector fsS[];


int main() {
  eps0 = 0.4;
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;
  muG = 1.e-3;
  L0 = 1.5*D0;
  DT = 1e-3;

  rhoS = 100.;
  rhoG = 1.;

  zeta_policy = ZETA_SWELLING;
  maxlevel = 7;
  init_grid(1 << maxlevel); 
  run();
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

event init(i = 0) {
  fraction (f, circle(x, y, 0.5*D0));

  foreach()
    porosity[] = eps0*f[];

#ifdef BALANCES
  //specify the output file name
  sprintf(mb.name, "balances-%d-%d", maxlevel, zeta_policy);
#endif

  char name[80];
  sprintf(name, "OutputFields-init");
  static FILE *fs = fopen(name, "w");
  output_field({u.x, u.y, omega, zeta, porosity, f}, fs);
  fflush(fs);
  fclose(fs);
}

//update the porosity field
event chemistry(i++) {  
  foreach()
    omega[] = rhoS/10.;
  foreach() {
    if (f[] > F_ERR) {
      porosity[] = porosity[]/f[];
      porosity[] += (omega[]*(1.-porosity[])*(1.-zeta[])/rhoS)*dt;
      porosity[] *= f[];
    }
  }
}

// event adapt (i++) {
//     adapt_wavelet_leave_interface({porosity, u.x, u.y}, {f}, 
//       (double[]){1e-2, 1e-2, 1e-2}, maxlevel, minlevel);
// }

event log (t+=0.1) {
  fprintf(stderr, "%g\n", t);
}

event stop (t = 10) {
  if (maxlevel == 7) {
    // Output fields for visualization
    char name[80];
    sprintf(name, "OutputFields-%d", zeta_policy);
    static FILE *fs = fopen(name, "w");

    output_field({u.x, u.y, omega, zeta, porosity, f}, fs);
    fflush(fs);
    fclose(fs);
  }
  write_balances();
}

// event movie (t += 0.1) {
//   if (maxlevel == 8) {
//     clear();
//     squares("porosity", min = 0, max = eps0);
//     draw_vof("f", lw=2);
//     save ("movie.mp4");
//   }
// }

/** 
#~~~gnuplot Mass conservation profiles at level 8
#reset
#set terminal svg size 1000, 700
#set multiplot layout 2,1 title "Mass conservation profiles, LEVEL 8 AXI" \
#  margins 0.1, 0.9, 0.1, 0.9 spacing 0.1, 0.1
#
## Calculate differences in total mass
#stats "balances-8-0" using 1:4 nooutput
#diff_00 = STATS_max_y - STATS_min_y
#
#stats "balances-8-1" using 1:4 nooutput
#diff_10 = STATS_max_y - STATS_min_y
#
## SHRINK
#set title "SHRINK"
#set xlabel "time"
#set ylabel "M/M_0"
#set key bottom right box width 1
#set xrange [0:10]
#set yrange [0:1.2]
#set grid
#plot "balances-8-0" u 1:2 w l lw 2 lc "red" t "sol mass", \
#     "balances-8-0" u 1:3 w l lw 2 lc "web-green" t "gas mass", \
#     "balances-8-0" u 1:4 w l lw 2 lc "web-blue" t sprintf("total mass (Δ=%.3f)", diff_00)
#
## SWELLING
#set title "SWELLING"
#plot "balances-8-1" u 1:2 w l lw 2 lc "red" t "sol mass", \
#     "balances-8-1" u 1:3 w l lw 2 lc "web-green" t "gas mass", \
#     "balances-8-1" u 1:4 w l lw 2 lc "web-blue" t sprintf("total mass (Δ=%.3f)", diff_10)
#
#unset multiplot
#~~~

~~~gnuplot Error vs maxlevel
reset
set terminal svg size 450, 450
#set terminal epslatex size 3.6, 3.6 color colortext
#set output "mass-conservation-error-maxlevel.tex"

# Initialize arrays to store differences
array diff_x0[3]  # For SHRINK cases
array diff_x1[3]  # For SWELLING cases
array diff_x3[3]  # For SMOOTH cases
analytical = 1.

# Loop through levels 5 to 7
do for [i=0:2] {
    level = i + 5
    
    # Calculate difference for SHRINK case (x-0)
    stats sprintf("balances-%d-0", level) u (last_y=$4) nooutput
    diff_x0[i+1] = abs(last_y - analytical)
    
    # Calculate difference for SWELLING case (x-1)
    stats sprintf("balances-%d-1", level) u (last_y=$4) nooutput
    diff_x1[i+1] = abs(last_y - analytical)

    # Calculate difference for SMOOTH case (x-3)
    stats sprintf("balances-%d-3", level) u (last_y=$4) nooutput
    diff_x3[i+1] = abs(last_y - analytical)
}

array x_levels[3] = [2**5, 2**6, 2**7]

set xlabel "maxlevel"
set ylabel "M/M_0"
set format y "10^{%T}"
set key top right box opaque width 0.8
set logscale x 2
set logscale y
set yr [0.001:1]
set xr [2**4:2**8]
set size square
set grid

# Define fitting functions
A0 = 10
B0 = -2
A1 = 10
B1 = -2
A3 = 10
B3 = -2

f0(x) = A0*x**(B0)
f1(x) = A1*x**(B1)
f3(x) = A3*x**(B3)
f4(x) = 0.7*x**(-1)


set print $SHRINK_DATA
do for [i=1:3] { 
    print sprintf("%g %g", x_levels[i], diff_x0[i])
}
set print

set print $SWELLING_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x1[i])
}
set print

set print $SMOOTH_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x3[i])
}
set print

# Perform curve fitting
fit f0(x) '$SHRINK_DATA' using 1:2 via A0, B0
fit f1(x) '$SWELLING_DATA' using 1:2 via A1, B1
fit f3(x) '$SMOOTH_DATA' using 1:2 via A3, B3

# Plot with equations in legend
plot  '$SHRINK_DATA' u 1:2 w p pt 2 ps 3 lc "dark-green" t "Shrink",\
      '$SWELLING_DATA' u 1:2 w p pt 8 ps 3 lc "black" t "Swelling",\
      '$SMOOTH_DATA' u 1:2 w p pt 4 ps 3 lc "blue" t "Smooth", \
      f4(x) lw 3 lc "red" title "first order"
      #f0(x) lw 3 lc "dark-green" title sprintf("fit: %.2f x^{%.2f}", A0, B0), \
      #f1(x) lw 3 lc "black" title sprintf("fit: %.2f x^{%.2f}", A1, B1), \
      #f3(x) lw 3 lc "blue" title sprintf("fit: %.2f x^{%.2f}", A3, B3), \

~~~

~~~gnuplot Error in solid mass
reset
set terminal svg size 450, 450
#set terminal epslatex size 3.6, 3.6 color colortext
#set output "solid-mass-error-maxlevel.tex"

analytical = exp(-1)

# Initialize arrays to store differences
array diff_x0[3]  # For SHRINK cases
array diff_x1[3]  # For SWELLING cases

# Loop through levels 5 to 7
do for [i=0:2] {
    level = i + 5
    
    # Calculate difference for SHRINK case (x-0)
    stats sprintf("balances-%d-0", level) u (last_y=$2) nooutput
    diff_x0[i+1] = abs(last_y - analytical)
    
    # Calculate difference for SWELLING case (x-1)
    stats sprintf("balances-%d-1", level) u (last_y=$2) nooutput
    diff_x1[i+1] = abs(last_y - analytical)

    # Calculate difference for SMOOTH case (x-3)
    stats sprintf("balances-%d-3", level) u (last_y=$2) nooutput
    diff_x3[i+1] = abs(last_y - analytical)
}

array x_levels[3] = [2**5, 2**6, 2**7]

A0 = 0.01
B0 = -1
A1 = 0.13
B1 = -2
A3 = 0.01
B3 = -1

f0(x) = A0*x**(B0)
f1(x) = A1*x**(B1)
f3(x) = A3*x**(B3)

set print $SHRINK_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x0[i])
}
set print
set print $SWELLING_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x1[i])
}
set print
set print $SMOOTH_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x3[i])
}
set print

fit f0(x) '$SHRINK_DATA' using 1:2 via A0, B0
fit f1(x) '$SWELLING_DATA' using 1:2 via A1, B1
fit f3(x) '$SMOOTH_DATA' using 1:2 via A3, B3

set key bottom right box width 1

set xlabel "maxlevel"
set ylabel "error"
set format y "10^{%T}"
set key top right box opaque
set logscale x 2
set logscale y
#set yr [0.001:1]
set xr [2**4:2**8]
set size square
set grid

plot '$SHRINK_DATA' u 1:2 w p pt 2 ps 3 lc "dark-green" t "Shrink",\
      '$SWELLING_DATA' u 1:2 w p pt 8 ps 3 lc "black" t "Swelling",\
      '$SMOOTH_DATA' u 1:2 w p pt 4 ps 3 lc "blue" t "Smooth", \
      f3(x) lw 3 lc "blue" title sprintf("fit: %.2f x^{%.2f}", A3, B3)
~~~

~~~gnuplot Error in gas mass
reset
set terminal svg size 450, 450
#set terminal epslatex size 3.6, 3.6 color colortext
#set output "gas-mass-error-maxlevel.tex"

analytical = 1-exp(-1)

# Initialize arrays to store differences
array diff_x0[3]  # For SHRINK cases
array diff_x1[3]  # For SWELLING cases
array diff_x3[3]  # For SMOOTH cases

# Loop through levels 5 to 7
do for [i=0:2] {
    level = i + 5
    
    # Calculate difference for SHRINK case (x-0)
    stats sprintf("balances-%d-0", level) u (last_y=$3) nooutput
    diff_x0[i+1] = abs(last_y - analytical)
    
    # Calculate difference for SWELLING case (x-1)
    stats sprintf("balances-%d-1", level) u (last_y=$3) nooutput
    diff_x1[i+1] = abs(last_y - analytical)

    # Calculate difference for SMOOTH case (x-3)
    stats sprintf("balances-%d-3", level) u (last_y=$3) nooutput
    diff_x3[i+1] = abs(last_y - analytical)
}

array x_levels[3] = [2**5, 2**6, 2**7]

set title "Gas mass conservation"
set key top right box opaque

set xlabel "maxlevel"
set ylabel "Error"
set format y "10^{%T}"
set logscale x 2
set logscale y
set yr [0.001:1]
set xr [2**4:2**8]
set size square
set grid

A0 = 0.01
B0 = -1
A1 = 0.13
B1 = -2
A3 = 0.01
B3 = -1

f0(x) = A0*x**(B0)
f1(x) = A1*x**(B1)
f3(x) = A3*x**(B3)

set print $SHRINK_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x0[i])
}
set print
set print $SWELLING_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x1[i])
}
set print
set print $SMOOTH_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x3[i])
}
set print
fit f0(x) '$SHRINK_DATA' using 1:2 via A0, B0
fit f1(x) '$SWELLING_DATA' using 1:2 via A1, B1
fit f3(x) '$SMOOTH_DATA' using 1:2 via A3, B3

plot '$SHRINK_DATA' u 1:2 w p pt 2 ps 3 lc "dark-green" t "Shrink",\
      '$SWELLING_DATA' u 1:2 w p pt 8 ps 3 lc "black" t "Swelling",\
      '$SMOOTH_DATA' u 1:2 w p pt 4 ps 3 lc "blue" t "Smooth", \
      f3(x) lw 3 lc "blue" title sprintf("fit: %.2f x^{%.2f}", A3, B3)

~~~
~~~gnuplot Maps for Shrink
reset
set terminal svg size 450, 450
set pm3d map interpolate 0,0
#load "/root/gnuplot-palettes/ylgn.pal"
load "/root/gnuplot-palettes/ylorrd.pal"

set multiplot

set xr [0:1]
set yr [0:1]
set size square
unset key
unset xtics 
unset ytics
unset colorbox
set cbrange [0:1]
splot "OutputFields-0" u 1:2:7

set contour base
set cntrparam levels discrete 0.5
set cntrparam bspline
set cntrlabel onecolor

unset surface
splot "OutputFields-0" u 1:2:8 lt 3 lc "black" lw 1.5, "OutputFields-init" u 1:2:8 lt 3 dt 2  lc "black" lw 1.5
unset multiplot
~~~

~~~gnuplot Maps for Swelling
reset
set terminal svg size 450, 450
set pm3d map interpolate 0,0
#load "/root/gnuplot-palettes/ylgn.pal"
load "/root/gnuplot-palettes/ylorrd.pal"

set multiplot

set xr [0:1]
set yr [0:1]
set size square
unset key
unset xtics 
unset ytics
unset colorbox
set cbrange [0:1]
splot "OutputFields-1" u 1:2:7

set contour base
set cntrparam levels discrete 0.5
set cntrparam bspline
set cntrlabel onecolor

unset surface
splot "OutputFields-1" u 1:2:8 lt 3 lc "black" lw 1.5
unset multiplot
~~~

~~~gnuplot Maps for Smooth
reset
set terminal svg size 450, 450
set pm3d map interpolate 0,0
#load "/root/gnuplot-palettes/ylgn.pal"
load "/root/gnuplot-palettes/ylorrd.pal"

set multiplot

set xr [0:1]
set yr [0:1]
set size square
unset key
unset xtics 
unset ytics
unset colorbox
set cbrange [0:0.8]
splot "OutputFields-3" u 1:2:7

set contour base
set cntrparam levels discrete 0.5
set cntrparam bspline
set cntrlabel onecolor

unset surface
splot "OutputFields-3" u 1:2:8 lt 3 lc "black" lw 1.5, "OutputFields-init" u 1:2:8 lt 3 dt 2  lc "black" lw 1.5
unset multiplot
~~~

~~~gnuplot Maps for initlal condition
reset
set terminal svg size 450, 450
set pm3d map interpolate 0,0
#load "/root/gnuplot-palettes/ylgn.pal"
load "/root/gnuplot-palettes/ylorrd.pal"

set multiplot

set xr [0:1]
set yr [0:1]
set size square
unset key
unset xtics 
unset ytics
unset colorbox
set cbrange [0:1]
splot "OutputFields-init" u 1:2:7

set contour base
set cntrparam levels discrete 0.5
set cntrparam bspline
set cntrlabel onecolor

unset surface
splot "OutputFields-init" u 1:2:8 lt 3 lc "black" lw 1.5
unset multiplot
~~~
**/