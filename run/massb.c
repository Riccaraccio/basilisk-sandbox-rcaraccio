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
  muG = 1.e-3;
  L0 = 1.5*D0;
  DT = 1e-2;

  rhoS = 100.;
  rhoG = 1.;

for (maxlevel=5; maxlevel<=7; maxlevel++) {
  for (zetamodel=0; zetamodel <= 1; zetamodel++) {
    for (amr=0; amr <=1; amr++) {
        fprintf(stderr, "Running maxlevel %d zetamodel %d, amr %d\n",
          maxlevel, zetamodel, amr);
        init_grid(1 << maxlevel);
        run();
    }
  }
}

// amr = 0;
// for (maxlevel=5; maxlevel<=7; maxlevel++) {
//   for (zetamodel=0; zetamodel <= 1; zetamodel++) {
//         fprintf(stderr, "Running maxlevel %d zetamodel %d, amr %d\n",
//           maxlevel, zetamodel, amr);
//         init_grid(1 << maxlevel);
//         run();
//     }
//   }

  // maxlevel = 6;
  // amr = 0;
  // zetamodel = 1;
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

#ifdef BALANCES
  //specify the output file name
  sprintf(mb.name, "balances-%d-%d-%d", maxlevel, zetamodel, amr);
#endif
}

//update the porosity field
event chemistry(i++) {  
  foreach()
    omega[] = 10.;
  foreach() {
    if (f[] > F_ERR) {
      porosity[] = porosity[]/f[];
      porosity[] += (omega[]*(1.-porosity[])*(1.-zeta[])/rhoS)*dt;
      porosity[] *= f[];
    }
  }
}

event adapt (i++) {
  if (amr > 0) {
  adapt_wavelet_leave_interface({porosity, u.x, u.y}, {f}, 
    (double[]){1e-2, 1e-3, 1e-3}, maxlevel, minlevel);
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
set multiplot layout 1,2 title "Mass errors at different grid refinments AXI"

# Initialize arrays to store differences
array diff_x00[4]  # For SHRINK cases
array diff_x10[4]  # For SWELLING cases
array diff_x01[4]  # For SHRINK cases AMR
array diff_x11[4]  # For SWELLING cases AMR
analytical = 1.

# Loop through levels 5 to 8
do for [i=0:3] {
    level = i + 5
    
    # Calculate difference for SHRINK case (x-0-0)
    stats sprintf("balances-%d-0-0", level) u (last_y=$4) nooutput
    diff_x00[i+1] = last_y - analytical
    
    # Calculate difference for SWELLING case (x-1-0)
    stats sprintf("balances-%d-1-0", level) u (last_y=$4) nooutput
    diff_x10[i+1] = last_y - analytical
    
    # Calculate difference for SWELLING case (x-0-1)
    stats sprintf("balances-%d-0-1", level) u (last_y=$4) nooutput
    diff_x01[i+1] = last_y - analytical
    
    # Calculate difference for SWELLING case (x-1-1)
    stats sprintf("balances-%d-1-1", level) u (last_y=$4) nooutput
    diff_x11[i+1] = last_y - analytical
}

array x_levels[4] = [2**5, 2**6, 2**7, 2**8]

set title "SHRINK"
set xlabel "maxlevel"
set ylabel "ΔM/M"
set key bottom left box width 1 spacing 1.5
set logscale x 2
set logscale y
set xr [2**4:2**9]
set size square
set grid

plot  x_levels u (x_levels[$1]):(diff_x00[$1]) w p pt 8 ps 2 t "SHRINK",\
      x_levels u (x_levels[$1]):(diff_x01[$1]) w p pt 8 ps 2 t "SHRINK AMR",\
      10*x**(-1) lw 2 title "1^{st} order", \
      30*x**(-2) lw 2 title "2^{nd} order"

set title "SWELLING"
plot  x_levels u (x_levels[$1]):(diff_x10[$1]) w p pt 8 ps 2 t "SWELLING",\
      x_levels u (x_levels[$1]):(diff_x11[$1]) w p pt 8 ps 2 t "SWELLING AMR",\
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

# Initialize arrays to store differences
array diff_x00[4]  # For SHRINK cases
array diff_x10[4]  # For SWELLING cases
array diff_x01[4]  # For SHRINK cases AMR
array diff_x11[4]  # For SWELLING cases AMR

# Loop through levels 5 to 8
do for [i=0:3] {
    level = i + 5
    
    # Calculate difference for SHRINK case (x-0-0)
    stats sprintf("balances-%d-0-0", level) u (last_y=$2) nooutput
    diff_x00[i+1] = abs(last_y - analytical)
    
    # Calculate difference for SWELLING case (x-1-0)
    stats sprintf("balances-%d-1-0", level) u (last_y=$2) nooutput
    diff_x10[i+1] = abs(last_y - analytical)
    
    # Calculate difference for SWELLING case (x-0-1)
    stats sprintf("balances-%d-0-1", level) u (last_y=$2) nooutput
    diff_x01[i+1] = abs(last_y - analytical)
    
    # Calculate difference for SWELLING case (x-1-1)
    stats sprintf("balances-%d-1-1", level) u (last_y=$2) nooutput
    diff_x11[i+1] = abs(last_y - analytical)
}

array x_levels[4] = [2**5, 2**6, 2**7, 2**8]

set title "Solid mass conservation"
set key bottom right box width 1

set xlabel "maxlevel"
set ylabel "error"
set key bottom left box width 1 spacing 1.5
set logscale x 2
set logscale y
set xr [2**4:2**9]
set size square
set grid

plot  x_levels u (x_levels[$1]):(diff_x00[$1]) w p pt 8 ps 2 t "SHRINK",\
      x_levels u (x_levels[$1]):(diff_x01[$1]) w p pt 8 ps 2 t "SHRINK AMR",\
      x_levels u (x_levels[$1]):(diff_x10[$1]) w p pt 8 ps 2 t "SWELLING",\
      x_levels u (x_levels[$1]):(diff_x11[$1]) w p pt 8 ps 2 t "SWELLING AMR",\
      0.01*x**(-1) lw 2 title "1^{st} order", \
      0.13*x**(-2) lw 2 title "2^{nd} order"
~~~

~~~gnuplot Error in gas mass
reset
#set terminal pdfcairo

rhoS = 100.
omega = 10.
analytical = 1-exp(-omega/rhoS*10)

# Initialize arrays to store differences
array diff_x00[4]  # For SHRINK cases
array diff_x10[4]  # For SWELLING cases
array diff_x01[4]  # For SHRINK cases AMR
array diff_x11[4]  # For SWELLING cases AMR

# Define analytical value (assuming it's already defined)
# analytical = your_value

# Loop through levels 5 to 8
do for [i=0:3] {
    level = i + 5
    
    # Calculate difference for SHRINK case (x-0-0)
    stats sprintf("balances-%d-0-0", level) u (last_y=$3) nooutput
    diff_x00[i+1] = last_y - analytical
    
    # Calculate difference for SWELLING case (x-1-0)
    stats sprintf("balances-%d-1-0", level) u (last_y=$3) nooutput
    diff_x10[i+1] = last_y - analytical
    
    # Calculate difference for SWELLING case (x-0-1)
    stats sprintf("balances-%d-0-1", level) u (last_y=$3) nooutput
    diff_x01[i+1] = last_y - analytical
    
    # Calculate difference for SWELLING case (x-1-1)
    stats sprintf("balances-%d-1-1", level) u (last_y=$3) nooutput
    diff_x11[i+1] = last_y - analytical
}

array x_levels[4] = [2**5, 2**6, 2**7, 2**8]

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

plot  x_levels u (x_levels[$1]):(diff_x00[$1]) w p pt 8 ps 2 t "SHRINK",\
      x_levels u (x_levels[$1]):(diff_x01[$1]) w p pt 8 ps 2 t "SHRINK AMR",\
      x_levels u (x_levels[$1]):(diff_x10[$1]) w p pt 8 ps 2 t "SWELLING",\
      x_levels u (x_levels[$1]):(diff_x11[$1]) w p pt 8 ps 2 t "SWELLING AMR",\
      0.5*x**(-1) lw 2 title "1^{st} order", \
      5*x**(-2) lw 2 title "2^{nd} order"
~~~
**/