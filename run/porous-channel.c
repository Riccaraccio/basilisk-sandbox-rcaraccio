/**
# Porous channel flow
This case simulates the flow in a 2D channel partially filled with a porous medium.
The porous medium occupies the lower half of the channel, while the upper half 
is free fluid. Velocity at the inlet is uniform and let to develop along 
the channel length. The simulation parameters are chosen to match the results
 reported by [Betchen et al. (2006)](https://doi.org/10.1080/10407780500430967).

## Simulation parameters
Particular attention must be given to the U0 value at the inlet.
As the original paper, we have to set U0 so that Re number is 1 in the free
fluid region when the profile is fully developed. This turns out to be U 
$\approx$ 1 or U0 = 1.17 for the current case. This is verified by computing
the average velocity in the free fluid region at the outlet, as done in the
`stop` event.
*/

int maxlevel = 10;        // Maximum refinement level
double H = 1.;            // Channel height
double U0 = 1.17  ;       // Inflow velocity
double eps0 = 0.7;        // Porosity
double tend = 1.;         // End time
double Re = 1.;           // Reynolds number

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "darcy.h"
#include "view.h"

/**
## Boundary conditions
+ left: inlet only in the free fluid region
+ right: outlet
+ top: wall
+ bottom: wall
*/
u.n[left] = dirichlet (U0*(1. - f[]));
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);
pf[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
pf[right] = dirichlet (0.);

u.n[top] = dirichlet (0.);
u.t[top] = dirichlet (0.);
p[top] = neumann (0.);
pf[top] = neumann (0.);

u.n[bottom] = dirichlet (0.);
u.t[bottom] = dirichlet (0.);
p[bottom] = neumann (0.);
pf[bottom] = neumann (0.);

scalar porosity[], f[];
double rhoG = 1., muG;

int main() {

  /**
  We set the fluid properties to match the desired Re number in the free fluid
  region.
  */
  muG = H*1.*rhoG/Re;
  const face vector muv[] = {muG, muG};
  mu = muv;

  Da = (coord) {1.e-2/sq(H), 1.e-2/sq(H)};

  /**
   Domain in 8H x 2H
  */
  size (8*H);
  dimensions (nx=8, ny=2);
  
  origin (0, -H);
  init_grid (1 << maxlevel);

  run();
}

/**
 We initialize the porous medium in the lower half of the channel.
*/
event init (i = 0) {
  fraction (f, -y);
  foreach() {
    porosity[] = f[]*eps0;
  }
}

/**
 We restrict the maximum CFL number to 0.1 to ensure stability.
*/

const double cfl_max = 0.1;
event stability (i++) {
  if (CFL > cfl_max)
    CFL = cfl_max;
}

/**
## Log event
We log the velocity profile at near the outflow region.
We also compute the average velocity in the free fluid region,
it should be close to 1.
*/

event stop (t = tend) {
  double step = 2*H/(1 << maxlevel);
  double x_interpolate = 0.99*8*H;
  int counter = 0;
  double avg_U = 0.;
  for (double y = -H; y<H; y+=step){
    double Y = y/H;
    double U = interpolate (u.x, x_interpolate, y)/1.;
    fprintf (stderr, "%g %g\n", U, Y);

    if (y > 0) {
      counter++;
      avg_U += U;
    }
  }
  fprintf (stderr, "# Avg u: %g\n", avg_U/counter);
}

/**
## Velocity profile plot
~~~gnuplot
reset
set terminal svg size 400,400
set output "velocity-profile.svg"

set xlabel "u/U"
set ylabel "y/H"
set grid
set size square
unset key

set xrange [0:1.5]
plot "log" u 1:2 w l lw 2 lc "blue", \
     "../../data/porouschannel/velocity-da-02" w p pt 4
~~~
*/
