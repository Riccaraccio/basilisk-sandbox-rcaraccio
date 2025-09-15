#ifndef DARCY
# define DARCY_INDEX 0
#endif

#define POROUS_ADVECTION 1
#define POROUS_MEDIA 1
#define AMR_ACTIVE 1

////INPUTS//////////////
int maxlevel = 10; 
double Re = 20;
double R0 = 0.5;
double U0 = 1.;
double epsi0 = 0.7;
double side_length = 20;
double tend = 2.;
////////////////////////

#ifdef POROUS_MEDIA
#include "navier-stokes/centered-phasechange.h"
#include "fractions.h"
#include "darcy.h"
#include "view.h"
#include "adapt_wavelet_leave_interface.h"

// BOUNDARIES: left and right open, top and bottom no-slip
u.n[left] = dirichlet (U0);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);
pf[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
pf[right] = dirichlet (0.);

scalar porosity[];
double rhoG, muG;

int main() {
  double DaList[] = {1.00E-05, 5.00E-05, 1.00E-04, 5.00E-04, 1.00E-03, 2.50E-03, 5.00E-03};
  muG = R0*2*U0/Re;
  const face vector muv[] = {muG, muG};
  mu = muv;
  rhoG = 1;
  Da = (coord) {DaList[DARCY_INDEX], DaList[DARCY_INDEX]};
  fprintf(stderr, "Da = %g\n", DaList[DARCY_INDEX]);

  L0 = side_length*R0*2;
  origin(-L0/2, 0);
  init_grid (1 << maxlevel);
  run();
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))
scalar f[];
event init (i = 0) {
  mask (y > L0/2 ? top : none);
  fraction (f, circle(x, y, R0));

  foreach()
    porosity[] = f[]*epsi0;
  
  foreach()
    u.x[] = f[] > 0 ? U0/10 : U0; 
}

#if AMR_ACTIVE
event adapt (i++) {
  adapt_wavelet_leave_interface ({u.x, u.y}, {f}, (double[]){1.e-3, 1.e-3}, maxlevel, 2, padding=2);
}
#endif

#else // !POROUS_MEDIA
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "adapt_wavelet_leave_interface.h"

// BOUNDARIES: left and right open, top and bottom no-slip
u.n[left] = dirichlet (U0);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

u.n[embed] = fabs(y) > R0 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > R0 ? neumann(0.) : dirichlet(0.);

face vector muv[];

int main() {
  mu = muv;

  L0 = side_length*R0*2;
  origin(-L0/2, 0);
  init_grid (1 << maxlevel);
  run();
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

event properties (i++){
  foreach_face()
    muv.x[] = fm.x[]*R0*2*U0/Re;
}

event init (i = 0) {
  mask (y > L0/2 ? top : none);
  solid (cs, fs, -circle(x, y, R0));
  foreach()
    u.x[] = cs[] > 1.- 1.e-10 ? U0 : 0;
}

# if   AMR_ACTIVE
event adapt (i++) {
  adapt_wavelet_leave_interface ({u.x, u.y}, {cs}, (double[]){1.e-3, 1.e-3}, maxlevel, 3, padding=1);
}
# endif // AMR_ACTIVE
#endif // !POROUS_MEDIA


scalar avcd[];
// Movie event
event movie (t=tend) {
  scalar omega[];
  vertex scalar stream[];
  foreach()
    omega[] = 0;

  vorticity (u, omega);

  // avcd[bottom] = dirichlet(0.);
  // avcd[top]    = ne(0.);
  // avcd[left]   = neumann(U0);
  // avcd[right]  = neumann(U0);

  avcd[top] = dirichlet(0);
  avcd[bottom] = dirichlet(U0*L0/2);

  poisson(avcd, omega);
  boundary({avcd});
  foreach_vertex()
    stream[] = (avcd[0,-1] + avcd[-1, -1]
              + avcd[] + avcd[-1])/4;

  view (quat = {0.000, 0.000, 0.000, 1.000},
        fov = 30, near = 0.01, far = 1000,
        tx = -0.036, ty = -0.029, tz = -0.182,
        width = 1920, height = 1080);
  draw_vof (c = "f");
  isoline (phi = "avcd", n = 400, min = 9.95, max = 11);
  save ("movie.ppm");

  // view (fov = 6, tx = -0.65, ty = -0.5,
  //       width = 1600, height = 800);
  // box ();
  // #ifdef POROUS_MEDIA
  //   draw_vof ("f", fc = {0.5,0.5,0.5});
  // #else  
  //   draw_vof ("cs", fc = {0.5,0.5,0.5});
  // #endif
  // isoline (phi = "psi", n = 80, min=0.4 , max=0.6);
  // save("streamlines.mp4");

  // clear();

  // view (fov = 6, tx = -0.65, ty = -0.5,
  //       width = 1600, height = 800);
  // box ();
  // #ifdef POROUS_MEDIA
  //   draw_vof ("f", fc = {0.5,0.5,0.5});
  // #else
  //   draw_vof ("cs", fc = {0.5,0.5,0.5});
  // #endif
  // squares (color = "u.x", spread = -1, linear = true);
  // save("velocity.mp4");
} 

//Print the current time on log file
event logfile (t += 1) {
  fprintf(stderr,"%g \n", t);
}

// Output data
void write_profile (void) {
  fprintf(stderr, "Writing wake length profile\n");
  char name[80];
  #ifdef POROUS_MEDIA
    sprintf (name, "WakeLength-Da-%d", DARCY_INDEX);
  #else
    sprintf (name, "WakeLength-Embed");
  #endif
  FILE * fp = fopen (name, "w");

  double step = L0/(1<<maxlevel);
  for (double x = 0; x<L0; x+=step){
    fprintf(fp, "%g %g\n", x, interpolate(u.x, x, 0));
  }

  fprintf (fp, "\n\n");
  fflush(fp);
  fclose(fp);
}

event stop (t = tend) {
  fprintf(stderr, "Reached final time\n");
  write_profile();
}

/**
~~~gnuplot wake length
reset
set terminal svg size 450,400
set output "wake-length.svg"

Wake_Embed = 0.923063153

set grid
set xlabel "Permeability [m^2]"
set ylabel "Normalized Wake length [-]"
set logscale x
set format x "10^{%T}"
set xrange [1e-6:1e-2]
set yrange [0:1.0]
set key bottom left box opaque

# Inline data for Wake_Da and Da
$WakeData << EOD
1.00E-05 0.878739855
5.00E-05 0.875961644
1.00E-04 0.872480282
5.00E-04 0.842615767
1.00E-03 0.799354374
2.50E-03 0.624881385
5.00E-03 0.000000000
EOD

plot Wake_Embed w l dt 2 lw 2 lc "black" t "Embed",\
     $WakeData u 1:2 w lp pt 4 lw 2 lc "dark-green" notitle,\
     "../../data/porouscylinder/yuData" u 1:2 w p pt 6 lc "blue" t "Yu et al. (2011)"

~~~
**/
