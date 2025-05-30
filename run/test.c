#define NO_ADVECTION_DIV 1

scalar fS[];
face vector fsS[];

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "const-prop.h"
#include "two-phase.h"
#include "shrinking.h"
#include "view.h"
#include "balances-interface.h"

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

int main() {
  eps0 = 0.4;
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;
  muG = 1.e-3;
  L0 = 1.5*D0;
  DT = 1e-2;

  rhoS = 100.;
  rhoG = 1.;

  zeta_policy = ZETA_SWELLING;
  // for (maxlevel=8; maxlevel<=8; maxlevel++) {
  //   init_grid(1 << maxlevel);
  //   run();
  // }

  init_grid(1 << maxlevel);
  run();
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))
double solid_mass0, radius0;

event init(i = 0) {
  fraction (f, circle(x, y, 0.5*D0));

  foreach()
    porosity[] = eps0*f[];

  solid_mass0 = 0.;
  foreach()
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f
  

  double gas_mass0 = 0.;
  foreach() {
    if (f[] > F_ERR) {
      gas_mass0 += porosity[]*rhoG*dv(); //(1-f)
    }
  }
}

//update the porosity field
event chemistry(i++) {  
  foreach()
    omega[] = rhoS/10;
  foreach() {
    if (f[] > F_ERR) {
      porosity[] = porosity[]/f[];
      porosity[] += (omega[]*(1.-porosity[])*(1.-zeta[])/rhoS)*dt;
      porosity[] *= f[];
    }
  }
}

event adapt (i++) {
    adapt_wavelet_leave_interface({porosity, u.x, u.y}, {f}, 
      (double[]){1e-2, 1e-2, 1e-2}, maxlevel, minlevel, 1);
}

event log (t+=0.1) {
  fprintf(stderr, "%g\n", t);

  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f
  
  double radius = sqrt (statsf(f).sum/pi)*2;

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g %g\n", t, solid_mass/solid_mass0, radius/(D0/2));
  fflush(fp);
}

// event movie (t += 0.1) {
//   if (maxlevel == 7) {
//     scalar eps[];
//     foreach()
//       eps[] = porosity[] + (1-f[]);

//     view(tx=-0.5, ty=-0.5);
//     clear();
//     cells();
//     squares("eps", min=eps0, max=1);
//     draw_vof("f", lw=2);
//     save ("movie.mp4");
//   }
// }

event stop (t = 10);

/** 
~~~gnuplot Mass profile
reset
set terminal svg size 400,350
set output "mass.svg"
set xlabel "t [s]"
set ylabel "M/M_0"
set key bottom left box width 1
set xrange [0:10]
set yrange [0:1.1]
set grid

analytical(x) = exp(-x/10)
set samples 20

plot  "OutputData-5" u 1:2 w l lw 2 lc "blue" t "level 5", \
      "OutputData-6" u 1:2 w l lw 2 lc "black" t "level 6", \
      "OutputData-7" u 1:2 w l lw 2 lc "dark-green" t "level 7", \
      analytical(x) w p pt 4 lc "black" t "analytical"
~~~

~~~gnuplot Radius profile
reset
set terminal svg size 400,350
set output "radius.svg"
set xlabel "t [s]"
set ylabel "R/R_0"
set key bottom left box width 1
set xrange [0:10]
set yrange [0:1.1]
set grid

analytical(x) = exp(-x/20) #if pure shrinking
#analytical(x) = 1 #if pure swelling
set samples 20 

plot  "OutputData-5" u 1:3 w l lw 2 lc "blue" t "level 5", \
      "OutputData-6" u 1:3 w l lw 2 lc "black" t "level 6", \
      "OutputData-7" u 1:3 w l lw 2 lc "dark-green" t "level 7", \
      analytical(x) w p pt 4 lc "black" t "analytical"
~~~
**/