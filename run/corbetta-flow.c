#define NO_ADVECTION_DIV 1
#define NO_EXPANSION 1
#define SOLVE_TEMPERATURE 1
#define CONST_DIFF 2.05e-5
#define FSOLVE_ABSTOL 1.e-3
#define RADIATION_INTERFACE 1

// #include "temperature-profile.h"
#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
#include "view.h"

double Uin = 16.4; //inlet velocity, 300 NL/min over 8.04e-4 m^2 area

u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
psi[left]    = neumann (0.);

u.n[top]      = dirichlet (0.);
u.t[top]      = dirichlet (0.);
p[top]        = neumann (0.);
psi[top]      = neumann (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = dirichlet (0.);

int maxlevel = 7; int minlevel = 2;
double D0 = 2.54e-2;
double solid_mass0 = 0.;

int main() {
  lambdaS = 0.1987;
  TS0 = 300.; TG0 = 750.;
  rhoS = 850;
  eps0 = 0.4;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  Da = 1e-10;
  zeta_policy = ZETA_SWELLING;
  
  L0 = 6*D0;
  origin (-L0/3, 0);

  DT = 1e-1;

  kinfolder = "biomass/dummy-solid";
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  mask (y > 5.25e-2/2 ? top : none); //5.25 cm diameter of tube
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 0.4807;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")] = 0.2611;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")] = 0.1325;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")] = 0.0957;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")] = 0.0214;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]  = 0.0086;

  foreach()
    porosity[] = eps0*f[];

  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  TG[left] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[left] = dirichlet (1.);
    } else {
      YG[left] = dirichlet (0.);
    }
  }

  foreach()
    u.x[] = f[] > F_ERR ? 0. : Uin; 
}

event output (t+=1) {
  fprintf (stderr, "%g\n", t);

  // char name[80];
  // sprintf(name, "OutputData-%d", maxlevel);
  // static FILE * fp = fopen (name, "w");

  // //log mass profile
  // double solid_mass = 0.;
  // foreach (reduction(+:solid_mass))
  //   solid_mass += (f[]-porosity[])*rhoS*dv();

  // //calculate radius
  // double radius = pow (3.*statsf(f).sum, 1./3.);

  // //save temperature profile
  // double Tcore  = interpolate (T, 0., 0.);
  // double Tr2    = interpolate (T, radius/2., 0.);
  // double Tsurf  = interpolate (T, radius, 0.);

  // #ifdef TEMPERATURE_PROFILE
  //   Tsurf = TemperatureProfile_GetT(t);
  // #endif

  // fprintf (fp, "%g %g %g %g %g\n", t, solid_mass/solid_mass0, Tcore, Tr2, Tsurf);
  // fflush(fp);
}

event adapt (i>10; i++) {
  adapt_wavelet_leave_interface ({T, u.x, u.y}, {f},
    (double[]){1.e-1, 1.e-1, 1.e-1}, maxlevel, minlevel, 1);
}

// event movie(t+=1) {
//   clear();
//   box();
//   view (ty=-0.5, tx=-0.5);
//   squares("T", max=TG0,  min=TS0);
//   draw_vof("f");
//   //cells();
//   save ("temperature.mp4");
//
//   clear();
//   box();
//   view (ty=-0.5, tx=-0.5);
//   squares("C6H10O5", max=1,  min=0);
//   draw_vof("f");
//   save ("LVG.mp4");
// }

// event movie(t+=5) {
//   clear();
//   box();
//   view (ty=-0.5, width = 1400.);
//   draw_vof("f", lw=2);
//   squares ("T", min=300, max=800, linear=true);
//   mirror ({1.,0.}) {
//     draw_vof ("f", lw=2);
//     squares ("C6H10O5_G+C6H10O5_S", min=0., max=0.3, linear=true);
//     // vectors ("u", scale=1);
//  }
//  save ("movie.mp4");

  // clear ();
  // box ();
  // view (ty = -0.5, width = 1400.);
  // draw_vof ("f", lw = 2);
  // scalar epsi[];
  // foreach()
  //   epsi[] = f[]>F_ERR ? porosity[]/f[] : 1.;
  // squares ("epsi", min = 0., max = 1., linear = true);
  // mirror ({1., 0.}) {
  //   draw_vof ("f", lw = 2);
  //   cells();
  //   vectors ("u", scale = 1e-2);
  // }
  // save("movie2.mp4");
// }

// #if DUMP
// int count = 0;
// event snapshots (t += 1) {
//   // we keep overwriting the last two snapshots
//   if (count == 0) {
//     dump ("snapshot-0");
//     count++;
//   } else {
//     dump ("snapshot-1");
//     count = 0;
//   }
// }
// #endif

// #if TRACE > 1
//   event profiling (i += 20) {
//   static FILE * fp = fopen ("profiling", "w");
//   trace_print (fp, 1);
// }
// #endif

event stop (t = 1000);

/** 
~~~gnuplot temperature profiles
reset
set xlabel "t [s]"
set ylabel "temperature [K]"
set key bottom right box width 1
set xrange [0:1000]
set yrange [300:800]
set grid

plot  "OutputData-6" u 1:3 w l lw 2 lc "red" t "Core", \
      "OutputData-6" u 1:4 w l lw 2 lc "web-green" t "R/2", \
      "OutputData-6" u 1:5 w l lw 2 lc "web-blue" t "Surface", \
      "data/corbetta-core.txt"  u 1:2 w p pt 7 lc "red" t "Corbetta core", \
      "data/corbetta-r2.txt"    u 1:2 w p pt 7 lc "web-green" t "Corbetta R/2", \
      "data/corbetta-surf.txt"  u 1:2 w p pt 7 lc "web-blue" t "Corbetta surface"
~~~

~~~gnuplot temperature profiles
reset
set xlabel "t [s]"
set ylabel "temperature [K]"
set key bottom right box width 1
set xrange [0:1000]
set yrange [300:800]
set grid

plot  "OutputData-7" u 1:3 w l lw 2 lc "red" t "Core", \
      "OutputData-7" u 1:4 w l lw 2 lc "web-green" t "R/2", \
      "OutputData-7" u 1:5 w l lw 2 lc "web-blue" t "Surface", \
      "data/biosmoke-core.txt"  u 1:2 w l lw 2 dt 2 lc "red" t "BioSMOKE core", \
      "data/biosmoke-r2.txt"    u 1:2 w l lw 2 dt 2 lc "web-green" t "BioSMOKE R/2", \
      "data/biosmoke-surf.txt"  u 1:2 w l lw 2 dt 2 lc "web-blue" t "BioSMOKE surface"

~~~
~~~gnuplot temperature profiles
reset
set xlabel "t [s]"
set ylabel "temperature [K]"
set key bottom right box width 1
set xrange [0:1000]
set yrange [300:800]
set grid

plot  "OutputData-7" u 1:3 w l lw 2 lc "red" t "Core", \
      "OutputData-7" u 1:4 w l lw 2 lc "web-green" t "R/2", \
      "OutputData-7" u 1:5 w l lw 2 lc "web-blue" t "Surface", \
      "data/corbetta-core.txt"  u 1:2 w p pt 7 lc "red" t "Corbetta core", \
      "data/corbetta-r2.txt"    u 1:2 w p pt 7 lc "web-green" t "Corbetta R/2", \
      "data/corbetta-surf.txt"  u 1:2 w p pt 7 lc "web-blue" t "Corbetta surface", \
      "data/biosmoke-core.txt"  u 1:2 w l lw 2 dt 2 lc "red" t "BioSMOKE core", \
      "data/biosmoke-r2.txt"    u 1:2 w l lw 2 dt 2 lc "web-green" t "BioSMOKE R/2", \
      "data/biosmoke-surf.txt"  u 1:2 w l lw 2 dt 2 lc "web-blue" t "BioSMOKE surface"
~~~
~~~gnuplot temperature profiles
reset
set xlabel "t [s]"
set ylabel "temperature [K]"
set key bottom right box width 1
set xrange [0:1000]
set yrange [300:800]
set grid

plot  "OutputData-7" u 1:3 w l lw 2 lc "red" t "Core", \
      "OutputData-7" u 1:4 w l lw 2 lc "web-green" t "R/2", \
      "OutputData-7" u 1:5 w l lw 2 lc "web-blue" t "Surface", \
      "data/corbetta-core.txt"  u 1:2 w p pt 7 lc "red" t "Corbetta core", \
      "data/corbetta-r2.txt"    u 1:2 w p pt 7 lc "web-green" t "Corbetta R/2", \
      "data/corbetta-surf.txt"  u 1:2 w p pt 7 lc "web-blue" t "Corbetta surface", \
      "data/temperature.ris"    u 1:2 w l lw 2 dt 2 lc "red" t "BioSMOKE core", \
      "data/temperature.ris"    u 1:3 w l lw 2 dt 2 lc "web-green" t "BioSMOKE R/2", \
      "data/temperature.ris"  u 1:11 w l lw 2 dt 2 lc "web-blue" t "BioSMOKE surface"
~~~
**/