#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MULTICOMPONENT 1

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "constant-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
//#include "flame.h"

const double Uin = 0.5; //inlet velocity
u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = neumann (0.);

double tend = 15;
int maxlevel = 9, minlevel = 2;
double solid_mass0 = 0.;
double D0 = 3e-3;

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

int main() {

  TS0 = 300.; TG0 = 1473.;
  rhoS = 1000;
  eps0 = 0.4;

  rhoG    = 0.31;      // air ~1100 K, 1 atm
  muG     = 4.5e-5;
  lambdaG = 0.08;      cpG = 1200.;
  lambdaS = 0.2;       cpS = 1500.;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_SWELLING;

  DT = 1e-2;

  kinfolder = "biomass/dummy-solid";
  shift_prod = true;

  L0 = 20*D0;
  origin (-L0/2, 0);

  emissivity = emissivity_constant;

  init_grid(1 << min (maxlevel, 8));
  refine (circle (x, y, 4.*D0) > 0. && level < maxlevel);

  run();
}

double r0;
event init (i= 0) {
  scalar f0[];

  fraction (f0, circle(x, y, 0.5 * D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  //gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  //gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  foreach()
    porosity[] = eps0*f0[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += f0[]*(1. - eps0)*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  TG[left] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      //YG[left] = dirichlet (0.765);
      //YG[top] = dirichlet (0.765);
      YG[left] = dirichlet (1.);
      YG[top] = dirichlet (1.);
    } 
    //else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
    //  YG[left] = dirichlet (0.235);
    //  YG[top] = dirichlet (0.235);
    //}
    else {
      YG[left] = dirichlet (0.);
      YG[top] = dirichlet (0.);
    }
  }

  if (restore (file = "last-snapshot", list = all)) {
    fprintf (stderr, "Restart file found!\n");
    restarted = true;
  } else {
    fprintf (stderr, "No restart file found, starting from scratch!\n");

    foreach() {
      f[] = f0[];
      porosity[] = eps0*f[];
    }
  }
}

// Calculates the H2O-density path-averaged temperature
double T_H2O_weigthed_average (double x_interp, int n_samples = 200, const double length = L0/2) {
  scalar YH2O = YGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
  //scalar XH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];

  double numerator = 0., denominator = 0.;
  coord pos, box[2] = {{x_interp, 0.}, {x_interp, length}}, nn = {1, n_samples};
  foreach_region (pos, box, nn, reduction(+:numerator) reduction(+:denominator)) {
    double yH2O_local = interpolate_linear (point, YH2O, pos.x, pos.y, pos.z);
    numerator   += yH2O_local;
    denominator += yH2O_local / interpolate_linear (point, T, pos.x, pos.y, pos.z);
  }

  if (denominator <= 0) // avoid division by 0
    return TG0;

  return numerator/denominator;
}

event output (t += 0.01) {

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, restarted ? "a" : "w");
  if (fp == NULL) {
    fprintf (stderr, "Error opening OutputData\n");
    exit(1);
  }

  // Interpolate Water vapor mass fraction
  double Tavg[3], sample_points[3] = {D0/2 + 1e-3, D0/2 + 2e-3, D0/2 + 4e-3};

  const double L_flame_exp = 20e-3/2;

  for (int ii = 0; ii < 3; ii++)
      Tavg[ii] = T_H2O_weigthed_average (sample_points[ii]);

  if (i == 0)
    fprintf (fp, "#t(1), Ms/Ms0(2), Tmax(3), Tavg_1mm(4), Tavg_2mm(5), Tavg_4mm(6)\n");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[] - porosity[])*rhoS*dv();

  fprintf (fp, "%g %g %g %g %g %g\n", t, solid_mass/solid_mass0, statsf(T).max, 
                                      Tavg[0], Tavg[1], Tavg[2]);

  fflush(fp);
}

#if TREE
event adapt (i++) {
  //scalar oxidiser = YGList_G[OpenSMOKE_IndexOfSpecies ("O2")];
  scalar oxidiser = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];

  //scalar zdiff[];
  //foreach()
  //  zdiff[] = zmix[] - zsto[];

  //adapt_wavelet_leave_interface ({T, oxidiser, zdiff}, {f},
  //  (double[]){5e0, 1.e-2, 1.e-2}, maxlevel, minlevel, 2);
  
  adapt_wavelet_leave_interface ({T, oxidiser}, {f}, (double[]){5e0, 1.e-2}, maxlevel, minlevel, 2);

  // Unrefine for outflow condition
  unrefine (x > L0*0.4);
}
#endif

//event movie (t += 0.1) {
//  clear();
//  view (theta=0, phi=0, psi=-pi/2., width = 1080, height = 1080);
//  squares ("T", min = 300, max = 2000, spread = -1, linear = true);
//  isoline ("T", val = statsf(T).max);
//  isoline ("zmix - zsto", lw = 1.5, lc = {1., 1., 1.});
//  draw_vof ("f", lw = 1.5);
//  mirror ({0, 1}) {
//    squares ("O2_G + O2_S", min = 0., max = 0.235, spread = -1, linear = true);
//    isoline ("zmix - zsto", lw = 1.5, lc = {1., 0., 0.});
//    draw_vof ("f", lw = 1.5);
//  }
//  save ("movie.mp4");
//}

event dump (t = 1; t += 1) {
  dump("last-snapshot");
}

event stop (t = tend);

/**
~~~gnuplot
~~~
**/
