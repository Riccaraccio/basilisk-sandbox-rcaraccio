#define NO_ADVECTION_DIV 1
#define FSOLVE_ABSTOL 1.e-3

//#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "var-prop.h"
#include "two-phase.h"
//#include "temperature-v.h"
#include "multicomponent_r.h"
#include "shrinking.h"
//#include "darcy.h"
#include "view.h"

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
psi[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
psi[right] = dirichlet (0.);

int maxlevel = 8; int minlevel = 2;
double D0 = 1;
double solid_mass0 = 0.;

//scalar omega[];

int main() {
  //lambdaS = 0.124069; lambdaG = 0.0295641;
  //cpS = 2244.92; cpG = 1041.52;
  //TS0 = 300.; TG0 = 1000.;
  rhoS = 1000; rhoG = 1.;
  muG = 2.02391e-5;

  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;
  L0 = 6*D0;
  DT = 0.1;

  kinfolder = "biomass/Solid-only-2003";
  init_grid(1 << maxlevel);
  run();

}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))
scalar TS[];
scalar T[];

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 1.;

  foreach()
    porosity[] = eps0*f[];

  foreach() {
    T[] = 650.;
    TS[] = 650.*f[];
  }
  zeta_policy = ZETA_SWELLING;

  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS; //Note: (1-e) = (1-ef)!= (1-e)f
}

event output (t+=1) {
  fprintf (stderr, "%g\n", t);

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double solid_mass = 0.;
  foreach (reduction(+:solid_mass)) {
    solid_mass += (f[]-porosity[])*rhoS;
  }

  fprintf (fp, "%g %g\n", t, solid_mass/solid_mass0);
  fflush(fp);
}

scalar totG[], totS[];
event end_timestep (i++) {
  foreach() {
    totG[] = 0.; 
    for (int ii=0; ii<NGS; ii++) {
      scalar YG = YGList[ii];
      totG[] += YG[];
    }

    totS[] = 0.; 
    for (int ii=0; ii<NSS; ii++) {
      scalar YS = YSList[ii];
      totS[] += YS[];
   }
  }
}

event bcs (i=0) {
  // TG[top] = dirichlet (TG0);
  // TG[right] = dirichlet (TG0);
  scalar inert = YGList[OpenSMOKE_IndexOfSpecies ("N2")];

  inert[top] = dirichlet (1.);
  inert[right] = dirichlet (1.);
}

//event phasechange (i++){
//  foreach()
//    omega[] = 0.; //fixed interface
//}

event adapt (i++) {
  scalar inert = YGList[OpenSMOKE_IndexOfSpecies ("N2")];
  scalar fuel = YSList[OpenSMOKE_IndexOfSolidSpecies ("CELL")];
  adapt_wavelet_leave_interface ({fuel,inert, T, u.x, u.y, ubf.x, ubf.y}, {f},
     (double[]){1.e-3, 1e-2,1.e0,1.e-1,1.e-1,1.e-1,1.e-1}, maxlevel, minlevel, 1);
}

event movie(t+=0.1) {
  clear();
  box();
  view (ty=-0.5, tx=-0.5);
  squares("C6H10O5");
  cells();
  save ("movie.mp4");
}

event stop (t = 10) {
  return 1;
}
