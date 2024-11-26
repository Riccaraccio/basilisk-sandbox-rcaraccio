#define NO_ADVECTION_DIV    1
#define FSOLVE_ABSTOL       1.e-3
#define SOLVE_TEMPERATURE   1
//#define EXPLICIT_REACTIONS  1
//#define EXPLICIT_DIFFUSION  1
//#define FIXED_INT_TEMP    1

#include "axi.h" 
#include "navier-stokes/centered-phasechange.h"
#include "prop.h"
#include "two-phase.h"
#include "temperature-vt.h"
#include "multicomponent.h"
#include "shrinking.h"
//#include "darcy.h"
#include "view.h"

u.n[top]      = neumann (0.);
u.t[top]      = neumann (0.);
p[top]        = dirichlet (0.);
psi[top]      = dirichlet (0.);
ubf.t[top]    = neumann (0.);
ubf.n[top]    = neumann (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = dirichlet (0.);
ubf.t[right]  = neumann (0.);
ubf.n[right]  = neumann (0.);

int maxlevel = 8; int minlevel = 2;
double D0 = 2*1.27e-2;
double solid_mass0 = 0.;

int main() {
  lambdaS = 0.1987; lambdaG = 0.076;
  cpS = 1600; cpG = 1167;
  TS0 = 300.; TG0 = 743.;
  rhoS = 850; rhoG = 0.674;
  muG = 3.53e-5;
  eps0 = 0.4;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  L0 = 3.5*D0;

#ifdef EXPLICIT_DIFFUSION
  fprintf(stderr, "Using EXPLICIT_DIFFUSION\n");
  DT = 1e-2;
#else
  fprintf(stderr, "Using IMPLICIT_DIFFUSION\n");
  DT = 1e-2;
#endif

  kinfolder = "biomass/Solid-only-2003";
  init_grid(1 << maxlevel);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 0.4807;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")] = 0.2611;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")] = 0.1325;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")] = 0.0957;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")] = 0.0214;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]  = 0.0086;

  foreach()
    porosity[] = eps0*f[];

  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  zeta_policy = ZETA_SWELLING;

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
  scalar inert = YGList[OpenSMOKE_IndexOfSpecies ("N2")];

  inert[top] = dirichlet (1.);
  inert[right] = dirichlet (1.);
}

event output (t+=1) {
  fprintf (stderr, "%g\n", t);

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();

  //save temperature profile
  double Tcore  = interpolate (T, 0., 0.);
  double Tr2    = interpolate (T, D0/4, 0.);
  double Tsurf  = interpolate (T, D0/2, 0.);

  fprintf (fp, "%g %g %g %g %g\n", t, solid_mass/solid_mass0, Tcore, Tr2, Tsurf);
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

event adapt (i++) {
  adapt_wavelet_leave_interface ({T, u.x, u.y}, {f},
     (double[]){1.e0, 1.e-1, 1.e-1}, maxlevel, minlevel, 1);
}

event movie(t+=1) {
  clear();
  box();
  view (ty=-0.5, tx=-0.5);
  squares("T", max=TG0,  min=TS0);
  draw_vof("f");
  //cells();
  save ("temperature.mp4");

  clear();
  box();
  view (ty=-0.5, tx=-0.5);
  squares("C6H10O5", max=1,  min=0);
  draw_vof("f");
  save ("LVG.mp4");
}

#if DUMP
int count = 0;
event snapshots (t += 1) {
  // we keep overwriting the last two snapshots
  if (count == 0) {
    dump ("snapshot-0");
    count++;
  } else {
    dump ("snapshot-1");
    count = 0;
  }
}
#endif

#if TRACE > 1
  event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif

event stop (t = 800);
