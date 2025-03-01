#define NO_ADVECTION_DIV 1
#define NO_EXPANSION 1
#define SOLVE_TEMPERATURE 1
#define CONST_DIFF 2.05e-5
#define FSOLVE_ABSTOL 1.e-3
// #define RADIATION_INTERFACE 1

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
// #include "darcy.h"
#include "view.h"

u.n[top]      = neumann (0.);
u.t[top]      = neumann (0.);
p[top]        = dirichlet (0.);
psi[top]      = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = dirichlet (0.);

u.n[left]    = neumann (0.);
u.t[left]    = neumann (0.);
p[left]      = dirichlet (0.);
psi[left]    = dirichlet (0.);

int maxlevel = 7; int minlevel = 2;
double D0 = 2e-2; //2cm
double H0 = 3e-2; //3cm
double tend = 800.; //800s

int main() {
  eps0 = 0.4;
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;
  // DT = 1e-2;

  lambdaS = 0.2;
  TS0 = 300.; TG0 = 723.;
  rhoS = 1200.;

  L0 = 6.1415e-2;
  DT = 1e-2;
  origin(-L0/2, 0);

  //0: SHRINK, 
  //1: SWELLING, 
  //2: SMOOTH, 
  //3: SHARP, 
  //4: LEVELSET
  zeta_policy = ZETA_SHRINK;
  kinfolder = "biomass/Solid-only-2407";
  init_grid(1 << maxlevel);
  run();
}

#define rect(x,y)(fabs(x) < 0.5*H0 && fabs(y) < 0.5*D0)

event init(i=0) {
  fraction (f, rect(x, y));
  mask (y > 0.5*L0 ? top : none);
  
  foreach()
    porosity[] = eps0*f[];
  
  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 0.4169;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("XYHW")] = 0.3147;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGO")] = 0.1039;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGH")] = 0.0595;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("LIGC")] = 0.0005;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TANN")] = 0.0616;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("TGL")]  = 0.0349;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]  = 0.0080;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")]= 0.0000;

  #ifdef SOLVE_TEMPERATURE
    TG[top] = dirichlet (TG0);
    TG[right] = dirichlet (TG0);
    TG[left] = dirichlet (TG0);
  #endif

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[top] = dirichlet (1.);
      YG[right] = dirichlet (1.);
      YG[left] = dirichlet (1.);
    } else {
      YG[top] = dirichlet (0.);
      YG[right] = dirichlet (0.);
      YG[left] = dirichlet (0.);
    }
  }
}

event adapt (i++) {
  adapt_wavelet_leave_interface ({T, u.x, u.y}, {f},
    (double[]){1.e0, 1.e-1, 1.e-1}, maxlevel, minlevel, 1);
}

event output (t += 0.1) {
  fprintf(stderr, "%g\n", t);
}

event movie (t += 0.1) {
  clear();
  view (tx = -0.5*L0, ty = 0);
  draw_vof ("f");
  squares("p", spread = -1);
  mirror ({0, 1}) {
    draw_vof("f");
    squares("(u.x^2 + u.y^2)^0.5", spread = -1);
  }
  save ("movie.mp4");
}

event stop (t = tend);