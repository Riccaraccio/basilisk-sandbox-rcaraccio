#define SOLVE_TEMPERATURE 1
#define NO_ADVECTION_DIV 1
#define TURN_OFF_HEAT_OF_REACTION 1

#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "balances.h"

u.n[top] = neumann(0.);
u.t[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);
psi[top] = dirichlet(0.);
TG[top] = dirichlet(TG0);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
psi[right] = dirichlet(0.);
TG[right] = dirichlet(TG0);

double D0 = 1.e-2, tend = 10.;

int maxlevel = 7, minlevel = 2;
int main() {
  TS0 = 700.; TG0 = 700.;
  rhoS = 1000;
  eps0 = 0.4;
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;
  DT = 1e-2;
  lambdaSmodel = L_LU;
  emissivity = emissivity_lu;
  zeta_policy = ZETA_SHRINK;
  shift_prod = true;
  kinfolder = "biomass/dummy-solid-gas";
  L0 = D0*5;
  init_grid (1 << maxlevel);
  run();
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))
event init (i = 0) {
  fraction (f, circle(x, y, 0.5*D0));
  foreach()
    porosity[] = eps0*f[];
  
  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 0.9; // 93.5% biomass
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]     = 0.1; // 0.6% ash

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      YG[top] = dirichlet (0.765);
      YG[right] = dirichlet (0.765);
    } else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      YG[top] = dirichlet (0.235);
      YG[right] = dirichlet (0.235);
    } else {
      YG[top] = dirichlet (0.);
      YG[right] = dirichlet (0.);
    }
  }

  divq_rad = opensmoke_optically_thin;
}

event adapt (i++) {
  scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];

  adapt_wavelet_leave_interface ({T, u.x, u.y, inert}, {f},
    (double[]){1.e-2, 1.e-2, 1.e-2, 1e-2}, maxlevel, minlevel, 1);
}

event stop (t = 10) {
#ifdef BALANCES
  write_balances();
#endif
}


