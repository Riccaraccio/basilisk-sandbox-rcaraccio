#include "vofToLs.h"
#include "velocity-potential.h"
#include "common-evaporation.h"

double eps0 = 0.5;
double rhoS = 100.;
double rhoG = 1.;
double muG = 1e-5;

extern scalar omega;
scalar porosity[];
scalar zeta[];
scalar levelset[];

typedef enum {
  ZETA_SHRINK = 0,
  ZETA_SWELLING,
  ZETA_SMOOTH,
  ZETA_SHARP,
  ZETA_LEVELSET
} zeta_types;

zeta_types zeta_policy = ZETA_SHRINK;

void set_zeta (zeta_types zeta_policy) {
  #if AXI
    double radius = pow (3.*statsf(f).sum, 1./3.);
  #else
    double radius = sqrt (statsf(f).sum/pi)*2;
  #endif

  switch (zeta_policy) {
    case 0: // ZETA_SHRINK
      foreach()
        zeta[] = 1.;
        break;
    case 1: // ZETA_SWELLING
      foreach()
        zeta[] = 0.;
      break;
    case 2: // ZETA_SMOOTH
      foreach()
          //zeta[] = 1 / (1 + exp(32*radius - 40*sqrt(sq(x)+sq(y)+sq(z))));
          zeta[] = 1 / (1 + exp(-(sqrt(sq(x)+sq(y)+sq(z))-radius/2)/pow(radius, 4./3.)));
      break;
    case 3: // ZETA_SHARP
      foreach()
        zeta[] = (sqrt(sq(x) + sq(y)) > radius*0.8) ? 1. : 0.;
      break;
    case 4: // ZETA_LEVELSET
      vof_to_ls (f, levelset, imax=10); //Solution is sensible to imax value
      foreach()
        zeta[] = levelset[] > 0.8*statsf(levelset).min ? 1. : 0.;
      break;
  }
}

event defaults (i=0) {
  set_zeta (zeta_policy);
}

event init (i=0) {
  f.tracers = list_append (f.tracers, porosity);
}

event reset_sources (i++);

event chemistry (i++);

event phasechange (i++) {

  foreach() {
    f[] = clamp(f[], 0.,1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    porosity[] = clamp(porosity[], 0., 1.);
    porosity[] = (porosity[] > F_ERR) ? porosity[] : 0.;
  }

  set_zeta (zeta_policy);

  mgpsf = project_sv (ubf, psi, alpha, dt, mgpsf.nrelax);

  foreach()
    gasSource[] = -omega[]*(f[]-porosity[])*(1/rhoG - 1/rhoS)*cm[]; // gas production, *(f-ef)
}

face vector ufsave[];
event vof(i++) {
  foreach_face()
    ufsave.x[] = uf.x[];

  foreach_face()
    uf.x[] = ubf.x[];
}

event tracer_diffusion (i++) {
  foreach_face()
    uf.x[] = ufsave.x[];
}

event stability (i++, last) {
  face vector us[];
  foreach_face()
    us.x[] = ubf.x[];
  dt = dtnext (timestep (us, dtmax));
}

event cleanup (t = end) {
  free(f.tracers), f.tracers = NULL;
}