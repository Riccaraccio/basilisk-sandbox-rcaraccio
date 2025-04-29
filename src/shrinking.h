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
scalar o[];

(const) face vector ccc = unityf;
extern double TG0, TS0;
extern scalar TS, T;
vector gTS[];
scalar modg[];
double o_max = 0.;
scalar ls[];

enum zeta_types {
  ZETA_SHRINK,
  ZETA_SWELLING,
  ZETA_SMOOTH,
  ZETA_SHARP,
  ZETA_LEVELSET,
  ZETA_REACTION
};

enum zeta_types zeta_policy;

void set_zeta (enum zeta_types zeta_policy) {
  #if AXI
    double radius = pow (3.*statsf(f).sum, 1./3.);
  #else
    double radius = sqrt (statsf(f).sum/pi)*2;
  #endif

  switch (zeta_policy) {
    case ZETA_SHRINK:
      foreach()
        zeta[] = 0.5;//temp
      break;

    case ZETA_SWELLING:
      foreach()
        zeta[] = 0.;
      break;

    case ZETA_SMOOTH:
      foreach()
          //zeta[] = 1 / (1 + exp(32*radius - 40*sqrt(sq(x)+sq(y)+sq(z))));
          zeta[] = 1 / (1 + exp(-(sqrt(sq(x)+sq(y)+sq(z))-radius/2)/pow(radius, 4./3.)));
      break;

    case ZETA_SHARP:
      foreach()
        zeta[] = (sqrt(sq(x) + sq(y)) > radius*0.8) ? 1. : 0.;
      break;

    case ZETA_LEVELSET: {
      // vof_to_ls (f, levelset, imax=5);
      // double lmin = statsf(levelset).min;
      // if (fabs(lmin) > F_ERR)
      //   foreach() {
      //     zeta[] = 1 - levelset[]/statsf(levelset).min;
      //     zeta[] = clamp(zeta[], 0., 1.);
      //   }
      
      // vof_to_ls (f, levelset, imax=60);
      break;
    }

    case ZETA_REACTION: {
      foreach()
        o[] = omega[]*f[];

      o_max = statsf(o).max;

      foreach() {
        zeta[] = o_max > F_ERR ? omega[]/o_max : 0.;
        zeta[] = clamp(zeta[], 0., 1.);
      }
      break;
    }
    default:
      fprintf (stderr, "Unknown Shrinking model\n");
      return;
  }
}

event defaults (i=0) {
  f.tracers = list_append (f.tracers, porosity);
  set_zeta (zeta_policy);
  #if AXI
    ccc = new face vector;
    foreach_face()
      ccc.x[] = fm.x[];
  #endif
}

event reset_sources (i++);

event chemistry (i++);

event phasechange (i++) {

  foreach() {
    f[] = clamp(f[], 0.,1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    porosity[] = clamp(porosity[], 0., 1.);
    porosity[] = (f[] > F_ERR) ? porosity[] : 0.;
  }

  set_zeta (zeta_policy);

  double delta = L0/(1 << maxlevel);
  mgpsf = project_sv (ubf, psi, ccc, delta, mgpsf.nrelax);

  foreach() {
    if (f[] > F_ERR) {
#ifdef VARPROP //*(f-ef)
      gasSource[] = -omega[]*(f[]-porosity[])*(1/rhoGv_S[] - 1/rhoSv[])*cm[];
#else
      gasSource[] = -omega[]*(f[]-porosity[])*(1/rhoG - 1/rhoS)*cm[];
#endif
    }
  }
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
  dt = dtnext (timestep (ubf, dtmax));
}

event cleanup (t = end) {
  free(f.tracers), f.tracers = NULL;
}
