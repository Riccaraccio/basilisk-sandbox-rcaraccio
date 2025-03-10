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

extern scalar TS;
vector gTS[];
scalar modg[];


enum zeta_types {
  ZETA_SHRINK,
  ZETA_SWELLING,
  ZETA_SMOOTH,
  ZETA_SHARP,
  ZETA_LEVELSET,
  ZETA_REACTION,
  ZETA_GRADIENT
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
        zeta[] = 1.;
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
      
      vof_to_ls (f, levelset, imax=60);
      foreach()
        zeta[] = 1.;

      break;
    }

    case ZETA_REACTION: {
      foreach ()
        o[] = f[] > 1. - F_ERR ? omega[] : 0.;
      double o_max = statsf(o).max;
      foreach () {
        zeta[] = o_max > F_ERR ? omega[] / o_max : 0.;
        zeta[] = clamp(zeta[], 0., 1.);
      }
      break;
    }
    case ZETA_GRADIENT: {
      #ifdef SOLVE_TEMPERATURE
      foreach()
        TS[] = f[] > F_ERR ? TS[] / f[] : 0.;

      gradients({TS}, {gTS});

      foreach() {
        modg[] = 0.;
        foreach_dimension()
            modg[] += sq(gTS.x[]);
        modg[] = sqrt(modg[]);
      }

      double delta = L0/(1<<7);
      foreach() {
        zeta[] = f[] > F_ERR ? tanh(modg[]/TS[]*delta) : 0.;
        // zeta[] = TS[] > 0. ? 1 - 1/((modg[]/TS[])+1) : 0.;
        TS[] = f[] > F_ERR ? TS[]*f[] : 0.;
      }
      #endif
      break;
    } 
    default:
      fprintf (stderr, "Unknown Shrinking model\n");
      return;
  }
}

event init (i=0) {
  f.tracers = list_append (f.tracers, porosity);
  set_zeta (zeta_policy);
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
  face vector us[];
  foreach_face()
    us.x[] = ubf.x[];
  dt = dtnext (timestep (us, dtmax));
}

event cleanup (t = end) {
  free(f.tracers), f.tracers = NULL;
}