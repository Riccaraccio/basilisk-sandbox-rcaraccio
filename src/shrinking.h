#include "vofToLs.h"
scalar mEvap[], * mEvapList = {mEvap};

double eps0 = 0.5;
double rhos = 100.;
double rhog = 1;

extern scalar omega;
extern scalar levelset;
scalar feps[], * f_tracers = NULL;
scalar zeta[];
scalar levelset[];

typedef enum {
  ZETA_SHRINK = 0,
  ZETA_SWELLING,
  ZETA_SMOOTH,
  ZETA_SHARP,
  ZETA_LEVELSET
} zeta_types;

zeta_types zeta_policy;

event defaults (i=0) {
  zeta_policy = ZETA_SHRINK;
}

event init (i=0) {
  foreach()
    feps[] = eps0*f[];
  f_tracers = f.tracers;
  f.tracers = list_append (f.tracers, feps);
}

scalar psi[];
face vector ubf[];
mgstats mgpsf;

trace
mgstats project_sv (face vector ubf, scalar psi,
    double dt = 1.,
    int nrelax = 4)
{
  scalar prod[];
  foreach()
    prod[] = (omega[]*f[]*zeta[]/rhos)/dt;

  mgstats mgp = poisson (psi, prod, tolerance = TOLERANCE/sq(dt),
      nrelax = nrelax);

  foreach_face()
    ubf.x[] = -dt*face_gradient_x (psi, 0);

  return mgp;
}

psi[right] = neumann (neumann_pressure(ghost));
psi[left]  = neumann (- neumann_pressure(0));

#if AXI
ubf.n[bottom] = 0.;
ubf.t[bottom] = dirichlet(0);
psi[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
psi[top]    = neumann (neumann_pressure(ghost));
psi[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
psi[front]  = neumann (neumann_pressure(ghost));
psi[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AX

event phasechange (i++) {

  // remove spurious values
  foreach() {
    f[] = clamp(f[], 0.,1.);
    feps[] = clamp(feps[], 0., 1.);
  }

  // compute radius
  double radius = sqrt(statsf(f).sum/pi)*2; //*2 TEMP

  //compute zeta
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
        zeta[] = 1 / (1 + exp(32*radius - 40*sqrt(sq(x)+sq(y))));
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

  //calculate gas source, epsi evolution and interface regression velocity
  mgpsf = project_sv (ubf, psi, dt, mgpsf.nrelax);

  foreach() {
  //drhodt[] = omega/rhos*f[]*zeta[]; // solid shrinking
  drhodt[] = -omega[]/rhog*f[]; // gas production
  feps[] += (feps[] > F_ERR) ? -dt*omega[]/rhos*feps[]*(1-zeta[]) : 0.; //epsi evolution
  }

}

face vector ufsave[];
event vof(i++) {
  foreach_face()
    ufsave.x[] = uf.x[];

  foreach_face()
    uf.x[] = ubf.x[];
}

event tracer_advection (i++) {
  foreach_face()
    uf.x[] = ufsave.x[];
}

event cleanup (t = end) {
  free (f.tracers), f.tracers = f_tracers;
}
