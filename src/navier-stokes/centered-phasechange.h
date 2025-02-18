#include "poisson.h"
extern double rhoG;

scalar gasSource[];
scalar drhodt[];

trace
mgstats project_sf (face vector uf, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] /= dt*Delta;
  }

  foreach() {
    div[] += gasSource[]/dt;
    div[] += drhodt[]/dt;
  }

  mgstats mgp = poisson (p, div, alpha,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  foreach_face()
    uf.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);

  return mgp;
}

#include "utils.h"
#include "bcg.h"

void advection_div (scalar * tracers, face vector u, double dt,
    scalar * src = NULL)
{
  scalar * psrc = src;
  if (!src)
    for (scalar s in tracers) {
      const scalar zero[] = 0.;
      src = list_append (src, zero);
    }
  assert (list_len (tracers) == list_len (src));

  scalar f, source;
  for (f,source in tracers,src) {
    face vector flux[];
    tracer_fluxes (f, u, flux, dt, source);
#if !EMBED
    foreach() {
#if NO_ADVECTION_DIV
      double fold = f[];
#endif
      foreach_dimension()
#if NO_ADVECTION_DIV
        f[] += dt*(flux.x[] - flux.x[1] + fold*(u.x[1] - u.x[]))/(Delta*cm[]);
#else
        f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
#endif
    }
#else // EMBED
    update_tracer (f, u, flux, dt);
#endif // EMBED
  }

  if (!psrc)
    free (src);
}

event defaults (i = 0) {
  foreach(){
    gasSource[] = 0.;
    drhodt[] = 0.;
  }
}

#define project(...) project_sf(__VA_ARGS__)
#define advection(...) advection_div(__VA_ARGS__)
#include "navier-stokes/centered.h"
#undef advection
#undef project
