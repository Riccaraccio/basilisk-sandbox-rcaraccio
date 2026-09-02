/**
# Centered phase change
This file extends the centered Navier-Stokes solver to account for phase change
in two-phase flows with porous media.
It modifies the projection method to include a gas source term in the continuity equation,
and it provides an advection scheme that can account for porosity variations.
*/

#include "poisson.h"
extern double rhoG;
extern scalar porosity, f;

scalar gas_source[];
scalar drhodt[];

/**
## Filter of the divergence source

`GAS_SOURCE_EXACT` enables two changes to the expansion source. The first one
lives in `chemistry.h`: the chemistry part of `drhodt` becomes the exact step
mean of the expansion, `ln(rho_start/rho_end)/dt`. The second one lives here: a
short diffusion filter acts on the sum `gas_source + drhodt` before the
projection. Compile with `-DGAS_SOURCE_EXACT=1` to select both.

Why a filter: in a stiff flame the reaction sheet is one cell thick, so the
source that reaches the Poisson solver is a delta on one cell. That cell hops
from one step to the next, and the velocity field hops with it. The heat that a
cell releases in one step diffuses over `sqrt(alpha*dt)` before the next step,
which is one to two cells in `fatehi-combustion.c`. A filter of that width puts
the expansion where the physics puts it anyway. A wider filter is not physical:
the velocity divergence no longer matches the density change of the cell.

The filter is explicit diffusion in `gas_source_filter_passes` passes. One
pass applies `D*tau = Delta_min^2/(4*dimension)`, with `Delta_min` the size of
the finest cell. At that level the centre keeps one half of its value and the
neighbours share the other half. The width grows as
`sigma^2 = passes*Delta_min^2/(2*dimension)`. In two dimensions 4 passes give
`sigma = Delta_min`. Coarse cells receive the same `D*tau`, so the physical
width is the same everywhere and the explicit step stays stable. Set
`gas_source_filter_passes = 0` to keep the exact source and switch the filter
off.

The filter conserves the total source. `gas_source` and `drhodt` carry the
metric `cm[]`, so the filter divides by `cm[]` first, diffuses the unweighted
rate with the face metric `fm.x`, and multiplies by `cm[]` after. The flux is
per unit face area, in the convention of `bcg.h`, so the face average that
Basilisk applies at a level jump keeps the flux the same on both sides. At the
axis of an axisymmetric case `fm.y` is zero, and at the domain boundary the
default symmetry condition gives a zero flux. `test/gas-source-filter.c` checks
the conservation on a tree with a level jump.

The source of `psi` in `velocity-potential.h` is not filtered. The interface
must move with the local consumption of the solid. */

#if GAS_SOURCE_EXACT
# ifndef GAS_SOURCE_FILTER_PASSES
#  define GAS_SOURCE_FILTER_PASSES 4
# endif
int gas_source_filter_passes = GAS_SOURCE_FILTER_PASSES;

static void filter_divergence_source (scalar s)
{
  scalar q[];
  foreach()
    q[] = (cm[] > 0.) ? s[]/cm[] : 0.;

  int maxl = 0;
  foreach (reduction(max:maxl))
    if (level > maxl)
      maxl = level;
  double dmin = L0/(1 << maxl);
  const double Dtau = sq(dmin)/(4.*dimension);

  for (int pass = 0; pass < gas_source_filter_passes; pass++) {
    face vector flux[];
    foreach_face()
      flux.x[] = fm.x[]*(q[] - q[-1])/Delta;
    foreach()
      if (cm[] > 0.) {
        double d = 0.;
        foreach_dimension()
          d += flux.x[1] - flux.x[];
        q[] += Dtau*d/(Delta*cm[]);
      }
  }

  foreach()
    s[] = q[]*cm[];
}
#endif

/**
## Projection method with gas source term
We modify the Projection method to account for the gas source term
in the continuity equation.
*/

trace
mgstats project_sf (face vector uf, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{

  /**
  The gas source term and the expansion term enter the divergence together.
  With `GAS_SOURCE_EXACT` the filter above acts on their sum. */

  scalar src[];
  foreach() {
    src[] = gas_source[];
#ifndef NO_EXPANSION
    src[] += drhodt[];
#endif
  }
#if GAS_SOURCE_EXACT
  if (gas_source_filter_passes > 0)
    filter_divergence_source (src);
#endif

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] /= dt*Delta;
    div[] += src[]/dt;
  }

#ifdef POROUS_ADVECTION
  scalar eps[];
  foreach()
    eps[] = porosity[] + (1. - f[]);

  face vector alpha_eff[];
  foreach_face()
    alpha_eff.x[] = face_value(eps, 0)*alpha.x[];
#else
  #define alpha_eff alpha
#endif

  mgstats mgp = poisson (p, div, alpha_eff,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  foreach_face()
    uf.x[] -= dt*alpha_eff.x[]*face_gradient_x (p, 0);

  return mgp;
}

/**
## Advection with non diverging velocity field
We provide an advection scheme that can account for a diverging velocity field.
*/

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


/**
## Default events and overrides
We set default values for the gas source and drhodt fields
*/

event defaults (i = 0) {
  foreach(){
    gas_source[] = 0.;
    drhodt[] = 0.;
  }
}

/** 
We set placeholder to set the correct order of the events 
*/

event set_dtmax (i++, last);
event stability (i++, last);
event reset_sources (i++, last);
event chemistry (i++, last);
event phasechange (i++, last);

/** 
We overwrite the project and advection events with the one defined above
*/
#define project(...) project_sf(__VA_ARGS__)
#define advection(...) advection_div(__VA_ARGS__)
#include "navier-stokes/centered.h"
#undef advection
#undef project


/**
# Porous media advection
We have the option to account for porous media advection by taking into 
account the porosity field in the advection term.
We set stokes=true to suppress the original advection term performed in the
centered.h file.
*/

#ifdef POROUS_ADVECTION
event defaults (i = 0) {
  stokes = true;
}

event advection_term (i++, last) {
  prediction();
  mgpf = project_sf (uf, pf, alpha, dt/2., mgpf.nrelax);

  /**
  porosity is a tracer field appended to f. Here we need the one-field form
  computed in the field 'eps'.
  */

  scalar eps[];
  foreach()
    eps[] = porosity[] + (1. - f[]);

  face vector ufn[];
  foreach_face() {
    double ef = face_value(eps, 0);
    ufn.x[] = uf.x[]/ef;
  }
  
  advection ((scalar *){u}, ufn, dt, (scalar *){g});
}

/** 
the stability event gets disable if stokes is set to true
since the solution is implicit. Therefore, we redefine the
stability event to set the timestep.
*/

event stability (i++, last) {
  dt = dtnext (timestep (uf, dtmax));
}
#endif
