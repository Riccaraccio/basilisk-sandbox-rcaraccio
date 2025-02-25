extern face vector ef;
#include "poisson.h"

trace
mgstats project_porous (face vector uf, scalar p,
            (const)face vector alpha = unityf,
            double dt = 1., int nrelax = 4)
{
  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] /= dt*Delta;
  }

  face vector t[];                                   //
  foreach_face()                                     //
      t.x[] = alpha.x[] * ef.x[];   //

  mgstats mgp = poisson(p, div, t,                                        //
                        tolerance = TOLERANCE / sq(dt), nrelax = nrelax); // OVERWITTEN

  foreach_face()
      uf.x[] -= dt * alpha.x[] * ef.x[] * face_gradient_x(p, 0); // OVERWITTEN

  return mgp;
}


// #include "utils.h"
// #include "bcg.h"
// trace 
// void advection_porous(scalar *tracers, face vector u, double dt,
//                             scalar *src = NULL)
// {
//   scalar *psrc = src;
//   if (!src)
//     for (scalar s in tracers)
//     {
//       const scalar zero[] = 0.;
//       src = list_append(src, zero);
//     }
//   assert(list_len(tracers) == list_len(src));

//   scalar f, source;
//   for (f, source in tracers, src)
//   {
//     face vector flux[];
//     tracer_fluxes(f, u, flux, dt, source);
// #if !EMBED
//     foreach ()
//       foreach_dimension() {
//           f[] += dt * (flux.x[]/ef.x[] - flux.x[1]/ef.x[1]) / (Delta * cm[]);
//       }
// #else  // EMBED
//     update_tracer(f, u, flux, dt);
// #endif // EMBED
//   }

//   if (!psrc)
//     free(src);
// }

// #define project(...) project_porous(__VA_ARGS__)
// #define centered_gradient_porous centered_gradient
#include "navier-stokes/centered.h"
// #undef project
// #undef centered_gradient_porous

// void centered_gradient (scalar p, vector g);

// void centered_gradient_porous (scalar p, vector g)
// {
//   face vector gf[];
//   foreach_face()
//     gf.x[] = fm.x[]*a.x[] - alpha.x[]*ef.x[]*(p[] - p[-1])/Delta;

//   trash ({g});
//   foreach()
//     foreach_dimension()
//       g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);
// }

event defaults(i=0) {
  stokes = true;
}

extern scalar porosity;
event advection_term(i++, last) {
  prediction();

  mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax);

  face vector uf_real[];
  foreach_face()
    uf_real.x[] = uf.x[] / ef.x[];
    // uf_real.x[] = uf.x[];
  
  advection ((scalar *){u}, uf_real, dt, (scalar *){g});
  
  // foreach()
  //   foreach_dimension()
  //     u.x[] = u.x[]*porosity[]; // OVERWITTEN
}

// #define advection(...) advection_porous(__VA_ARGS__)
// #undef advection