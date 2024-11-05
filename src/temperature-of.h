#include "intgrad.h"
#include "fracface.h"
#include "diffusion.h"
#include "common-evaporation.h"

#define SOLVE_TEMPERATURE

extern scalar porosity;
extern double rhoS, rhoG;

double lambdaS = 1.; double lambdaG = 1.;
double cpS = 1.; double cpG = 1.;
double TS0 = 300.; double TG0 = 600.;

scalar T[];
scalar sGT[];
face vector lambdaf[];

face vector fsS[];

#ifndef RADIATION_INTERFACE
# define RADIATION_INTERFACE 0.
#endif

double divq_rad_int (double TInti, double Tbulk = 300., double alphacorr = 1.) {
  return alphacorr*5.670373e-8*(pow(Tbulk, 4.) - pow(TInti, 4.));
}

event defaults (i=0){
  sGT.nodump = true;

  T.refine = refine_linear;
  T.restriction  = restriction_volume_average;
  T.dirty = true;
}

event init (i=0) {
  foreach()
    T[] = TS0*f[] + TG0*(1-f[]);
}

event reset_sources (i++) {
  foreach()
    sGT[] = 0.;
}

#include "bcg.h"
extern face vector ufsave;
event tracer_advection (i++) {

  face_fraction (f, fsS);
  foreach()
    porosity[] = f[] > F_ERR ? porosity[]/f[] : 0.;

  face vector darcyv[];
  foreach_face() {
    darcyv.x[] = fsS.x[] > F_ERR ? ufsave.x[]*porosity[]: ufsave.x[];
  }

  //Advection of T using the darcy velocity
  advection ({T}, darcyv, dt);

  foreach()
    porosity[] *= f[];
}

event tracer_diffusion (i++) {

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
  }

  face_fraction (f, fsS);

  //Compute radiation source terms
  foreach () {
    if (f[] > F_ERR && f[] < 1-F_ERR) {

      coord n = facet_normal (point, f, fsS), p;
      double alpha = plane_alpha (f[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

#ifdef AXI
      sGT[] += divq_rad_int(T[], TG0, RADIATION_INTERFACE)*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      sGT[] += divq_rad_int(T[], TG0, RADIATION_INTERFACE)*area/Delta*cm[];
#endif
    }
  }

  scalar theta[];

#if TREE
  theta.refine = theta.prolongation = fraction_refine;
  theta.dirty = true;
#endif

  foreach_face()
    // lambdaf.x[] = (lambda1v[]*fsS.x[] + (1-fsS.x[])*lambda2v[])*fm.x[];
    lambdaf.x[] = (face_value(lambda1v, 0)*fsS.x[] + (1-fsS.x[])*face_value(lambda2v, 0))*fm.x[];

  foreach()
    theta[] = cm[]*max(rhocp1v[]*f[] + (1-f[])*rhocp2v[], F_ERR);

    diffusion (T, dt, D=lambdaf, r=sGT, theta=theta);
}
