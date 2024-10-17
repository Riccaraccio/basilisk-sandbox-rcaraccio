#include "intgrad.h"
#include "fracface.h"
#include "diffusion.h"
#include "common-evaporation.h"

extern scalar porosity;
extern double lambda1, lambda2, cp1, cp2;
extern double TS0, TG0;

scalar T[];
scalar sT[];
face vector lambdaf[];
face vector fs[];

scalar rhocpeff[], lambdaeff[];

#ifndef RADIATION_INTERFACE
# define RADIATION_INTERFACE 0.
#endif

double divq_rad_int (double TInti, double Tbulk = 300., double alphacorr = 1.) {
  return alphacorr*5.670373e-8*(pow(Tbulk, 4.) - pow(TInti, 4.));
}

event init (i=0) {
  foreach()
    T[] = TS0*f[] + TG0*(1-f[]);
}

event reset_sources (i++) {
  foreach()
    sT[] = 0.;
}

#include "bcg.h" // TO BE CHANGED
extern face vector ufsave;
event tracer_advection (i++) {

 //Advection of T using ug !!WARN!! should reconstruct the darcy velocity
 //advection ({T}, ufsave, dt);

}

event tracer_diffusion (i++) {

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    porosity[] = f[] > F_ERR ? porosity[]/f[] : 1.;
  }

  face_fraction (f, fs);

  foreach() {
    rhocpeff[] = (f[]-porosity[]*f[])*rho1*cp1 + (1-f[]+porosity[]*f[])*rho2*cp2;
    lambdaeff[] = (f[]-porosity[]*f[])*lambda1 + (1-f[]+porosity[]*f[])*lambda2;
 }

  //Compute source terms
  foreach ()
    if (f[] > F_ERR && f[] < 1-F_ERR) {

      coord n = facet_normal (point, f, fs), p;
      double alpha = plane_alpha (f[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      //rho and cp should be weighted by porosity
#ifdef AXI
      sT[] += divq_rad_int(T[], TG0, RADIATION_INTERFACE)/rho2/cp2*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      sT[] += divq_rad_int(T[], TG0, RADIATION_INTERFACE)/rho2/cp2*area/Delta*cm[];
#endif
    }

  scalar theta[];

#if TREE
  theta.refine = theta.prolongation = fraction_refine;
  theta.dirty = true;
#endif


  foreach_face(){
    lambdaf.x[] = face_value(lambdaeff, 0)/ face_value(rhocpeff, 0) *fm.x[];
  }

  foreach()
    theta[] = cm[];

  diffusion (T, dt, D=lambdaf, r=sT, theta=theta);
}
