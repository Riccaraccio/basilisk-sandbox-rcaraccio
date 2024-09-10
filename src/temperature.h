#include "intgrad.h"
#include "fracface.h"
#include "diffusion.h"
#include "common-evaporation.h"

extern double lambda1, lambda2, dhR, cp1, cp2;
extern double TS0, TG0;
bool success;

scalar T[], TInt[];
scalar TS, TG;
scalar sST[], sGT[], sSTimp[], sGTimp[];
face vector lambda1f[], lambda2f[];

scalar fG[], fS[];
face vector fsS[], fsG[];
scalar f0[];

event defaults (i=0) {

  fS.nodump = true;
  fG.nodump =true;

  TS = new scalar;
  TG = new scalar;

  TS.inverse = false;
  TG.inverse = true;

  sST.nodump = true;
  sGT.nodump = true;
  sGTimp.nodump = true;
  sSTimp.nodump = true;

  f.tracers = list_append (f.tracers, TS);
  f.tracers = list_append (f.tracers, TG);

  TS.refine = refine_linear;
  TS.restriction  = restriction_volume_average;
  TS.dirty = true;

  TG.refine = refine_linear;
  TG.restriction  = restriction_volume_average;
  TG.dirty = true;

}

event cleanup (t=end) {
  delete ({TS,TG});
}

event init (i=0) {
  foreach() {
    TS[] = TS0*f[];
    TG[] = TG0*(1. - f[]);
    T[]  = TS[] + TG[];
  }
}

event reset_sources (i++) {
  foreach(){
    sST[] = 0.;
    sGT[] = 0.;
    sSTimp[] = 0.;
    sGTimp[] = 0.;
  }
}

event phasechange (i++) {

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    f0[] = f[];
    fS[] = f[]; fG[] = 1. - f[];

    TS[] = f[] > F_ERR ? TS[]/f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]/(1. - f[]) : 0.;
  }

  //Compute face gradients
  face_fraction (fS, fsS);
  face_fraction (fG, fsG);

  //Assign interface temperature
  foreach() {
    TInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      TInt[] = avg_neighbor (point, TS, f);
    fprintf(stderr, "fs = %g\t fsS = %g\n", fS[], fsS.x[]);
  }

  //Compute source terms 
  fprintf(stderr, "i = %d\n", i);
  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord n = facet_normal (point, fS, fsS), p;
      double alpha = plane_alpha (fS[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      double bc = TInt[];
      double Gtrgrad = ebmgrad (point, TG, fS, fG, fsS, fsG, true, bc, &success);
      double Strgrad = ebmgrad (point, TG, fS, fG, fsS, fsG, false, bc, &success);
      fprintf(stderr, "%d\n", success);
      fprintf(stderr, "%g\n",Strgrad);
      double Sheatflux = lambda1*Strgrad;
      double Gheatflux = lambda2*Gtrgrad;
#ifdef AXI
      sST[] += Sheatflux/rho1/cp1*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sGT[] += Gheatflux/rho2/cp2*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      sST[] += Sheatflux/rho1/cp1*area/Delta*cm[];
      sGT[] += Gheatflux/rho2/cp2*area/Delta*cm[];
#endif
    }
  }
  foreach() {
    TS[] *= f[];
    TG[] *= (1. - f[]);
  }
}

event tracer_advection (i++);

event tracer_diffusion (i++) {

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    fS[] = f[]; fG[] = 1. - f[];
    TS[] = f[] > F_ERR ? TS[]/f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]/(1. - f[]) : 0.;
  }

  //Compute face gradients
  face_fraction (fS, fsS);
  face_fraction (fG, fsG);

  scalar theta1[], theta2[];

  foreach_face() {
    lambda1f.x[] = lambda1/rho1/cp1*fsS.x[]*fm.x[];
    lambda2f.x[] = lambda2/rho2/cp2*fsG.x[]*fm.x[];
  }

  foreach() {
    theta1[] = cm[]*max(fS[], F_ERR);
    theta2[] = cm[]*max(fG[], F_ERR);
  }

  diffusion (TS, dt, D=lambda1f, r=sST, theta=theta1);
  diffusion (TG, dt, D=lambda2f, r=sGT, theta=theta2);

  foreach() {
    TS[] *= f[];
    TG[] *= (1. - f[]); //TODO should check for f[]?
  }
  foreach() 
    T[] = TS[] + TG[];
}

