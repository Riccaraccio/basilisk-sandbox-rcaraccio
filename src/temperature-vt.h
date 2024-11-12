#include "intgrad.h"
#include "fracface.h"
#include "diffusion.h"
#include "common-evaporation.h"
#include "int-temperature-v.h"

#define SOLVE_TEMPERATURE

extern scalar porosity;
extern double rhoS, rhoG;

double lambdaS = 1.; double lambdaG = 1.;
double cpS = 1.; double cpG = 1.;
double TS0 = 300.; double TG0 = 600.;
bool success;

scalar T[], TInt[];
scalar TS, TG;
scalar sST[], sGT[];
face vector lambda1f[], lambda2f[];

scalar fG[], fS[];
face vector fsS[], fsG[];
scalar f0[];
scalar fu[]; //dummy tracer

event defaults (i=0) {

  fS.nodump = true;
  fG.nodump =true;

  TS = new scalar;
  TG = new scalar;

  TS.inverse = false;
  TG.inverse = true;

  sST.nodump = true;
  sGT.nodump = true;

  fu.tracers = NULL;

  fu.tracers = list_append (fu.tracers, TS);
  fu.tracers = list_append (fu.tracers, TG);

  TS.refine = refine_linear;
  TS.restriction  = restriction_volume_average;
  TS.dirty = true;

  TG.refine = refine_linear;
  TG.restriction  = restriction_volume_average;
  TG.dirty = true;
}

event init (i=0) {

  foreach()
    fu[] = f[];

#if TREE
  {
    fu.refine = fu.prolongation = fraction_refine;
    fu.dirty = true;
    scalar* tracers = fu.tracers;
    for (scalar t in tracers) {
      t.restriction = restriction_volume_average;
      t.refine = t.prolongation = vof_concentration_refine;
      t.dirty = true;
      t.c = fu;
    }
  }
#endif
  {
    scalar* tracers = fu.tracers;
    for (scalar t in tracers)
      t.depends = list_add (t.depends, fu);
  }

  foreach() {
    TS[] = TS0*f[];
    TG[] = TG0*(1. - f[]);
    T[]  = TS[] + TG[];
    TInt[] = f[]<1-F_ERR && f[] > F_ERR ? TS0 : 0.;
  }
}

event cleanup (t = end) {
  delete ({TS,TG});
  delete (fu.tracers), free(fu.tracers), fu.tracers = NULL;
}

event reset_sources (i++) {
  foreach() {
    sST[] = 0.;
    sGT[] = 0.;
  }
}

extern face vector ufsave;
event tracer_advection (i++) {
  foreach() 
    fu[] = f[];

  //recover pure form
  foreach() {
  //   TS[] *= f[] > F_ERR ? 1./f[] : 0.;
  //   TG[] *= f[] < 1-F_ERR ? 1./(1-f[]) : 0.;
    porosity[] *= f[] > F_ERR ? 1./f[] : 0.;
  }

  //calculate darcy velocity
  face vector darcyv[];
  foreach_face() {
    double epsif = face_value(porosity, 0);
    darcyv.x[] = ufsave.x[]*epsif;
  }

  foreach_face()
    uf.x[] = ufsave.x[];

  vof_advection ({fu}, i);

  // //advect both temperature fields
  // advection_div ({TS}, ufsave, dt);
  // advection_div ({TG}, ufsave, dt);

  //remove values
  foreach() {
  //   TS[] *= f[];
  //   TG[] *= (1-f[]);
    porosity[] *= f[];
  }
}

event tracer_diffusion (i++) {

  foreach() {
   f[] = clamp (f[], 0., 1.);
   f[] = (f[] > F_ERR) ? f[] : 0.;
   fS[] = f[]; fG[] = 1. - f[];
   TS[] = fu[] > F_ERR ? TS[]/fu[] : 0.;
   TG[] = ((1. - fu[]) > F_ERR) ? TG[]/(1. - fu[]) : 0.;
  }

  //Compute face gradients
  face_fraction (fS, fsS);
  face_fraction (fG, fsG);

  //Assign interface temperature first guess
  foreach() {
    TInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      TInt[] = avg_neighbor (point, TS, f);
  }

#ifdef FIXED_INT_TEMP //Force interface temperature = TG0
  foreach()
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      TInt[] = TG0;
#else //default: solve interface balance
  ijc_CoupledTemperature();
#endif

  //Compute source terms 
  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord n = facet_normal (point, fS, fsS), p;
      double alpha = plane_alpha (fS[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      double bc = TInt[];
      double Gtrgrad = ebmgrad (point, TG, fS, fG, fsS, fsG, true, bc, &success);
      double Strgrad = ebmgrad (point, TS, fS, fG, fsS, fsG, false, bc, &success);

      double Sheatflux = lambda1v[]*Strgrad;
      double Gheatflux = lambda2v[]*Gtrgrad;

#ifdef AXI
      sST[] += Sheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sGT[] += Gheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      sST[] += Sheatflux*area/Delta*cm[];
      sGT[] += Gheatflux*area/Delta*cm[];
#endif
    }
  }

  scalar theta1[], theta2[];

#if TREE
  theta1.refine = theta1.prolongation = fraction_refine;
  theta2.refine = theta2.prolongation = fraction_refine;
  theta1.dirty = true;
  theta2.dirty = true;
#endif

  foreach_face() {
    lambda1f.x[] = face_value(lambda1v, 0)*fsS.x[]*fm.x[];
    lambda2f.x[] = face_value(lambda2v, 0)*fsG.x[]*fm.x[];
  }

  foreach() {
    theta1[] = cm[]*max(fS[]*rhocp1v[], F_ERR);
    theta2[] = cm[]*max(fG[]*rhocp2v[], F_ERR);
  }

  diffusion (TS, dt, D=lambda1f, r=sST, theta=theta1);
  diffusion (TG, dt, D=lambda2f, r=sGT, theta=theta2);

  foreach() {
    TS[] *= f[];
    TG[] *= (1. - f[]);
    T[] = TS[] + TG[];
  }
}

//event stability (i++, last) {
//  // does not work, dt found is 100x smaller than actual dtmin
//  double dtmax = DT;
//  foreach_face (reduction(min:dtmax))
//    if (lambda1f.x[] > 0 || lambda2f.x[] > 0){
// 
//      dt = Delta*Delta/(2*dimension)/max(lambda1f.x[],lambda2f.x[]);
// 
//      if (dt < dtmax) dtmax = dt;
//    }
//  dt = dtnext (dtmax);
//}
