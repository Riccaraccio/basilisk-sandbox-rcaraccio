#include "intgrad.h"
#include "fracface.h"
#include "diffusion.h"
#include "common-evaporation.h"
#include "int-temperature.h"

extern double lambda1, lambda2, cp1, cp2;
extern double TS0, TG0;
bool success;

mgstats diffstats;
scalar T[], TInt[];
scalar TS, TG;
scalar sST[], sGT[];
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

  f.tracers = list_append (f.tracers, TS);
  f.tracers = list_append (f.tracers, TG);

  TS.refine = refine_linear;
  TS.restriction  = restriction_volume_average;
  TS.dirty = true;

  TG.refine = refine_linear;
  TG.restriction  = restriction_volume_average;
  TG.dirty = true;

}

event init (i=0) {
  foreach() {
    TS[] = TS0*f[];
    TG[] = TG0*(1. - f[]);
    T[]  = TS[] + TG[];
    TInt[] = f[]<1-F_ERR && f[] > F_ERR ? TS0 : 0.;
  }

 //foreach_face() {
 //  lambda1f.x[] = lambda1/rho1/cp1; //first guess needed for stability event
 //  lambda2f.x[] = lambda2/rho2/cp2;
 //}
}

event cleanup (t = end) {
  delete ({TS,TG});
}

event reset_sources (i++) {
  foreach(){
    sST[] = 0.;
    sGT[] = 0.;
  }
}

#include "bcg.h"
extern face vector ufsave;
event tracer_advection (i++) {
 //Reconstrucy T
 foreach()
   T[] = TS[] + TG[];

 //Advection of T using ug !!WARN!! should reconstruct the darcy velocity
 advection ({T}, ufsave, dt);

 //Reconstruct TG, TS is kept
 foreach()
    TG[] =T[]*(1-f[]);
}

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

 //Assign interface temperature first guess
  foreach() {
    TInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      TInt[] = avg_neighbor (point, TS, f);
  }

//  //Force interface temperature 
//  foreach() { 
//    if (f[] > F_ERR && f[] < 1.-F_ERR)
//      TInt[] = TG0;
//  }

  //Solve interface balance
  ijc_CoupledTemperature();

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

      double Sheatflux = lambda1*Strgrad;
      double Gheatflux = lambda2*Gtrgrad;
#ifdef AXI
      sST[] += Sheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sGT[] += Gheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
      //sST[] += Sheatflux/rho1/cp1*area*(y + p.y*Delta)/(Delta*y)*cm[];
      //sGT[] += Gheatflux/rho2/cp2*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      sST[] += Sheatflux*area/Delta*cm[]; // add dhR
      sGT[] += Gheatflux*area/Delta*cm[];
      //sST[] += Sheatflux/rho1/cp1*area/Delta*cm[]; // add dhR
      //sGT[] += Gheatflux/rho2/cp2*area/Delta*cm[];
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
    lambda1f.x[] = lambda1*fsS.x[]*fm.x[];
    lambda2f.x[] = lambda2*fsG.x[]*fm.x[];
    //lambda1f.x[] = lambda1/rho1/cp1*fsS.x[]*fm.x[];
    //lambda2f.x[] = lambda2/rho2/cp2*fsG.x[]*fm.x[];
 }

  foreach() {
    theta1[] = cm[]*max(fS[]*rho1*cp1, F_ERR);
    theta2[] = cm[]*max(fG[]*rho2*cp2, F_ERR);
    //theta1[] = cm[]*max(fS[], F_ERR);
    //theta2[] = cm[]*max(fG[], F_ERR);
  }

  diffstats = diffusion (TS, dt, D=lambda1f, r=sST, theta=theta1);
  diffusion (TG, dt, D=lambda2f, r=sGT, theta=theta2);
  fprintf(stderr, "Number of iterations: %d\n", diffstats.i);
  fprintf(stderr, "Maximum residual before iterations: %f\n", diffstats.resb);
  fprintf(stderr, "Maximum residual after iterations: %f\n", diffstats.resa);
  fprintf(stderr, "Sum of r.h.s.: %f\n", diffstats.sum);
  fprintf(stderr, "Number of relaxations: %d\n", diffstats.nrelax);
  fprintf(stderr, "Minimum level of the multigrid hierarchy: %d\n", diffstats.minlevel);

  foreach() {
    TS[] *= f[];
    TG[] *= (1. - f[]);
  }
  foreach()
    T[] = TS[] + TG[];
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
