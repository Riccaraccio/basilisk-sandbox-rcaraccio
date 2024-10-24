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

vector gT[], grc[];


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
face vector darcyv[];

event tracer_advection (i++) {
  //Reconstrucy T
  foreach()
    T[] = TS[] + TG[];

  foreach_face() {
    double ff = face_value (f, 0);
    darcyv.x[] = ff>F_ERR ? ufsave.x[]*face_value(porosity, 0)/ff : ufsave.x[];
  }

  //Advection of T using the darcy velocity
  advection ({T}, darcyv, dt);

  //Reconstruct TG, TS is kept
  foreach()
    TG[] = T[]*(1-f[]);
}

void gradients_f (scalar c, scalar* s, vector* g) {
  assert (list_len(s) == vectors_len(g));
  foreach() {
    scalar a; vector v;
    for (a,v in s,g) {
      foreach_dimension() {
        if (c[1] && c[-1]) //centered
            v.x[] = (a[1]-a[-1])/(2.*Delta);
        else
            v.x[] = (a[]-a[-1])/Delta*ceil(c[-1]) + (a[1]-a[])/Delta*ceil(c[1]);

        if (c[] == 0) v.x[] = 0.; //handle edgecase
      }
    }
  }
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

#ifdef FIXED_INT_TEMP
  //Force interface temperature 
  foreach() { 
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      TInt[] = TG0;
  }
#else //default

  //Solve interface balance
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

      double Sheatflux = lambda1v[]*Strgrad; //fprintf(stderr,"sST[] = %g\n", Sheatflux);
      double Gheatflux = lambda2v[]*Gtrgrad;

#ifdef AXI
# ifdef SHIFT
      sST[] += Sheatflux/rhocp1v[]*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sGT[] += Gheatflux/rhocp2v[]*area*(y + p.y*Delta)/(Delta*y)*cm[];
# else
      sST[] += Sheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sGT[] += Gheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
# endif
#else
# ifdef SHIFT
      sST[] += Sheatflux/rhocp1v[]*area/Delta*cm[];
      sGT[] += Gheatflux/rhocp2v[]*area/Delta*cm[];
# else
      sST[] += Sheatflux*area/Delta*cm[];
      sGT[] += Gheatflux*area/Delta*cm[];
# endif
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
#ifdef SHIFT
    lambda1f.x[] = face_value(lambda1v, 0)/face_value(rhocp1v, 0)*fsS.x[]*fm.x[];
    lambda2f.x[] = face_value(lambda2v, 0)/face_value(rhocp2v, 0)*fsG.x[]*fm.x[];
#else
    lambda1f.x[] = face_value(lambda1v, 0)*fsS.x[]*fm.x[];
    lambda2f.x[] = face_value(lambda2v, 0)*fsG.x[]*fm.x[];
#endif
  }

  foreach() {
#ifdef SHIFT
    theta1[] = cm[]*max(fS[], F_ERR);
    theta2[] = cm[]*max(fG[], F_ERR);
#else
    theta1[] = cm[]*max(fS[]*rhocp1v[], F_ERR);
    theta2[] = cm[]*max(fG[]*rhocp2v[], F_ERR);
#endif
  }

#ifdef SHIFT
  scalar t[];

  foreach()
    t[] = fS[]>F_ERR ? 1/rhocp1v[] : 0.;
  gradients_f (fS, {TS, t}, {gT, grc});
  foreach()
    foreach_dimension()
      sST[] -= fS[] > F_ERR ? lambda1v[]*gT.x[]*grc.x[] : 0.;

  trash({gT, grc, t});

  foreach()
    t[] = fG[]>F_ERR ? 1/rhocp2v[] : 0.;
  gradients_f (fG, {TG, t}, {gT, grc});
  foreach()
    foreach_dimension()
      sGT[] -= fG[] > F_ERR ? lambda2v[]*gT.x[]*grc.x[] : 0.;

#endif

  diffusion (TS, dt, D=lambda1f, r=sST, theta=theta1);
  diffusion (TG, dt, D=lambda2f, r=sGT, theta=theta2);

  foreach() {
    TS[] *= f[];
    TG[] *= (1. - f[]);
  }
  foreach()
    T[] = TS[] + TG[];
}

//event properties (i++) { //to be moved to var-prop module
//  foreach() {
//    lambda1v[] = f[] > F_ERR ? pavg (porosity[]/f[], lambdaG, lambdaS) : lambdaG;
//    lambda2v[] = lambdaG;
//
//    rho1v[] = f[] > F_ERR ? pavg (porosity[]/f[], rhoG, rhoS) : rhoG;
//    rho2v[] = rhoG;
//
//    rhocp1v[] =  f[] > F_ERR ? pavg (porosity[]/f[], rhoG*cpG, rhoS*cpS) : rhoG*cpG;
//    rhocp2v[] =  rhoG*cpG;
//  }
//
//}

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
