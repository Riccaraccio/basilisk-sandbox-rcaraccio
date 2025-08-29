/**
# Darcy flow

We implement the acceleration term due to flow in porous media.

## Implementation
We linearize the Forchheimer term and apply the correction 
before the viscous term is computed.
*/

#ifndef F_ERR
  #define F_ERR 1e-10
#endif

extern scalar porosity;
extern scalar f;

#ifdef VARPROP
extern scalar rhoGv_S, muGv_S;
#else
extern double rhoG, muG;
#endif

coord Da = {1e-10, 1e-10};

event viscous_term (i++) {
  foreach() {
    if (f[] > F_ERR) {
      double e = porosity[]/f[];
      double F = 1.75/sqrt(150*pow (e, 3));
      double Umag = sqrt(sq(u.x[]) + sq(u.y[]));

      double muGh, rhoGh;
      #ifdef VARPROP
      muGh = muGv_S[];
      rhoGh = rhoGv_S[];
      #else
      muGh = muG;
      rhoGh = rhoG;
      #endif

      foreach_dimension() {
        double A = (1./rhoGh)*(muGh*e/Da.x);
        double B = F*e/sqrt(Da.x)*Umag;

        // u.x[] = u.x[]/(1. + dt*(A+B)*f[]); // Implicit v1 does not work
        u.x[] = u.x[]*exp(-(A+B)*dt*f[]);     // Implicit v2 does not work, but better than v1
        // u.x[] = u.x[]*(1. - dt*(A+B)*f[]); // Explicit v1 works
        // u.x[] = u.x[]*(1. - dt/2*(A+B)*f[])/(1. + dt/2*(A+B)*f[]); // Crank-Nicolson in between, does not work
      }
    }
  }
}


// Old approach, works but explicit
// event defaults (i = 0) {
//   if (is_constant(a.x)) {
//     a = new face vector;
//     foreach_face() {
//       a.x[] = 0.;
//       dimensional (a.x[] == Delta/sq(DT));
//     }
//   }
// }

// event acceleration (i++){
//   face vector av = a;
//   foreach_face() {
//     double ff = face_value(f, 0);
//     if (ff > F_ERR) {
//       double ef = face_value(porosity, 0);
//       double F  = 1.75/pow (150*pow (ef, 3), 0.5);

//       // Darcy contribution, weighted by the face fraction of the interface
//       av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (muG*ef/Da.x) *uf.x[] *ff; 

//       // Forcheimer contribution
//       av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (F*ef*rhoG/pow(Da.x,0.5)) *fabs(uf.x[])*uf.x[] *ff;
//     }
//   }
// }