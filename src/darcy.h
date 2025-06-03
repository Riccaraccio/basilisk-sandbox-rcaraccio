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
  correction (dt);
  foreach() {
    if (f[] > F_ERR) {
      double e = porosity[]/f[];
      double F = 1.75/pow (150*pow (e, 3), 0.5);
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
        double A = (1./rhoGh)*(muGh*e/Da.x)*f[]; 
        double B = (F*e/pow(Da.x,0.5))*Umag*f[];
        u.x[] /= (1. + (A+B)*dt);
      }
    }
  }
  correction (-dt);
}