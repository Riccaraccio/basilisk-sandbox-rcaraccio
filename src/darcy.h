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
extern double rhoG, muG;
coord Da = {5e-3, 5e-3};

event viscous_term (i++) {
  correction (dt);
  foreach() {
    if (f[] > F_ERR) {
      double e = porosity[]/f[];
      double F = 1.75/pow (150*pow (e, 3), 0.5);
      double Umag = sqrt(sq(u.x[]) + sq(u.y[]));

      foreach_dimension() {
        double A = alpha.x[]/(fm.x[] + SEPS)*(muG*e/Da.x)*f[]; 
        double B = alpha.x[]/(fm.x[] + SEPS)*(F*e*rhoG/pow(Da.x,0.5))*Umag*f[];

        u.x[] /= (1. + (A+B)*dt);
      }
    }
  }
  correction (-dt);
}