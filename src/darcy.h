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
double Da = 5e-3; //to be chaged to coord Da

event viscous_term (i++) {
  correction (dt);
  // foreach() {
  //   if (f[] > F_ERR) {
  //     double e = porosity[]/f[];
  //     double F = 1.75/pow (150*pow (e, 3), 0.5);

  //     foreach_dimension() {
  //       double A = alpha.x[]/(fm.x[] + SEPS)*(muG*e/Da)*f[]; //note: porosity[] = e[]*f[] here
  //       double B = alpha.x[]/(fm.x[] + SEPS)*(F*e*rhoG/pow(Da,0.5))*fabs(u.x[])*f[];

  //       u.x[] /= (1. + (A+B)*dt);
  //     }
  //   }
  // }
  foreach() {
    if (f[] > F_ERR) {
      double e = porosity[]/f[];
      double F = 1.75/pow (150*pow (e, 3), 0.5);

      foreach_dimension() {
        double A = alpha.x[]/(fm.x[] + SEPS)*(muG*e/Da)*f[];
        double B = alpha.x[]/(fm.x[] + SEPS)*(F*e*rhoG/pow(Da,0.5))*f[];

        u.x[] = (-(1+A*dt) + sqrt(sq(1+A*dt) + 4*fabs(u.x[])*B*dt))/(2*B*dt)*sign(u.x[]);
      }
    }
  }
  correction (-dt);
}