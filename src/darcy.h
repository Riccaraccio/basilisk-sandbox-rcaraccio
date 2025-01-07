/**
# Darcy flow

We implement the acceleration term due to flow in porous media.
*/

#include "fracface.h"

#ifndef F_ERR
  #define F_ERR 1e-10
#endif

extern scalar porosity;
extern double rhoG, muG;

event defaults (i = 0) {
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() {
      a.x[] = 0.;
      dimensional (a.x[] == Delta/sq(DT));
    }
  }
}

double Da = 5e-3; //to be chaged to coord Da

event acceleration (i++){
  face vector av = a;
  //face_fraction(f,fS);
  foreach_face() {
    if (f[]>F_ERR) {
      double ef = face_value (porosity, 0);
      double F  = 1.75/pow (150*pow (ef, 3), 0.5);

      // Darcy contribution
      av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (muG*ef/Da) *uf.x[];

      // Forcheimer contribution
      av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (F*pow(ef,2)*rhoG/pow(Da,0.5)) *fabs(uf.x[])*uf.x[];
    }
  }
}
