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
double Da = 5e-3; //to be chaged to coord Da

event defaults (i = 0) {
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() {
      a.x[] = 0.;
      dimensional (a.x[] == Delta/sq(DT));
    }
  }
}

face vector ef;
event acceleration (i++){
  face vector av = a;
  face vector ff[];
  face_fraction(f,ff);
  foreach_face() {
    if (ff.x[] > F_ERR) {
      ef.x[] = 1. - (1-porosity[]/f[])*ff.x[];
      double F  = 1.75/pow (150*pow (ef.x[], 3), 0.5);

      // Darcy contribution, weighted by the face fraction of the interface
      av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (muG*ef.x[]/Da) *uf.x[] *ff.x[]; 

      // Forcheimer contribution
      av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (F*ef.x[]*rhoG/pow(Da,0.5)) *fabs(uf.x[])*uf.x[] *ff.x[];
    }
  }
}
