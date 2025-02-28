/**
# Darcy flow

We implement the acceleration term due to flow in porous media.
*/

#ifndef F_ERR
  #define F_ERR 1e-10
#endif

extern scalar porosity;
extern scalar f;
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

event acceleration (i++){
  face vector av = a;
  foreach_face() {
    double ff = face_value(f, 0);
    if (ff > F_ERR) {
      double ef = face_value(porosity, 0);
      double F  = 1.75/pow (150*pow (ef, 3), 0.5);

      // Darcy contribution, weighted by the face fraction of the interface
      av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (muG*ef/Da) *uf.x[] *ff; 

      // Forcheimer contribution
      av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (F*ef*rhoG/pow(Da,0.5)) *fabs(uf.x[])*uf.x[] *ff;
    }
  }
}
