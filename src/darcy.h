/**
# Darcy flow

We implement the acceleration term due to flow in porous media.
*/

event defaults (i = 0) {  
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() {
      a.x[] = 0.;
      dimensional (a.x[] == Delta/sq(DT));
    }
  }
}

double Da = 200; //to be chaged to coord Da

event acceleration (i++){
    face vector av = a;
    foreach_face(){
        double ff = face_value (f, 0);
        av.x[] -= alpha.x[]/(fm.x[] + SEPS)*(mu.x[]*Da*ff*uf.x[]);
    }
}
