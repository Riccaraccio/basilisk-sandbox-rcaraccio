/**
# Darcy flow

We solve the Darcy equation for the pressure field $p$ in a porous medium
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

double Da = 0;
event acceleration (i++){
    face vector av = a;
    foreach_face(){
        double ff = face_value (f, 0);
        av.x[] += (mu2*Da*ff*uf.x[]/rho2);
    }
}
