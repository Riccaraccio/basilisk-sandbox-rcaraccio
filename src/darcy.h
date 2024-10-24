/**
# Darcy flow

We implement the acceleration term due to flow in porous media.
*/

// extern scalar epsi;
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
    foreach_face(){
      if (f[] > 1e-6)
      {
        //double epsif = face_value(epsi, 0);
        double epsif = face_value(porosity, 0);
        double Cff = 1.75/pow(150*pow(1-epsif,3), 0.5);

       av.x[] -= alpha.x[]/(fm.x[] + SEPS)*(muG*uf.x[]*(1-epsif)/Da); // Darcy term 
       av.x[] -= alpha.x[]/(fm.x[] + SEPS)*(rhoG*pow((1-epsif),2)* Cff* fabs(uf.x[])* uf.x[]/ pow(Da,0.5) ); // Forcheimer term
      }
    }
}
