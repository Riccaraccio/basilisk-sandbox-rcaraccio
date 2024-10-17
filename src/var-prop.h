extern face vector alphav;
extern scalar rhov;

extern double rhoG;
extern double muG;
extern double mu1, mu2, rho1, rho2;

event properties (i++) {

  foreach_face() {
    alphav.x[] = fm.x[]/rhoG;
    {
      face vector muv = mu;
      muv.x[] = muG*fm.x[];
    }
  }

  foreach()
    rhov[] = rhoG*cm[];
}
