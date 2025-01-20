#ifndef F_ERR
# define F_ERR 1.e-10
#endif

#define pavg(p,vg,vs) (p*vg+(1-p)*vs)

scalar lambda1v[], lambda2v[];
scalar rho1v[], rho2v[];
scalar rhocp1v[], rhocp2v[];

extern face vector alphav;
extern scalar rhov;
extern scalar f;
extern scalar porosity;

extern double rhoG, rhoS;
extern double lambdaG, lambdaS;
extern double cpG, cpS;
extern double muG;

event properties (i++) {
  //Navier-Stokes equation are solved using the gas properties on the whole domain
  foreach_face() {
    alphav.x[] = fm.x[]/rhoG;
    {
      face vector muv = mu;
      muv.x[] = muG*fm.x[];
    }
  }

  foreach()
    rhov[] = rhoG*cm[];

  //Update other fields properties
  foreach() {

    rho1v[] = f[] > F_ERR ? pavg (porosity[]/f[], rhoG, rhoS) : rhoG;
    rho2v[] = rhoG;

#ifdef SOLVE_TEMPERATURE
    lambda1v[] = f[] > F_ERR ? pavg (porosity[]/f[], lambdaG, lambdaS) : lambdaG;
    lambda2v[] = lambdaG;

    rhocp1v[] =  f[] > F_ERR ? pavg (porosity[]/f[], rhoG*cpG, rhoS*cpS) : rhoG*cpG;
    rhocp2v[] =  rhoG*cpG;
#endif
  }
}
