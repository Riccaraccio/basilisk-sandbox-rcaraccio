extern face vector alphav;
extern scalar rhov;
extern scalar f;
extern scalar porosity;

extern double rhoG, rhoS;
extern double lambdaG, lambdaS;
extern double cpG, cpS;
extern double muG;

scalar lambda1v[], lambda2v[];
scalar rhocp1v[], rhocp2v[];
scalar rho1v[], rho2v[];
scalar cpSv[], cpGv[];

//NOTEtoSELF: later introduce rhoSv, rhoGv, cpSv ...

#define F_ERR 1.e-10
#define pavg(p,vg,vs) (p*vg+(1-p)*vs)

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

//#ifdef SOLVE_TEMPERATURE
  foreach() {
    lambda1v[] = f[] > F_ERR ? pavg (porosity[]/f[], lambdaG, lambdaS) : lambdaG;
    lambda2v[] = lambdaG;

    rho1v[] = f[] > F_ERR ? pavg (porosity[]/f[], rhoG, rhoS) : rhoG;
    rho2v[] = rhoG;

    cpSv[] = cpS;
    cpGv[] = cpG;

    rhocp1v[] =  f[] > F_ERR ? pavg (porosity[]/f[], rhoG*cpG, rhoS*cpS) : rhoG*cpG;
    rhocp2v[] =  rhoG*cpG;
  }
//#endif
}
