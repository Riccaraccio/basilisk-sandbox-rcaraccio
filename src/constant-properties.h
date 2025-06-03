#ifndef F_ERR
# define F_ERR 1.e-10
#endif

#define pavg(p,vg,vs) (p*vg+(1-p)*vs)

extern double rhoG, rhoS;
extern double muG;
extern face vector alphav;
extern scalar rhov;
extern scalar f;
extern scalar porosity;

#ifdef SOLVE_TEMPERATURE
extern vector lambda1v, lambda2v;
lambda2v.n[bottom] = neumann(0.); 
lambda2v.t[bottom] = neumann(0.);
lambda2v.n[left] = neumann(0.);
lambda2v.t[left] = neumann(0.);
lambda2v.n[right] = neumann(0.);
lambda2v.t[right] = neumann(0.);
lambda2v.n[top] = neumann(0.);
lambda2v.t[top] = neumann(0.);

lambda1v.n[bottom] = neumann(0.);
lambda1v.t[bottom] = neumann(0.);
lambda1v.n[left] = neumann(0.);
lambda1v.t[left] = neumann(0.);
lambda1v.n[right] = neumann(0.);
lambda1v.t[right] = neumann(0.);
lambda1v.n[top] = neumann(0.);
lambda1v.t[top] = neumann(0.);

extern double lambdaG, lambdaS;
extern double cpG, cpS;
#endif
#ifdef MULTICOMPONENT
extern unsigned int NGS;
extern scalar* DmixGList_S;
extern scalar* DmixGList_G;
#endif

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

#ifdef SOLVE_TEMPERATURE
  foreach() {
    foreach_dimension() {
      lambda1v.x[] = f[] > F_ERR ? pavg (porosity[]/f[], lambdaG, lambdaS) : 0.;
      lambda2v.x[] = f[] < 1.-F_ERR ? lambdaG : 0.;
    }
  }
#endif

#ifdef MULTICOMPONENT
  double Dmixv =  2.05e-5; //Diff of CO in N2 at 500K, 1 atm

  foreach() {
    // set the same for all species
    for (int jj = 0; jj < NGS; jj++) {
      scalar Dmix2 = DmixGList_G[jj];
      Dmix2[] = f[] < 1. - F_ERR ? Dmixv : 0.;
    }
    for (int jj = 0; jj < NGS; jj++) {
      scalar Dmix2 = DmixGList_S[jj];
      Dmix2[] = f[] > F_ERR ? Dmixv * pow(porosity[], 3. / 2.) : 0.; // effect of solid, to be revised
    }
  }
#endif
}
