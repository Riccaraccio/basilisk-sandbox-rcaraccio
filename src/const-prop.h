#ifndef F_ERR
# define F_ERR 1.e-10
#endif

#define pavg(p,vg,vs) (p*vg+(1-p)*vs)

scalar rho1v[], rho2v[];
extern vector lambda1v, lambda2v;
extern face vector alphav;
extern scalar rhov;
extern scalar f;
extern scalar porosity;

extern double rhoG, rhoS;
extern double lambdaG, lambdaS;
extern double cpG, cpS;
extern double muG;
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

  //Update other fields properties
  foreach() {

    rho1v[] = f[] > F_ERR ? pavg (porosity[]/f[], rhoG, rhoS) : rhoG;
    rho2v[] = rhoG;

#ifdef SOLVE_TEMPERATURE
    foreach_dimension() {
      lambda1v.x[] = f[] > F_ERR ? pavg (porosity[]/f[], lambdaG, lambdaS) : lambdaG;
      lambda2v.x[] = lambdaG;
    }
#endif
  }

#ifdef MULTICOMPONENT
  double Dmixv =  2.05e-5; //Diff of CO in N2 at 500K, 1 atm

  foreach() {
    if (f[] < 1.-F_ERR) {
    //set the same for all species
      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix2 = DmixGList_G[jj];
        Dmix2[] = Dmixv;
      }
    } 
    if (f[] > F_ERR) {
      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix2 = DmixGList_S[jj];
        Dmix2[] = Dmixv*pow(porosity[], 3./2.); //effect of solid, to be revised
      }
    }
  }
#endif
}
