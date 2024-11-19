#define VARPROP

//Phases propeties
// (const) scalar lambda1v = zeroc, lambda2v = zeroc;
// (const) scalar rhocp1v = zeroc, rhocp2v = zeroc;
// (const) scalar rho1v = zeroc, rho2v = zeroc;
// (const) scalar muGv = zeroc;
//
// //Gas and solid properties
// (const) scalar cpSv = zeroc, cpGv = zeroc;
// (const) scalar rhoGv = zeroc, rhoSv = zeroc;
// (const) scalar lambdaSv = zeroc, lambdaGv = zeroc;

(const) scalar lambda1v;
(const) scalar lambda2v;
(const) scalar rhocp1v;
(const) scalar rhocp2v;
(const) scalar rho2v;
(const) scalar rho1v;
(const) scalar muGv;

//Gas and solid properties
(const) scalar cpSv;
(const) scalar cpGv;
(const) scalar rhoSv;
(const) scalar rhoGv;
(const) scalar lambdaSv;
(const) scalar lambdaGv;
#ifndef F_ERR
# define F_ERR 1.e-10
#endif

#define pavg(p,vg,vs) (p*vg+(1-p)*vs)

typedef struct {
  double T, P;
  double* x;
} ThermoState;

typedef struct {
  // Mixture properties
  double (* rhov) (void *);
  double (* muv) (void *);
  double (* lambdav) (void *);
  double (* cpv) (void *);
  // Species properties
  double (* diff) (void *, int);
  double (* cps)  (void *, int); //note: maybe unused
} ThermoProps;

extern face vector alphav;
extern scalar rhov;
extern scalar f;
extern scalar porosity;

// extern double rhoG, rhoS;
// extern double lambdaG, lambdaS;
// extern double cpG, cpS;
// extern double muG;

event defaults (i=0) {
  lambda1v[] = 0.; lambda2v[] = 0.;
  rhocp1v[] = 0.; rhocp2v[] = 0.;
  rho1v[] = 0.; rho2v[] = 0.;

  muGv[] = 0.,
  cpSv[] = 0.; cpGv[] = 0.;
  rhoSv[] = 0.; rhoGv[] = 0.;
  lambdaSv[] = 0.; lambdaGv[] = 0.;

}

event properties (i++) {
  //NS equation are solved using the gas properties every were
  foreach_face() {
    alphav.x[] = fm.x[]/rhoGv[];
    {
      face vector muv = mu;
      muv.x[] = muGv[]*fm.x[];
    }
  }

  foreach()
    rhov[] = rhoGv[]*cm[];

  //Update other fields properties
  foreach() {
    lambda1v[] = f[] > F_ERR ? pavg (porosity[]/f[], lambdaGv[], lambdaSv[]) : lambdaGv[];
    lambda2v[] = lambdaGv[];

    rho1v[] = f[] > F_ERR ? pavg (porosity[]/f[], rhoGv[], rhoSv[]) : rhoGv[];
    rho2v[] = rhoGv[];

    rhocp1v[] =  f[] > F_ERR ? 0. : rhoGv[]*cpGv[];
    // rhocp1v[] =  f[] > F_ERR ? pavg (porosity[]/f[], rhoGv[]*cpGv[], rhoSv[]*cpSv[]) : rhoGv[]*cpGv[];
    rhocp2v[] =  rhoGv[]*cpGv[];

  }
}

//Useful functions from EDOcip

int check_thermostate (ThermoState * ts, int NS) {
  double sum = 0.;
  for (int jj=0; jj<NS; jj++)
    sum += ts->x[jj];

  int T_ok = (ts->T > 180. && ts->T < 4000.) ? true : false;
  int P_ok = (ts->P > 1e3 && ts->P < 1e7) ? true : false;
  int X_ok = (sum > 1.-1.e-3 && sum < 1.+1.e-3) ? true : false;

  return T_ok*P_ok*X_ok;
}

void print_thermostate (ThermoState * ts, int NS, FILE * fp = stdout) {
  fprintf (fp, "Temperature = %g - Pressure = %g\n", ts->T, ts->P);
  for (int jj=0; jj<NS; jj++)
    fprintf (fp, "  Composition[%d] = %g\n", jj, ts->x[jj]);
  fprintf (fp, "\n");
}

void print_thermoprop (ThermoProps * tp, ThermoState * ts, int NS, FILE * fp = stdout) {
  if (tp->rhov)
    fprintf (fp, "rhov = %g\n", tp->rhov (ts));
  if (tp->muv)
    fprintf (fp, "muv = %g\n", tp->muv (ts));
  if (tp->lambdav)
    fprintf (fp, "lambdav = %g\n", tp->lambdav (ts));
  if (tp->cpv)
    fprintf (fp, "cpv = %g\n", tp->cpv (ts));
  if (tp->diff) {
    for (int jj=0; jj<NS; jj++)
      fprintf (fp, "diff[%d] = %g\n", jj, tp->diff (ts, jj));
  }
  if (tp->cps) {
    for (int jj=0; jj<NS; jj++)
      fprintf (fp, "cps[%d] = %g\n", jj, tp->cps (ts, jj));
  }
}
