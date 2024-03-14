/*
Header for porous media

We solve the porosity and calculate the velocity field to move the interface
*/
extern double eps0;
extern double rhop1, rhop2;
extern double omega;

scalar * mEvapList = NULL;
scalar eps[], feps[], * f_tracers = NULL;
scalar zeta[];

typedef enum {
  ZETA_SHRINK = 0,
  ZETA_SWELLING,
  ZETA_SMOOTH,
  ZETA_SHARP
} zeta_types;

zeta_types zeta_policy;

event defaults (i = 0) {
  zeta_policy = ZETA_SHARP;
}

event init (i = 0) {
  foreach() {
    feps[] = eps0*f[];
    eps[] = eps0;
  }
  f_tracers = f.tracers;
  f.tracers = list_append (f.tracers, feps);
}

event cleanup (t = end)  {
  free (f.tracers), f.tracers = f_tracers;
}

event phasechange (i++) {

  // Remove spurious values
  foreach() {
    f[] = clamp(f[], 0., 1.);
    feps[] = clamp(feps[], 0., 1.);
  }
  double radius = sqrt(statsf(f).sum/pi);
  // Compute zeta
  switch (zeta_policy) {
    case 0: // ZETA_SHRINK
      foreach() 
        zeta[] = 1.;
      break;
    case 1: // ZETA_SWELLING
      foreach() 
        zeta[] = 0.;
      break;
    case 2: // ZETA_SMOOTH
      foreach() 
        zeta[] = 1 / (1 + exp(32*radius - 40*sqrt(sq(x)+sq(y))));
      break;
    case 3: // ZETA_SHARP
      foreach()
        zeta[] = (sqrt(sq(x) + sq(y)) > radius*0.8) ? 1. : 0.;
      break;
  }

  foreach(){
    drhodt[] = omega/rhop1*f[]*zeta[];
    feps[] += (feps[] > F_ERR) ? -dt*omega/rhop1*feps[]*(1-zeta[]) : 0.;
  }
}