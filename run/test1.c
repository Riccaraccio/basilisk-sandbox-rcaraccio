#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};
double Reynolds = 160.;
int maxlevel = 9;
face vector muv[];

int main()
{
  L0 = 8. [1];
  origin (-0.5, -L0/2.);
  N = 512;
  mu = muv;

  run(); 
}

double D = 0.125, U0 = 1.;

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*D*U0/Reynolds;
}

u.n[left]  = dirichlet(U0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);


u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);

event init (t = 0)
{
  solid (cs, fs, intersection (intersection (0.5 - y, 0.5 + y),
			       sqrt(sq(x) + sq(y)) - D/2.));
  
  foreach()
    u.x[] = cs[] ? U0 : 0.;
}
event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-2,3e-2,3e-2,3e-2}, maxlevel, 4);
}

event stop (t = 2);