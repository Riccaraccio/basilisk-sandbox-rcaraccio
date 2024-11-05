#include "embed.h"
#include "navier-stokes/centered.h"

double Reynolds = 20.;
int maxlevel = 10;
face vector muv[];

int main()
{
  L0 = 5;
  mu = muv;
  init_grid (1<< maxlevel);
  run(); 
}

double D = 1., U0 = 1;

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*D*U0/Reynolds;
}

u.n[left]  = dirichlet(U0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[embed] = fabs(y-L0/2) > D/2 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y-L0/2) > D/2 ? neumann(0.) : dirichlet(0.);

event init (t = 0)
{
  solid (cs, fs, sqrt(sq(x - L0/2) + sq(y - L0/2)) - D/2.);
  foreach()
    u.x[] = cs[] ? U0 : 0.;
}

event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-2,1e-2,1e-2}, maxlevel, 6);
}

event logfile (t+=1) {
  fprintf (stderr, "%g\n", t);
}

event stop (t=200) {
  char name[80];
  sprintf (name, "WakeLength-%d", maxlevel);
  FILE * fp = fopen (name, "w");

  double step = L0/(1<<maxlevel);
  for (double x = 0; x<L0; x+=step){
    fprintf(fp, "%g %g\n", x, interpolate(u.x, x, L0/2));
  }

  fprintf (fp, "\n\n");
  fflush(fp);

  return 1;
}
