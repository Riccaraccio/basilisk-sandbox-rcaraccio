#include "embed.h"
#include "navier-stokes/centered.h"

double Reynolds = 20.;
int maxlevel = 7;
face vector muv[];

double D = 1., U0 = 0.1;

int main()
{
  L0 = 10*D;
  mu = muv;
  init_grid (1<<maxlevel);
  run(); 
}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*D*U0/Reynolds;
}

u.n[left]  = dirichlet(U0);
p[left]    = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);

u.n[embed] = fabs(y-L0/2) > D/2 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y-L0/2) > D/2 ? neumann(0.) : dirichlet(0.);

event init (t = 0) {
  solid (cs, fs, -(sqrt(sq(x - 0.5*L0) + sq(y - 0.5*L0)) - D/2.));
  //fraction(f, -(sqrt(sq(x - 0.5*L0) + sq(y - 0.5*L0)) - D/2.));
  foreach()
   u.x[] = cs[] ? U0 : 0.;
}

event logprofile (t+=1) {
  fprintf (stderr, "%g\n", t);
}

event stop (t=100) {
  char name[80];
  sprintf (name, "WakeLength");
  FILE * fp = fopen (name, "w");

  double step = L0/(1<<maxlevel);
  for (double x = 0; x<L0; x+=step){
    fprintf(fp, "%g %g\n", x, interpolate(u.x, x, L0/2));
  }

  fprintf (fp, "\n\n");
  fflush(fp);

  return 1;
}
