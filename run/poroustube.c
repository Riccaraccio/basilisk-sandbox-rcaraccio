#include "navier-stokes/centered.h"
#include "fractions.h"
// #define DARCY 1
#define F_ERR 1e-8

double uin = 0.8;
u.n[left] = dirichlet (uin);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);
pf[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
pf[right] = dirichlet (0.);

int maxlevel = 8;

double eps0 = 0.5;
scalar eps[];
face vector epsf[];

int main () {
  const face vector muv[] = {1e-1, 1e-1};
  mu = muv;

  stokes = true;
  init_grid (1 << maxlevel);
  run();
}

scalar f[];
vector v[];
scalar p_prime[];

event init (i = 0) {

  fraction (f, x > 0.4 && x < 0.6);

  foreach()
    eps[] = f[]*eps0 + (1. - f[]);

  foreach_face()
    epsf.x[] = face_value (eps, 0);

  foreach()
    u.x[] = (1. - f[])*uin;
}

event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}

event endtimestep (i++) {
  boundary({u});
  foreach()
    foreach_dimension()
      v.x[] = u.x[]/eps[];

  foreach()
    p_prime[] = p[]/eps[];
}

event advection_term (i++,last) {
  prediction();
  mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax);

  face vector ufn[];
  foreach_face()
    ufn.x[] = uf.x[]/epsf.x[];
  
  advection ((scalar *){u}, ufn, dt, (scalar *){g});
}

#if DARCY
double Da = 5e-3;
event defaults (i = 0) {
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() {
      a.x[] = 0.;
      dimensional (a.x[] == Delta/sq(DT));
    }
  }
}

face vector ef;
event acceleration (i++){
  face vector av = a;
  foreach_face() {
    double ff = face_value(f,0);
    if (ff > F_ERR) {
      double F  = 1.75/pow (150*pow (epsf.x[], 3), 0.5);

      // Darcy contribution
      av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (mu.x[]*epsf.x[]/Da) *uf.x[] *ff; 

      // Forcheimer contribution
      av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (F*epsf.x[]*rho[]/pow(Da,0.5)) *fabs(uf.x[])*uf.x[] *ff;
    }
  }
}
#endif

event logprofile (i = 100) {
  char name[80];
  sprintf (name, "Output");
  static FILE * fp = fopen (name, "w");

  for (double x = 0.; x < L0; x += L0/(1<<maxlevel))
    fprintf (fp, "%g %g %g %g %g\n", x, interpolate (u.x, x, 0.5), interpolate (v.x, x, 0.5),  interpolate (p, x, 0.5),  interpolate (p_prime, x, 0.5));
  fflush (fp);
}

event stop (i = 100);

/** 
~~~gnuplot velocity profile
reset
set xlabel "x"
set ylabel "velocity"
set yrange [0:2.5]
set key top right box width 1
set samples 40

u_in = 0.8
epsilon = 0.5
analytical_u(x) = u_in
analytical_v(x) = (x < 0.4) ? u_in : (x > 0.6) ? u_in : u_in/epsilon
plot "Output" u 1:2 w l lw 2 lc "red" t "u", \
     "Output" u 1:3 w l lw 2 lc "web-green" t "v", \
      analytical_u(x) w p pt 4 lc "red" t "analytical u", \
      analytical_v(x) w p pt 4 lc "web-green" t "analytical v"
~~~

~~~gnuplot pressure profile
reset
set xlabel "x"
set ylabel "pressure"
set yrange [-2.5:0.5]
set samples 40
set key bottom right box width 1

u_in = 0.8
epsilon = 0.5

analytical_p(x) = (x < 0.4) ? 0 : (x > 0.6) ? 0 : u_in*u_in*(1-1/epsilon)
analytical_p_prime(x) = (x < 0.4) ? 0 : (x > 0.6) ? 0 : u_in*u_in*(1-1/epsilon)/epsilon

plot "Output" u 1:4 w l lw 2 lc "red" t "p", \
     "Output" u 1:5 w l lw 2 lc "web-green" t "p'", \
      analytical_p(x) w p pt 4 lc "red" t "analytical p", \
      analytical_p_prime(x) w p pt 4 lc "web-green" t "analytical p'"
~~~
*/