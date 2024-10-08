#include "navier-stokes/centered.h"
#include "two-phase.h"

scalar psi[];
face vector ubf[];
mgstats mgpsf;
scalar omega[];
scalar zeta[];

psi[bottom]  = neumann (- neumann_pressure(0));
psi[left]  = neumann (- neumann_pressure(0));
psi[top]  = dirichlet (0);
psi[right]  = dirichlet (0);

trace
mgstats project_sv (face vector ubf, scalar psi,
    (const) face vector alpha = unityf,
    double dt = 1.,
    int nrelax = 4)
{
  scalar prod[];
  foreach()
    prod[] =(omega[]*zeta[]*f[]/100)/dt;

  mgstats mgp = poisson (psi, prod, alpha,
      tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  foreach_face()
    ubf.x[] = -dt*alpha.x[]*face_gradient_x (psi, 0);

  return mgp;
}

int main () {
  init_grid (1<<5);
  DT = 1e-1;
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))
event init (i=0) {
  fraction (f, circle (x, y, 0.5*L0));
  foreach()
    omega[] = 1.;
}

event phasechange (i++){
  double radius = sqrt(statsf(f).sum/pi);
  foreach()
    zeta[] = 1 / (1 + exp(32*radius - 40*sqrt(sq(x)+sq(y))));
  mgpsf = project_sv (ubf, psi, alpha, dt, mgpsf.nrelax);
}

face vector uf_save[];
event vof (i++) {
  foreach_face()
    uf_save.x[] = uf.x[];
  foreach_face()
    uf.x[] = ubf.x[];
}

event tracer_advection (i++) {
  foreach_face()
    uf.x[] = uf_save.x[];
}

event stability (i++) {
  dt = dtnext (timestep (ubf, dtmax));
}

event logprofile (i++) {
  fprintf(stderr, "========== i=%d ==========\n",i);
  fprintf(stderr, "Number of iterations: %d\n", mgpsf.i);
  fprintf(stderr, "Maximum residual before iterations: %f\n", mgpsf.resb);
  fprintf(stderr, "Maximum residual after iterations: %f\n", mgpsf.resa);
  fprintf(stderr, "Sum of r.h.s.: %f\n", mgpsf.sum);
  fprintf(stderr, "Number of relaxations: %d\n", mgpsf.nrelax);
  fprintf(stderr, "Minimum level of the multigrid hierarchy: %d\n", mgpsf.minlevel);
}
event stop (t = 100) {
  return 1;
}
