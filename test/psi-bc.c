/**
# Boundary conditions of the velocity potential

The test checks *psi_dirichlet_faces()* of `velocity-potential.h`. With the
default conditions, no boundary face has a Dirichlet condition. After the test
sets a Dirichlet condition on the right boundary, the number of faces is the
number of cells along that boundary. The test then solves for a disc of
shrinkage. All the volume must leave through the right boundary, because the
other boundaries have a zero normal velocity. */

#include "grid/quadtree.h"
#include "run.h"
#include "fractions.h"

#define F_ERR 1.e-10

/* The module reads these fields of the phase change stack. */
scalar omega[], zeta[], f[], fS[];
double rhoS = 1.;
face vector fsS[];
void shift_field (scalar fts, scalar f, int dir) {}

#include "velocity-potential.h"

int main() {
  init_grid (1 << 6);
  run();
}

event init (i = 0) {
  refine (sq(x - 0.5) + sq(y - 0.5) < sq(0.3) && level < 7);

  int n0 = psi_dirichlet_faces (psi);

  psi[right] = dirichlet (0.);
  int n1 = psi_dirichlet_faces (psi);
  int nright = 0;
  foreach_boundary (right, reduction(+:nright))
    nright++;

  fprintf (stderr, "Dirichlet faces: default %d, with psi[right] %d, "
           "faces on right %d\n", n0, n1, nright);
  assert (n0 == 0);
  assert (n1 == nright && n1 > 0);

  fraction (f, sq(0.2) - sq(x - 0.5) - sq(y - 0.5));
  foreach() {
    omega[] = 1.;
    zeta[] = 1.;
  }
  mgpsf = project_sv (ubf, psi);

  double vol = 0.;
  foreach (reduction(+:vol))
    vol += prod[]*sq(Delta);
  double out = 0.;
  foreach_boundary (right, reduction(+:out))
    out += ubf.x[1]*Delta;
  double leak = 0.;
  foreach_boundary (left, reduction(+:leak))
    leak += fabs (ubf.x[])*Delta;

  fprintf (stderr, "source %.10e outflow right %.10e (rel %+.3e) "
           "flux left %.3e, cycles %d, residual %.3e\n",
           vol, out, out/vol + 1., leak, mgpsf.i, mgpsf.resa);
  assert (fabs (out/vol + 1.) < 1e-4);
  assert (leak == 0.);
}
