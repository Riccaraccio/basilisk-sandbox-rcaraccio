#include "grid/multigrid.h"
#include "run.h"
#include "timestep.h"
#include "diffusion.h"

scalar T[];
mgstats dstats;

double TR = 300.;
double TL = 1000.;

T[left] = dirichlet (TL);
T[top] = neumann (0.);
T[bottom] = neumann (0.);

int main () {
  init_grid (1<<5);
  TOLERANCE = 1e-3;
  DT = 10.;
  L0 = 1e-2;
  run();
}
double dtmax;
event init (i=0) {
  dtmax = DT;
  event ("stability");

  foreach()
    T[] = 300.;
}

event stability (i++,last) {
  dt = dtnext (DT);
}

event integration (i++) {

  //(const) face vector lambda; lambda[] = {1e-2, 1e-2};
  //(const) scalar s; s[] = 0.;
  //(const) scalar theta; theta[] = 1e6; //this does NOT work

  face vector lambda[];
  scalar s[];
  scalar theta[];

  foreach() {
    theta[] = 1e6;
    s[] = 0.;
  }

  foreach_face()
    lambda.x[] = 1e-1;

  dstats = diffusion (T, dt, D=lambda, r=s, theta=theta);
  fprintf (stderr, "n iter = %d\n", dstats.i);

}

event stop (t = 1000) {
  return 1;
}
