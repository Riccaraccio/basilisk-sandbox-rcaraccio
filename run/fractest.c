#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "fracface.h"

int main () {
  init_grid (1<<3);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

face vector fsS[];
face vector ff[];

scalar porosity[];
face vector pff[], pfS[];

double D0 = 1;
event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  foreach()
    porosity[] = f[]*0.5;

  face_fraction (f, fsS);

  foreach_face() {
    ff.x[] = face_value(f, 0);
    pff.x[] = face_value(porosity, 0);
    pfS.x[] = fsS.x[] ? porosity[]/fsS.x[] : 0.;
  }
}

event stop(i=10) {

  foreach_face()
    fprintf(stderr, "%g %g %g\n", fsS.x[], ff.x[], (fsS.x[] - ff.x[]));

  return 1;
}
