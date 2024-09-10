#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "temperature.h"
#include "view.h"
double TS0 = 300.;
double TG0 = 500.;

double lambda1 = 1e-2;
double cp1 =1.;

double lambda2 = 1e-1;
double cp2 =1.;

int main() {
  DT = 1e-2;
  init_grid(1 << 5);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5));
}

event movie (t += 0.1) {
  clear();
  view (tx=-0.5, ty=-0.5);
  draw_vof ("f", lw=1.5);
  box();
  squares ("T", linear = true, min=TS0, max=TG0);
  save ("movie.mp4");
}


event stop (t = 10);

