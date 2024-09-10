#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "temperature.h"

double TS0 = 25.;
double TG0 = 300.;

double lambda1 = 1.;
double cp1 =1.;

double lambda2 = 1.;
double cp2 =1.;

int main() {
  init_grid(1 << 6);
  run();
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5));
}

event stop (t = 100);

