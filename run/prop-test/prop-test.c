#include "navier-stokes/centered-phasechange.h"
#include "var-prop.h"
//#include "opensmoke-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent.h"
#include "temperature-vt.h"

//scalar omega[];

int main() {
  init_grid (1<<4);
  run();
}

double D0 = 0.5;
#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init(i=0) {
  fraction (f, circle (x, y, 0.5*D0));

  // gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  // sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")] = 1.;

  foreach()
    porosity[] = eps0*f[];
  
/*   update_properties_initial(); */
}

event stop (i=2);


