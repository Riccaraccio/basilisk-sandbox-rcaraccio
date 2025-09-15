#include "vofi.h"
#pragma autolink -L$HOME/Vofi/build/lib -lvofi

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "view.h"

static double xc = 0., yc = 0., Rc = 0.3;

static double sphere (vofi_creal p[dimension]) {
  return sq(p[0] - xc) + sq(p[1] - yc) - sq(Rc);
}

static void vofi (scalar c, int maxlevel) {
  double fh = vofi_Get_fh (sphere, NULL, 1./(1 << maxlevel), dimension, 0);
  foreach() {
    vofi_creal p[3] = {x - Delta/2., y - Delta/2., z - Delta/2.};
    c[] = vofi_Get_cc (sphere, p, Delta, fh, dimension);
  }
}

int maxlevel = 7;

int main (void) {
  init_grid (1 << maxlevel);
  run();
}

event init (i = 0) {
  vofi (f, maxlevel);

  clear();
  view (tx = -0.5, ty = -0.5);
  squares ("f", spread = -1);
  draw_vof ("f", lw = 2.5);
  save ("picture.png");
}
