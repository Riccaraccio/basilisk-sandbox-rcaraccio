#include "run.h"
#include "view.h"
//#include "display.h"

const int maxlevel = 4;

int main () {
  init_grid (1<<maxlevel);
  run();
}

scalar w[];
event init (i = 0) {
  int idx = 0;
  foreach() {
    w[] = idx;
    idx++;
  }

  clear();
  view (tx=-0.5, ty=-0.5);
  box();
  cells();
  labels ("w");
  squares ("w", spread=-1);
  save ("default_order.png");
}

#if 0
  fputc ('\n', stderr);
  display_url (stderr);
  fputc ('\n', stderr);

  display ("box();");
  while (1) {
    if (display_poll (-1))
      display_update (INT_MAX);
  }
  display_destroy();
#endif
