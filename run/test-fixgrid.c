/**
# The base case on a fixed grid

This case is `test.c` with no mesh adaptation. It tests one question: does
the adaptation of the mesh cause the one-sample spikes of the velocity?

In `test-fbest` the point probes 0 and 2 spike on the same samples. At a
spike the residual of the projection is 10 to 20 times its median, and the
grid changes by about 3 times more cells than at a normal sample. The
release does not change. If the adapt causes the spikes, this case has no
spikes. If the spikes stay, the adapt is not the cause.

The case restores the snapshot of `test-base` at t = 15 s, as the restart
rungs `test-s*` of `run/Makefile`. After the restore it refines the region
below to `maxlevel` one time, and it never adapts again. Outside the region
the grid of the snapshot stays as it is. The region holds the particle, the
nine point probes and the near wake.

Compare with `test-sprobe`: the same restart, the same probe of the
projection (`projstep.dat`), and the normal adapt.

Caution: this file is a real source, not one of the `test-*.c` links to
`test.c` that `run/Makefile` makes. Two rules keep it: an explicit rule
`test-fixgrid.c: ;` in `run/Makefile` stops the pattern rule from replacing
it with a link, and an exception in `.gitignore` keeps it tracked. Do not
remove either rule, and do not add `test-fixgrid` to `TLAD_CASES`, because
`make test-unlink` removes those files.

Cost: the region holds about 18 000 cells at level 10, against about 8 500
cells in the whole grid of `test-base`. Expect about 2 to 2.5 times the cost
of `test-sref` for each step. The CFL limit comes from the finest cells, so
`dt` does not change. */

#define FIXED_GRID 1

/**
The region, in metres. The particle has its centre at the origin and a
diameter `D0` of 8 mm. The flow goes in the +x direction. The defaults give
x from -12 mm to +24 mm and y from 0 to 12 mm. Set the three flags to change
the region. */

#ifndef FIXED_GRID_XMIN
# define FIXED_GRID_XMIN (-1.5*D0)
#endif

#ifndef FIXED_GRID_XMAX
# define FIXED_GRID_XMAX (3.*D0)
#endif

#ifndef FIXED_GRID_YMAX
# define FIXED_GRID_YMAX (1.5*D0)
#endif

#define FIXED_GRID_REGION \
  (x > FIXED_GRID_XMIN && x < FIXED_GRID_XMAX && y < FIXED_GRID_YMAX)

#include "test.c"
