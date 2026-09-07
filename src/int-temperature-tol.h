/**
# The tolerance of the two temperature solves

`TOLERANCE` is one number for the whole run, but every solver compares it
against a residual in its own units. `run/test.c` sets `TOLERANCE = 1e-5`
because the projection needs it there. The temperature solves inherit it, and
for them it means something quite different.

`diffusion()` solves

    theta*(T^{n+1} - T^n)/dt = div(D grad T^{n+1}) + src

and `poisson.h` compares `TOLERANCE` against the largest cell residual of
that equation, which carries the units of the right side, `theta/dt*T`. So the
temperature error that `TOLERANCE` really asks for is

    eps_K = TOLERANCE*dt/theta

## How wrong the inherited value is

Measured in `test-fullrobin`, the gas solve at the step where it stalls:

    res: 2.0e-05  sum: -1.09e+13  nrelax: 50  tolerance: 1e-05

`sum` is the plain sum of the right side over about 8000 cells, so the mean
right side is about 1.4e9 and `theta/dt` is about 1.4e6. The inherited
tolerance therefore asks the gas temperature to

    eps_K = 1e-5 / 1.4e6 = 7e-12 K

Seven picokelvin. The solve reaches a relative residual of 1.4e-14, which is
about 60 times the double precision epsilon, and it cannot go further.

The cost is not only the wasted digits. `poisson.h:195` raises `nrelax`
whenever the residual is still above the tolerance, so a tolerance that can
never be met drives `nrelax` from 4 up to 50 and beyond, and the solver then
runs the full `NITERMAX` of 100 iterations. One such solve costs about a
thousand normal ones. It happens exactly when `dt` has already collapsed,
which turns a slow patch of a run into a dead one. `test-fullrobin` and
`test-fullpicard` both died that way.

## The repair

Scale the tolerance of each solve to its own residual:

    TOLERANCE = max (TOLERANCE, max(theta)/dt * INT_TEMP_TOL_K)

`INT_TEMP_TOL_K` is the temperature error to ask for, in kelvin. The default
is 1e-6 K, which is six orders below anything the physics can use and still
about five orders looser than the inherited value.

Two properties matter.

- `max(theta)` is read from the `theta1` and `theta2` fields that the solve is
  about to use, so it needs no property constant. An earlier version used
  `rhoS*cpS`, which had no gas counterpart: `run/test.c` never sets `cpG`, so
  `rhoG*cpG` is meaningless, and under `VARPROP` the gas heat capacity is the
  field `cpGv_G`.
- The `max` with the old value means this can only LOOSEN the tolerance, never
  tighten it. No other solver can be made stricter by accident, and a
  configuration whose physical scale falls below 1e-5 keeps 1e-5.

Caution: read `max(theta1)` and `max(theta2)` BEFORE the first `diffusion()`
call. `diffusion()` overwrites `theta` in place with `theta*(-1/dt)`.

Caution: with `INT_TEMP_PICARD` the outer loop cannot converge below what the
linear solve delivers. Keep `INT_TEMP_TOL_K` well under `INT_TEMP_PICARD_TOL`,
which is 1e-2 K by default, and check `rel_max` after any change.

`INT_TEMP_TOL` is the switch, and it is 1 by default. Set it to 0 to restore
the plain inherited `TOLERANCE`. Keep the switch and the value separate: the
preprocessor cannot compare a floating point value in an `#if`.

## How to read `tsolve.dat`

    #t(1) dt(2) tolS(3) iS(4) nrelaxS(5) resaS(6) tolG(7) iG(8) nrelaxG(9) resaG(10)

- `iS`, `iG` the multigrid cycles each solve used. `run/test.c` sets
  `NITERMIN = 2`, so 2 is the floor and 2 means the tolerance never bound.
  A value of 100 is `NITERMAX`: that solve failed.
- `nrelaxS`, `nrelaxG` the relaxation sweeps per cycle. This is the real cost.
  It starts at 4 and `poisson.h` raises it while the residual stays above the
  tolerance. Anything above about 10 means the solve is struggling.
- `resaS`, `resaG` the largest cell residual left. Compare it against the
  tolerance of the same line, not against 1e-5.

With `INT_TEMP_PICARD` the line holds the LAST pass of the step.
*/

#ifndef INT_TEMP_TOL_K
# define INT_TEMP_TOL_K 1e-6
#endif

double ITT_tolS, ITT_tolG;
double ITT_iS, ITT_iG, ITT_nrelaxS, ITT_nrelaxG, ITT_resaS, ITT_resaG;

event tsolve_output (i++, last) {
  static FILE * fp = NULL;
  if (!fp) {
    fp = fopen ("tsolve.dat", restarted ? "a" : "w");
    if (!fp) {
      fprintf (stderr, "Error opening tsolve.dat\n");
      return 0;
    }
    if (!restarted)
      fprintf (fp, "#t(1) dt(2) tolS(3) iS(4) nrelaxS(5) resaS(6)"
                   " tolG(7) iG(8) nrelaxG(9) resaG(10)\n");
  }

  fprintf (fp, "%g %g %g %g %g %g %g %g %g %g\n",
           t, dt, ITT_tolS, ITT_iS, ITT_nrelaxS, ITT_resaS,
           ITT_tolG, ITT_iG, ITT_nrelaxG, ITT_resaG);
  fflush (fp);

  return 0;
}
