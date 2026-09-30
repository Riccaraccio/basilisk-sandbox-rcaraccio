/**
# Conservation of *shift_field()*

The test fills a field in the cut cells of a solid region and shifts it to the
pure solid cells. The sum of the field times the cell volume must not change.
The solid region has a thin strip (thinner than one cell), a small cap, and a
disc on the bottom boundary. A level jump crosses the interface. The test runs
the old and the new version of the function. The old version loses the value of
a cut cell with no pure solid neighbour. It also loses or doubles the parts
that go to halo cells and to ghost cells. */

#include "grid/quadtree.h"
#include "fractions.h"

/* common-phasechange.h reads the global f and two OpenSMOKE functions. The
test does not use the species, so the stubs return no species. */
scalar f[];
int OpenSMOKE_IndexOfSolidSpeciesWithoutError (const char * s) { return -1; }
int OpenSMOKE_IndexOfSolidSpecies (const char * s) { return -1; }

#include "common-phasechange.h"

/* The old version, copied from commit 93ed0a0 */
void shift_field_old (scalar fts, scalar f, int dir) {
  scalar avg[];
  foreach() {
    avg[] = 0.;
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      int count = 0;
      foreach_neighbor (1)
        if (dir == 1 ? f[] > 1.-F_ERR : f[] < F_ERR)
          count ++;
      avg[] = count;
    }
  }
  scalar sf0[];
  foreach() {
    sf0[] = fts[];
    if (f[] > F_ERR && f[] < 1. - F_ERR)
      fts[] = 0.;
  }
  foreach() {
    if (dir == 1 ? f[] > 1.-F_ERR : f[] < F_ERR) {
      double val = 0.;
      foreach_neighbor (1)
        if (f[] > F_ERR && f[] < 1. - F_ERR && avg[] > 0)
          val += sf0[]/avg[];
      fts[] += val;
    }
  }
}

scalar prod[];

double total (scalar s) {
  double tot = 0.;
  foreach (reduction(+:tot))
    tot += s[]*sq(Delta);
  return tot;
}

void fill (void) {
  foreach()
    prod[] = (f[] > F_ERR && f[] < 1. - F_ERR) ? 1. + x + 2.*y : 0.;
}

int main() {
  origin (-0.5, 0.);
  init_grid (1 << 6);
  refine (x > 0.02 && level < 7);   // the level jump at x = 0.02 crosses the interface

  fraction (f, max (max (sq(0.3) - sq(x) - sq(y),
                         0.004 - fabs (y - 0.6)),
                    sq(0.03) - sq(x + 0.3) - sq(y - 0.85)));

  for (int dir = 1; dir >= 0; dir--) {
    fill();
    double t0 = total (prod);
    shift_field_old (prod, f, dir);
    double told = total (prod);

    fill();
    shift_field (prod, f, dir);
    double tnew = total (prod);

    int nkeep = 0;
    foreach (reduction(+:nkeep))
      if (f[] > F_ERR && f[] < 1. - F_ERR && prod[] != 0.)
        nkeep++;

    fprintf (stderr, "dir %d before %.15e old %.15e (rel %+.3e) "
             "new %.15e (rel %+.3e) cut cells that keep the value %d\n",
             dir, t0, told, told/t0 - 1., tnew, tnew/t0 - 1., nkeep);
    assert (fabs (tnew/t0 - 1.) < 1e-12);
  }
}
