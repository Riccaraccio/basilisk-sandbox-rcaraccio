/**
# Conservation of the divergence-source filter

`navier-stokes/centered-phasechange.h` filters the divergence source with a
short explicit diffusion. This test checks that the filter conserves the
total source on an axisymmetric tree with a level jump.

The source is two small discs of about two fine cells. One disc sits on the
level jump, so that fluxes cross coarse-fine faces. The other disc touches the
axis, where `fm.y` is zero. Both carry the metric `cm[]`, as `gas_source`
does. The volume integral of the source is the cell sum of
`gas_source*Delta^2`, up to the factor `2*pi`. The test fails if the filter
changes it by more than `1e-12` in relative terms, or if one pass does not
reduce the maximum of the unweighted rate. */

#include "grid/quadtree.h"
#include "axi.h"

double rhoG = 1.;
scalar f[], porosity[];

#include "navier-stokes/centered-phasechange.h"

static double total (scalar s) {
  double tot = 0.;
  foreach (reduction(+:tot))
    tot += s[]*sq(Delta);
  return tot;
}

/**
The maximum of the unweighted rate in the fine disc (`x > 0.5`) or in the
coarse disc at the axis (`x < 0.5`). */

static double maximum (scalar s, bool fine) {
  double m = 0.;
  foreach (reduction(max:m))
    if ((x > 0.5) == fine && cm[] > 0. && fabs (s[]/cm[]) > m)
      m = fabs (s[]/cm[]);
  return m;
}

int main() {
  L0 = 1.;
  origin (0., 0.);
  init_grid (1 << 5);
  refine (sq(x - 0.5) + sq(y - 0.3) < sq(0.15) && level < 8);

  foreach() {
    double disc1 = (sq(x - 0.648) + sq(y - 0.3) < sq(0.005)) ? 1. : 0.;
    double disc2 = (sq(x - 0.2) + sq(y) < sq(0.04)) ? 1. : 0.;
    gas_source[] = -cm[]*(disc1 + disc2);
  }

  double tot0 = total (gas_source);
  double max0 = maximum (gas_source, true), cmax0 = maximum (gas_source, false);

  int status = 0;
  for (int passes = 1; passes <= 8; passes *= 2) {
    scalar s[];
    foreach()
      s[] = gas_source[];
    gas_source_filter_passes = passes;
    filter_divergence_source (s);
    double tot1 = total (s);
    double max1 = maximum (s, true), cmax1 = maximum (s, false);
    double err = fabs (tot1/tot0 - 1.);
    fprintf (stderr, "passes %d total %.12e ratio-1 %.3e"
             " fine max %.4f -> %.4f coarse max %.4f -> %.4f\n",
             passes, tot1, err, max0, max1, cmax0, cmax1);
    if (err > 1e-12 || max1 >= max0 || cmax1 >= cmax0)
      status = 1;
  }
  fprintf (stderr, status ? "FAIL\n" : "PASS\n");
  return status;
}
