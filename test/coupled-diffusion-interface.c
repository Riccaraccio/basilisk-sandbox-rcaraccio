/**
# The coupled diffusion next to an interface

The test reproduces the setting of the species solves of
`multicomponent-varprop.h` in a small problem. A circle of radius 0.2 cuts
a quadtree, and the volume fraction `f` and the face fractions `fs` come
from the same distance function. The phase is either the inside of the
circle (as the pore gas, `PHASE = 0`) or the outside (as the external gas,
`PHASE = 1`).

As in `tracer_diffusion`:

* the species hold the value of the phase, and $Y = 0$ where the phase
  fraction is below `F_ERR`;
* $\theta = cm \max(\rho \varepsilon f_{ph}, F_{ERR})$ and
  $\alpha = \rho f_s f_m$;
* the species are VOF tracers, with `restriction_volume_average`.

The three species have MW 2, 28 and 162, and D 1, 0.2 and 0.05. The test
runs the molar and the Fick forms, with upwind face values. Each run must
converge, keep $\sum Y = 1$ in the full and the cut cells of the phase,
keep the mass of each species, and keep the species positive. */

#include "grid/quadtree.h"
#include "fractions.h"
#include "coupled-diffusion.h"
#include "utils.h"

#ifndef PHASE
# define PHASE 0
#endif

#define NS 3
#define F_ERR 1e-10
#define EPS 0.5

scalar f[], Ya[], Yb[], Yc[], Da[], Db[], Dc[];
scalar thetac[];
face vector fs[], alpha[];
double MWs[NS] = {2., 28., 162.};

static double phase (Point point)
{
  return PHASE ? 1. - f[] : f[];
}

static void initial (scalar * YL)
{
  scalar a = YL[0], b = YL[1], c = YL[2];
  foreach() {
    if (phase (point) > F_ERR) {
      double s = 0.5*(1. + tanh ((x - 0.5)/0.05));
      a[] = 0.05 + 0.6*s;
      b[] = 0.85 - 0.8*s;
      c[] = 1. - a[] - b[];
    }
    else
      a[] = b[] = c[] = 0.;
  }
}

static void check (scalar * YL, double * mass0, const char * name,
                   mgstats st)
{
  double efull = 0., ecut = 0., ymin = HUGE, dm = 0.;
  int ncut = 0, nsliver = 0;
  foreach (reduction(max:efull) reduction(max:ecut) reduction(min:ymin)
           reduction(+:ncut) reduction(+:nsliver)) {
    double ph = phase (point);
    if (ph > F_ERR) {
      double s = 0.;
      for (scalar Y in YL)
        s += Y[], ymin = min (ymin, Y[]);
      if (ph < 1. - F_ERR)
        ecut = max (ecut, fabs (s - 1.)), ncut++, nsliver += (ph < 1e-3);
      else
        efull = max (efull, fabs (s - 1.));
    }
  }
  for (int i = 0; i < NS; i++) {
    scalar Y = YL[i];
    double m = 0.;
    foreach (reduction(+:m))
      m += thetac[]*Y[]*sq(Delta);
    dm = max (dm, fabs (m - mass0[i])/mass0[i]);
  }
  fprintf (stderr, "%s: cycles %d, res %.2e, masked faces %d, cut cells %d "
           "(slivers %d), |sum Y - 1| full %.2e cut %.2e, min Y %.3e, "
           "rel d(mass) %.2e\n", name, st.i, st.resa, cd_nmasked, ncut,
           nsliver, efull, ecut, ymin, dm);
  assert (st.resa <= 1e-6);
  assert (efull < 1e-9 && ecut < 1e-9 && dm < 1e-10 && ymin >= 0.);
}

int main()
{
  init_grid (32);
  double R = 0.2 + 1e-3;
  refine (fabs (sqrt (sq(x - 0.5) + sq(y - 0.5)) - R) < 0.05 && level < 7);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(R) - sq(x - 0.5) - sq(y - 0.5);
  fractions (phi, f, fs);

  scalar * YL = {Ya, Yb, Yc}, * DL = {Da, Db, Dc};
  for (scalar Y in YL)
    set_restriction (Y, restriction_volume_average);

  foreach() {
    thetac[] = cm[]*max (EPS*phase (point), F_ERR);
    Da[] = 1., Db[] = 0.2, Dc[] = 0.05;
  }
  foreach_face()
    alpha.x[] = (PHASE ? 1. - fs.x[] : fs.x[])*fm.x[];

  int nslivers = 0;
  foreach (reduction(+:nslivers))
    nslivers += (phase (point) > F_ERR && phase (point) < 1e-3);
  fprintf (stderr, "phase %s, cells with 0 < phase fraction < 1e-3: %d\n",
           PHASE ? "outside" : "inside", nslivers);

  for (int molar = 1; molar >= 0; molar--) {
    initial (YL);
    double mass0[NS];
    for (int i = 0; i < NS; i++) {
      scalar Y = YL[i];
      double m = 0.;
      foreach (reduction(+:m))
        m += thetac[]*Y[]*sq(Delta);
      mass0[i] = m;
    }
    for (int s = 0; s < 3; s++) {
      mgstats st = coupled_diffusion (YL, 1e-3, DL, alpha, thetac,
                                      MWs = molar ? MWs : NULL,
                                      tolerance = 1e-8);
      char name[80];
      sprintf (name, "%s, step %d", molar ? "molar" : "Fick ", s + 1);
      check (YL, mass0, name, st);
    }
  }
}
