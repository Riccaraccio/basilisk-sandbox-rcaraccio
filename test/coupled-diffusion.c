/**
# The coupled diffusion of three species in 1-D

The test checks `coupled_diffusion()` of `coupled-diffusion.h` in a closed
box. Species 0 is light ($D = 1$), species 1 is the carrier ($D = 0.2$), and
species 2 is heavy ($D = 0.05$) with a uniform initial value. Only the Fick
correction moves species 2.

The test does five checks:

1. With equal diffusivities, the correction is zero. The result must be
   equal to `diffusion()` for each species.
2. With central face values, the result must converge with first order in
   `dt` to an explicit reference with a very small step. The reference uses
   the same fluxes, without the linearization.
3. Check 2 for the molar form, where the reference uses the exact $X$.
4. Each run must keep $\sum Y = 1$ and the mass of each species.
5. With upwind face values and a step at a corrective Courant number near 8,
   the solver must converge and the species must stay positive.

Check 5 does not use central face values. At this step, the central system
loses its diagonal dominance, and the block Gauss-Seidel diverges. */

#include "grid/multigrid1D.h"
#include "coupled-diffusion.h"
#include "diffusion.h"
#include "utils.h"

#define NS 3
#define TEND 0.01

scalar Ya[], Yb[], Yc[], Da[], Db[], Dc[];
scalar thetac[]; // not theta: utils.h has a global theta
face vector alpha[];
double MWs[NS] = {2., 28., 100.};

static void initial (scalar * YL)
{
  scalar a = YL[0], b = YL[1], c = YL[2];
  foreach() {
    double s = 0.5*(1. + tanh ((x - 0.5)/0.05));
    a[] = 0.05 + 0.8*s;
    b[] = 0.85 - 0.8*s;
    c[] = 0.1;
  }
}

static void set_diffusivities (scalar * DL, double d0, double d1, double d2)
{
  scalar a = DL[0], b = DL[1], c = DL[2];
  foreach()
    a[] = d0, b[] = d1, c[] = d2;
}

/**
One explicit step of the corrected system, with the central face value of
$Y_i$ and, for the molar form, the exact $X$. */

static void explicit_step (scalar * YL, scalar * DL, double * mw, double dt)
{
  scalar * XL = list_clone (YL), * dYL = list_clone (YL);
  scalar MW[];
  foreach() {
    double s = 0.;
    for (int k = 0; k < NS; k++) {
      scalar Y = YL[k];
      s += mw ? Y[]/mw[k] : 0.;
    }
    MW[] = s > 0. ? 1./s : 0.;
    for (int k = 0; k < NS; k++) {
      scalar Y = YL[k], X = XL[k];
      X[] = mw ? Y[]*MW[]/mw[k] : Y[];
    }
  }

  face vector U[], g[];
  foreach_face() {
    U.x[] = 0.;
    for (int k = 0; k < NS; k++) {
      scalar X = XL[k], D = DL[k];
      double A = alpha.x[]*face_value (D, 0);
      if (mw)
        A *= mw[k]/face_value (MW, 0);
      U.x[] += A*(X[] - X[-1])/Delta;
    }
  }
  for (int i = 0; i < NS; i++) {
    scalar X = XL[i], Y = YL[i], D = DL[i], dY = dYL[i];
    foreach_face() {
      double A = alpha.x[]*face_value (D, 0);
      if (mw)
        A *= mw[i]/face_value (MW, 0);
      g.x[] = A*(X[] - X[-1])/Delta - face_value (Y, 0)*U.x[];
    }
    foreach() {
      dY[] = 0.;
      foreach_dimension()
        dY[] += dt*(g.x[1] - g.x[])/(Delta*thetac[]);
    }
  }
  for (int i = 0; i < NS; i++) {
    scalar Y = YL[i], dY = dYL[i];
    foreach()
      Y[] += dY[];
  }
  delete (XL), free (XL);
  delete (dYL), free (dYL);
}

static double max_diff (scalar * al, scalar * bl)
{
  double m = 0.;
  for (int i = 0; i < NS; i++) {
    scalar a = al[i], b = bl[i];
    foreach (reduction(max:m))
      m = max (m, fabs (a[] - b[]));
  }
  return m;
}

/**
The largest $|\sum Y - 1|$, the smallest $Y$, and the largest change of the
mass of one species. */

static void check_state (scalar * YL, double * mass0,
                         double * sumerr, double * ymin, double * dmass)
{
  double se = 0., ym = HUGE;
  foreach (reduction(max:se) reduction(min:ym)) {
    double s = 0.;
    for (scalar Y in YL) {
      s += Y[];
      ym = min (ym, Y[]);
    }
    se = max (se, fabs (s - 1.));
  }
  *sumerr = se, *ymin = ym, *dmass = 0.;
  for (int i = 0; i < NS; i++) {
    scalar Y = YL[i];
    double m = statsf (Y).sum;
    *dmass = max (*dmass, fabs (m - mass0[i]));
  }
}

static void species_mass (scalar * YL, double * mass)
{
  for (int i = 0; i < NS; i++) {
    scalar Y = YL[i];
    mass[i] = statsf (Y).sum;
  }
}

/**
The implicit run to `TEND` with `nstep` steps. It returns the largest number
of multigrid cycles of one step. */

static int implicit_run (scalar * YL, scalar * DL, double * mw, bool upwind,
                         int nstep, double * ymin_run)
{
  int imax = 0;
  *ymin_run = HUGE;
  initial (YL);
  for (int s = 0; s < nstep; s++) {
    mgstats st = coupled_diffusion (YL, TEND/nstep, DL, alpha, thetac,
                                    MWs = mw, upwind = upwind,
                                    tolerance = 1e-9);
    assert (st.resa <= 1e-9);
    imax = max (imax, st.i);
    double ym = HUGE;
    foreach (reduction(min:ym))
      for (scalar Y in YL)
        ym = min (ym, Y[]);
    *ymin_run = min (*ymin_run, ym);
  }
  return imax;
}

int main()
{
  init_grid (128);
  foreach()
    thetac[] = cm[];
  foreach_face()
    alpha.x[] = fm.x[];

  scalar * YL = {Ya, Yb, Yc}, * DL = {Da, Db, Dc};
  scalar * RL = list_clone (YL);
  double mass0[NS], sumerr, ymin, dmass;
  initial (YL);
  species_mass (YL, mass0);

  /**
  ## 1. Equal diffusivities */

  set_diffusivities (DL, 0.3, 0.3, 0.3);
  for (int up = 0; up <= 1; up++) {
    initial (RL);
    for (int s = 0; s < 10; s++)
      for (int i = 0; i < NS; i++) {
        scalar R = RL[i], th[];
        face vector Df[];
        foreach()
          th[] = thetac[];
        foreach_face()
          Df.x[] = alpha.x[]*0.3;
        diffusion (R, TEND/10, D = Df, theta = th, tolerance = 1e-9);
      }
    double ymin_run;
    implicit_run (YL, DL, NULL, up, 10, &ymin_run);
    double d = max_diff (YL, RL);
    fprintf (stderr, "equal D, %s: max |coupled - diffusion()| = %.3e\n",
             up ? "upwind" : "central", d);
    assert (d < 1e-9);
  }

  /**
  ## 2. and 3. Convergence in `dt` against the explicit reference */

  set_diffusivities (DL, 1., 0.2, 0.05);
  for (int molar = 0; molar <= 1; molar++) {
    double * mw = molar ? MWs : NULL;
    initial (RL);
    double h = 1./128, dtr = 0.05*sq(h);
    int nref = (int) ceil (TEND/dtr);
    for (int s = 0; s < nref; s++)
      explicit_step (RL, DL, mw, TEND/nref);

    double err[4];
    int nsteps[4] = {5, 10, 20, 40};
    for (int j = 0; j < 4; j++) {
      double ymin_run;
      int imax = implicit_run (YL, DL, mw, false, nsteps[j], &ymin_run);
      err[j] = max_diff (YL, RL);
      check_state (YL, mass0, &sumerr, &ymin, &dmass);
      fprintf (stderr, "%s, central, dt = T/%-2d: error %.4e, "
               "|sum Y - 1| %.2e, d(mass) %.2e, min Y %.4e, cycles %d\n",
               molar ? "molar" : "Fick ", nsteps[j], err[j], sumerr, dmass,
               ymin_run, imax);
      assert (sumerr < 1e-9 && dmass < 1e-10);
    }
    double order = log (err[2]/err[3])/log (2.);
    fprintf (stderr, "%s: observed order %.3f (T/20 to T/40)\n",
             molar ? "molar" : "Fick ", order);
    assert (order > 0.8 && order < 1.2);
  }

  /**
  ## 5. A large step with upwind face values

  The corrective flux of the initial state sets the Courant number of the
  explicit scheme. */

  {
    initial (YL);
    face vector un[];
    foreach_face() {
      un.x[] = 0.;
      for (int k = 0; k < NS; k++) {
        scalar Y = YL[k], D = DL[k];
        un.x[] += alpha.x[]*face_value (D, 0)*(Y[] - Y[-1])/Delta;
      }
    }
    double umax = 0.;
    foreach_face (reduction(max:umax))
      umax = max (umax, fabs (un.x[]));
    double dt = TEND, cfl = umax*dt*128;
    for (int molar = 0; molar <= 1; molar++) {
      initial (YL);
      mgstats st = coupled_diffusion (YL, dt, DL, alpha, thetac,
                                      MWs = molar ? MWs : NULL,
                                      upwind = true, tolerance = 1e-9);
      check_state (YL, mass0, &sumerr, &ymin, &dmass);
      fprintf (stderr, "%s, upwind, one step at corrective Courant %.2f: "
               "cycles %d, res %.2e, |sum Y - 1| %.2e, d(mass) %.2e, "
               "min Y %.4e\n", molar ? "molar" : "Fick ",
               cfl, st.i, st.resa, sumerr, dmass, ymin);
      assert (st.resa <= 1e-9 && sumerr < 1e-9 && dmass < 1e-10);
      assert (ymin >= 0.);
    }
  }

  delete (RL), free (RL);
}
