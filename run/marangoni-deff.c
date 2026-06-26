/**
# Marangoni → effective-diffusivity correction (single-case baseline)

Transient surface-driven uptake/heating of a drop: the interface holds the
tracer at `C = 1` (Dirichlet), the interior starts uniform at `C = 0`, and the
front propagates inward by diffusion + (optionally) a steady Marangoni
recirculation. The diffusivity multiplier `D_eff/D` is read off from how much
faster the drop fills when the toroidal cell is present.

This file is the parameter-clean, validated baseline derived from the
`recirculation.c` prototype. Dimensionless control is `(Ma, Sc, Ca)`; physical
units use `rho_d = mu_d = 1` and a fixed drop radius `R0`, so:

  D      = lambda = 1/Sc            (Sc = mu_d/(rho_d D))
  Dsigma = Ma/(R0*Sc)              (Ma = Dsigma*R0/(mu_d D))
  sigma0 = Dsigma/Ca               (Ca = Dsigma/sigma0, kept << 1 and fixed)
  U      = Dsigma                  (= Dsigma/mu_d), Re = Ma/Sc = Dsigma*R0

The drop is phase 1 (`f = 1` inside the circle). Matched properties this step:
`rho1 = rho2 = 1`, `mu1 = mu2 = 1`. The sigma-driver is decoupled (a fixed
imposed field, not coupled to the tracer) so the flow is steady.
*/

#define F_ERR 1e-10

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "integral.h"
#include "gradients.h"
#include "diffusion.h"
#include "view.h"

scalar sigmav[], tr[];

int maxlevel = 6;
double R0 = 0.3, tbc = 1.;

/** Dimensionless control (overridable from the command line). */
double Ma = 20., Sc = 1., Ca = 0.05;
bool marangoni = true;

/** Derived quantities (set in main). */
double lambda, dsigma, sigma0, U, tEnd;

int main (int argc, char ** argv) {

  /** `Ma Sc Ca marangoni(0/1) maxlevel`. */
  if (argc > 1) Ma        = atof (argv[1]);
  if (argc > 2) Sc        = atof (argv[2]);
  if (argc > 3) Ca        = atof (argv[3]);
  if (argc > 4) marangoni = atoi (argv[4]);
  if (argc > 5) maxlevel  = atoi (argv[5]);

  /** Dimensionless mapping (rho_d = mu_d = 1). `sigma0` is taken from the
  case's Dsigma so that a baseline run (marangoni off) keeps the *same*
  capillary stiffness — only the tangential driver is switched off. */
  lambda = 1./Sc;
  dsigma = Ma/(R0*Sc);
  sigma0 = dsigma/Ca;
  U      = dsigma;
  tEnd   = R0*R0/lambda;            // ~1 diffusive time

  rho1 = rho2 = 1.;
  mu1  = mu2  = 1.;

  d.sigmaf = sigmav;
  f.tracers = {tr};

  init_grid (1 << maxlevel);
  assert (R0 < L0);

  fprintf (stderr, "# Ma=%g Sc=%g Ca=%g marangoni=%d level=%d "
           "| D=%g Dsigma=%g sigma0=%g U=%g Re=%g tEnd=%g\n",
           Ma, Sc, Ca, marangoni, maxlevel,
           lambda, dsigma, sigma0, U, Ma/Sc, tEnd);
  fprintf (stderr, "# 1:t 2:tau 3:Cbar 4:1-Cbar 5:ln(1-Cbar) "
           "6:maxu 7:maxu/U 8:ke 9:Ccenter\n");

  run();
}

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

event init (i = 0) {
  foreach()
    d[] = circle (x, y, R0);

  fraction (f, circle (x, y, R0));

  /** Uniform interior IC: C = 0 (matches the analytic Dirichlet series). */
  foreach()
    tr[] = 0.;
}

/** Decoupled surface-tension driver. The gradient is a function of `y`
(radial), the only form compatible with the `x = 0` symmetry plane of the
quarter domain. `marangoni = false` keeps `sigma0` (drop stays spherical) but
removes the tangential driver -> pure-diffusion / spurious-current baseline. */
event properties (i++) {
  foreach()
    sigmav[] = sigma0 + (marangoni ? dsigma*(y/R0) : 0.);
}

scalar sexp[];
event tracer_diffusion (i++) {
  foreach() {
    f[] = (f[] > F_ERR) ? clamp (f[], 0., 1.) : 0.;
    tr[] = (f[] > F_ERR) ? tr[]/f[] : 0.;
  }

  scalar theta[];
  foreach()
    theta[] = max (cm[]*f[], 1e-10);

  face vector fs[];
  face_fraction (f, fs);

  face vector D[];
  foreach_face()
    D.x[] = lambda*fs.x[]*fm.x[];

  foreach()
    sexp[] = 0.;

  foreach()
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      coord m = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], m);
      double area = plane_area_center (m, alpha, &p);
#if AXI
      double dirac = area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      double dirac = area/Delta*cm[];
#endif
      sexp[] = lambda*plic_gradient (point, tr, f, fs, tbc, true, NULL)*dirac;
    }

  diffusion (tr, dt, D = D, theta = theta, r = sexp);

  foreach()
    tr[] *= f[];
}

/**
## Diagnostics

Uptake history, spurious-current / steadiness monitor, and the mean-decay data
used to validate against the Dirichlet sphere eigenvalue `lambda*R = pi`:
`1 - Cbar(tau) ~ (6/pi^2) exp(-pi^2 tau)`, `tau = t*D/R^2`. */
event logfile (i += 5) {
  double sum = 0, vol = 0, maxu = 0, ke = 0;
  foreach (reduction(+:sum) reduction(+:vol) reduction(max:maxu)
           reduction(+:ke)) {
    if (f[] > F_ERR) {
      double C = tr[]/f[];
      sum += C*dv();
      vol += dv();
      double umag = sqrt (sq(u.x[]) + sq(u.y[]));
      if (umag > maxu) maxu = umag;
      ke += 0.5*(sq(u.x[]) + sq(u.y[]))*dv();
    }
  }
  double Cbar = sum/vol;
  double tau  = t*lambda/(R0*R0);
  double Cc   = interpolate (tr, 0, 0)/max (interpolate (f, 0, 0), F_ERR);
  double lg   = (1. - Cbar > 0.) ? log (1. - Cbar) : nodata;

  // t  tau  Cbar  1-Cbar  ln(1-Cbar)  maxu  maxu/U  ke  Ccenter
  fprintf (stderr, "%g %g %g %g %g %g %g %g %g\n",
           t, tau, Cbar, 1. - Cbar, lg, maxu, maxu/U, ke, Cc);
}

/** Radial mean profile C(r) at successive tau -> feeds the later 1D-solver fit.
Scheduled by a tau threshold inside an `i++` event: Basilisk only accepts a
compile-time constant as the `t += ` period, so a runtime `tEnd/5` never arms. */
event profiles (i++) {
  static double tau_next = 0.;
  double tau = t*lambda/(R0*R0);
  if (tau < tau_next)
    return 0;
  tau_next += 0.2;

  char name[80];
  sprintf (name, "profile-tau%.2f.dat", tau);
  FILE * fp = fopen (name, "w");
  int nb = 24;
  for (int b = 0; b < nb; b++) {
    double ra = b*R0/nb, rb = (b + 1.)*R0/nb, s = 0., c = 0.;
    foreach (reduction(+:s) reduction(+:c))
      if (f[] > F_ERR) {
        double r = sqrt (x*x + y*y);
        if (r >= ra && r < rb) {
          s += tr[]/f[]*dv();
          c += dv();
        }
      }
    if (c > 0.)
      fprintf (fp, "%g %g\n", 0.5*(ra + rb)/R0, s/c);
  }
  fclose (fp);
}

event movie (i++) {
  static double tau_next = 0.;
  double tau = t*lambda/(R0*R0);
  if (tau < tau_next)
    return 0;
  tau_next += 0.05;

  scalar C[];
  foreach()
    C[] = (f[] > F_ERR) ? tr[]/f[] : 0.;

  view (fov = 18, tx = -0.15, ty = -0.15, width = 600, height = 600);
  clear();
  draw_vof ("f", lw = 2);
  squares ("C", min = 0., max = 1., linear = true);
  save ("uptake.mp4");

  clear();
  draw_vof ("f", lw = 2);
  squares ("u.x", linear = true);
  save ("velocity.mp4");
}

event end (t = tEnd);

/**
## Results (default run: Ma=20, Sc=1, Ca=0.05, Marangoni on)

Surface-driven uptake: the volume-mean and the drop-center concentration both
climb from 0 towards 1 as the front propagates inward.

~~~gnuplot Surface-driven uptake history
reset
set terminal svg size 500,400
set output 'uptake.svg'
set xlabel 'tau = t D / R^2'
set ylabel 'concentration'
set grid
set key bottom right
plot 'log' u 2:3 w l lw 2 t 'volume mean C', \
     'log' u 2:9 w l lw 2 t 'center C'
~~~

Pure-diffusion baseline validation: in the exponential regime the mean decays
as the Dirichlet sphere eigenvalue, ln(1 - C) = ln(6/pi^2) - pi^2 tau.

~~~gnuplot Mean-decay vs Dirichlet eigenvalue (slope = -pi^2)
reset
set terminal svg size 500,400
set output 'decay.svg'
set xlabel 'tau = t D / R^2'
set ylabel 'ln(1 - mean C)'
set grid
set xrange [0:0.6]
set key bottom left
plot 'log' u 2:5 w l lw 2 t 'VOF', \
     log(6./pi**2) - pi*pi*x w l dt 2 lc rgb 'black' \
     t 'analytic: ln(6/pi^2) - pi^2 tau'
~~~

Steady toroidal cell: the peak velocity (normalised by U = Dsigma/mu_d) rises
then plateaus, confirming a single steady recirculation.

~~~gnuplot Peak velocity vs time (steadiness)
reset
set terminal svg size 500,400
set output 'steadiness.svg'
set xlabel 'tau = t D / R^2'
set ylabel 'max|u| / U'
set grid
plot 'log' u 2:7 w l lw 2 notitle
~~~

Radial mean profiles C(r) at successive times -- the quantity matched against a
1D radial solver to calibrate D_eff.

~~~gnuplot Radial mean profiles at successive times
reset
set terminal svg size 500,400
set output 'profiles.svg'
set xlabel 'r / R0'
set ylabel 'mean C(r)'
set grid
set key top left
plot for [f in system('ls profile-tau*.dat')] f u 1:2 w lp t f
~~~
*/
