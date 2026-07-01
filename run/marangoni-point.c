/**
# Marangoni → D_eff: single-point solver for cluster job arrays

One `(Ma, Sc, Ca, marangoni, level)` point per invocation, driven from the
command line, so a scheduler job array can fan the `(Ma,Sc)` grid across nodes.
Physics is identical to the validated `marangoni-deff.c` / `marangoni-sweep.c`
(Dirichlet surface-driven uptake + decoupled steady Marangoni cell); this file
just strips the in-binary loop and the `view.h` movies (no GL on compute nodes),
adds MPI-safe master-only file writes, and an early-stop at saturation so
high-enhancement / high-Sc points don't integrate a full (long) diffusive time.

  ./marangoni-point Ma Sc Ca marangoni(0/1) cells_per_bl [taumax]

`cells_per_bl` is the number of grid points to place across the interfacial
concentration boundary layer; `maxlevel` is derived from it (see main()) rather
than being set directly. With this nondimensionalization Pe = U*R0/D = Ma, so
delta_c ~ R0/sqrt(Ma) and level = ceil(log2(L0*cells_per_bl*sqrt(Ma)/R0)),
clamped to [level_min, level_max].

Outputs into the cwd: `uptake-Ma%g-Sc%g.dat` (or `baseline-Sc%g.dat`),
`profile-...-tau%.2f.dat`, `angnu-...dat`. Aggregate + fit with marangoni-fit.py.

Build (serial):  make marangoni-point.tst
Build (MPI):     CC='mpicc -D_MPI=1' make marangoni-point.tst   (then mpirun)
or directly:     qcc -D_MPI=1 -O2 -disable-dimensions -I../src \
                   -I$OPENSMOKE_INTERFACE/src \
                   -I$HOME/basilisk/basilisk-sandbox-ecipriano/src \
                   marangoni-point.c -o marangoni-point -lm

Movies (opt-in): compile with -DMOVIE and link the headless software GL
backend (fb_tiny -> no OpenGL/OSMesa needed); needs ffmpeg in PATH at run time.
  qcc -O2 -disable-dimensions -DMOVIE [-DMOVIE_NFRAMES=50] \
      -I$HOME/basilisk/basilisk-sandbox-ecipriano/src \
      marangoni-point.c -o marangoni-point \
      -L$BASILISK/gl -lglutils -lfb_tiny -lm
Writes `movie-<tag>.mp4` into the cwd: a fixed MOVIE_NFRAMES snapshots per run
regardless of (Ma,Sc) -- under MOVIE the saturation early-stop is disabled so
every case integrates the full taumax window. Movie runs are manual-only (the
sweep script builds without -DMOVIE); launch the points you want by hand.
*/

#define F_ERR 1e-10

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "integral.h"
#include "gradients.h"
#include "diffusion.h"
#include "adapt_wavelet_leave_interface.h"

#ifdef MOVIE
#include "view.h"
#ifndef MOVIE_NFRAMES
#define MOVIE_NFRAMES 50        // snapshots per run, independent of setup (override -D)
#endif
double movie_dt = 1.;           // frame interval in t; real value set in main() once
                                // tEnd is known. MUST be a positive placeholder here:
                                // Basilisk reads it at event registration (before main)
                                // to enable the periodic event; 0 would disable it.
#endif

scalar sigmav[], tr[], sexp[];

int runlevel = 7;                  // derived in main() from cells_bl (see below)
int level_min = 6, level_max = 13; // floor / cap on the derived level
double cells_bl = 2.;              // target grid points across the concentration BL
double R0 = 0.3, tbc = 1.;
double Ma = 100., Sc = 1., Ca = 0.05, taumax = 1.0;
bool marangoni = true;
double lambda, dsigma, sigma0, U, tEnd;

FILE * fp_uptake = NULL, * fp_angnu = NULL;
char tag[64];
double tau_next;

int main (int argc, char ** argv) {
  if (argc > 1) Ma        = atof (argv[1]);
  if (argc > 2) Sc        = atof (argv[2]);
  if (argc > 3) Ca        = atof (argv[3]);
  if (argc > 4) marangoni = atoi (argv[4]);
  if (argc > 5) cells_bl  = atof (argv[5]);
  if (argc > 6) taumax    = atof (argv[6]);

  lambda = 1./Sc;
  dsigma = Ma/(R0*Sc);
  sigma0 = dsigma/Ca;
  U      = dsigma;
  tEnd   = taumax*R0*R0/lambda;
#ifdef MOVIE
  movie_dt = tEnd/MOVIE_NFRAMES;         // fixed snapshot count over the whole run
#endif

  /* Resolve the interfacial concentration boundary layer. With this
     nondimensionalization Pe = U*R0/D = Ma, so delta_c ~ R0/sqrt(Ma). To put
     cells_bl grid points across it we need Delta = delta_c/cells_bl, i.e.
     level = ceil(log2(L0/Delta)) = ceil(log2(L0*cells_bl/delta_c)), clamped to
     [level_min, level_max]. delta_c uses the (thinner) plain Ma^-1/2 estimate,
     so the Sc^1/6 high-Sc correction never causes under-resolution. */
  double delta_c = R0/sqrt(Ma);
  double lvl = ceil (log2 (L0*cells_bl/delta_c));
  runlevel = (int) clamp (lvl, (double) level_min, (double) level_max);
  double cells_actual = delta_c*(1 << runlevel)/L0;

  rho1 = rho2 = 1.;
  mu1  = mu2  = 1.;
  d.sigmaf = sigmav;
  f.tracers = {tr};

  init_grid (1 << runlevel);

  if (pid() == 0) {
    fprintf (stderr, "# Ma=%g Sc=%g Ca=%g marangoni=%d level=%d taumax=%g | "
             "D=%g Dsigma=%g sigma0=%g U=%g Re=%g tEnd=%g | "
             "delta_c=%g cells_bl(target=%g actual=%g)\n",
             Ma, Sc, Ca, marangoni, runlevel, taumax,
             lambda, dsigma, sigma0, U, Ma/Sc, tEnd,
             delta_c, cells_bl, cells_actual);
    if (cells_actual < cells_bl - 0.5)
      fprintf (stderr, "# WARNING: BL under-resolved (level capped at %d): "
               "%g < %g cells across delta_c\n",
               level_max, cells_actual, cells_bl);
  }

  run();
}

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

event init (i = 0) {
  foreach()
    d[] = circle (x, y, R0);
  fraction (f, circle (x, y, R0));
  foreach()
    tr[] = 0.;

  tau_next = 0.;
  if (pid() == 0) {                         // master writes files
    char name[120];
    if (marangoni) {
      sprintf (tag, "Ma%g-Sc%g", Ma, Sc);
      sprintf (name, "uptake-%s.dat", tag); fp_uptake = fopen (name, "w");
      sprintf (name, "angnu-%s.dat", tag);  fp_angnu  = fopen (name, "w");
    }
    else {
      sprintf (tag, "Sc%g", Sc);
      sprintf (name, "baseline-Sc%g.dat", Sc); fp_uptake = fopen (name, "w");
      fp_angnu = NULL;
    }
    fprintf (fp_uptake, "# tau Cbar 1-Cbar maxu/U ke Ccore "
             "| Ma=%g Sc=%g Ca=%g level=%d\n", Ma, Sc, Ca, runlevel);
  }
  else
    sprintf (tag, marangoni ? "Ma%g-Sc%g" : "Sc%g",
             marangoni ? Ma : Sc, Sc);       // tag still needed for names below
}

event properties (i++) {
  foreach()
    sigmav[] = sigma0 + (marangoni ? dsigma*(y/R0) : 0.);
}

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

event adapt (i++) {
  adapt_wavelet_leave_interface ({tr}, {f},
    (double[]){1e-1}, runlevel, padding=2);
}

/** Uptake history + early-stop on saturation. The core concentration is a
reduction over the inner r < 0.15 R0 (MPI-safe; avoids interpolate at a point). */
event logfile (i += 5) {
  double sum = 0, vol = 0, maxu = 0, ke = 0, cs = 0, cv = 0;
  foreach (reduction(+:sum) reduction(+:vol) reduction(max:maxu)
           reduction(+:ke) reduction(+:cs) reduction(+:cv)) {
    if (f[] > F_ERR) {
      double C = tr[]/f[];
      sum += C*dv(); vol += dv();
      double u2 = sq(u.x[]) + sq(u.y[]);
      if (sqrt(u2) > maxu) maxu = sqrt(u2);
      ke += 0.5*u2*dv();
      if (sqrt(x*x + y*y) < 0.15*R0) { cs += C*dv(); cv += dv(); }
    }
  }
  double Cbar = sum/vol;
  double tau  = t*lambda/(R0*R0);
  double Ccore = cv > 0. ? cs/cv : 0.;

  if (pid() == 0 && fp_uptake)
    fprintf (fp_uptake, "%g %g %g %g %g %g\n",
             tau, Cbar, 1. - Cbar, maxu/U, ke, Ccore);

#ifndef MOVIE
  if (Cbar > 0.995 && tau > 0.02)           // saturated -> stop (collective)
    return 1;                               // (disabled under MOVIE: run the full
#endif                                      //  window so every case yields MOVIE_NFRAMES)
}

event profiles (i += 20) {
  if (!marangoni)
    return 0;
  double tau = t*lambda/(R0*R0);
  if (tau < tau_next)
    return 0;
  tau_next += 0.2;

  FILE * fp = NULL;
  if (pid() == 0) {
    char name[120];
    sprintf (name, "profile-%s-tau%.2f.dat", tag, tau);
    fp = fopen (name, "w");
  }
  int nb = 24;
  double gnum = 0., gden = 0.;
  for (int b = 0; b < nb; b++) {
    double ra = b*R0/nb, rb = (b + 1.)*R0/nb, s = 0., c = 0.;
    foreach (reduction(+:s) reduction(+:c))
      if (f[] > F_ERR) {
        double r = sqrt (x*x + y*y);
        if (r >= ra && r < rb) { s += tr[]/f[]*dv(); c += dv(); }
      }
    if (c <= 0.)
      continue;
    double meanb = s/c, ss = 0.;
    foreach (reduction(+:ss))
      if (f[] > F_ERR) {
        double r = sqrt (x*x + y*y);
        if (r >= ra && r < rb) { double C = tr[]/f[]; ss += sq(C - meanb)*dv(); }
      }
    if (pid() == 0 && fp)
      fprintf (fp, "%g %g %g\n", 0.5*(ra + rb)/R0, meanb, sqrt (ss/c));
    gnum += ss; gden += c;
  }
  if (pid() == 0 && fp) fclose (fp);
  if (pid() == 0 && fp_angnu)
    fprintf (fp_angnu, "%g %g\n", tau, gden > 0. ? sqrt (gnum/gden) : 0.);
}

event end (t = tEnd) {
  if (pid() == 0) {
    if (fp_uptake) { fclose (fp_uptake); fp_uptake = NULL; }
    if (fp_angnu)  { fclose (fp_angnu);  fp_angnu  = NULL; }
  }
}

#ifdef MOVIE
/** A fixed number of snapshots (MOVIE_NFRAMES) per run, independent of setup:
the interval is tEnd/MOVIE_NFRAMES, and the saturation early-stop is disabled
under MOVIE (see logfile) so every case integrates the full window. The period
must be a global set in main() -- Basilisk captures a `t += <expr>` increment at
event registration, before main() assigns tEnd, so an inline expression would
evaluate with tEnd=0 and never fire (the placeholder init enables the event).
Rendering uses the software fb_tiny backend. */
event movie (t += movie_dt) {
  clear();
  view(tx=-0.20, ty=-0.20, fov=8);
  draw_vof ("f", lw=2);
  squares ("tr", min=0, max=1);
  char name[120];
  sprintf (name, "movie-%s.mp4", tag);
  save(name);
}
#endif
