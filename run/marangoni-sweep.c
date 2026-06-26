/**
# Marangoni → effective-diffusivity correction: (Ma, Sc) sweep + fit

In-binary parameter sweep over `(Ma, Sc)` at fixed `Ca`, reusing the physics
validated in `marangoni-deff.c` (Dirichlet surface-driven uptake of a passive
tracer in a drop, with a decoupled steady Marangoni cell). For each point we
record the volume-mean uptake `C̄(τ)`; a 1-parameter fit against the *exact*
analytic Dirichlet series then gives `κ = D_eff/D`, and a global fit yields the
correlation `D_eff/D = 1 + C(Sc)·Ma^n`.

Dimensionless mapping (rho_d = mu_d = 1, R0 fixed) — identical to marangoni-deff:

  D = 1/Sc,  Dsigma = Ma/(R0*Sc),  sigma0 = Dsigma/Ca,  U = Dsigma,  Re = Ma/Sc.

Sweep idiom follows `porous-cylinder.c`: a parameter list looped in main() with
repeated `run()` calls. All runs write per-case files into the case run-dir; the
`pythonplot` block at the bottom does the per-point and global fits.
*/

#define F_ERR 1e-10

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "integral.h"
#include "gradients.h"
#include "diffusion.h"

scalar sigmav[], tr[], sexp[];

int runlevel = 6;
double R0 = 0.3, tbc = 1.;

/** Control + derived (set per run by setmap). */
double Ma, Sc, Ca = 0.05;
bool marangoni = true;
double lambda, dsigma, sigma0, U, tEnd;

/** Per-run output state. */
FILE * fp_uptake = NULL, * fp_angnu = NULL;
char tag[64];
double tau_next;

/** Sweep grid (pilot; widen once the pipeline is validated). */
double MaList[] = {20., 100., 500.};
double ScList[] = {1., 10.};

void setmap (void) {
  lambda = 1./Sc;
  dsigma = Ma/(R0*Sc);
  sigma0 = dsigma/Ca;
  U      = dsigma;
  tEnd   = R0*R0/lambda;            // ~1 diffusive time (tau = 1)
}

int main (void) {

  rho1 = rho2 = 1.;
  mu1  = mu2  = 1.;
  d.sigmaf = sigmav;
  f.tracers = {tr};

  init_grid (1 << runlevel);

  size_t nMa = sizeof(MaList)/sizeof(MaList[0]);
  size_t nSc = sizeof(ScList)/sizeof(ScList[0]);

  for (size_t j = 0; j < nSc; j++) {
    Sc = ScList[j];

    /** One baseline per Sc (no Marangoni) for the discrete eigenvalue. Uses the
    smallest Ma so sigma0 -- hence the capillary dt -- stays cheap. */
    marangoni = false; Ma = MaList[0]; setmap();
    fprintf (stderr, "# RUN baseline Sc=%g (D=%g tEnd=%g)\n", Sc, lambda, tEnd);
    run();

    /** Marangoni points. */
    marangoni = true;
    for (size_t k = 0; k < nMa; k++) {
      Ma = MaList[k]; setmap();
      fprintf (stderr, "# RUN Ma=%g Sc=%g Ca=%g | Dsigma=%g sigma0=%g U=%g "
               "Re=%g tEnd=%g\n", Ma, Sc, Ca, dsigma, sigma0, U, Ma/Sc, tEnd);
      run();
    }
  }
}

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

event init (i = 0) {
  foreach()
    d[] = circle (x, y, R0);
  fraction (f, circle (x, y, R0));
  foreach()
    tr[] = 0.;

  /** Reset per-run state (statics persist across run() calls otherwise). */
  tau_next = 0.;
  char name[120];
  if (marangoni) {
    sprintf (tag, "Ma%g-Sc%g", Ma, Sc);
    sprintf (name, "uptake-%s.dat", tag); fp_uptake = fopen (name, "w");
    sprintf (name, "angnu-%s.dat", tag);  fp_angnu  = fopen (name, "w");
    fprintf (fp_uptake, "# tau Cbar 1-Cbar maxu/U ke Ccenter "
             "| Ma=%g Sc=%g Ca=%g level=%d\n", Ma, Sc, Ca, runlevel);
  }
  else {
    sprintf (tag, "Sc%g", Sc);
    sprintf (name, "baseline-Sc%g.dat", Sc); fp_uptake = fopen (name, "w");
    fp_angnu = NULL;
    fprintf (fp_uptake, "# tau Cbar 1-Cbar maxu/U ke Ccenter "
             "| baseline Sc=%g level=%d\n", Sc, runlevel);
  }
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

/** Uptake history -> per-case file (one fit point per Marangoni run). */
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
  fprintf (fp_uptake, "%g %g %g %g %g %g\n",
           tau, Cbar, 1. - Cbar, maxu/U, ke, Cc);
}

/** Radial profiles (mean + std per bin) and a global angular-nonuniformity
metric -> quantifies how non-radial the field is (i.e. whether a single radial
D_eff is even adequate). Marangoni runs only. */
event profiles (i++) {
  if (!marangoni)
    return 0;
  double tau = t*lambda/(R0*R0);
  if (tau < tau_next)
    return 0;
  tau_next += 0.2;

  char name[120];
  sprintf (name, "profile-%s-tau%.2f.dat", tag, tau);
  FILE * fp = fopen (name, "w");
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
    fprintf (fp, "%g %g %g\n", 0.5*(ra + rb)/R0, meanb, sqrt (ss/c));
    gnum += ss; gden += c;
  }
  fclose (fp);
  if (fp_angnu)
    fprintf (fp_angnu, "%g %g\n", tau, gden > 0. ? sqrt (gnum/gden) : 0.);
}

event end (t = tEnd) {
  if (fp_uptake) { fclose (fp_uptake); fp_uptake = NULL; }
  if (fp_angnu)  { fclose (fp_angnu);  fp_angnu  = NULL; }
}

/**
## Fit: per-point D_eff/D and the global correlation

The 1D radial Dirichlet-uptake reference is exact:
`C̄(τ;κ) = 1 − (6/π²) Σ_n (1/n²) exp(−n²π²·κ·τ)`, with `κ = D_eff/D` rescaling τ.
We fit κ per point (dominant-mode slope, refined on the full series), then fit
`D_eff/D = 1 + C(Sc)·Ma^n` per Sc.

~~~pythonplot D_eff/D vs Ma with fitted correlation
import glob, re, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PI2 = np.pi**2

def load(fn):
    d = np.loadtxt(fn)
    return d[:,0], d[:,1], d[:,3]   # tau, Cbar, maxu/U

def header_vals(fn):
    with open(fn) as f:
        h = f.readline()
    ca = re.search(r'Ca=([0-9.eE+-]+)', h)
    lv = re.search(r'level=([0-9]+)', h)
    return (float(ca.group(1)) if ca else float('nan'),
            int(lv.group(1)) if lv else 0)

def series_mean(tau, kappa, N=60):
    n = np.arange(1, N+1)[:,None]
    return 1. - (6./PI2)*np.sum(np.exp(-n**2*PI2*kappa*tau[None,:])/n**2, axis=0)

def tail_slope(tau, Cbar, lo=0.6, hi=0.99):
    # slope of ln(1-Cbar) vs tau in the SINGLE-MODE tail: above the late-time
    # discretization floor, and past the multi-mode VOF startup transient that
    # contaminates the early [0.1,0.9] window. This is the eigenvalue regime.
    m = (Cbar > lo) & (Cbar < hi) & (1.-Cbar > 0)
    if m.sum() < 4:
        return np.nan
    return np.polyfit(tau[m], np.log(1.-Cbar[m]), 1)[0]

def series_residual(tau, Cbar, kappa, lo=0.05, hi=0.99):
    # RMS misfit of the best single-kappa 1D series over the whole uptake; large
    # -> a single D_eff cannot capture the (two-stage) mixing.
    m = (Cbar > lo) & (Cbar < hi)
    return np.sqrt(np.mean((series_mean(tau[m], kappa) - Cbar[m])**2))

# Discrete baseline tail slope per Sc (~ -pi^2). Normalising the Marangoni slope
# by it cancels the grid discretization AND the startup transient shared by both
# runs, so the ratio is D_eff/D directly.
base = {}
for fn in glob.glob('baseline-Sc*.dat'):
    sc = float(re.search(r'Sc([0-9.eE+-]+)\.dat', fn).group(1))
    tau, Cbar, _ = load(fn)
    base[sc] = tail_slope(tau, Cbar)
    print(f"# baseline Sc={sc}: tail slope={base[sc]:.4f} "
          f"(analytic -pi^2={-PI2:.4f})")

rows = []
for fn in sorted(glob.glob('uptake-Ma*-Sc*.dat')):
    g = re.search(r'Ma([0-9.eE+-]+)-Sc([0-9.eE+-]+)\.dat', fn)
    Ma, Sc = float(g.group(1)), float(g.group(2))
    tau, Cbar, maxu = load(fn)
    sm = tail_slope(tau, Cbar)
    kappa = sm/base.get(Sc, -PI2)              # = D_eff/D
    resid = series_residual(tau, Cbar, kappa)
    # steadiness: maxu/U plateau over the last 30% of the record
    tail = maxu[int(0.7*len(maxu)):]
    steady = 1 if (tail.mean() > 0 and tail.std()/tail.mean() < 0.05) else 0
    ca, lvl = header_vals(fn)
    rows.append((Ma, Sc, ca, lvl, kappa, -sm/PI2, tail.mean(), steady, resid))

rows.sort()
with open('deff.dat','w') as f:
    f.write("# Ma Sc Ca level Deff kappa_raw maxuU steady residual\n")
    for r in rows:
        f.write("%g %g %g %d %.4f %.4f %.4f %d %.2e\n" % r)
        print("Ma=%g Sc=%g  D_eff/D=%.3f  steady=%d  resid=%.1e"
              % (r[0], r[1], r[4], r[7], r[8]))

# global fit kappa-1 = C*Ma^n per Sc, plus the plot
data = np.array([[r[0],r[1],r[4]] for r in rows])
plt.figure(figsize=(5.5,4.2))
for sc in sorted(set(data[:,1])):
    s = data[data[:,1]==sc]
    Ma, k = s[:,0], s[:,2]
    plt.plot(Ma, k, 'o', label=f'Sc={sc:g}')
    pos = k>1.0001
    if pos.sum() >= 2:
        n, logC = np.polyfit(np.log(Ma[pos]), np.log(k[pos]-1.), 1)
        C = np.exp(logC)
        xx = np.logspace(np.log10(Ma.min()), np.log10(Ma.max()), 50)
        plt.plot(xx, 1.+C*xx**n, '-', lw=1,
                 label=f'  1+{C:.2e}·Ma^{n:.2f}')
        print(f"# global fit Sc={sc:g}: D_eff/D = 1 + {C:.3e} * Ma^{n:.3f}")
plt.xscale('log'); plt.xlabel('Ma'); plt.ylabel('D_eff / D')
plt.grid(True, which='both', ls=':'); plt.legend(fontsize=8)
plt.tight_layout(); plt.savefig('deff_vs_Ma.svg')

# fit-quality / validity: residual & angular-nonuniformity vs Ma
plt.figure(figsize=(5.5,4.2))
for fn in sorted(glob.glob('angnu-Ma*-Sc*.dat')):
    g = re.search(r'Ma([0-9.eE+-]+)-Sc([0-9.eE+-]+)\.dat', fn)
    try:
        a = np.loadtxt(fn)
    except Exception:
        continue
    if a.ndim == 2 and len(a):
        plt.plot(a[:,0], a[:,1], '-o', ms=3, label=f"Ma{g.group(1)} Sc{g.group(2)}")
plt.xlabel('tau'); plt.ylabel('angular nonuniformity  rms(C - <C>_r)')
plt.grid(True, ls=':'); plt.legend(fontsize=7)
plt.tight_layout(); plt.savefig('deff_quality.svg')
~~~
*/
