#!/usr/bin/env python3
"""
Aggregate + fit the Marangoni (Ma,Sc) sweep produced by marangoni-point.c.

Reads all `uptake-Ma*-Sc*.dat` and `baseline-Sc*.dat` in a results directory,
extracts kappa = D_eff/D per point as the SINGLE-MODE late-tail decay-rate ratio
(normalised by the discrete baseline -> cancels grid discretization and the VOF
startup transient), then fits the GLOBAL Lorentzian/Hill correlation
    g = Sc/(Sc+s0);  A = A0 - A1*g;  Mc = M0 - M1*g
    D_eff/D = 1 + A / (1 + (Mc/Ma)^2)
one parameter set (A0,A1,M0,M1,s0) across ALL Sc.  Low-Ma limit -> A*(Ma/Mc)^2
(the closed-streamline Pe^2 rise); saturates at 1+A.

The "1D radial solver with adjustable D" is exact for this Dirichlet-uptake
problem:  C_bar(tau; kappa) = 1 - (6/pi^2) sum_n exp(-n^2 pi^2 kappa tau)/n^2,
so kappa just rescales tau; the tail slope of ln(1-C_bar) is -kappa*pi^2.

Usage:  python3 marangoni-fit.py [results_dir]   (default: current dir)
Writes: deff.dat, deff_vs_Ma.svg, deff_quality.svg
"""
import sys, os, glob, re, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# quiet the expected mid-sweep noise (an in-flight file is empty/being written)
warnings.filterwarnings('ignore', message='.*input contained no data.*')

PI2 = np.pi**2
D = sys.argv[1] if len(sys.argv) > 1 else '.'

def load(fn):
    # ndmin=2 + try/except so a file still being written (header only, a single
    # row, or a torn last line) is skipped, not fatal -> safe to run mid-sweep.
    try:
        a = np.loadtxt(fn, ndmin=2)
    except Exception:
        return None
    if a.shape[0] < 4 or a.shape[1] < 4:
        return None
    return a[:, 0], a[:, 1], a[:, 3]          # tau, Cbar, maxu/U

def header_vals(fn):
    with open(fn) as f:
        h = f.readline()
    ca = re.search(r'Ca=([0-9.eE+-]+)', h)
    lv = re.search(r'level=([0-9]+)', h)
    return (float(ca.group(1)) if ca else float('nan'),
            int(lv.group(1)) if lv else 0)

def series_mean(tau, kappa, N=60):
    n = np.arange(1, N + 1)[:, None]
    return 1. - (6./PI2)*np.sum(np.exp(-n**2*PI2*kappa*tau[None, :])/n**2, axis=0)

def tail_slope(tau, Cbar, lo=0.6, hi=0.99):
    m = (Cbar > lo) & (Cbar < hi) & (1. - Cbar > 0)
    if m.sum() < 4:
        return np.nan
    return np.polyfit(tau[m], np.log(1. - Cbar[m]), 1)[0]

def series_residual(tau, Cbar, kappa, lo=0.05, hi=0.99):
    m = (Cbar > lo) & (Cbar < hi)
    if m.sum() < 4:
        return np.nan
    return np.sqrt(np.mean((series_mean(tau[m], kappa) - Cbar[m])**2))

def deff_model(p, Ma, Sc):
    """Global Lorentzian/Hill correction (one param set for all Sc):
         g = Sc/(Sc+s0);  A = A0 - A1*g;  Mc = M0 - M1*g
         D_eff/D = 1 + A / (1 + (Mc/Ma)^2)."""
    A0, A1, M0, M1, s0 = p
    g = Sc/(Sc + s0)
    return 1. + (A0 - A1*g)/(1. + ((M0 - M1*g)/Ma)**2)

def _nelder_mead(f, x0, step=0.05, it=8000):
    """Minimal numpy-only Nelder-Mead (no scipy dependency)."""
    n = len(x0)
    simp = np.array([x0.copy() for _ in range(n + 1)], float)
    for i in range(n):
        simp[i + 1][i] = x0[i]*(1. + step) if x0[i] else step
    fv = np.array([f(s) for s in simp])
    for _ in range(it):
        o = np.argsort(fv); simp, fv = simp[o], fv[o]
        c = simp[:-1].mean(0)
        xr = c + (c - simp[-1]); fr = f(xr)
        if fr < fv[0]:
            xe = c + 2.*(c - simp[-1]); fe = f(xe)
            simp[-1], fv[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < fv[-2]:
            simp[-1], fv[-1] = xr, fr
        else:
            xc = c + 0.5*(simp[-1] - c); fc = f(xc)
            if fc < fv[-1]:
                simp[-1], fv[-1] = xc, fc
            else:
                simp[1:] = simp[0] + 0.5*(simp[1:] - simp[0])
                fv[1:] = np.array([f(s) for s in simp[1:]])
    j = np.argmin(fv)
    return simp[j], fv[j]

def fit_lorentzian(Ma, Sc, k, p0=(4.03, 0.98, 768., 321., 7.75)):
    """GLOBAL fit of deff_model over all Sc at once.  numpy-only.
    Returns (params, rms, n_used) or None. Seeded at the campaign values."""
    m = np.isfinite(Ma) & np.isfinite(Sc) & np.isfinite(k)
    if m.sum() < 5:                                    # 5 free params
        return None
    Mm, Sm, km = Ma[m], Sc[m], k[m]
    def cost(p):
        if p[4] <= 0:                                  # s0 must be positive
            return 1e9
        return np.sqrt(np.mean((deff_model(p, Mm, Sm) - km)**2))
    p, rms = _nelder_mead(cost, np.array(p0, float))
    return p, rms, int(m.sum())

# --- discrete baseline tail slope (per Sc; should be ~ -pi^2) -----------------
base = {}
for fn in glob.glob(os.path.join(D, 'baseline-Sc*.dat')):
    sc = float(re.search(r'Sc([0-9.eE+-]+)\.dat', fn).group(1))
    r = load(fn)
    if r is None:
        continue
    tau, Cbar, _ = r
    base[sc] = tail_slope(tau, Cbar)
    print(f"# baseline Sc={sc:g}: tail slope={base[sc]:.4f} (analytic {-PI2:.4f})")
# fall back to a single baseline (Sc-independent in tau) if a given Sc is missing
base_any = np.mean([v for v in base.values() if np.isfinite(v)]) if base else -PI2

# --- per-point kappa ----------------------------------------------------------
rows = []
for fn in sorted(glob.glob(os.path.join(D, 'uptake-Ma*-Sc*.dat'))):
    g = re.search(r'Ma([0-9.eE+-]+)-Sc([0-9.eE+-]+)\.dat', fn)
    Ma, Sc = float(g.group(1)), float(g.group(2))
    r = load(fn)
    if r is None:
        continue
    tau, Cbar, maxu = r
    sm = tail_slope(tau, Cbar)
    sb = base.get(Sc, base_any)
    kappa = sm/sb if np.isfinite(sm) else np.nan
    resid = series_residual(tau, Cbar, kappa) if np.isfinite(kappa) else np.nan
    tail = maxu[int(0.7*len(maxu)):]
    steady = 1 if (tail.mean() > 0 and tail.std()/tail.mean() < 0.05) else 0
    ca, lvl = header_vals(fn)
    rows.append((Ma, Sc, ca, lvl, kappa, -sm/PI2, tail.mean(), steady, resid))

rows.sort()
with open(os.path.join(D, 'deff.dat'), 'w') as f:
    f.write("# Ma Sc Ca level Deff kappa_raw maxuU steady residual\n")
    for r in rows:
        f.write("%g %g %g %d %.4f %.4f %.4f %d %.2e\n" % r)
        print("Ma=%-7g Sc=%-5g D_eff/D=%.3f  steady=%d  resid=%.1e"
              % (r[0], r[1], r[4], r[7], r[8]))

# --- GLOBAL Lorentzian fit (one param set for all Sc), overlaid per Sc --------
if not rows:
    print("# no completed uptake-*.dat yet (sweep still warming up) -- "
          "rerun once some points finish.")
    sys.exit(0)
data = np.array([[r[0], r[1], r[4], r[7]] for r in rows])  # Ma Sc kappa steady
good = (data[:, 3] == 1) & np.isfinite(data[:, 2])         # steady points only
fit = fit_lorentzian(data[good, 0], data[good, 1], data[good, 2])
if fit:
    p, rms, nused = fit
    A0, A1, M0, M1, s0 = p
    print(f"# GLOBAL fit (n={nused}, rms={rms:.4f}):")
    print(f"#   g = Sc/(Sc+{s0:.3f});  A = {A0:.3f} - {A1:.3f}*g;  "
          f"Mc = {M0:.1f} - {M1:.1f}*g")
    print(f"#   D_eff/D = 1 + A/(1 + (Mc/Ma)^2)")
    for sc in sorted(set(data[good, 1])):
        m = data[good, 1] == sc
        r = np.sqrt(np.mean((deff_model(p, data[good][m, 0], sc)
                             - data[good][m, 2])**2))
        print(f"#   Sc={sc:<5g} rms={r:.4f}")

plt.figure(figsize=(6, 4.4))
for sc in sorted(set(data[:, 1])):
    s = data[data[:, 1] == sc]
    Ma, k, st = s[:, 0], s[:, 2], s[:, 3]
    mk = ['o' if x else 'x' for x in st]
    for i in range(len(Ma)):
        plt.plot(Ma[i], k[i], mk[i], color=f'C{int(np.log10(sc))%10}')
    plt.plot([], [], 'o', color=f'C{int(np.log10(sc))%10}', label=f'Sc={sc:g}')
    if fit:                                            # overlay the global model
        xx = np.logspace(np.log10(Ma.min()), np.log10(Ma.max()), 80)
        plt.plot(xx, deff_model(p, xx, sc), '-', lw=1,
                 color=f'C{int(np.log10(sc))%10}')
plt.xscale('log'); plt.yscale('log')
plt.xlabel('Ma'); plt.ylabel('D_eff / D')
plt.title('o = steady,  x = unsteady (excluded from fit)', fontsize=9)
plt.grid(True, which='both', ls=':'); plt.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(D, 'deff_vs_Ma.svg'))

# --- fit-quality / single-D_eff adequacy: residual & angular nonuniformity ----
plt.figure(figsize=(6, 4.4))
nang = 0
for fn in sorted(glob.glob(os.path.join(D, 'angnu-Ma*-Sc*.dat'))):
    g = re.search(r'Ma([0-9.eE+-]+)-Sc([0-9.eE+-]+)\.dat', fn)
    try:
        a = np.loadtxt(fn, ndmin=2)
    except Exception:
        continue
    if len(a):
        plt.plot(a[:, 0], a[:, 1], '-o', ms=3,
                 label=f"Ma{g.group(1)} Sc{g.group(2)}")
        nang += 1
plt.xlabel('tau'); plt.ylabel('angular nonuniformity  rms(C - <C>_r)')
plt.grid(True, ls=':')
if nang:
    plt.legend(fontsize=7, ncol=2)
plt.tight_layout(); plt.savefig(os.path.join(D, 'deff_quality.svg'))
print("# wrote deff.dat, deff_vs_Ma.svg, deff_quality.svg in", os.path.abspath(D))
