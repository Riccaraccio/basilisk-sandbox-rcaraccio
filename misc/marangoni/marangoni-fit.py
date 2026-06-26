#!/usr/bin/env python3
"""
Aggregate + fit the Marangoni (Ma,Sc) sweep produced by marangoni-point.c.

Reads all `uptake-Ma*-Sc*.dat` and `baseline-Sc*.dat` in a results directory,
extracts kappa = D_eff/D per point as the SINGLE-MODE late-tail decay-rate ratio
(normalised by the discrete baseline -> cancels grid discretization and the VOF
startup transient), then fits the correlation D_eff/D = 1 + C(Sc)*Ma^n per Sc.

The "1D radial solver with adjustable D" is exact for this Dirichlet-uptake
problem:  C_bar(tau; kappa) = 1 - (6/pi^2) sum_n exp(-n^2 pi^2 kappa tau)/n^2,
so kappa just rescales tau; the tail slope of ln(1-C_bar) is -kappa*pi^2.

Usage:  python3 marangoni-fit.py [results_dir]   (default: current dir)
Writes: deff.dat, deff_vs_Ma.svg, deff_quality.svg
"""
import sys, os, glob, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PI2 = np.pi**2
D = sys.argv[1] if len(sys.argv) > 1 else '.'

def load(fn):
    a = np.loadtxt(fn)
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

# --- discrete baseline tail slope (per Sc; should be ~ -pi^2) -----------------
base = {}
for fn in glob.glob(os.path.join(D, 'baseline-Sc*.dat')):
    sc = float(re.search(r'Sc([0-9.eE+-]+)\.dat', fn).group(1))
    tau, Cbar, _ = load(fn)
    base[sc] = tail_slope(tau, Cbar)
    print(f"# baseline Sc={sc:g}: tail slope={base[sc]:.4f} (analytic {-PI2:.4f})")
# fall back to a single baseline (Sc-independent in tau) if a given Sc is missing
base_any = np.mean([v for v in base.values() if np.isfinite(v)]) if base else -PI2

# --- per-point kappa ----------------------------------------------------------
rows = []
for fn in sorted(glob.glob(os.path.join(D, 'uptake-Ma*-Sc*.dat'))):
    g = re.search(r'Ma([0-9.eE+-]+)-Sc([0-9.eE+-]+)\.dat', fn)
    Ma, Sc = float(g.group(1)), float(g.group(2))
    tau, Cbar, maxu = load(fn)
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

# --- global fit D_eff/D = 1 + C(Sc) Ma^n, per Sc ------------------------------
data = np.array([[r[0], r[1], r[4], r[7]] for r in rows])  # Ma Sc kappa steady
plt.figure(figsize=(6, 4.4))
for sc in sorted(set(data[:, 1])):
    s = data[data[:, 1] == sc]
    Ma, k, st = s[:, 0], s[:, 2], s[:, 3]
    mk = ['o' if x else 'x' for x in st]
    for i in range(len(Ma)):
        plt.plot(Ma[i], k[i], mk[i], color=f'C{int(np.log10(sc))%10}')
    plt.plot([], [], 'o', color=f'C{int(np.log10(sc))%10}', label=f'Sc={sc:g}')
    # fit only steady points with measurable enhancement
    good = (st == 1) & (k > 1.0001) & np.isfinite(k)
    if good.sum() >= 2:
        n, logC = np.polyfit(np.log(Ma[good]), np.log(k[good] - 1.), 1)
        C = np.exp(logC)
        xx = np.logspace(np.log10(Ma.min()), np.log10(Ma.max()), 80)
        plt.plot(xx, 1. + C*xx**n, '-', lw=1, color=f'C{int(np.log10(sc))%10}')
        print(f"# fit Sc={sc:g}: D_eff/D = 1 + {C:.3e} * Ma^{n:.3f}")
plt.xscale('log'); plt.yscale('log')
plt.xlabel('Ma'); plt.ylabel('D_eff / D')
plt.title('o = steady,  x = unsteady (excluded from fit)', fontsize=9)
plt.grid(True, which='both', ls=':'); plt.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(D, 'deff_vs_Ma.svg'))

# --- fit-quality / single-D_eff adequacy: residual & angular nonuniformity ----
plt.figure(figsize=(6, 4.4))
for fn in sorted(glob.glob(os.path.join(D, 'angnu-Ma*-Sc*.dat'))):
    g = re.search(r'Ma([0-9.eE+-]+)-Sc([0-9.eE+-]+)\.dat', fn)
    try:
        a = np.loadtxt(fn)
    except Exception:
        continue
    if a.ndim == 2 and len(a):
        plt.plot(a[:, 0], a[:, 1], '-o', ms=3,
                 label=f"Ma{g.group(1)} Sc{g.group(2)}")
plt.xlabel('tau'); plt.ylabel('angular nonuniformity  rms(C - <C>_r)')
plt.grid(True, ls=':'); plt.legend(fontsize=7, ncol=2)
plt.tight_layout(); plt.savefig(os.path.join(D, 'deff_quality.svg'))
print("# wrote deff.dat, deff_vs_Ma.svg, deff_quality.svg in", os.path.abspath(D))
