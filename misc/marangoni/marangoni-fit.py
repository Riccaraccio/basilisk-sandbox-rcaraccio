#!/usr/bin/env python3
"""
End-to-end Marangoni D_eff pipeline: raw (Ma,Sc) sweep -> deff.dat -> fit + plots.

Ma is the MARANGONI number (not Mach); in this nondimensionalization Ma == Pe.

Stage 1 (aggregate, from the cluster marangoni-fit.py):
  Read every `uptake-Ma*-Sc*.dat` / `baseline-Sc*.dat` in a results dir and get
  kappa = D_eff/D per point as the SINGLE-MODE late-tail decay-rate ratio,
  normalised by the discrete baseline (cancels grid discretization + the VOF
  startup transient). The 1D radial solver is exact for this Dirichlet uptake:
    C_bar(tau; kappa) = 1 - (6/pi^2) sum_n exp(-n^2 pi^2 kappa tau)/n^2,
  so kappa just rescales tau and the tail slope of ln(1-C_bar) is -kappa*pi^2.
  Writes deff.dat. (If no raw files are present, an existing deff.dat is read.)

Stage 2 (fit + plot, from deff_fit.py):
  Global monotonic-saturating model, one parameter set across all Sc; the high-Ma
  rollover is treated as a numerical/under-converged artifact (no decline term):
    g       = Sc / (Sc + s0)            saturating Schmidt factor in [0,1]
    A       = a0 + a1*g                 plateau amplitude (= max enhancement)
    Mc      = b0 + b1*g                 onset Marangoni/Peclet (Pe == Ma here)
    D_eff/D = 1 + A / (1 + (Mc/Ma)**2)
  Algebraically identical to a tanh/Hill sigmoid in log(Ma) with slope p=2
  (0.5*(1+tanh(ln(Ma/Mc))) == Ma**2/(Ma**2+Mc**2)); the Ma**2 low-Ma limb is the
  weak-convection (Rhines-Young) Pe**2 scaling. Fit uses steady points only.

Usage:  python3 marangoni-fit.py [results_dir_or_deff.dat]   (default: .)
Outputs (next to the data): deff.dat, deff_global_fit.png
"""
import sys
import os
import glob
import re
import warnings
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PI2 = np.pi ** 2
warnings.filterwarnings("ignore", message=".*input contained no data.*")

ARG = sys.argv[1] if len(sys.argv) > 1 else "."


# ===================== stage 1: raw sweep -> kappa ===========================
def load(fn):
    # ndmin=2 + try/except so a file still being written (header only, one row,
    # or a torn last line) is skipped, not fatal -> safe to run mid-sweep.
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
    ca = re.search(r"Ca=([0-9.eE+-]+)", h)
    lv = re.search(r"level=([0-9]+)", h)
    return (float(ca.group(1)) if ca else float("nan"),
            int(lv.group(1)) if lv else 0)


def series_mean(tau, kappa, N=60):
    n = np.arange(1, N + 1)[:, None]
    return 1. - (6. / PI2) * np.sum(np.exp(-n ** 2 * PI2 * kappa * tau[None, :]) / n ** 2, axis=0)


def tail_slope(tau, Cbar, lo=0.6, hi=0.99):
    m = (Cbar > lo) & (Cbar < hi) & (1. - Cbar > 0)
    if m.sum() < 4:
        return np.nan
    return np.polyfit(tau[m], np.log(1. - Cbar[m]), 1)[0]


def series_residual(tau, Cbar, kappa, lo=0.05, hi=0.99):
    m = (Cbar > lo) & (Cbar < hi)
    if m.sum() < 4:
        return np.nan
    return np.sqrt(np.mean((series_mean(tau[m], kappa) - Cbar[m]) ** 2))


def full_kappa(tau, Cbar, lo=0.05, hi=0.99):
    """WHOLE-PROFILE kappa: minimise the exact-series misfit in log(1-Cbar) over
    the FULL curve [lo,hi] (vs tail_slope's narrow late window) by golden section,
    so it differs from tail_slope ONLY in the window -> the tail-vs-whole spread
    isolates single-D_eff breakdown (a single kappa that fits the late tail but
    not the early/middle = two-timescale mixing)."""
    m = (Cbar > lo) & (Cbar < hi) & (1. - Cbar > 0)
    if m.sum() < 4:
        return np.nan
    t, lc = tau[m], np.log(1. - Cbar[m])
    def cost(k):
        return np.mean((np.log(1. - series_mean(t, k)) - lc) ** 2)
    a, b = 0.3, 9.0                                    # kappa bracket
    gr = (np.sqrt(5.) - 1.) / 2.
    x1, x2 = b - gr * (b - a), a + gr * (b - a)
    f1, f2 = cost(x1), cost(x2)
    for _ in range(50):
        if f1 < f2:
            b, x2, f2 = x2, x1, f1; x1 = b - gr * (b - a); f1 = cost(x1)
        else:
            a, x1, f1 = x1, x2, f2; x2 = a + gr * (b - a); f2 = cost(x2)
    return 0.5 * (a + b)


def aggregate(D):
    """Read the raw sweep in dir D, write deff.dat, return the row tuples."""
    base = {}
    for fn in glob.glob(os.path.join(D, "baseline-Sc*.dat")):
        sc = float(re.search(r"Sc([0-9.eE+-]+)\.dat", fn).group(1))
        r = load(fn)
        if r is None:
            continue
        tau, Cbar, _ = r
        base[sc] = tail_slope(tau, Cbar)
        print(f"# baseline Sc={sc:g}: tail slope={base[sc]:.4f} (analytic {-PI2:.4f})")
    base_any = np.mean([v for v in base.values() if np.isfinite(v)]) if base else -PI2

    rows = []
    for fn in sorted(glob.glob(os.path.join(D, "uptake-Ma*-Sc*.dat"))):
        g = re.search(r"Ma([0-9.eE+-]+)-Sc([0-9.eE+-]+)\.dat", fn)
        Ma, Sc = float(g.group(1)), float(g.group(2))
        r = load(fn)
        if r is None:
            continue
        tau, Cbar, maxu = r
        sm = tail_slope(tau, Cbar)
        sb = base.get(Sc, base_any)
        kappa = sm / sb if np.isfinite(sm) else np.nan
        resid = series_residual(tau, Cbar, kappa) if np.isfinite(kappa) else np.nan
        kfull = full_kappa(tau, Cbar)
        deff_whole = kfull * (-PI2) / sb if (np.isfinite(kfull) and np.isfinite(sb)
                                             and sb != 0) else np.nan
        div = (deff_whole / kappa - 1.) if (np.isfinite(deff_whole)
                                            and np.isfinite(kappa) and kappa != 0) else np.nan
        tail = maxu[int(0.7 * len(maxu)):]
        steady = 1 if (tail.mean() > 0 and tail.std() / tail.mean() < 0.05) else 0
        ca, lvl = header_vals(fn)
        rows.append((Ma, Sc, ca, lvl, kappa, -sm / PI2, tail.mean(), steady, resid,
                     deff_whole, div))
    rows.sort()

    out = os.path.join(D, "deff.dat")
    with open(out, "w") as f:
        f.write("# Ma Sc Ca level Deff kappa_raw maxuU steady residual Deff_whole div\n")
        for r in rows:
            f.write("%g %g %g %d %.4f %.4f %.4f %d %.2e %.4f %.3f\n" % r)
            print("Ma=%-7g Sc=%-5g D_eff/D=%.3f (whole %.3f, tail-vs-whole %+.0f%%)"
                  "  steady=%d  resid=%.1e"
                  % (r[0], r[1], r[4], r[9], 100 * r[10], r[7], r[8]))
    print("# wrote", os.path.abspath(out))
    return rows


# ===================== stage 2: model + fit ==================================
def model(X, a0, a1, b0, b1, s0):
    """D_eff/D as a function of (Ma, Sc) and the 5 constants."""
    Ma, Sc = X
    g = Sc / (Sc + s0)
    A = a0 + a1 * g
    Mc = b0 + b1 * g
    return 1.0 + A / (1.0 + (Mc / Ma) ** 2)


def get_data(arg):
    """Resolve `arg` to (Ma, Sc, Deff, steady, datadir). If it points at raw
    sweep output, aggregate it first; otherwise read an existing deff.dat."""
    if os.path.isdir(arg):
        datadir = arg
        if glob.glob(os.path.join(arg, "uptake-Ma*-Sc*.dat")):
            rows = aggregate(arg)
            if not rows:
                sys.exit("# no completed uptake-*.dat yet -- rerun once points finish.")
            d = np.array([[r[0], r[1], r[4], r[7]] for r in rows])
            return d[:, 0], d[:, 1], d[:, 2], d[:, 3].astype(int), datadir
        deff = os.path.join(arg, "deff.dat")
    else:
        deff, datadir = arg, os.path.dirname(arg) or "."
    if not os.path.exists(deff):
        sys.exit(f"# no raw sweep and no deff.dat at {arg!r}")
    a = np.loadtxt(deff)                  # Ma Sc Ca lvl Deff ... steady ...
    return a[:, 0], a[:, 1], a[:, 4], a[:, 7].astype(int), datadir


def main():
    Ma, Sc, D, steady, datadir = get_data(ARG)

    # ---- fit (steady, finite points only) -------------------------------
    fitm = (steady == 1) & np.isfinite(D)
    p0 = [4.0, -1.0, 770.0, -320.0, 8.0]
    bounds = ([0, -10, 1, -5000, 1e-2], [50, 10, 1e5, 5000, 1e3])
    x, _ = curve_fit(model, (Ma[fitm], Sc[fitm]), D[fitm],
                     p0=p0, bounds=bounds, maxfev=400000)
    r = D[fitm] - model((Ma[fitm], Sc[fitm]), *x)
    rmse = np.sqrt(np.mean(r ** 2))
    r2 = 1 - np.sum(r ** 2) / np.sum((D[fitm] - D[fitm].mean()) ** 2)

    a0, a1, b0, b1, s0 = x
    print("Fitted constants x = [a0, a1, b0, b1, s0]  (fit on %d steady points):" % fitm.sum())
    print("  ", np.array2string(x, precision=4, max_line_width=200))
    print(f"  g  = Sc/(Sc + {s0:.3g})")
    print(f"  A  = {a0:.3g} + ({a1:.3g})*g")
    print(f"  Mc = {b0:.3g} + ({b1:.3g})*g")
    print(f"  D_eff/D = 1 + A/(1 + (Mc/Ma)**2)")
    print(f"  R2 = {r2:.5f}   RMSE = {rmse:.4f}")
    print("Per-Sc RMSE:")
    for s in sorted(set(Sc)):
        m = (Sc == s) & fitm
        if not m.any():
            continue
        rr = D[m] - model((Ma[m], Sc[m]), *x)
        print(f"  Sc={int(s):<4d} RMSE={np.sqrt(np.mean(rr**2)):.4f}  max|r|={np.max(np.abs(rr)):.3f}")

    scs = sorted(set(Sc))
    cols = plt.cm.viridis(np.linspace(0, 0.85, len(scs)))
    Mg = np.logspace(np.log10(Ma.min() * 0.8), np.log10(Ma.max() * 1.1), 300)

    def scatter(ax, s, c, sz):
        """filled o = steady (in fit), x = unsteady (excluded)."""
        m = Sc == s
        st, un = m & (steady == 1), m & (steady == 0)
        ax.scatter(Ma[st], D[st], color=c, s=sz, zorder=3, edgecolor="k", linewidth=.4)
        ax.scatter(Ma[un], D[un], color=c, s=sz, zorder=3, marker="x")

    # ---- plot 1: global fit ---------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for s, c in zip(scs, cols):
        scatter(ax, s, c, 42)
        ax.plot(Mg, model((Mg, np.full_like(Mg, s)), *x), color=c, lw=1.8, label=f"Sc={int(s)}")
    ax.set_xscale("log")
    ax.set_xlabel("Ma")
    ax.set_ylabel(r"$D_{eff}/D$")
    ax.set_title("Global fit of D_eff/D (Ma, Sc)   (o = steady, x = unsteady)")
    ax.legend(title="Sc", frameon=False, loc="lower right")
    ax.grid(alpha=.3, which="both")
    ann = (f"g  = Sc/(Sc+{s0:.3g})\nA  = {a0:.3g} {a1:+.3g} g\n"
           f"Mc = {b0:.3g} {b1:+.3g} g\nD_eff/D = 1 + A/(1+(Mc/Ma)**2)\n"
           f"R2={r2:.3f}  RMSE={rmse:.3f}")
    ax.text(0.02, 0.97, ann, transform=ax.transAxes, va="top", fontsize=7.5,
            family="monospace", bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=.85))
    fig.tight_layout()
    fig.savefig(os.path.join(datadir, "deff_global_fit.png"), dpi=150)
    print("wrote deff_global_fit.png")


if __name__ == "__main__":
    main()
