#!/usr/bin/env python3
"""
Compare the profiles of test/diffusion-paths.c on a grid ladder.

Usage: diffusion-paths.py <dir> [levels]

<dir> holds r-c<path>-L<level>/prof, with path 0 (the old path) and 1
(COUPLED_DIFFUSION). The finest level is the reference of each path. The
script averages a fine profile onto the cells of a coarse level (2 cells to
1), and gives the mean of |difference| over x for each quantity:

1. the error of each level against the finest level of the same path, and
   the observed order between two levels;
2. the difference of the two paths at the same level;
3. the error of each level against the finest level of the OTHER path.
"""
import sys, os
import numpy as np

QUANT = ["H2", "TAR", "CO2", "N2", "u"]


def load(d, c, L):
    p = os.path.join(d, f"r-c{c}-L{L}", "prof")
    if not os.path.exists(p):
        return None
    head = open(p).readline().split()[1:]
    a = np.loadtxt(p)
    a = a[np.argsort(a[:, 0])]
    return {k: a[:, i] for i, k in enumerate(head)}


def restrict(v, n):
    """Average v onto n cells (2^k fine cells to 1)."""
    r = len(v)//n
    return v.reshape(n, r).mean(1)


def err(a, b):
    """Mean |a - b| for each quantity, b restricted onto a."""
    n = len(a["x"])
    return [np.mean(np.abs(a[q] - restrict(b[q], n))) for q in QUANT]


def main():
    d = sys.argv[1]
    levels = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 \
        else list(range(5, 11))
    P = {c: {L: load(d, c, L) for L in levels} for c in (0, 1)}
    for c in (0, 1):
        levels_c = [L for L in levels if P[c][L] is not None]
        P[c] = {L: P[c][L] for L in levels_c}
    fine = {c: max(P[c]) for c in (0, 1)}
    fmt = lambda v: " ".join(f"{x:10.3e}" for x in v)
    hdr = "level " + " ".join(f"{q:>10s}" for q in QUANT)

    for c in (0, 1):
        R = P[c][fine[c]]
        print(f"\n1. path {c}: error against its own level {fine[c]}")
        print(hdr)
        prev = None
        for L in sorted(P[c]):
            if L == fine[c]:
                continue
            e = err(P[c][L], R)
            line = f"{L:5d} " + fmt(e)
            if prev is not None:
                line += "   order " + " ".join(
                    f"{np.log2(p/x):5.2f}" if x > 0 and p > 0 else "  nan"
                    for p, x in zip(prev, e))
            print(line)
            prev = e

    print("\n2. |path 0 - path 1| at the same level")
    print(hdr)
    for L in sorted(set(P[0]) & set(P[1])):
        print(f"{L:5d} " + fmt(err(P[0][L], P[1][L])))

    for c in (0, 1):
        o = 1 - c
        print(f"\n3. path {c}: error against level {fine[o]} of path {o}")
        print(hdr)
        for L in sorted(P[c]):
            print(f"{L:5d} " + fmt(err(P[c][L], P[o][fine[o]])))

    print("\nscale: max of each quantity at the finest level of path 0")
    print("      " + fmt([np.max(np.abs(P[0][fine[0]][q])) for q in QUANT]))


if __name__ == "__main__":
    main()
