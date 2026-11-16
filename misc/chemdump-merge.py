#!/usr/bin/env python3
"""Join the dump files of all the MPI ranks (`CHEMISTRY_DUMP`, chemistry.h).

    python3 chemdump-merge.py <dir> [<out-dir>]

The script finds each set `chemdump-<iter>-<half>-<pid>.cells`, `.meta` and
`.end` in <dir>. It writes one file `chemdump-<iter>-<half>.cells`, `.meta`
and `.end` for each set into <out-dir> (default: <dir>). The cells keep the
order of the ranks. The `.cells` and `.end` files have the format of
`export.py`, so `make_reference --cells` of cuSMOKE reads them; it finds the
`.meta` and `.end` files next to the `.cells` file. A dump made before the
`.end` file existed has none: the script then writes none.

The script then prints, for each set, the number of cells, the number of
cells that the gate skipped, the CPU time of the solves and the spread of
the evaluations of the right-hand side, and the cells whose solve
diverged (non-finite end state).
"""
import glob
import os
import re
import struct
import sys

import numpy as np


def read_cells(path):
    with open(path, "rb") as fp:
        n, nsp = struct.unpack("<qi", fp.read(12))
        names = fp.read(32 * nsp)
        rows = np.fromfile(fp, "<f8")
    if rows.size != n * (4 + nsp):
        sys.exit(f"{path}: {rows.size} doubles, expected {n * (4 + nsp)}")
    return nsp, names, rows.reshape(n, 4 + nsp)


def main():
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    src = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) == 3 else src
    os.makedirs(out, exist_ok=True)

    sets = {}
    for path in glob.glob(os.path.join(src, "chemdump-*-*-*.cells")):
        m = re.search(r"chemdump-(\d+)-(\d)-(\d+)\.cells$", path)
        if m:
            sets.setdefault((int(m[1]), int(m[2])), []).append((int(m[3]), path))
    if not sets:
        sys.exit(f"no chemdump-<iter>-<half>-<pid>.cells in {src}")

    for (it, half), files in sorted(sets.items()):
        files.sort()
        nsp0, names0, cells, ends, metas = None, None, [], [], []
        for pid, path in files:
            nsp, names, rows = read_cells(path)
            if names0 is None:
                nsp0, names0 = nsp, names
            elif names != names0:
                sys.exit(f"{path}: the species differ from rank {files[0][0]}")
            meta = np.fromfile(path[:-len("cells")] + "meta", "<f8").reshape(-1, 4)
            if len(meta) != len(rows):
                sys.exit(f"{path}: {len(rows)} cells, {len(meta)} rows of meta")
            cells.append(rows)
            metas.append(meta)
            end_path = path[:-len("cells")] + "end"
            if os.path.exists(end_path):
                _, end_names, end_rows = read_cells(end_path)
                if end_names != names or end_rows.shape != rows.shape:
                    sys.exit(f"{end_path}: not the same cells as {path}")
                if not np.array_equal(end_rows[:, :3], rows[:, :3]):
                    sys.exit(f"{end_path}: x, y or w differ from {path}")
                ends.append(end_rows)
        if ends and len(ends) != len(files):
            sys.exit(f"iter {it} half {half}: {len(ends)} of {len(files)} ranks have a .end file")
        cells = np.concatenate(cells)
        metas = np.concatenate(metas)

        base = os.path.join(out, f"chemdump-{it}-{half}")
        with open(base + ".cells", "wb") as fp:
            fp.write(struct.pack("<qi", len(cells), nsp0))
            fp.write(names0)
            cells.astype("<f8").tofile(fp)
        metas.astype("<f8").tofile(base + ".meta")
        if ends:
            ends = np.concatenate(ends)
            with open(base + ".end", "wb") as fp:
                fp.write(struct.pack("<qi", len(ends), nsp0))
                fp.write(names0)
                ends.astype("<f8").tofile(fp)

        frozen = metas[:, 1] > 0.5
        cpu, nrhs = metas[~frozen, 2], metas[~frozen, 3]
        print(f"iter {it} half {half}: {len(cells)} cells from {len(files)} ranks, "
              f"{frozen.sum()} skipped by the gate")
        if len(cpu):
            q = np.percentile(nrhs, [50, 90, 99, 100])
            print(f"  solve CPU time: total {cpu.sum():.2f} s, median {1e3*np.median(cpu):.2f} ms, "
                  f"max {1e3*cpu.max():.2f} ms")
            print(f"  RHS evaluations: median {q[0]:.0f}, p90 {q[1]:.0f}, p99 {q[2]:.0f}, max {q[3]:.0f}")
        if len(ends):
            bad = (~np.isfinite(ends[:, 3:])).any(axis=1).sum()
            print(f"  end states: {bad} cells with a non-finite value (solve diverged)")
        print(f"  -> {base}.cells, {base}.meta" + (f", {base}.end" if len(ends) else ""))


if __name__ == "__main__":
    main()
