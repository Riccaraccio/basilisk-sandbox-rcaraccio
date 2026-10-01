/**
# Binning/Cell-Agglomeration Technique

We wish to accelerate the direct integration of the reactive step by grouping
together cells with similar thermodynamic states in *bins*.  Finally, the bins
are integrated as if they were a single computational cell, and the results are
then re-mapped into the cell values in a conservative manner.  We follow the
algorithm proposed by [Goldin et al. 2009](#goldin2009cell) (a copy is in
`papers/Goldin2009.md`).

This module is generic: it knows fields, targets and weights, not species
or enthalpy. `chemistry.h` uses it for the gas-phase reactions, and it does
the reaction step and the remap of section 2 of the paper itself.

## Usage

~~~literatec
BinTable * table = binning (fields, targets, tol, w, mask);

foreach_bin (table) {
  // integrate the state bin->phi[0..nfields-1]
}
binning_remap (table, fields);
binning_cleanup (table);
~~~

## The algorithm (section 3 of the paper)

1. Each target `k` is normalized over the cells of the mask, with the min
   and max of all MPI ranks (Eq. 7):
   `phi_hat = (phi - phi_min)/(phi_max - phi_min)`.
2. Each normalized target is cut into `L = 1/tol` equal intervals. The index
   of a cell in target `l` is `i_l = floor (phi_hat*L)`, in [0, L], as in
   Eq. 8. The cells exactly at the maximum (`phi_hat = 1`) thus get the
   index L and a bin of their own. Keep this: in an ignition the hottest
   cell ignites first. In `test/opensmoke-binning.c` the peak error of the
   temperature is 7e-4 K with the index L, and 757 K with a clamp to L-1
   (the hottest cell then goes into the bin below). The tolerance is
   relative to the current range of each target, and it is the same for
   all targets.
3. The cells with the same tuple `(i_0, ..., i_{M-1})` make one bin. Eq. 8
   of the paper writes the tuple as one integer `sum_l i_l L^l`. That
   integer does not fit in 64 bits for M = 11 and L = 100, so a hash table
   (Eq. 9) stores the tuple itself. Two different tuples never make one
   bin.
4. The state of a bin is the mean of `fields` over its cells, with the weight
   `w*dV` (Eq. 3). If the weights of a bin sum to 0, the mean uses the
   volume.

The table keeps the bin number of each cell in the field `table->id` (-1 for
the cells out of the mask). The mean and the remap are then one `foreach`
each.

Each MPI rank makes its own bins. Two cells on different ranks never make
one bin, also with the same tuple. The normalization uses the min and max of
all ranks, so a tuple means the same state on every rank.
`binning_mpi_report()` measures what one global bin map would save.

Caution: build the table, use it and clean it in the same step. The field
`table->id` does not follow a change of the grid (`adapt`). */

/**
`BINNING` makes `chemistry.h` use the binning sweep for the gas-phase
reactions. A case must include this file before `chemistry.h`, thus before
`multicomponent-varprop.h`. If the order is wrong, `chemistry.h` builds the
per-cell sweep and the case gets no binning, so stop the build. */

#ifdef CHEMISTRY_H
# error "Include binning.h before chemistry.h (thus before multicomponent-varprop.h)."
#endif

#define BINNING

#include <stdint.h>

/**
## Bin and bin table

`phi` and `phi0` are the mean of each field in the bin, after and before the
integration. `key` is the tuple of the bin (one index for each target). */

typedef struct {
  int nfields;
  long ncells;
  double mass, volume;
  double * phi, * phi0;
  int64_t * key;
} Bin;

typedef struct {
  int ntargets, nfields;
  int64_t L;
  int nbins, capacity;
  long nmask; // the cells of the mask on this rank
  Bin * bins;
  scalar id;
} BinTable;

/**
## The hash table of the tuples

Open addressing with linear probing. The size is a power of 2 with at least
two slots for each key, so a probe sequence is short.

The hash is 64-bit FNV-1a over the indices (offset basis
14695981039346656037, prime 1099511628211), then the final mix of
MurmurHash3, so that the low bits that select the slot depend on all bits.
The keys and the hash use `int64_t` and `uint64_t`, not `long`: `long` has
32 bits on 64-bit Windows, and `L` can be 1e12. */

static uint64_t binning_key_hash (const int64_t * key, int n)
{
  uint64_t h = UINT64_C(14695981039346656037);
  for (int l = 0; l < n; l++) {
    h ^= (uint64_t) key[l];
    h *= UINT64_C(1099511628211);
  }
  h ^= h >> 33;
  h *= UINT64_C(0xff51afd7ed558ccd);
  h ^= h >> 33;
  return h;
}

static int binning_slots_size (long nkeys)
{
  int size = 16;
  while (size < 2*nkeys + 2)
    size *= 2;
  return size;
}

/**
Return the slot of `key` in `slots` (the slot holds an index into `keys`,
or -1 if it is empty). */

static int binning_slot (const int * slots, int size, const int64_t * keys,
                         int n, const int64_t * key)
{
  int s = binning_key_hash (key, n) & (size - 1);
  while (slots[s] >= 0) {
    const int64_t * k = keys + (int64_t) slots[s]*n;
    bool same = true;
    for (int l = 0; l < n && same; l++)
      if (k[l] != key[l])
        same = false;
    if (same)
      return s;
    s = (s + 1) & (size - 1);
  }
  return s;
}

static Bin * binning_add_bin (BinTable * table, const int64_t * key)
{
  if (table->nbins == table->capacity) {
    table->capacity = table->capacity ? 2*table->capacity : 64;
    table->bins = realloc (table->bins, table->capacity*sizeof (Bin));
  }
  Bin * bin = table->bins + table->nbins++;
  bin->nfields = table->nfields;
  bin->ncells = 0;
  bin->mass = bin->volume = 0.;
  bin->phi = calloc (2*table->nfields, sizeof (double));
  bin->phi0 = bin->phi + table->nfields;
  bin->key = malloc (table->ntargets*sizeof (int64_t));
  for (int l = 0; l < table->ntargets; l++)
    bin->key[l] = key[l];
  return bin;
}

/**
## Iterators

`foreach_bin` loops over the bins of the table, with the index `i` and the
pointer `bin`. `foreach_bin_field` loops over the fields of a bin, with the
index `j` and the values `phi` and `phi0`. */

macro foreach_bin (BinTable * table) {
  for (int i = 0; i < table->nbins; i++) {
    Bin * bin = table->bins + i; NOT_UNUSED (bin);
    {...}
  }
}

macro foreach_bin_field (Bin * bin) {
  for (int j = 0; j < bin->nfields; j++) {
    double phi = bin->phi[j], phi0 = bin->phi0[j];
    {...}
    bin->phi[j] = phi, bin->phi0[j] = phi0;
  }
}

/**
## Program interface

`binning()` builds the table: the bin of each cell of `mask`, and the mean
of `fields` in each bin. `tol` is the tolerance of the paper (`L = 1/tol`
intervals for each target, with `tol` at least 1e-12). */

trace
BinTable * binning (scalar * fields, scalar * targets, double tol,
                    (const) scalar w = unity, (const) scalar mask = unity)
{
  int M = list_len (targets);
  int nf = list_len (fields);

  BinTable * table = calloc (1, sizeof (BinTable));
  table->ntargets = M;
  table->nfields = nf;
  table->L = (int64_t) (1./max (tol, 1e-12) + 0.5);
  if (table->L < 1)
    table->L = 1;
  int64_t L = table->L;

  scalar id = new scalar;
  table->id = id;

  /**
  The min and max of each target over the mask, on all ranks (Eq. 7). */

  double tmin[M + 1], tmax[M + 1];
  int l = 0;
  for (scalar t in targets) {
    double mn = HUGE, mx = - HUGE;
    foreach (reduction(min:mn) reduction(max:mx))
      if (mask[]) {
        if (t[] < mn) mn = t[];
        if (t[] > mx) mx = t[];
      }
    tmin[l] = mn, tmax[l] = mx;
    l++;
  }

  /**
  The bin of each cell (Eq. 8 and 9). */

  long nmask = 0;
  foreach (serial)
    if (mask[])
      nmask++;

  int size = binning_slots_size (nmask);
  int * slots = malloc (size*sizeof (int));
  for (int s = 0; s < size; s++)
    slots[s] = -1;
  int64_t * keys = malloc ((nmask > 0 ? nmask : 1)*(M > 0 ? M : 1)*sizeof (int64_t));

  foreach (serial) {
    id[] = -1;
    if (mask[]) {
      int64_t key[M + 1];
      int m = 0;
      for (scalar t in targets) {
        double range = tmax[m] - tmin[m];
        int64_t il = 0;
        if (range > 0. && isfinite (t[]))
          il = (int64_t) floor ((t[] - tmin[m])/range*L);
        key[m] = il < 0 ? 0 : il > L ? L : il;
        m++;
      }
      int s = binning_slot (slots, size, keys, M, key);
      if (slots[s] < 0) {
        slots[s] = table->nbins;
        for (int k = 0; k < M; k++)
          keys[(int64_t) table->nbins*M + k] = key[k];
        binning_add_bin (table, key);
      }
      id[] = slots[s];
      table->bins[slots[s]].ncells++;
    }
  }
  free (slots);
  free (keys);
  table->nmask = nmask;

  /**
  The mean of each field (Eq. 3). `phi0` holds the volume sum first, for a
  bin whose weights sum to 0. */

  foreach (serial) {
    int b = id[];
    if (b >= 0) {
      Bin * bin = table->bins + b;
      double v = dv(), m = w[]*v;
      bin->mass += m;
      bin->volume += v;
      int j = 0;
      for (scalar s in fields) {
        bin->phi[j] += s[]*m;
        bin->phi0[j] += s[]*v;
        j++;
      }
    }
  }
  foreach_bin (table)
    for (int j = 0; j < nf; j++) {
      bin->phi[j] = bin->mass > 0. ? bin->phi[j]/bin->mass :
        bin->volume > 0. ? bin->phi0[j]/bin->volume : 0.;
      bin->phi0[j] = bin->phi[j];
    }

  return table;
}

/**
The simple remap of the paper (Eq. 5): each cell gets the change of its bin.
`chemistry.h` does not use it: it scales the change to keep the cells
realizable, and it remaps the energy as enthalpy. */

void binning_remap (BinTable * table, scalar * fields)
{
  scalar id = table->id;
  foreach() {
    int b = id[];
    if (b >= 0) {
      Bin * bin = table->bins + b;
      int j = 0;
      for (scalar s in fields) {
        s[] += bin->phi[j] - bin->phi0[j];
        j++;
      }
    }
  }
}

void binning_cleanup (BinTable * table)
{
  for (int i = 0; i < table->nbins; i++) {
    free (table->bins[i].phi);
    free (table->bins[i].key);
  }
  free (table->bins);
  scalar id = table->id;
  delete ({id});
  free (table);
}

/**
### Post-Processing

The statistics of the table of this rank. */

typedef struct {
  // total number of (active/non-empty) bins
  size_t nactive;
  // average number of cells per bin
  size_t navg;
  // maximum number of cells in a bin
  size_t nmax;
  // minimum number of cells in a bin
  size_t nmin;
  // number of cells in the mask
  size_t ncells;
} bstats;

bstats binning_stats (const BinTable * table)
{
  bstats s = {0};
  s.nmin = table->nbins > 0 ? (size_t) 1e9 : 0;
  size_t nreal = 0;
  for (int i = 0; i < table->nbins; i++) {
    size_t n = table->bins[i].ncells;
    nreal += n;
    s.nmax = max (s.nmax, n);
    s.nmin = min (s.nmin, n);
  }
  s.nactive = table->nbins;
  s.ncells = nreal;
  s.navg = s.nactive > 0 ? nreal/s.nactive : 0;
  return s;
}

void binning_ids (const BinTable * table, scalar ids)
{
  scalar id = table->id;
  foreach()
    ids[] = id[] >= 0 ? id[] : nodata;
}

/**
## Local or global bins: the MPI diagnostic

Each rank makes its own bins. One global bin map would integrate each tuple
one time, and it could give each rank an equal share of the bins. This
function measures what that would save, before anyone writes the global map.

`cost[i]` is the time of the integration of the bin `i` of this rank. The
function gathers the tuples and the costs on rank 0, and writes to `fp` (on
rank 0 only):

    t half nbins_sum nkeys Tlocal_max Tlocal_mean Tglobal_ideal

* `nbins_sum`: the bins of all ranks (the solves with local bins);
* `nkeys`: the different tuples of all ranks (the solves with a global map);
* `Tlocal_max`, `Tlocal_mean`: the max and the mean over the ranks of the
  sum of `cost`. The step waits for the max;
* `Tglobal_ideal`: the sum over the different tuples of the mean cost of the
  tuple, divided by the number of ranks. This is the time of a global map
  with a perfect share and no communication.

The gather costs communication, so call it only for a measurement. */

@if _MPI
void binning_mpi_report (const BinTable * table, const double * cost,
                         double t, int half, FILE * fp)
{
  int np = npe(), M = table->ntargets, nb = table->nbins;
  double tlocal = 0.;
  for (int i = 0; i < nb; i++)
    tlocal += cost[i];

  int * counts = NULL;
  if (pid() == 0)
    counts = malloc (np*sizeof (int));
  MPI_Gather (&nb, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int64_t * keys = malloc ((nb > 0 ? nb : 1)*(M > 0 ? M : 1)*sizeof (int64_t));
  for (int i = 0; i < nb; i++)
    for (int l = 0; l < M; l++)
      keys[(int64_t) i*M + l] = table->bins[i].key[l];

  int * kcounts = NULL, * kdispl = NULL, * ccounts = NULL, * cdispl = NULL;
  int64_t * allkeys = NULL;
  double * allcost = NULL;
  long ntot = 0;
  if (pid() == 0) {
    kcounts = malloc (np*sizeof (int)), kdispl = malloc (np*sizeof (int));
    ccounts = malloc (np*sizeof (int)), cdispl = malloc (np*sizeof (int));
    for (int r = 0; r < np; r++) {
      ccounts[r] = counts[r], cdispl[r] = ntot;
      kcounts[r] = counts[r]*M, kdispl[r] = ntot*M;
      ntot += counts[r];
    }
    allkeys = malloc ((ntot > 0 ? ntot : 1)*(M > 0 ? M : 1)*sizeof (int64_t));
    allcost = malloc ((ntot > 0 ? ntot : 1)*sizeof (double));
  }
  MPI_Gatherv (keys, nb*M, MPI_INT64_T, allkeys, kcounts, kdispl, MPI_INT64_T,
               0, MPI_COMM_WORLD);
  MPI_Gatherv ((void *) cost, nb, MPI_DOUBLE, allcost, ccounts, cdispl,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double tmax = tlocal, tsum = tlocal;
  MPI_Reduce (&tlocal, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (&tlocal, &tsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (pid() == 0) {
    int size = binning_slots_size (ntot);
    int * slots = malloc (size*sizeof (int));
    for (int s = 0; s < size; s++)
      slots[s] = -1;
    int64_t * ukeys = malloc ((ntot > 0 ? ntot : 1)*(M > 0 ? M : 1)*sizeof (int64_t));
    double * ucost = calloc (ntot > 0 ? ntot : 1, sizeof (double));
    int * ucount = calloc (ntot > 0 ? ntot : 1, sizeof (int));
    int nkeys = 0;
    for (long i = 0; i < ntot; i++) {
      const int64_t * key = allkeys + i*M;
      int s = binning_slot (slots, size, ukeys, M, key);
      if (slots[s] < 0) {
        slots[s] = nkeys;
        for (int l = 0; l < M; l++)
          ukeys[(int64_t) nkeys*M + l] = key[l];
        nkeys++;
      }
      ucost[slots[s]] += allcost[i];
      ucount[slots[s]]++;
    }
    double tglobal = 0.;
    for (int k = 0; k < nkeys; k++)
      tglobal += ucost[k]/ucount[k];
    tglobal /= np;
    fprintf (fp, "%g %d %ld %d %g %g %g\n", t, half, ntot, nkeys,
             tmax, tsum/np, tglobal);
    fflush (fp);
    free (slots), free (ukeys), free (ucost), free (ucount);
    free (counts), free (kcounts), free (kdispl), free (ccounts), free (cdispl);
    free (allkeys), free (allcost);
  }
  free (keys);
}
@endif

/**
## References

~~~bib
@article{goldin2009cell,
  title={A cell agglomeration algorithm for accelerating detailed chemistry in CFD},
  author={Goldin, Graham M and Ren, Zhuyin and Zahirovic, Selma},
  journal={Combustion Theory and Modelling},
  volume={13},
  number={4},
  pages={721--739},
  year={2009},
  publisher={Taylor \& Francis}
}
~~~
*/
