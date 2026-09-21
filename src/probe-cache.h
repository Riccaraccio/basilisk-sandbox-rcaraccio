/**
# Cache of the quantities that two probes share

Some probe events of a case ask for the same global quantity at the same
instant. `run/test.c` asks two times for `statsf(T)` and two times for the
integral of the discrete divergence of `uf`. Each request costs one sweep of
the grid.

The two functions below hold the value and the step index of the last
evaluation. A second request in the same step returns the value that the
first request computed, and it costs nothing. The value is therefore the
same in every file that prints it, to the last bit.

The step counter of these functions is the global `iter`. The name `i` is
visible inside an event body only.

## How to use it

Call these functions from the probe events only. Every probe event of
`run/test.c` runs before the `stability` event, so all of them read the same
state. A call from a later event of the same step returns the value of the
start of the step, which is not the value of that later state.

`probe_div_uf_set()` gives the cache a value that the caller computed in a
loop of its own. Use it when the loop needs the per-cell divergence for
another column as well.

Caution: the integral of `probe_div_uf()` uses the same weight as
`gas_source`. The discrete divergence carries `cm[]`, so the integral uses
`sq(Delta)` and gives a value per radian in an axisymmetric case. */

static int probe_statsT_i = -1;
static stats probe_statsT_val;

static stats probe_stats_T (void)
{
  if (probe_statsT_i != iter) {
    probe_statsT_val = statsf (T);
    probe_statsT_i = iter;
  }
  return probe_statsT_val;
}

static int probe_div_i = -1;
static double probe_div_val = 0.;

static void probe_div_uf_set (double value)
{
  probe_div_val = value;
  probe_div_i = iter;
}

static double probe_div_uf (void)
{
  if (probe_div_i != iter) {
    double q = 0.;
    foreach (reduction(+:q)) {
      double d = 0.;
      foreach_dimension()
        d += uf.x[1] - uf.x[];
      d /= Delta;
      q += d*sq(Delta);
    }
    probe_div_uf_set (q);
  }
  return probe_div_val;
}
