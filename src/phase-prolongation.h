/**
# The restriction and prolongation of phase values

The phase-change modules divide their VOF tracers (the gas species, `TS`,
`TG`) by the phase fraction before the solves, so that the fields hold the
value of one phase. A field of that kind needs its own restriction and
prolongation:

* The restriction of the VOF tracers, `restriction_volume_average()`,
  averages over all the children. Next to the interface some children hold
  no phase, and their value is 0. The parent then is not a value of the
  phase: for the species, $\sum Y < 1$ there.
* The bilinear prolongation reads the parent and its neighbours. A
  neighbour without the phase, or a parent of the kind above, gives a ghost
  cell with $\sum Y \ne 1$. In a fatehi restart at level 11, 1 to 3 faces
  at some steps saw $\sum Y$ off by up to 0.125, and the coupled species
  solve then kept $|\sum Y - 1|$ up to 2.9e-2 in a gas cell.

`restriction_phase()` averages only the children that hold the phase, with
the weight of the phase volume `cm*phase`. `refine_phase_bilinear()` keeps
the weights of `bilinear()` but drops the coarse cells without the phase
and divides by the sum of the weights that remain. The weights depend only
on the phase fraction, so they are the same for all the species: if the
coarse values have $\sum Y = 1$, the children have $\sum Y = 1$. Where all
the cells hold the phase, the two functions are the plain average and the
bilinear prolongation.

The phase fraction of a field `s` is `s.c` (or `1 - s.c` if `s.inverse`),
the attributes of the VOF tracers. Include this file after `vof.h`. A cell
holds the phase when its phase fraction is above `F_ERR`. */

#ifndef F_ERR
# define F_ERR 1e-10
#endif

static inline double phase_fraction (double c, bool inverse)
{
  return inverse ? 1. - c : c;
}

static void restriction_phase (Point point, scalar s)
{
  scalar c = s.c;
  double sw = 0., sv = 0.;
  foreach_child() {
    double ph = phase_fraction (c[], s.inverse);
    if (ph > F_ERR) {
      double w = cm[]*ph;
      sw += w, sv += w*s[];
    }
  }
  s[] = sw > 0. ? sv/sw : 0.;
}

static void refine_phase_bilinear (Point point, scalar s)
{
  scalar c = s.c;
  foreach_child() {
#if dimension == 1
    double w[2] = {3., 1.};
    double v[2] = {coarse(s), coarse(s,child.x)};
    double p[2] = {coarse(c), coarse(c,child.x)};
    int n = 2;
#elif dimension == 2
    double w[4] = {9., 3., 3., 1.};
    double v[4] = {coarse(s), coarse(s,child.x), coarse(s,0,child.y),
                   coarse(s,child.x,child.y)};
    double p[4] = {coarse(c), coarse(c,child.x), coarse(c,0,child.y),
                   coarse(c,child.x,child.y)};
    int n = 4;
#else // dimension == 3
    double w[8] = {27., 9., 9., 9., 3., 3., 3., 1.};
    double v[8] = {coarse(s), coarse(s,child.x), coarse(s,0,child.y),
                   coarse(s,0,0,child.z), coarse(s,child.x,child.y),
                   coarse(s,child.x,0,child.z), coarse(s,0,child.y,child.z),
                   coarse(s,child.x,child.y,child.z)};
    double p[8] = {coarse(c), coarse(c,child.x), coarse(c,0,child.y),
                   coarse(c,0,0,child.z), coarse(c,child.x,child.y),
                   coarse(c,child.x,0,child.z), coarse(c,0,child.y,child.z),
                   coarse(c,child.x,child.y,child.z)};
    int n = 8;
#endif
    double sw = 0., sv = 0.;
    for (int k = 0; k < n; k++)
      if (phase_fraction (p[k], s.inverse) > F_ERR)
        sw += w[k], sv += w[k]*v[k];
    s[] = sw > 0. ? sv/sw : 0.;
  }
}

/**
`phase_attributes (list, keepr, keepp, true)` keeps the restriction and the
prolongation of each field of `list` in `keepr` and `keepp`, and sets the
two functions above. `phase_attributes (list, keepr, keepp, false)` puts
the kept ones back. The arrays must hold `list_len (list)` entries.

Caution: every field of `list` must be a VOF tracer, with its attribute `c`
set. The default of `c` is the field of index 0, not "no field", so the
function cannot test it.

The coarse levels of the phase fraction must be up to date, so the
function restricts each phase fraction once when it sets the phase
functions. */

void phase_attributes (scalar * list,
                       void (** keepr)(Point, scalar),
                       void (** keepp)(Point, scalar),
                       bool phase)
{
  int k = 0, last = -1;
  for (scalar s in list) {
    if (phase) {
      keepr[k] = s.restriction, keepp[k] = s.prolongation;
      scalar c = s.c;
      if (c.i != last) {
        restriction ({c});
        last = c.i;
      }
      set_restriction (s, restriction_phase);
      set_prolongation (s, refine_phase_bilinear);
    }
    else {
      set_restriction (s, keepr[k]);
      set_prolongation (s, keepp[k]);
    }
    k++;
  }
}
