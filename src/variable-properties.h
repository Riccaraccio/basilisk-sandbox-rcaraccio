/**
# Variable properties
This file defines the variable properties for the gas and solid phases, as well as the functions to compute them. 
*/


#define VARPROP

#ifndef F_ERR
# define F_ERR 1e-10
#endif

/**
Variable properties fields. 
Given that we consider posous media, we define separate properties for the gas phase inside the pores and for the solid matrix.
The suffixe "_S" refers to internal gas  and "_G" to the surrounding gas.
 */
scalar rhoGv_G[], rhoGv_S[], rhoSv[];
scalar muGv_G[], muGv_S[];
scalar lambdaGv_G[], lambdaGv_S[], lambdaSv[];
scalar cpGv_G[], cpGv_S[], cpSv[];

/**
## Phase aware prolongation and restriction

`update_properties()` resets every field above to 0 and fills each one only
where its gate passes. A 0 therefore does not mean "a small density". It
means "this phase is not present in this cell".

The default operators of Basilisk do not know that. `refine_bilinear` mixes
the 0 of a cell of the other phase into a new cell of this phase, and
`restriction_average` divides by the full count of the children even when
some of them hold no value.

`PHASE_AWARE_PROPERTIES` at 1 gives the two operators below to those
fields. They use the weights of the Basilisk operators, but they take only
the donors that hold a value. A new cell then takes the value of the same
phase from its neighbours. If no donor holds a value, the result stays 0.

Caution: the default is 0. These operators do not repair the stop of
2026-09-03. That stop was measured, and the parent of the failing cell held
`rhoGv_G` of 0.5476 while three of its four donors were positive. The 0 of
that cell did not come from the prolongation. See the section below. Turn
these operators on only after you measure what they change: they alter the
grid transfer of every property field, and the ladder campaign measures
differences of 0.5 to 1.5 %.

These operators suit a quantity that is positive where its phase is present
and 0 where it is not. Do not give them to a field that can be 0 or negative
for a physical reason. */

#ifndef PHASE_AWARE_PROPERTIES
# define PHASE_AWARE_PROPERTIES 0
#endif

#if TREE

static inline double bilinear_phase (Point point, scalar s)
{
#if dimension == 1
  double w[2] = {3., 1.};
  double v[2] = {coarse(s), coarse(s,child.x)};
#elif dimension == 2
  double w[4] = {9., 3., 3., 1.};
  double v[4] = {coarse(s), coarse(s,child.x), coarse(s,0,child.y),
                 coarse(s,child.x,child.y)};
#else
  double w[8] = {27., 9., 9., 9., 3., 3., 3., 1.};
  double v[8] = {coarse(s), coarse(s,child.x), coarse(s,0,child.y),
                 coarse(s,0,0,child.z), coarse(s,child.x,child.y),
                 coarse(s,child.x,0,child.z), coarse(s,0,child.y,child.z),
                 coarse(s,child.x,child.y,child.z)};
#endif
  int n = sizeof(v)/sizeof(v[0]);
  double sum = 0., wsum = 0.;
  for (int k = 0; k < n; k++)
    if (v[k] > 0.) {
      sum += w[k]*v[k];
      wsum += w[k];
    }
  return wsum > 0. ? sum/wsum : 0.;
}

static inline void refine_phase (Point point, scalar s)
{
  foreach_child()
    foreach_blockf (s)
      s[] = bilinear_phase (point, s);
}

static inline void restriction_phase (Point point, scalar s)
{
  foreach_blockf (s) {
    double sum = 0.;
    int n = 0;
    foreach_child()
      if (s[] > 0.) {
        sum += s[];
        n++;
      }
    s[] = n ? sum/n : 0.;
  }
}

/**
Give the operators to every field that `update_properties()` fills per
phase. `set_restriction()` keeps the coarse levels consistent with the
prolongation: the coarse donors that `bilinear_phase()` reads come from
`restriction_phase()`. */

event defaults (i = 0) {
#if PHASE_AWARE_PROPERTIES
  for (scalar s in {rhoGv_G, rhoGv_S, rhoSv,
                    muGv_G, muGv_S,
                    lambdaGv_G, lambdaGv_S, lambdaSv,
                    cpGv_G, cpGv_S, cpSv}) {
    s.refine = s.prolongation = refine_phase;
    s.restriction = restriction_phase;
  }
#endif
#ifdef PROPERTIES_VERBOSE
  fprintf (stderr, "variable-properties: rhoGv_G.refine=%s prolongation=%s"
                   " restriction=%s\n",
           rhoGv_G.refine == refine_phase ? "phase" : "OTHER",
           rhoGv_G.prolongation == refine_phase ? "phase" : "OTHER",
           rhoGv_G.restriction == restriction_phase ? "phase" : "OTHER");
#endif
}

#endif // TREE

typedef struct {
  double T, P;
  double * x;
} ThermoState;

typedef struct {
  // Mixture properties
  double (* rhov)     (void *);
  double (* muv)      (void *);
  double (* lambdav)  (void *);
  double (* cpv)      (void *);
  // Species properties
  void   (* diff)     (void *, double *);
  double (* cps)      (void *, int);
  void   (* cpvs)     (void *, double *);
} ThermoProps;

#define aavg(f,v1,v2) (clamp(f,0.,1.)*(v1 - v2) + v2)
#define havg(f,v1,v2) (1./(clamp(f,0,1)*(1./(v1) - 1./(v2)) + 1./(v2)))

extern scalar f;
extern face vector alphav;
extern scalar rhov;
#ifdef FILTERED
extern scalar sf;
#else
# define sf f
#endif

/**
We overwrite the properties for the Navier-Stokes solver with the variable properties.
*/

/**
Caution: `rhov[]` holds `cm[]*rhomix`, so it already carries the metric. Do
not give it metric-aware tree operators. Commit `f44028d` added

    rhov.refine = refine_linear;
    set_restriction (rhov, restriction_volume_average);

in a `defaults` event, and from that commit no case of this sandbox started
from scratch: every one stopped at step 1 with SIGFPE at `viscosity.h:167`,
which computes `dt/rho[]`. Both operators apply the metric a second time to a
field that already has it. `restriction_volume_average` sums `cm[]*s[]` over
the children and divides by the `cm[]` of the parent. Near the axis `cm` is
`y`, and `y` changes by a factor of 3 between sibling cells, so the coarse
value becomes non-positive. The block below is commented out for this reason.
Do not restore it. `two-phase-generic.h` gives the same field the default
operators, in the same axisymmetric cases. Keep that choice.

The symptom hides on the leaves. This event rewrites every leaf from
`rhoGv_G`, `rhoGv_S` and `f` on each timestep, so a `foreach()` loop finds
`rhov` positive everywhere. The measurement on the failing run gave a minimum
of 9.8e-5 and no cell at or below zero, out of 65536. But `viscosity.h`
solves on a multigrid with `foreach_level_or_leaf`, so it reads the coarse
levels too, and a leaf loop never visits those. Examine the coarse levels
before you conclude that the density is good.

This event carried three guards until 2026-09-03: a fallback that read the
density back from `rhov[]/cm[]`, a test before `1./rhomix`, and a test before
the write of `rhov[]`. They are gone, and they are not needed. `rhomix` stays
positive in every configuration of this sandbox. A guard on `rhov[]` also
does not work: it stops a good value from becoming 0, and it cannot repair a
cell whose stored value is already bad. That was measured on the from-scratch
start, and the crash did not change.

If a later change does make `rhomix` reach 0, the run stops at `1./rhomix`
in the loop below, because Basilisk arms the floating point traps. Fix the
cause in `update_properties()`, which sets both densities to 0 and refills
them only where its gates pass. Do not add a guard here that hides it.

Caution: do not identify the operation from the signal code alone.
`fsolve-gsl.h` masks the traps around the GSL solve and the masked operations
set sticky status flags. On x86 the SSE unit tests each instruction on its
own, thus a stale flag raises no signal there, and only the x87 unit defers a
trap in this way. That header now calls `feclearexcept (FE_ALL_EXCEPT)`
before it arms the traps again, thus the doubt is removed. Read the signal
code together with the line that `llvm-addr2line` gives, and remember that
`-O2` moves the line of an instruction.

On 2026-09-03 four runs of the ladder case `run/test.c` stopped at the line
`1./rhomix` below, with the signal code 7 that the caution above describes.
`llvm-addr2line` gave the same three frames in all four:

~~~literatec
properties_0   src/variable-properties.h:105
event          basilisk/src/grid/events.h:266
adapt          basilisk/src/navier-stokes/centered.h:461
~~~

The third frame is the important one. `centered.h` has

~~~literatec
event adapt (i++,last) {
  event ("properties");
}
~~~

so this event runs a second time in each step, after `adapt_wavelet` changes
the grid. `update_properties()` does not run again there. It runs in the
`tracer_diffusion` slot, which comes earlier in the step. The cells that
adapt creates therefore hold prolongated values of `rhoGv_G` and `rhoGv_S`,
not computed ones.

Those two fields are 0 outside their own phase, because `update_properties()`
resets them and refills each one only where its gate passes. The default
prolongation of a scalar is `refine_bilinear`, which interpolates across that
step to 0 at every interface cell. `rhomix` is a sum of the two, so the two
terms can cancel and give 0.

The block below reports the cell before the trap. It does not repair it.
Read the report, then fix the cause. Do not add a guard here.

Caution: `i` and `t` in a report of this event can read 0 and 0.
`centered.h` calls the event with `event ("properties")`, and `events.h`
gives the action dummy values for the step and the time. Take the time of
the stop from the last line of the output files, not from the report.

## The cause of the stop of 2026-09-03

The run of `run/test.c` with `MOISTURE=1` stops at t = 5.942, at the line
`1./rhomix` below. It stops in serial and with MPI, with the interior pin
and without it, from a restart and from the start. `x/i $pc` gives
`divsd %xmm0,%xmm2` and MXCSR holds `[ DE ZE PE ... ]`, thus a division by
zero. The local variable `rhomix` reads 0.

This event is not the cause. It is the place where the run stops.

`GASGATE_DEBUG` in `multicomponent-properties.h` reports every cell that
holds gas and does not pass the test `TG[] > 0.` of `update_properties()`.
The report of the failing run:

~~~literatec
t         x           y          f          TG
5.94      -0.00273437 0.00273438 0.99984    -0.0145
5.9405    -0.00273437 0.00273438 0.99984    -0.333
5.941     -0.00273437 0.00273438 0.99984    -1.099
5.9415    -0.00273437 0.00273438 0.99984    -2.793
5.942     -0.00273437 0.00273438 0.99984    -6.321
5.94217   -0.00273437 0.00273438 0.99984    -17.71    (spreads to 3 cells)
5.94228   -0.00289062 0.00289062 0.        -334.13    (the cell that stops it)
~~~

One cell holds the gas temperature. It has `f` of 0.99984, thus a gas
fraction of 1.6e-4. Its `TG` becomes negative and grows by a factor of about
3 at each step. After 7 steps it reaches the neighbours, and then a cell
with `f` of 0.

`update_properties()` is correct here. It refuses to give a density to a gas
at a negative temperature, thus it leaves `rhoGv_G` at the reset value of 0.
A cell with `f` of 0 and `rhoGv_G` of 0 gives `rhomix` of 0, and this event
divides by it.

The cause is the energy equation of the gas in a cell that holds almost no
gas. `TG` is a tracer of the form `TG*(1-f)`, and `theta2` carries the same
factor. At a gas fraction of 1.6e-4 the equation is close to singular and
the solution runs away. The threshold that admits such a cell to the solve
is `F_ERR`, which is 1e-10. That is a threshold for a volume fraction, not
for the conditioning of an energy equation.

Do not repair this here. Repair it in the energy equation: give the gas
temperature its own threshold, well above `F_ERR`, and give a cell below
that threshold the temperature of the solid.

The block below reports the cell before the trap. It does not repair it.
Read the report, then fix the cause. Do not add a guard here.

Caution: `i` and `t` in a report of this event can read 0 and 0.
`centered.h` calls the event with `event ("properties")`, and `events.h`
gives the action dummy values for the step and the time. Take the time of
the stop from the last line of the output files, not from the report.

## The cause of the stop of 2026-09-03

The run of `run/test.c` with `MOISTURE=1` stops at t = 5.94, at the line
`1./rhomix` below. It stops in serial and with MPI, with the interior pin
and without it, from a restart and from the start. `x/i $pc` gives
`divsd %xmm0,%xmm2` and MXCSR holds `[ DE ZE PE ... ]`, thus a division by
zero. The local variable `rhomix` reads 0.

`PROPERTIES_TRAP` at 0 gives the cell and its 3 by 3 neighbourhood:

~~~literatec
x -0.00289062  y 0.00289062  Delta 0.00015625  level 10
f 0  porosity 0  rhoGv_G 0  rhoGv_S 0

di dj      f     porosity   rhoGv_G    rhoGv_S
-1 -1      0        0        0.3946       0
 0 -1   0.4753    0.0950       0        0.5367
 0  0      0        0           0          0     <- the cell
 1 -1   0.9998    0.2045       0        0.5605
 1  0   0.4757    0.0951       0        0.5430
 1  1      0        0        1.2735       0
~~~

The cell sits at r = 4.088e-3, one half of a cell outside a particle of
radius 4e-3. Every neighbour with `f` of 0 holds a good `rhoGv_G`. Every
neighbour that holds the interface has `rhoGv_G` of 0 and a good `rhoGv_S`,
because `update_properties()` resets both fields and fills each one only
where its gate passes. The cell has `f` of 0, thus it is gas, but it holds
the 0 of the solid side.

The stop is in the second call of this event, the one that `centered.h`
makes after `adapt_wavelet`. `update_properties()` runs in the
`tracer_diffusion` slot, thus it does not run again there. The cells that
the adaptation makes hold prolongated values. `f` uses `fraction_refine`,
which is geometric, and `rhoGv_G` uses `refine_bilinear`. The two rules do
not agree on the side of the interface that a new cell is on: `f` gives 0,
which means gas, and `rhoGv_G` gives 0, which means solid. `rhomix` is then
0 and the division stops the run.

Note: `refine_bilinear` is a convex combination, thus it cannot make a
negative value. It can carry a 0 without any difficulty. An earlier note in
this file rejected the whole prolongation path for the first reason. That
was too broad.

The repair belongs in the adaptation, not here. Fill the property fields of
the new cells before this event reads them, or give those fields a
prolongation that takes only the donors of the same phase. A guard here
cannot know the phase of the cell, because both densities are 0. */

#ifndef PROPERTIES_REPORT_MAX
# define PROPERTIES_REPORT_MAX 50
#endif

/**
'RHOMIX_MIN' is the smallest density that this sandbox accepts. The gas at
1000 K and 1 atm gives about 0.35 kg/m3, and the solid gives about 1550. A
value below 1e-3 has no physical meaning.

The first form of the test below was 'rhomix > 0.'. It has two holes, and the
run of 2026-09-03 fell through them. Basilisk arms only 'FE_DIVBYZERO' and
'FE_INVALID' (see 'grid/config.h'). It does not arm 'FE_OVERFLOW'. Thus:

* A 'rhomix' of 1e-320 passes 'rhomix > 0.'. The division gives 'inf' and
  raises no signal. The cell then holds an enormous 'alphav', the velocity
  grows, the timestep collapses, and the run stops later, at a place that
  has no relation to the cause.
* A 'rhomix' of NaN makes the test itself raise 'FE_INVALID', because C maps
  '>' to a signalling comparison. The run stops on the line of the test,
  before it writes the report.

'rhomix_is_bad()' below uses no floating point operation on a value that can
be NaN. It reads the exponent field and the sign bit as an integer. */

#include <stdint.h>

#ifndef RHOMIX_MIN
# define RHOMIX_MIN 1.e-3
#endif

#ifndef PROPERTIES_SCAN
# define PROPERTIES_SCAN 1
#endif

/**
`PROPERTIES_TRAP` at 0 masks the traps for the length of this event and
reports every cell that gives a bad density, after the loops. Use it when
the run stops on `1./rhomix` and the scan above reports nothing. The report
holds the density that the scan read and the density that the division used.
If the two differ, a step between the two loops changed the fields.

The build of the campaign keeps the default of 1. */

#ifndef PROPERTIES_TRAP
# define PROPERTIES_TRAP 1
#endif

static inline bool not_finite (double x) {
  union { double d; uint64_t u; } v = {.d = x};
  return (v.u & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL;
}

static inline bool rhomix_is_bad (double x) {
  union { double d; uint64_t u; } v = {.d = x};
  if ((v.u & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL)
    return true;                       // inf or NaN
  if (v.u >> 63)
    return true;                       // negative, -0. included
  return v.d < RHOMIX_MIN;             // x is finite and positive here
}

event properties (i++) {

  /**
  ## The scan

  This pass holds no division and no product of a field with a field, thus
  the compiler has nothing to move into it. The first form of this report
  sat in the loop below, on the line above `1./rhomix`. It never wrote a
  file, and `llvm-addr2line` still gave the line of the division. At `-O2`
  GCC schedules the division across the branch of the report, thus the trap
  came before `fflush()` reached the disk. A separate pass removes that
  whole class of failure.

  Read the result this way:

  * The file exists: the density is the cause. The columns give the cell and
    both densities.
  * The run stops and the file does not exist: the density is good and the
    cause is elsewhere in this event. Look at the products in the two loops
    below, and remember that `fm.x[]` is 0 on the axis, thus `0.*inf` is an
    invalid operation.

  Set `PROPERTIES_SCAN` to 0 to remove the pass. Keep it at 1 while the
  cause of the stop of 2026-09-03 is open. */

#if PROPERTIES_SCAN
  foreach() {

    /**
    Test the inputs first. `rhoGv_S[]` of `inf` with `f[]` of 0 gives an
    invalid operation in the sum itself, before any test of the sum. */

    bool bad = (not_finite (rhoGv_G[]) || not_finite (rhoGv_S[]) ||
                not_finite (muGv_G[])  || not_finite (muGv_S[])  ||
                not_finite (f[]));
    double rhomix = bad ? 0. : rhoGv_G[]*(1.-f[]) + rhoGv_S[]*f[];

    if (bad || rhomix_is_bad (rhomix)) {
      static FILE * fp = NULL;
      static int nrep = 0;
      if (nrep < PROPERTIES_REPORT_MAX) {
        if (!fp) {
          char name[80];
          snprintf (name, sizeof(name), "rhomix-dbg-%d.dat", pid());
          fp = fopen (name, "w");
          fprintf (fp, "#i t dt x y Delta level f porosity rhomix"
                       " rhoGv_G rhoGv_S muGv_G muGv_S badinput\n");
        }
        fprintf (fp, "%d %g %g %g %g %g %d %.17g %.17g %.17g %.17g %.17g"
                     " %.17g %.17g %d\n",
                 i, t, dt, x, y, Delta, level, f[], porosity[], rhomix,
                 rhoGv_G[], rhoGv_S[], muGv_G[], muGv_S[], (int) bad);
        fflush (fp);
        nrep++;
      }
    }
  }
#endif

#if !PROPERTIES_TRAP
  scalar rhopre[], rhouse[];
  foreach()
    rhopre[] = rhoGv_G[]*(1.-f[]) + rhoGv_S[]*f[];
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
#endif

  scalar alphacenter[], mucenter[];
  foreach() {
    double rhomix = rhoGv_G[]*(1.-f[]) + rhoGv_S[]*f[];
#if !PROPERTIES_TRAP
    rhouse[] = rhomix;
#endif
    alphacenter[] = 1./rhomix;
    mucenter[] = (muGv_G[]*(1.-f[]) + muGv_S[]*f[]);
    rhov[] = cm[]*rhomix;
  }

  foreach_face() {
    alphav.x[] = fm.x[]*face_value(alphacenter, 0);
    {
      face vector muv = mu;
      muv.x[] = fm.x[]*face_value(mucenter, 0);
    }
  }

#if !PROPERTIES_TRAP
  feclearexcept (FE_ALL_EXCEPT);
  enable_fpe (FE_DIVBYZERO|FE_INVALID);

  int nbad = 0;
  foreach (reduction(+:nbad))
    if (rhomix_is_bad (rhouse[]))
      nbad++;

  if (nbad) {
    FILE * fp = fopen ("rhomix-post.dat", "w");
    fprintf (fp, "#i t dt x y Delta level f porosity rhopre rhouse"
                 " rhoGv_G rhoGv_S muGv_G muGv_S\n");
    foreach()
      if (rhomix_is_bad (rhouse[])) {
        fprintf (fp, "%d %g %g %g %g %g %d %.17g %.17g %.17g %.17g %.17g"
                     " %.17g %.17g %.17g\n",
                 i, t, dt, x, y, Delta, level, f[], porosity[],
                 rhopre[], rhouse[], rhoGv_G[], rhoGv_S[], muGv_G[], muGv_S[]);

        /**
        The 3 by 3 neighbourhood of the cell. It says whether the cell took
        the 0 from the solid side, where `update_properties()` leaves
        `rhoGv_G` at the reset value. */

        fprintf (fp, "# parent: level %d f %.17g porosity %.17g"
                     " rhoGv_G %.17g rhoGv_S %.17g\n",
                 level - 1, coarse(f), coarse(porosity),
                 coarse(rhoGv_G), coarse(rhoGv_S));
        fprintf (fp, "# parent neighbours: dk dl f rhoGv_G rhoGv_S\n");
        fprintf (fp, "#  %d 0 %.17g %.17g %.17g\n", child.x,
                 coarse(f,child.x), coarse(rhoGv_G,child.x),
                 coarse(rhoGv_S,child.x));
        fprintf (fp, "#  0 %d %.17g %.17g %.17g\n", child.y,
                 coarse(f,0,child.y), coarse(rhoGv_G,0,child.y),
                 coarse(rhoGv_S,0,child.y));
        fprintf (fp, "#  %d %d %.17g %.17g %.17g\n", child.x, child.y,
                 coarse(f,child.x,child.y), coarse(rhoGv_G,child.x,child.y),
                 coarse(rhoGv_S,child.x,child.y));
        fprintf (fp, "# neighbourhood: di dj f porosity rhoGv_G rhoGv_S\n");
        for (int di = -1; di <= 1; di++)
          for (int dj = -1; dj <= 1; dj++)
            fprintf (fp, "#  %d %d %.17g %.17g %.17g %.17g\n",
                     di, dj, f[di,dj], porosity[di,dj],
                     rhoGv_G[di,dj], rhoGv_S[di,dj]);
      }
    fclose (fp);
    fprintf (stderr, "PROPERTIES_TRAP: %d bad cells at i=%d t=%g,"
                     " written to rhomix-post.dat\n", nbad, i, t);
    fflush (stderr);
    exit (1);
  }
#endif
}

/**
## Useful functions

We define functions that are useful for variable properties
simulations.
*/

/**
### *check_termostate()*: check that the thermodynamic state is
reasonable. */

int check_thermostate (ThermoState * ts, int NS) {
  double sum = 0.;
  for (int jj=0; jj<NS; jj++)
    sum += ts->x[jj];

  int T_ok = (ts->T > 180. && ts->T < 4000.) ? true : false;
  int P_ok = (ts->P > 1e3 && ts->P < 1e7) ? true : false;
  int X_ok = (sum > 1.-1.e-3 && sum < 1.+1.e-3) ? true : false;

  return T_ok*P_ok*X_ok;
}

/**
### *print_thermostate()*: print the thermodynamic state of the mixture.
*/

void print_thermostate (ThermoState * ts, int NS, FILE * fp = stdout) {
  fprintf (fp, "Temperature = %g - Pressure = %g\n", ts->T, ts->P);
  for (int jj=0; jj<NS; jj++)
    fprintf (fp, "  Composition[%d] = %g\n", jj, ts->x[jj]);
  fprintf (fp, "\n");
}

/**
### *gasprop_thermal_expansion()*: Thermal expansion coefficient of an ideal gas
*/

double gasprop_thermal_expansion (ThermoState * ts) {
  return ts->T > 0. ? 1./ts->T : 0.;
}
