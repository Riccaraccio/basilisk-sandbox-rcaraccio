/**
# The interface conductance, and how it keeps the heat

`INT_TEMP_ROBIN` puts a conductance `K` on the diagonal of each temperature
solve. The block that builds `K` is in `multicomponent-varprop.h`. This file
holds only the counters and the event that writes them.

## Why the plain conductance loses heat

`ebmgrad` builds the interface gradient from the NEIGHBOURS of a cut cell, so
the source that heats the cell does not answer to the temperature of that
cell. In a sliver the heat capacity `theta` goes to zero with the phase
fraction, so one step delivers far more heat than the cell can hold. The
exchange number

    S = dt*lambda*h*aov/theta

reached 100 in a measured run, and the gas temperature of one cut cell went
through zero.

The conductance damps that. It sets `K = max (0, A/SMAX - theta/dt)`, so when
`K` fires the diagonal becomes exactly `A/SMAX` and the change of the step is

    dT_capped / dT_correct = (theta/dt) / (A/SMAX) = SMAX / S

With `SMAX = 1` and `S = 100` the cell keeps one per cent of the heat that the
interface gave it.

**The other 99 per cent does not come back.** Sum the discrete equation over
the domain: the added term `-K(T^{n+1} - T^n)` cancels only when the field
stops changing. While the particle heats, the conductance is a heat sink in
proportion to the heating rate. A measured run lost 20 to 25 K of surface
temperature and 25 per cent of the mass loss rate.

## The repair: give the heat back

Keep the withheld heat in a field and add it to the source of the next step:

    (theta/dt + K)(T^{n+1} - T^n) = div(D grad T^{n+1}) + q + D^n
    D^{n+1} = K*(T^{n+1} - T^n)

`D` is a rate, in the units of the source `q`. Two properties follow.

1. **The steady heating rate is correct.** Put `T^{n+1} - T^n = dT` for every
   step. Then `(theta/dt + K) dT = q + K dT`, so `(theta/dt) dT = q`. The `K`
   terms cancel and the cell heats at the rate the interface demands. The
   damping now costs a lag, not energy.
2. **The scheme conserves energy.** Over `N` steps the stored energy is
   `theta (T^N - T^0) = sum (dt*q) - K*dt*dT^{N-1}`. The only term that does
   not cancel is the heat withheld on the LAST step. `D` is a one-step carry,
   not a running sum, so the energy in flight is bounded by one step and goes
   to zero when the field settles.

The debt channel has gain `K/(theta/dt + K) = 1 - SMAX/S`, which is below one,
so it cannot amplify. It relaxes the cell toward the correct rate over about
`S/SMAX` steps.

Because `D` is a one-step carry, the loss of it costs one step. So the fields
are `nodump`, and a restart or a cell that leaves the interface loses at most
one step of heat.

`INT_TEMP_ROBIN_DEBT` selects the behaviour. It is 1 by default, which keeps
the heat. Set it to 0 to reproduce the lossy scheme for a comparison; do not
run physics that way.

## Measured: the debt is necessary, and it is not the cure

Five runs of the `full` ladder, matched, at t = 0.3 s. `Tsurf` against the
same case built without `INT_TEMP_ROBIN`:

| build | Tsurf | ur at 1.5 mm | cells clamped |
|---|---|---|---|
| no conductance | 378.266 | -0.0122 | - |
| `SMAX = 1e30` (K never fires) | 378.266 (+0.000) | -0.0122 | 0 |
| `SMAX = 20`, debt on | 378.266 (+0.000) | -0.0122 | 2 |
| `SMAX = 1`, debt off | 376.108 (-2.158) | **+0.0066** | 63 |
| `SMAX = 1`, debt on | 376.161 (-2.105) | **+0.0066** | 63 |

The debt gives back the heat, and it recovers 0.053 K of a 2.158 K error.
So the withheld heat was real but it was never the thing that broke the
answer. What breaks the answer is the clamp itself: at `SMAX = 1` the cut
cell may change by only `1/S` of what the physics asks, and no bookkeeping of
the withheld energy restores a response the cell was not allowed to make.

**`SMAX` is the control that matters.** At 20 the conductance fires on 2 cells
instead of 63 and the answer is identical to the run without it. Keep the debt
anyway: it costs six fields and one `foreach`, and it removes a heat sink that
would otherwise grow without bound whenever the clamp does fire.

## How to read `robin-debt.dat`

    #t(1) dt(2) ncap(3) debt_max(4) debt_l1(5) dTcap_max(6)

- `ncap` the interface cells whose conductance fired, that is `S > SMAX`.
  When this is zero the flag changes nothing at all.
- `debt_max` the largest withheld rate, in the units of the heat source.
- `debt_l1` the withheld energy rate over the whole domain. It must stay
  bounded. A value that grows without limit means the cells never catch up,
  and `INT_TEMP_ROBIN_SMAX` is too small for the step.
- `dTcap_max` the largest `|T^{n+1} - T^n|` over the clamped cells. This is
  the change the conductance actually allowed.

Caution: this file shows that the heat is carried, not that the answer is
right. The proof is a comparison of `Tsurf` and `mdot` against the same case
built without `INT_TEMP_ROBIN`. The deficit must collapse.
*/

double ITR_debt_max, ITR_debt_l1, ITR_dTcap_max, ITR_ncap;

event robin_debt_output (i++, last) {
  static FILE * fp = NULL;
  if (!fp) {
    fp = fopen ("robin-debt.dat", restarted ? "a" : "w");
    if (!fp) {
      fprintf (stderr, "Error opening robin-debt.dat\n");
      return 0;
    }
    if (!restarted)
      fprintf (fp, "#t(1) dt(2) ncap(3) debt_max(4) debt_l1(5) dTcap_max(6)\n");
  }

  fprintf (fp, "%g %g %g %g %g %g\n",
           t, dt, ITR_ncap, ITR_debt_max, ITR_debt_l1, ITR_dTcap_max);
  fflush (fp);

  return 0;
}
