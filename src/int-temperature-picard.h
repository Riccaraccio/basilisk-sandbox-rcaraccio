/**
# The record of the Picard loop on the interface temperature

The loop lives in `multicomponent-varprop.h`, in the `tracer_diffusion` event.
This file holds only the counters that the loop fills and the event that
writes them. It changes no field.

The case turns the loop on with `-DINT_TEMP_PICARD=1`. Three more flags tune
it, and each has a default:

- `INT_TEMP_PICARD_MAXITER` (5) the largest number of passes. A value of 0
  gives the code without the loop, exactly. Use it as the inertness control.
- `INT_TEMP_PICARD_TOL` (1e-2 K) the loop stops when the largest change of
  `TInt` in one pass is below this.
- `INT_TEMP_PICARD_OMEGA` (1, no under-relaxation) the loop applies
  `TInt = w*TInt_new + (1-w)*TInt_old`. Lower it if `picard.dat` shows that
  `dTInt` does not fall from pass to pass.

## How to read `picard.dat`

    #t(1) dt(2) niter(3) dTInt(4) dTInt0(5) nint(6) nskip(7) rel_max(8)

- `niter` the passes that the step used. A value equal to `MAXITER` means the
  step did not reach the tolerance.
- `dTInt` the largest change of `TInt` over the LAST pass, in kelvin.
- `dTInt0` the same, over the FIRST pass. The contraction rate of the map is
  `(dTInt/dTInt0)^(1/(niter-1))`. A rate well below 1 means the map contracts
  and the loop is sound. A rate near or above 1 means it does not, and the
  answer is `INT_TEMP_PICARD_OMEGA`, not more passes.

  Caution: compare the two within one step. A ratio taken between two steps
  measures how the flow changed, not how the loop converged.
- `nskip` the interface cells that `ijc_CoupledTemperature()` did not solve,
  because `TS` or `TG` was not positive there. Such a cell keeps its old
  `TInt` and adds nothing to `dTInt`.

  Caution: a step whose `nskip` is not zero has cells that the loop never
  touched. Do not call that step converged, whatever `dTInt` says.
- `rel_max` the normalised interface flux imbalance of
  `int-temperature-probe.h`, or -1 when `INT_TEMP_PROBE` is off. This is the
  quantity the loop exists to reduce. Compare it against a build with
  `MAXITER = 0` on the same case.
*/

/**
The loop writes these. `event picard_output` reads them. */

double ITP_niter, ITP_dTInt, ITP_dTInt0, ITP_nskip;

/**
`int-temperature-probe.h` declares `ITP_nint` too, and it writes the same
count. Declare it here only when that file is absent. */

#if !INT_TEMP_PROBE
double ITP_nint;
#endif

/**
One line per step, at the end of the step. */

event picard_output (i++, last) {
  static FILE * fp = NULL;
  if (!fp) {
    fp = fopen ("picard.dat", restarted ? "a" : "w");
    if (!fp) {
      fprintf (stderr, "Error opening picard.dat\n");
      return 0;
    }
    if (!restarted)
      fprintf (fp, "#t(1) dt(2) niter(3) dTInt(4) dTInt0(5) nint(6)"
                   " nskip(7) rel_max(8)\n");
  }

  double relmax = -1.;
#if INT_TEMP_PROBE
  relmax = ITP_rel_max;
#endif

  fprintf (fp, "%g %g %g %g %g %g %g %g\n",
           t, dt, ITP_niter, ITP_dTInt, ITP_dTInt0,
           ITP_nint, ITP_nskip, relmax);
  fflush (fp);

  return 0;
}
