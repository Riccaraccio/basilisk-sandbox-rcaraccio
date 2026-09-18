/**
# Check of the thermal expansion of the pore gas

This header measures the defect of section 2 of
`~/discretization-report/discretization-consistency.md`. The thermal part of
`divu1` in `update_divergence()` (`multicomponent-properties.h`) is

$$
d_{\text{code}} = \frac{Q}{T\,(\rho c)_{\text{eff}}},
\qquad
(\rho c)_{\text{eff}} = \varepsilon\rho_g c_{p,g} + (1-\varepsilon)\rho_s c_s ,
$$

with $Q$ = `DTDtS`. Only the pore gas expands, so the correct term is
$d_{\text{fix}} = \varepsilon\,d_{\text{code}}$. The header also measures the
expansion of the pore gas from the temperature that the solver computed:

$$
d_{\text{meas}} = \varepsilon\,\frac{T^{n+1} - T^{n}}{\Delta t\,T^{n}}\,c_m .
$$

It integrates the three terms over the cells of pure solid ($f = 1$ at the
start and at the end of the step). In these cells section 3 (the double
weight of $f$) has no effect, so the ratios show section 2 alone:

- `r_code` = $\int d_{\text{code}} / \int d_{\text{meas}}$. Without the fix it
  is about $1/\varepsilon$ (5 at the start).
- `r_fix` = $\int d_{\text{fix}} / \int d_{\text{meas}}$. It must be about 1.

The code approximates $DT/Dt$ by $Q/(\rho c)_{\text{eff}}$. $T^{n+1} - T^{n}$
also contains the advection of `TS` by `u_prime`. So `r_fix` is not exactly 1.
The difference between `r_code` and `r_fix` is the defect.

## Usage

Include the header after `multicomponent-varprop.h` and `darcy.h`:

~~~literatec
#include "divt-check.h"
~~~

Then restart the case from `last-snapshot`. The header writes
`divt-check.dat` in the working directory. Set `DIVT_EVERY` to write every
N steps. The fields of the header are not dumped, so a dump of this build
restores in a build without the header.

Columns: 1 `i`, 2 `t`, 3 `dt`, 4 number of cells, 5 volume of the cells,
6 mean porosity (volume weight), 7 porosity that the heat rate weights
(`I_fix/I_code`), 8 `I_code`, 9 `I_fix`, 10 `I_meas`, 11 `I_drhodt`,
12 `r_code`, 13 `r_fix`.

`I_drhodt` is $\int -$`drhodt` over the same cells. It contains the species
term of `divu1` too. With `TURN_OFF_REACTIONS` that term holds only the
diffusion, and `I_drhodt` is close to `I_code`.

All the integrals carry `cm` and are per radian in `AXI`, as `drhodt` is. */

#ifndef DIVT_EVERY
# define DIVT_EVERY 1
#endif

scalar divt_TS0[], divt_f0[];

event defaults (i = 0) {
  divt_TS0.nodump = true;
  divt_f0.nodump = true;
}

/**
The `stability` event of this header registers after the other headers, so
it runs first in its chain. No module changes `TS` or `f` between `adapt` of
the previous step and this event. `TS` and `f` are in tracer form here, and
in a cell with $f = 1$ the tracer form equals the intrinsic value. */

event stability (i++) {
  if (i % DIVT_EVERY == 0)
    foreach() {
      divt_TS0[] = TS[];
      divt_f0[] = f[];
    }
}

/**
The function reads `DTDtS`, `drhodt`, the properties and `TS` of the current
step. Call it after the projection and before `adapt`. The accumulators then
hold the values that `update_divergence()` used. `reset_sources` sets them to
zero at the start of the next step. */

void divt_check (FILE * fp) {
  double I_code = 0., I_fix = 0., I_meas = 0., I_drhodt = 0.;
  double vol = 0., epsvol = 0.;
  int n = 0;

  foreach (reduction(+:I_code) reduction(+:I_fix) reduction(+:I_meas)
           reduction(+:I_drhodt) reduction(+:vol) reduction(+:epsvol)
           reduction(+:n)) {
    if (f[] > 1. - F_ERR && divt_f0[] > 1. - F_ERR && divt_TS0[] > 0.) {
      double dV = pow (Delta, dimension);
      double eps = porosity[];    // tracer form, and f = 1
      double rc = rhoGv_S[]*cpGv_S[]*eps + rhoSv[]*cpSv[]*(1. - eps);
      double T0 = divt_TS0[];
      double dcode = rc > 0. ? DTDtS[]/(T0*rc) : 0.;

      I_code   += dcode*dV;
      I_fix    += eps*dcode*dV;
      I_meas   += eps*cm[]*(TS[] - T0)/(dt*T0)*dV;
      I_drhodt += -drhodt[]*dV;
      vol      += cm[]*dV;
      epsvol   += eps*cm[]*dV;
      n++;
    }
  }

  if (pid() == 0) {
    double nan = NAN;   // not 0./0.: that division trips the FPE traps
    fprintf (fp, "%d %g %g %d %g %g %g %g %g %g %g %g %g\n",
             iter, t, dt, n, vol,
             vol > 0. ? epsvol/vol : nan,
             I_code != 0. ? I_fix/I_code : nan,
             I_code, I_fix, I_meas, I_drhodt,
             I_meas != 0. ? I_code/I_meas : nan,
             I_meas != 0. ? I_fix/I_meas : nan);
    fflush (fp);
  }
}

event end_timestep (i++) {
  if (i % DIVT_EVERY == 0) {
    static FILE * fp = NULL;
    if (!fp && pid() == 0) {
      fp = fopen ("divt-check.dat", "a");
      fprintf (fp, "# i t dt ncell vol eps_mean eps_heat"
               " I_code I_fix I_meas I_drhodt r_code r_fix\n");
    }
    divt_check (fp);
  }
}
