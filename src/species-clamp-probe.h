/**
# Probe TL-5: the clamp of the gas species after the implicit solves

Item TL-5 of `~/discretization-report/time-level-review.md`, section 3.5 and
the end of section 4(b). This header changes no field in the run. It writes
one file.

## The question

The `tracer_diffusion` event of `multicomponent-varprop.h` solves each gas
species with an implicit `diffusion()`. The interface flux of a cut cell goes
into the solve as the explicit source `sSexp` (pore side) or `sGexp` (gas
side). The code builds that source from the state after the advection and
from `YGInt` of that state. Nothing bounds the exchange number

    S_j = dt*rho*D_j*h*a/theta

where `h` is the slope of the interface gradient against the interface
value, `a` the interface area per unit volume and `theta` the capacity of the
solve. A large `S_j` lets the explicit source push `Y_j` out of `[0,1]`.
`check_and_correct_fractions()` then sets each value below `Y_CUT` to 0 and
scales the sum back to 1. That changes the mass of each species, and of the
pore gas of the cell.

The probe measures three things on one step:

1. How far `Y` leaves `[0,1]` after the solves and before the clamp.
2. The mass that the clamp adds or removes.
3. The size of the exchange number, and the spread of `MWmixG_S`, in the cut
   cells.

The event calls `check_and_correct_fractions()` at its start too, so every
`Y` is in `[0,1]` before the solves. The only other step that changes `Y`
between that clamp and the probe is the corrective flux of `FICK_CORRECTED`.
With `FICK_CORRECTED` off, the solves alone make the values that the probe
counts.

## The two sides and the two regions

* Gas side `G`: the list `YGList_G` and the capacity `theta2`. The phase mass
  of a cell is `rhoG*(1-f)*dv()`.
* Pore side `S`: the list `YGList_S` and the capacity `theta1`. The phase mass
  of a cell is `rhoGS*porosity*dv()`. `porosity` is in tracer form here, so it
  is `eps*f`.

`rhoG` and `rhoGS` are `rhoGv_G` and `rhoGv_S` with `VARPROP`, and the
constant `rhoG` without it. Each side has two regions:

* `C` the cut cells, `F_ERR < f < 1-F_ERR`.
* `R` the rest: the other cells where the phase exists. For the gas side that
  is `f <= F_ERR`, for the pore side `f >= 1-F_ERR`. The clamp sets every `Y`
  of a cell without the phase to 0, and the recovery of the tracer form has
  already set it to 0 there, so such a cell adds nothing.

## What the probe measures

The value of a species in a cell is the intrinsic value `y = Y/(1-f)` (gas)
or `y = Y/f` (pore), computed as the clamp computes it. The probe repeats the
arithmetic of the clamp to get the value after it, and it does not write it.

* `nneg` the number of cells of the region with at least one `y_j < 0`.
* `nover` the number of cells of the region with at least one `y_j > 1`.
* `ymin`, `ymax` the smallest and the largest `y_j` of the region, over all
  species. `ymin < 0` or `ymax > 1` is the worst violation. 0 when the region
  has no cell.
* `dM` the signed mass that the clamp adds to the region, summed over the
  species, in kg. The mass of species `j` in a cell is `rhoG*Y_j*dv()` (gas)
  or `rhoGS*(porosity/f)*Y_j*dv()` (pore), with `Y_j` in tracer form. A
  positive value is a gain.
* `dMabs` the same sum of `|change of the mass of species j|`. A clamp that
  moves mass from one species to another gives `dM` near 0 and a large
  `dMabs`.
* `MC` the mass of the phase in the cut cells, in kg (definitions above).
* `rel` = `dMabsC/MC`, the relative mass change of the step in the cut cells.
* `Smax` the largest exchange number `S_j` of the cut cells, over all
  species. The probe builds it as the Robin conductance of
  `int-temperature.h` builds the thermal one:

      gas:  S_j = dt*rhoG*D_j*h*aov/max(fG*rhoG, F_ERR)
      pore: S_j = dt*rhoGS*D_j*h*aov/max(rhoGS*porosity, F_ERR)

  `h = |ebmgrad(1) - ebmgrad(0)|` is the exact slope of the interface
  gradient of that side against the interface value. `aov` is the factor of
  the source build: `area/Delta`, or `area*(y + p.y*Delta)/(Delta*y)` in an
  axisymmetric case. `D_j` is the cell value of `DmixGList_G` or
  `DmixGList_S`. The capacity is the one of the solve, without `cm`, because
  the source carries `cm` as well. The probe computes `Smax` before the
  solves.
* `dYsrc` the largest change of `y_j` that the explicit interface source
  alone gives in one step in a cut cell, `dt*|sexp_j|/theta`. It is read
  before the solves, because `diffusion()` overwrites the source. A value
  above the local `y_j` can push `y_j` below 0.
* `MWmin`, `MWmax`, `MWmean` the spread of `MWmixG_S` over the cut cells,
  before the solves, in kg/kmol. The value 0 of a cell that the properties
  did not fill counts in the spread.
* `nC` the number of cut cells.

In an axisymmetric case `dv()` carries `cm`, so each mass is a mass per
radian. The ratio `rel` does not depend on it.

## How to read `speciesclamp.dat`

    #t(1) i(2) dt(3)
     nnegGC(4) noverGC(5) yminGC(6) ymaxGC(7) dMGC(8) dMabsGC(9) MGC(10) relGC(11)
     nnegGR(12) noverGR(13) yminGR(14) ymaxGR(15) dMGR(16) dMabsGR(17)
     SmaxG(18) dYsrcG(19)
     nnegSC(20) noverSC(21) yminSC(22) ymaxSC(23) dMSC(24) dMabsSC(25) MSC(26) relSC(27)
     nnegSR(28) noverSR(29) yminSR(30) ymaxSR(31) dMSR(32) dMabsSR(33)
     SmaxS(34) dYsrcS(35)
     MWmin(36) MWmax(37) MWmean(38) nC(39)

The file holds one line for each output time. The line covers the step that
starts at that time. Columns 18, 19, 34 to 39 come from the state before the
solves, the other columns from the state after the solves and before the
clamp.

Caution: the mass of one step is small. Multiply `dM` by the number of steps
per second (1/`dt`) to compare it with a rate, for example with `mdot`.

Caution: `Smax` uses the Fick form. With `MOLAR_DIFFUSION` the source uses
the gradient of the mole fraction, and `S_j` changes by the factor
`MW_j/MWmix*dX_j/dY_j`, which is of the order of 1. The `FICK_CORRECTED`
term is not in `Smax` either.

## Flags

* `SPECIES_CLAMP_PROBE` 1 turns the probe on. The default is 0. Set the flag
  before `multicomponent-varprop.h`, which includes this header and calls its
  two functions.
* `SPECIES_CLAMP_PROBE_FILE` the name of the file.

## Cost

No field. On a step that the probe covers: one sweep of the cut cells before
the solves, with four calls of `ebmgrad()` per cut cell, and one sweep of the
grid after the solves. The probe covers the step that starts at an output
time, one step in about 20 at `DT = 5e-4`. On the other steps it tests one
flag. */

#ifndef SPECIES_CLAMP_PROBE
# define SPECIES_CLAMP_PROBE 0
#endif

#if SPECIES_CLAMP_PROBE

#include <time.h>

#ifndef SPECIES_CLAMP_PROBE_FILE
# define SPECIES_CLAMP_PROBE_FILE "speciesclamp.dat"
#endif

static bool scp_armed = false, scp_want = false;
static double scp_cpu = 0.;
static int scp_ncall = 0;
static double scp_SmaxG = 0., scp_SmaxS = 0., scp_dYG = 0., scp_dYS = 0.;
static double scp_MWmin = 0., scp_MWmax = 0., scp_MWmean = 0., scp_nC = 0.;

/**
The arm. An event with a time condition and a new name runs at the start of
the step, before `stability`. It asks for the step that starts at the output
time. The interval is the one of the output events of `run/test.c`, so the
event adds no new event time and does not change `dt`. */

event species_clamp_arm (t += 0.01) {
  scp_want = true;
}

/**
Call this after the interface sources `sSexp` and `sGexp` are complete and
before the first species solve. */

static void species_clamp_presolve (void)
{
  scp_armed = (dt > 0.) && scp_want;
  scp_want = false;
  if (!scp_armed)
    return;
  clock_t c0 = clock();

  double SmaxG = 0., SmaxS = 0., dYG = 0., dYS = 0.;
  double MWmin = HUGE, MWmax = -HUGE, MWsum = 0., nC = 0.;
  scalar YS0 = YGList_S[0], YG0 = YGList_G[0];

  foreach (reduction(max:SmaxG) reduction(max:SmaxS)
           reduction(max:dYG) reduction(max:dYS)
           reduction(min:MWmin) reduction(max:MWmax)
           reduction(+:MWsum) reduction(+:nC)) {
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      bool success = false;
      coord n = interface_source_normal (point, fS, fsS), p;
      double alpha = plane_alpha (fS[], n);
      double area = plane_area_center (n, alpha, &p);
#ifdef AXI
      double aov = area*(y + p.y*Delta)/(Delta*y);
#else
      double aov = area/Delta;
#endif
      double hS = fabs (ebmgrad (point, YS0, fS, fG, fsS, fsG, false, 1., &success)
                        - ebmgrad (point, YS0, fS, fG, fsS, fsG, false, 0., &success));
      double hG = fabs (ebmgrad (point, YG0, fS, fG, fsS, fsG, true, 1., &success)
                        - ebmgrad (point, YG0, fS, fG, fsS, fsG, true, 0., &success));
#ifdef VARPROP
      double rS = rhoGv_S[], rG = rhoGv_G[];
#else
      double rS = rhoG, rG = rhoG;
#endif
      double capS = max (rS*porosity[], F_ERR);
      double capG = max (fG[]*rG, F_ERR);

      for (int jj = 0; jj < NGS; jj++) {
        scalar DS = DmixGList_S[jj], DG = DmixGList_G[jj];
        scalar sS = sSexpList[jj], sG = sGexpList[jj];
        SmaxS = max (SmaxS, dt*rS*DS[]*hS*aov/capS);
        SmaxG = max (SmaxG, dt*rG*DG[]*hG*aov/capG);
        dYS = max (dYS, dt*fabs (sS[])/(cm[]*capS));
        dYG = max (dYG, dt*fabs (sG[])/(cm[]*capG));
      }

      MWmin = min (MWmin, MWmixG_S[]);
      MWmax = max (MWmax, MWmixG_S[]);
      MWsum += MWmixG_S[];
      nC += 1.;
    }
  }

  scp_SmaxG = SmaxG; scp_SmaxS = SmaxS;
  scp_dYG = dYG; scp_dYS = dYS;
  scp_nC = nC;
  scp_MWmin = (nC > 0.) ? MWmin : 0.;
  scp_MWmax = (nC > 0.) ? MWmax : 0.;
  scp_MWmean = (nC > 0.) ? MWsum/nC : 0.;
  scp_cpu += (double)(clock() - c0)/CLOCKS_PER_SEC;
}

/**
The effect of the clamp on one side of one cell. `YList` holds the tracer
form, `denom` is `f` or `1-f`, and `wm` converts a tracer value into a mass
per unit volume. The arithmetic is that of `check_and_correct_fractions()`,
step by step, so the value after the clamp is the same to the last bit. */

static void scp_cell (Point point, scalar * YList, double denom, double wm,
                      bool * neg, bool * over, double * ymin, double * ymax,
                      double * dm, double * dmabs)
{
  double inv_denom = 1./denom;
  double sum = 0.;
  *neg = false; *over = false;
  *dm = 0.; *dmabs = 0.;
  for (int jj = 0; jj < NGS; jj++) {
    scalar Y = YList[jj];
    double val = Y[]*inv_denom;
    if (val < 0.) *neg = true;
    if (val > 1.) *over = true;
    if (val < *ymin) *ymin = val;
    if (val > *ymax) *ymax = val;
    sum += (val < Y_CUT) ? 0. : val;
  }
  double scale = (sum > F_ERR) ? denom/sum : 0.;
  for (int jj = 0; jj < NGS; jj++) {
    scalar Y = YList[jj];
    double val = Y[]*inv_denom;
    double ynew = (sum > F_ERR) ? ((val < Y_CUT) ? 0. : val)*scale : 0.;
    double d = wm*(ynew - Y[])*dv();
    *dm += d;
    *dmabs += fabs (d);
  }
}

/**
Call this after the species solves and the recovery of the tracer form, and
before `check_and_correct_fractions()`. */

static void species_clamp_postsolve (void)
{
  if (!scp_armed)
    return;
  scp_armed = false;
  clock_t c0 = clock();

  double nnegGC = 0., noverGC = 0., yminGC = HUGE, ymaxGC = -HUGE;
  double dMGC = 0., dMabsGC = 0., MGC = 0.;
  double nnegGR = 0., noverGR = 0., yminGR = HUGE, ymaxGR = -HUGE;
  double dMGR = 0., dMabsGR = 0.;
  double nnegSC = 0., noverSC = 0., yminSC = HUGE, ymaxSC = -HUGE;
  double dMSC = 0., dMabsSC = 0., MSC = 0.;
  double nnegSR = 0., noverSR = 0., yminSR = HUGE, ymaxSR = -HUGE;
  double dMSR = 0., dMabsSR = 0.;

  foreach (reduction(+:nnegGC) reduction(+:noverGC)
           reduction(min:yminGC) reduction(max:ymaxGC)
           reduction(+:dMGC) reduction(+:dMabsGC) reduction(+:MGC)
           reduction(+:nnegGR) reduction(+:noverGR)
           reduction(min:yminGR) reduction(max:ymaxGR)
           reduction(+:dMGR) reduction(+:dMabsGR)
           reduction(+:nnegSC) reduction(+:noverSC)
           reduction(min:yminSC) reduction(max:ymaxSC)
           reduction(+:dMSC) reduction(+:dMabsSC) reduction(+:MSC)
           reduction(+:nnegSR) reduction(+:noverSR)
           reduction(min:yminSR) reduction(max:ymaxSR)
           reduction(+:dMSR) reduction(+:dMabsSR)) {
#ifdef VARPROP
    double rS = rhoGv_S[], rG = rhoGv_G[];
#else
    double rS = rhoG, rG = rhoG;
#endif
    bool cut = (f[] > F_ERR && f[] < 1. - F_ERR);
    bool neg, over;
    double dm, dmabs;

    /* gas side */
    double dG = 1. - f[];
    if (dG >= F_ERR) {
      double lo = HUGE, hi = -HUGE;
      scp_cell (point, YGList_G, dG, rG, &neg, &over, &lo, &hi, &dm, &dmabs);
      if (cut) {
        nnegGC += neg; noverGC += over;
        yminGC = min (yminGC, lo); ymaxGC = max (ymaxGC, hi);
        dMGC += dm; dMabsGC += dmabs;
        MGC += rG*dG*dv();
      }
      else {
        nnegGR += neg; noverGR += over;
        yminGR = min (yminGR, lo); ymaxGR = max (ymaxGR, hi);
        dMGR += dm; dMabsGR += dmabs;
      }
    }

    /* pore side */
    double dS = f[];
    if (dS >= F_ERR) {
      double lo = HUGE, hi = -HUGE;
      scp_cell (point, YGList_S, dS, rS*porosity[]/dS,
                &neg, &over, &lo, &hi, &dm, &dmabs);
      if (cut) {
        nnegSC += neg; noverSC += over;
        yminSC = min (yminSC, lo); ymaxSC = max (ymaxSC, hi);
        dMSC += dm; dMabsSC += dmabs;
        MSC += rS*porosity[]*dv();
      }
      else {
        nnegSR += neg; noverSR += over;
        yminSR = min (yminSR, lo); ymaxSR = max (ymaxSR, hi);
        dMSR += dm; dMabsSR += dmabs;
      }
    }
  }

#define SCP_Z(v) ((v) == HUGE || (v) == -HUGE ? 0. : (v))

  if (pid() == 0) {
    static FILE * fp = NULL;
    if (!fp) {
      fp = fopen (SPECIES_CLAMP_PROBE_FILE, restarted ? "a" : "w");
      if (fp == NULL) {
        fprintf (stderr, "Error opening %s\n", SPECIES_CLAMP_PROBE_FILE);
        exit (1);
      }
      if (!restarted)
        fprintf (fp, "#t(1) i(2) dt(3)"
                 " nnegGC(4) noverGC(5) yminGC(6) ymaxGC(7) dMGC(8)"
                 " dMabsGC(9) MGC(10) relGC(11)"
                 " nnegGR(12) noverGR(13) yminGR(14) ymaxGR(15) dMGR(16)"
                 " dMabsGR(17) SmaxG(18) dYsrcG(19)"
                 " nnegSC(20) noverSC(21) yminSC(22) ymaxSC(23) dMSC(24)"
                 " dMabsSC(25) MSC(26) relSC(27)"
                 " nnegSR(28) noverSR(29) yminSR(30) ymaxSR(31) dMSR(32)"
                 " dMabsSR(33) SmaxS(34) dYsrcS(35)"
                 " MWmin(36) MWmax(37) MWmean(38) nC(39)\n");
    }
    fprintf (fp, "%g %d %g"
             " %g %g %g %g %g %g %g %g"
             " %g %g %g %g %g %g %g %g"
             " %g %g %g %g %g %g %g %g"
             " %g %g %g %g %g %g %g %g"
             " %g %g %g %g\n",
             t, iter, dt,
             nnegGC, noverGC, SCP_Z(yminGC), SCP_Z(ymaxGC), dMGC, dMabsGC,
             MGC, MGC > 0. ? dMabsGC/MGC : 0.,
             nnegGR, noverGR, SCP_Z(yminGR), SCP_Z(ymaxGR), dMGR, dMabsGR,
             scp_SmaxG, scp_dYG,
             nnegSC, noverSC, SCP_Z(yminSC), SCP_Z(ymaxSC), dMSC, dMabsSC,
             MSC, MSC > 0. ? dMabsSC/MSC : 0.,
             nnegSR, noverSR, SCP_Z(yminSR), SCP_Z(ymaxSR), dMSR, dMabsSR,
             scp_SmaxS, scp_dYS,
             scp_MWmin, scp_MWmax, scp_MWmean, scp_nC);
    fflush (fp);
  }

#undef SCP_Z

  scp_cpu += (double)(clock() - c0)/CLOCKS_PER_SEC;
  scp_ncall++;
}

/**
The cost of the probe, in the log at the end of the run. */

event species_clamp_cost (t = end) {
  if (pid() == 0)
    fprintf (stderr, "# species-clamp-probe: %g s CPU in %d steps\n",
             scp_cpu, scp_ncall);
}

#endif // SPECIES_CLAMP_PROBE
