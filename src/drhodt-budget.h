/**
# Probe TL-2: the explicit fluxes of `drhodt` against the implicit solves

Item TL-2 of `~/discretization-report/time-level-review.md`, section 3.2.
This header changes nothing in the run. It writes one file.

## The question

`update_divergence()` (`multicomponent-properties.h`) runs before the two
temperature solves. It builds the transport part of `drhodt` from EXPLICIT
fluxes of the state after the advection, which the review calls `T**`:

    DTDtG_transport = sum_d (lambda2f.x[1]*grad T** - lambda2f.x[]*grad T**)/Delta + sGT

The solve that follows is implicit. `diffusion()` gives

    theta2*(T^{n+1} - T**)/dt = div (lambda2f grad T^{n+1}) + beta*T^{n+1} + sGT

So the heat that the solve really moves is `theta2*(T^{n+1} - T**)/dt`, and
the projection enforces an expansion that is built from another number. For
one Fourier mode the ratio of the two forms is `1 + alpha*dt*k^2`. At the
flame `alpha*dt/Delta^2` is about 9, so the local error can be as large as
the expansion of the chemistry. Its integral over the domain is near zero,
so the columns `X` and `Qdiv` of `expansion.dat` do not show it.

The probe measures the difference of the two forms, cell by cell, and gives
its size against `|drhodt|`.

## What the probe measures

The probe converts the difference of the two forms into the units of
`drhodt`, with the same factors that `update_divergence()` uses:

    dd = -[ w_S*eps/(TS***(rhoGv_S*cpGv_S*eps + rhoSv*cpSv*(1-eps)))*(implS - expS)
          + w_G/(TG***rhoGv_G*cpGv_G)*(implG - expG) ]

with `eps = porosity/f`. `w_S` and `w_G` are the weights of
`update_divergence()`: with `DRHODT_CELL_AVERAGE` they are 1 where the phase
exists and 0 elsewhere, and without it they are `f` and `1-f`. `impl` is the
implicit form `theta*(T^{n+1} - T**)/dt` and `exp` is the explicit form of
the code. The divisors use `T**`, because `update_divergence()` runs at that
state.

The probe gives, for the whole domain and for three regions:

* `n` the number of cells of the region,
* `L1d` the integral of `|dd|` over the region,
* `L1r` the integral of `|drhodt|` over the same region,
* `maxd` the largest `|dd|` of the region,
* `maxq` the largest `|dd|/|drhodt|` of one cell of the region. The ratio
  covers the cells whose `|drhodt|` is above `DRHODT_BUDGET_QFLOOR` times the
  largest `|drhodt|` of the domain, 1e-3 by default. Without that floor a cell
  with no expansion at all gives a ratio of any size.

`L1d/L1r` is the relative size of the error in the L1 norm. The prediction of
the review is that it grows with `dt`, and that it reaches the size of the
chemistry part at the flame.

The three regions are:

* `H` the hot gas: the cells with `1-f > F_ERR` and an intrinsic gas
  temperature above `DRHODT_BUDGET_THOT`, 1500 K by default.
* `C` the cut cells, `F_ERR < f < 1-F_ERR`. A cut cell books all of its
  interface heat itself, and the review expects the largest local error
  there.
* `I` the first `DRHODT_BUDGET_NLAYER` layers of pure gas cells outside the
  cut cells, 3 layers by default. The probe builds the layer index with one
  dilation of the cut cells for each layer.

## How to read `drhodtbudget.dat`

    #t(1) i(2) dt(3) L1d(4) L1r(5) maxd(6) maxr(7)
     nH(8) L1dH(9) L1rH(10) maxdH(11) maxqH(12)
     nC(13) L1dC(14) L1rC(15) maxdC(16) maxqC(17)
     nI(18) L1dI(19) L1rI(20) maxdI(21) maxqI(22)

Columns 4 to 7 cover the whole domain: the two L1 norms, the largest `|dd|`
and the largest `|drhodt|`. The file holds one line for each output time,
and the line covers the step that starts at that time.

Read `L1dH/L1rH` and `L1dI/L1rI` against `dt`. Both must fall with `dt`. If
they do not fall, the error is not the time level of the fluxes.

Caution: the probe reads the state of the last pass of the loop of
`INT_TEMP_PICARD`. With `INT_TEMP_VOFBC` the explicit form carries the
diagonal part `beta*T**`, as `update_divergence()` does.

Caution: the probe uses `lambda1f` and `lambda2f`, which the solves use.
`update_divergence()` builds the same product from `lambda1v`, `fsS` and
`fm`, in another order. The two differ by one bit at most.

Caution: the species part of `drhodt` has the same defect and the probe does
NOT measure it. That part needs the state of every species before the
species solves, thus `2*NGS` more fields.

## The meaning under `DRHODT_IMPLICIT`

With `DRHODT_IMPLICIT` the explicit form is gone from `drhodt`. The
temperature part of `drhodt` is then the implicit form
`theta*(T^{n+1} - T**)/dt` that `multicomponent-varprop.h` computes (the
fields `dri_*`). The probe compares that form against the operator of the
solve, evaluated again at the new state:

    ref = div (lambdaf grad T^{n+1}) + r + beta*T^{n+1} + (interface term)

`r` is the source of the solve (`sST`, `sGT`), which the probe keeps before
the solve. `beta` is `betaST`, `betaGT` under `INT_TEMP_ROBIN`. The interface
term is the one of `plic_flux()` under `INT_TEMP_VOFBC`, with the same sign
as in the residual of `poisson.h`. The equation of the solve says that
`carried = ref`, so `dd` is the residual of the solve plus any error of the
fix: a wrong weight, a wrong divisor, or a wrong time level. `L1d/L1r` must
then be small, of the order of the tolerance of the solve. `L1r` reads the
`drhodt` that the projection receives, thus it includes the temperature part
that the fix adds after the loop of `INT_TEMP_PICARD`.

Without `DRHODT_IMPLICIT` the probe is the same as before, bit for bit.

## Flags

* `DRHODT_BUDGET` 1 turns the probe on. The default is 0. The flag must be
  set before `multicomponent-varprop.h`, which includes this header and
  calls its two functions.
* `DRHODT_BUDGET_EVERY_STEP` 1 covers every step instead of the steps of the
  output times. The default is 0. Caution: this multiplies the cost of the
  probe by about 20.
* `DRHODT_BUDGET_THOT` the temperature of the hot region, 1500 K.
* `DRHODT_BUDGET_NLAYER` the number of gas layers of the region `I`, 3.
* `DRHODT_BUDGET_QFLOOR` the floor of the denominator of `maxq`, 1e-3 of the
  largest `|drhodt|` of the domain.
* `DRHODT_BUDGET_FILE` the name of the file.

## Cost

Six fields in each cell. On a step that the probe covers: one sweep of the
faces, three sweeps of the cells, and one dilation of the layer index for
each layer. The probe covers the step that starts at an output time, one
step in about 20 at `DT = 5e-4`. On the other steps it tests one flag. */

#ifndef DRHODT_BUDGET_EVERY_STEP
# define DRHODT_BUDGET_EVERY_STEP 0
#endif

#ifndef DRHODT_BUDGET_THOT
# define DRHODT_BUDGET_THOT 1500.
#endif

#ifndef DRHODT_BUDGET_NLAYER
# define DRHODT_BUDGET_NLAYER 3
#endif

#ifndef DRHODT_BUDGET_QFLOOR
# define DRHODT_BUDGET_QFLOOR 1e-3
#endif

#ifndef DRHODT_BUDGET_FILE
# define DRHODT_BUDGET_FILE "drhodtbudget.dat"
#endif

/**
`dbg_TSpre` and `dbg_TGpre` hold the state before the solves. `dbg_th1` and
`dbg_th2` hold the two heat capacities, because `diffusion()` overwrites
`theta` in place. `dbg_exS` and `dbg_exG` hold the explicit transport form of
`update_divergence()`. */

#if DRI_ON
scalar dbg_exS[], dbg_exG[];
# if INT_TEMP_ROBIN
scalar dbg_bS[], dbg_bG[];
# endif
#else
scalar dbg_TSpre[], dbg_TGpre[], dbg_th1[], dbg_th2[], dbg_exS[], dbg_exG[];
#endif

#include <time.h>

static bool dbg_armed = false, dbg_want = false;
static double dbg_cpu = 0.;
static int dbg_ncall = 0;

/**
The arm. An event with a time condition and a new name runs at the start of
the step, before `stability`. It asks for the step that starts at the output
time. The interval is the one of the output events of `run/test.c`, so the
event adds no new event time and does not change `dt`. */

event drhodt_budget_arm (t += 0.01) {
  dbg_want = true;
}

event defaults (i = 0) {
#if DRI_ON
  dbg_exS.nodump = true; dbg_exG.nodump = true;
# if INT_TEMP_ROBIN
  dbg_bS.nodump = true; dbg_bG.nodump = true;
# endif
#else
  dbg_TSpre.nodump = true; dbg_TGpre.nodump = true;
  dbg_th1.nodump = true; dbg_th2.nodump = true;
  dbg_exS.nodump = true; dbg_exG.nodump = true;
#endif
}

/**
Call this after `theta1` and `theta2` are complete and before the first
solve. */

static void drhodt_budget_presolve (scalar theta1, scalar theta2)
{
  dbg_armed = (dt > 0.) && (dbg_want || DRHODT_BUDGET_EVERY_STEP);
  dbg_want = false;
  if (!dbg_armed)
    return;
  clock_t c0 = clock();

#if DRI_ON

  /**
  Keep the sources of the solves: `diffusion()` overwrites them. */

  foreach() {
    dbg_exS[] = sST[];
    dbg_exG[] = sGT[];
# if INT_TEMP_ROBIN
    dbg_bS[] = betaST[];
    dbg_bG[] = betaGT[];
# endif
  }
#else
  face vector qS[], qG[];
  foreach_face() {
    qS.x[] = lambda1f.x[]*face_gradient_x (TS, 0);
    qG.x[] = lambda2f.x[]*face_gradient_x (TG, 0);
  }

  foreach() {
    double dS = 0., dG = 0.;
    foreach_dimension() {
      dS += (qS.x[1] - qS.x[])/Delta;
      dG += (qG.x[1] - qG.x[])/Delta;
    }
    dbg_exS[] = dS + sST[];
    dbg_exG[] = dG + sGT[];
#if INT_TEMP_VOFBC
    dbg_exS[] += betaST[]*TS[];
    dbg_exG[] += betaGT[]*TG[];
#endif
    dbg_TSpre[] = TS[];
    dbg_TGpre[] = TG[];
    dbg_th1[] = theta1[];
    dbg_th2[] = theta2[];
  }
#endif
  dbg_cpu += (double)(clock() - c0)/CLOCKS_PER_SEC;
}

/**
Call this after the two temperature solves, and before any line that changes
`TS` or `TG` again. */

static void drhodt_budget_postsolve (void)
{
  if (!dbg_armed)
    return;
  dbg_armed = false;
  clock_t c0 = clock();

  /**
  The layer index. It is 1 in the cut cells, then 2, 3, ... outward in the
  pure gas cells. One dilation of the 3x3 stencil adds one layer. */

  /**
  Caution: declare the temporary field `lnew` here, outside the loop. `qcc`
  does not free a field that a `for` body declares, so a declaration inside
  the loop leaks one field for each layer at each call. The leak grew the
  field count from 195 to 3234 in 1012 steps, and each new field
  reallocates every cell of the grid: the run became two times slower. The
  `var` count of the last line of the log shows such a leak. */

  scalar lay[], lnew[];
  foreach()
    lay[] = (f[] > F_ERR && f[] < 1. - F_ERR) ? 1. : 0.;

  for (int k = 0; k < DRHODT_BUDGET_NLAYER; k++) {
    foreach() {
      lnew[] = lay[];
      if (lay[] == 0. && f[] <= F_ERR) {
        double m = 0.;
        foreach_neighbor (1)
          if (lay[] > m)
            m = lay[];
        if (m > 0.)
          lnew[] = m + 1.;
      }
    }
    foreach()
      lay[] = lnew[];
@if _MPI
    boundary ({lay});
@endif
  }

  /**
  The scale of `drhodt` of this step. It gives the floor of the denominator of
  `maxq`. */

#if DRI_ON

  /**
  The operator of the solve at the new state. It goes into `dbg_exS` and
  `dbg_exG`, on top of the sources. */

  face vector gS[], gG[];
  foreach_face() {
    gS.x[] = lambda1f.x[]*face_gradient_x (TS, 0);
    gG.x[] = lambda2f.x[]*face_gradient_x (TG, 0);
  }
  foreach() {
    foreach_dimension() {
      dbg_exS[] += (gS.x[1] - gS.x[])/Delta;
      dbg_exG[] += (gG.x[1] - gG.x[])/Delta;
    }
# if INT_TEMP_ROBIN
    dbg_exS[] += dbg_bS[]*TS[];
    dbg_exG[] += dbg_bG[]*TG[];
# endif
  }
# if INT_TEMP_VOFBC
  plicbc_phase (fS, fsS);
  foreach() {
    double c, e = plic_flux (point, TS, lambda1f, &c);
    dbg_exS[] += e*TS[] - c;
  }
  plicbc_phase (fG, fsG);
  foreach() {
    double c, e = plic_flux (point, TG, lambda2f, &c);
    dbg_exG[] += e*TG[] - c;
  }
# endif

  /**
  The `drhodt` that the projection receives. The temperature part `dri_cT`
  goes into `drhodt` only after the loop of `INT_TEMP_PICARD`. */

  scalar drfin[];
  foreach()
    drfin[] = drhodt[] + dri_cT[];
  stats sdr = statsf (drfin);
#else
  stats sdr = statsf (drhodt);
#endif
  double qfloor = DRHODT_BUDGET_QFLOOR*max (fabs (sdr.min), fabs (sdr.max));

  double L1d = 0., L1r = 0., maxd = 0., maxr = 0.;
  double nH = 0., L1dH = 0., L1rH = 0., maxdH = 0., maxqH = 0.;
  double nC = 0., L1dC = 0., L1rC = 0., maxdC = 0., maxqC = 0.;
  double nI = 0., L1dI = 0., L1rI = 0., maxdI = 0., maxqI = 0.;

  foreach (reduction(+:L1d) reduction(+:L1r)
           reduction(max:maxd) reduction(max:maxr)
           reduction(+:nH) reduction(+:L1dH) reduction(+:L1rH)
           reduction(max:maxdH) reduction(max:maxqH)
           reduction(+:nC) reduction(+:L1dC) reduction(+:L1rC)
           reduction(max:maxdC) reduction(max:maxqC)
           reduction(+:nI) reduction(+:L1dI) reduction(+:L1rI)
           reduction(max:maxdI) reduction(max:maxqI)) {

#if DRI_ON
    double implS = dri_th1[]*(TS[] - dri_TS[])/dt;
    double implG = dri_th2[]*(TG[] - dri_TG[])/dt;
    double TSpre = dri_TS[], TGpre = dri_TG[], drf = drfin[];
#else
    double implS = dbg_th1[]*(TS[] - dbg_TSpre[])/dt;
    double implG = dbg_th2[]*(TG[] - dbg_TGpre[])/dt;
    double TSpre = dbg_TSpre[], TGpre = dbg_TGpre[], drf = drhodt[];
#endif
    double difS = implS - dbg_exS[], difG = implG - dbg_exG[];

    double eps = f[] > F_ERR ? porosity[]/f[] : 0.;
    double dens = TSpre*(rhoGv_S[]*cpGv_S[]*eps
                               + rhoSv[]*cpSv[]*(1. - eps));
    double cS = (TSpre*rhoGv_S[]*cpGv_S[] > 0. && dens > 0.) ?
      eps/dens : 0.;
    double dG = TGpre*rhoGv_G[]*cpGv_G[];
    double cG = (dG > 0.) ? 1./dG : 0.;

#if DRHODT_CELL_AVERAGE
    double wS = (f[] > F_ERR) ? 1. : 0., wG = (f[] < 1. - F_ERR) ? 1. : 0.;
#else
    double wS = f[], wG = 1. - f[];
#endif

    double dd = -(wS*cS*difS + wG*cG*difG);
    double ad = fabs (dd), ar = fabs (drf);
    double q = (ar > qfloor) ? ad/ar : 0.;

    L1d += ad*dv();
    L1r += ar*dv();
    maxd = max (maxd, ad);
    maxr = max (maxr, ar);

    if (1. - f[] > F_ERR && TGpre > DRHODT_BUDGET_THOT) {
      nH += 1.; L1dH += ad*dv(); L1rH += ar*dv();
      maxdH = max (maxdH, ad); maxqH = max (maxqH, q);
    }
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      nC += 1.; L1dC += ad*dv(); L1rC += ar*dv();
      maxdC = max (maxdC, ad); maxqC = max (maxqC, q);
    }
    if (lay[] >= 2. && lay[] <= 1. + DRHODT_BUDGET_NLAYER) {
      nI += 1.; L1dI += ad*dv(); L1rI += ar*dv();
      maxdI = max (maxdI, ad); maxqI = max (maxqI, q);
    }
  }

  if (pid() == 0) {
    static FILE * fp = NULL;
    if (!fp) {
      fp = fopen (DRHODT_BUDGET_FILE, restarted ? "a" : "w");
      if (fp == NULL) {
        fprintf (stderr, "Error opening %s\n", DRHODT_BUDGET_FILE);
        exit (1);
      }
      if (!restarted)
        fprintf (fp, "#t(1) i(2) dt(3) L1d(4) L1r(5) maxd(6) maxr(7)"
                     " nH(8) L1dH(9) L1rH(10) maxdH(11) maxqH(12)"
                     " nC(13) L1dC(14) L1rC(15) maxdC(16) maxqC(17)"
                     " nI(18) L1dI(19) L1rI(20) maxdI(21) maxqI(22)\n");
    }
    fprintf (fp, "%g %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
                 " %g %g %g %g %g\n",
             t, iter, dt, L1d, L1r, maxd, maxr,
             nH, L1dH, L1rH, maxdH, maxqH,
             nC, L1dC, L1rC, maxdC, maxqC,
             nI, L1dI, L1rI, maxdI, maxqI);
    fflush (fp);
  }
  dbg_cpu += (double)(clock() - c0)/CLOCKS_PER_SEC;
  dbg_ncall++;
}

/**
The cost of the probe, in the log at the end of the run. */

event drhodt_budget_cost (t = end) {
  if (pid() == 0)
    fprintf (stderr, "# drhodt-budget: %g s CPU in %d steps\n",
             dbg_cpu, dbg_ncall);
}
