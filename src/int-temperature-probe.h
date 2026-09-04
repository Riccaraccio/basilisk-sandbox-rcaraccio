/**
# Probe of the interface temperature balance

This header measures the error of the **partitioned** coupling between the
solid temperature `TS` and the gas temperature `TG`. It changes no field of
the solution. Compile it in with `-DINT_TEMP_PROBE=1`.

## What the probe measures

`multicomponent-varprop.h` couples the two temperatures in three steps:

1. `ijc_CoupledTemperature()` finds `TInt` from the flux balance, with the
   values of `TS` and `TG` at the start of the step.
2. The code freezes the two fluxes into the source fields `sST` and `sGT`.
3. `diffusion()` solves `TS` and `TG` separately, each with its frozen source.

The balance of step 1 is true at the start of the step. The probe evaluates
the same balance again after step 3, with the new fields. The probe reports
two independent quantities.

**The residual `res`.** The imbalance of the interface energy budget, in
W/m2:

    res = -q_rad(TInt) + lambda1*grad(TS_new) + lambda2*grad(TG_new)

`ijc_CoupledTemperature()` makes this quantity zero for the old fields. A
large value at the end of the step shows that the lag of `TInt` matters.

**The exchange number `S`.** The ratio of the heat that the frozen source
delivers in one step to the heat capacity of the cell that receives it:

    S = dt*lambda*h*(area/Delta) / theta

`h` is the derivative of the interface gradient with respect to the interface
temperature, and `theta` is the heat capacity of the phase in the cell.
`S > 1` means that one application of the source overshoots the capacity of
the cell. The gas side is the side at risk: `theta2` carries the gas fraction
`1 - f`, and the interface area does not, so `S` grows as `1/(1 - f)`.

## How to read the output

The probe writes `intres.dat`, one line per step. Read it like this:

- `Srel_max` stays under about 0.3: the frozen source is safe. The interface
  coupling is not the origin of the oscillation. Stop.
- `Srel_max` is O(1) or more, and `epsg_at_Smax` is small: the source
  overwhelms thin gas cells. This is the sliver mechanism.
- `rel_max` is large and it oscillates at the frequency of the campaign: the
  lag of `TInt` feeds the mode.

## Exact derivative of the interface gradient

`concentration_gradient_x()` in `intgrad.h` is **affine** in the boundary
value `bc`. Every branch has the form `(bc - v)/(d*Delta)`. So two calls,
with `bc = 0` and `bc = 1`, give the derivative with no truncation error:

    h = grad(1) - grad(0)

Do not replace this with a finite difference. It is exact.
*/

#ifndef INT_TEMP_PROBE
# define INT_TEMP_PROBE 1
#endif

extern scalar fS, fG;
extern face vector fsS, fsG;
extern vector lambda1v, lambda2v;
extern scalar TInt, TS, TG;
extern scalar * YSList;
extern bool restarted;

/**
Diagnostic fields. They are live scalars, so a movie or a dump can show
where the trouble sits. */

scalar resint[];      // interface energy residual, W/m2
scalar relint[];      // the same, divided by the largest of the three fluxes
scalar Sint_G[];      // exchange number of the gas side
scalar Sint_S[];      // exchange number of the solid side

/**
Reduced quantities of the last step. `int_temperature_probe()` fills them
and `event intres_output` writes them. */

double ITP_res_max, ITP_res_l2, ITP_rel_max, ITP_rel_l2;
double ITP_SG_max, ITP_SS_max, ITP_epsg_at_Smax, ITP_epsg_at_resmax;
double ITP_epsg_min, ITP_nint, ITP_nS1, ITP_nrel10, ITP_TG_min;
double ITP_QS, ITP_QG, ITP_Qrad;

void int_temperature_probe (void)
{
  int ash_index = OpenSMOKE_IndexOfSolidSpeciesWithoutError ("ASH");
  scalar YASH = {-1};
  if (ash_index >= 0)
    YASH = YSList[ash_index];

  foreach() {
    resint[] = 0.; relint[] = 0.;
    Sint_G[] = 0.; Sint_S[] = 0.;
  }

  double res_max = 0., res_l2 = 0., rel_max = 0., rel_l2 = 0.;
  double SG_max = 0., SS_max = 0., epsg_min = 1.;
  double nint = 0., nS1 = 0., nrel10 = 0.;
  double QS = 0., QG = 0., Qrad = 0.;

  foreach (reduction(max:res_max) reduction(+:res_l2)
           reduction(max:rel_max) reduction(+:rel_l2)
           reduction(max:SG_max)  reduction(max:SS_max)
           reduction(min:epsg_min)
           reduction(+:nint) reduction(+:nS1) reduction(+:nrel10)
           reduction(+:QS) reduction(+:QG) reduction(+:Qrad)) {

    if (f[] > F_ERR && f[] < 1. - F_ERR && TS[] > 0. && TG[] > 0.) {

    coord n = facet_normal (point, fS, fsS), pc;
    double alpha = plane_alpha (fS[], n);
    double area = plane_area_center (n, alpha, &pc);
    normalize (&n);

    /**
    The same normal-weighted blend of the conductivity that `EqTemperature`
    and the source block use. Keep the three sites identical. */

    n.x = fabs(n.x); n.y = fabs(n.y);
    double wsum = n.x + n.y;
    double lambda1vh = (wsum > 0.) ?
      n.x/wsum*lambda1v.x[] + n.y/wsum*lambda1v.y[] : 0.;
    double lambda2vh = (wsum > 0.) ?
      n.x/wsum*lambda2v.x[] + n.y/wsum*lambda2v.y[] : 0.;

    bool success = false;
    double gS = ebmgrad (point, TS, fS, fG, fsS, fsG, false, TInt[], &success);
    double gG = ebmgrad (point, TG, fS, fG, fsS, fsG, true,  TInt[], &success);

    /**
    The exact derivative. `ebmgrad` is affine in the interface value. */

    double hS = ebmgrad (point, TS, fS, fG, fsS, fsG, false, 1., &success)
              - ebmgrad (point, TS, fS, fG, fsS, fsG, false, 0., &success);
    double hG = ebmgrad (point, TG, fS, fG, fsS, fsG, true,  1., &success)
              - ebmgrad (point, TG, fS, fG, fsS, fsG, true,  0., &success);

    double char_fraction = calculate_char_fraction (point, YSList, f);
    double ash_fraction = (ash_index >= 0) ? YASH[]/f[] : 0.;
    double emis = emissivity (char_fraction, ash_fraction);
    double qrad = divq_rad_int (TInt[], RADIATION_TEMP, emis);

    double qS = lambda1vh*gS, qG = lambda2vh*gG;
    double res = -qrad + qS + qG;

    double qscale = max (max (fabs(qS), fabs(qG)), max (fabs(qrad), 1e-30));
    double rel = fabs(res)/qscale;

    /**
    The area of the interface per unit volume, with the metric of the source
    block of `multicomponent-varprop.h`. Keep the two identical. */

#ifdef AXI
    double aov = area*(y + pc.y*Delta)/(Delta*y)*cm[];
#else
    double aov = area/Delta*cm[];
#endif

    /**
    The heat capacity of each phase in this cell, as `theta1` and `theta2` of
    the diffusion call. */

    double theta1vh, theta2vh;
#ifdef VARPROP
    theta1vh = fS[] > F_ERR ?
      porosity[]/fS[]*rhoGv_S[]*cpGv_S[] + (1. - porosity[]/fS[])*rhoSv[]*cpSv[] : 0.;
    theta2vh = rhoGv_G[]*cpGv_G[];
#else
    theta1vh = fS[] > F_ERR ?
      porosity[]/fS[]*rhoG*cpG + (1. - porosity[]/fS[])*rhoS*cpS : 0.;
    theta2vh = rhoG*cpG;
#endif
    double th1 = cm[]*fS[]*theta1vh;
    double th2 = cm[]*fG[]*theta2vh;

    double SS = (th1 > 0.) ? dt*fabs(lambda1vh*hS)*aov/th1 : nodata;
    double SG = (th2 > 0.) ? dt*fabs(lambda2vh*hG)*aov/th2 : nodata;

    resint[] = res; relint[] = rel;
    Sint_S[] = SS; Sint_G[] = SG;

    res_max = max (res_max, fabs(res));
    rel_max = max (rel_max, rel);
    res_l2 += sq(res); rel_l2 += sq(rel);
    if (SS != nodata) SS_max = max (SS_max, SS);
    if (SG != nodata) SG_max = max (SG_max, SG);
    epsg_min = min (epsg_min, fG[]);
    nint  += 1.;
    nS1   += (SG != nodata && SG > 1.) ? 1. : 0.;
    nrel10 += (rel > 0.1) ? 1. : 0.;

    /**
    The total power that crosses the interface, per radian in axi. It gives
    an independent check: `QS + QG - Qrad` must be small. */

    QS   += qS*aov*dv();
    QG   += qG*aov*dv();
    Qrad += qrad*aov*dv();
    }
  }

  /**
  The gas fraction of the two worst cells. A second pass is the simplest way
  to get an argument of the maximum under MPI. */

  double epsg_at_Smax = nodata, epsg_at_resmax = nodata;
  if (nint > 0.) {
    epsg_at_Smax = -1.; epsg_at_resmax = -1.;
    foreach (reduction(max:epsg_at_Smax) reduction(max:epsg_at_resmax)) {
      if (f[] > F_ERR && f[] < 1. - F_ERR) {
        if (SG_max > 0. && Sint_G[] >= SG_max*(1. - 1e-12))
          epsg_at_Smax = max (epsg_at_Smax, fG[]);
        if (res_max > 0. && fabs(resint[]) >= res_max*(1. - 1e-12))
          epsg_at_resmax = max (epsg_at_resmax, fG[]);
      }
    }
  }

  /**
  The smallest gas temperature of the whole domain. It turns negative several
  steps before the run stops with a floating point exception, so it is the
  early warning of the runaway. */

  double TG_min = nodata;
  foreach (reduction(min:TG_min))
    if (fG[] > F_ERR)
      TG_min = min (TG_min, TG[]);

  ITP_res_max = res_max;
  ITP_rel_max = rel_max;
  ITP_res_l2  = (nint > 0.) ? sqrt (res_l2/nint) : 0.;
  ITP_rel_l2  = (nint > 0.) ? sqrt (rel_l2/nint) : 0.;
  ITP_SG_max  = SG_max;
  ITP_SS_max  = SS_max;
  ITP_epsg_at_Smax   = epsg_at_Smax;
  ITP_epsg_at_resmax = epsg_at_resmax;
  ITP_epsg_min = epsg_min;
  ITP_nint = nint; ITP_nS1 = nS1; ITP_nrel10 = nrel10;
  ITP_TG_min = TG_min;
  ITP_QS = QS; ITP_QG = QG; ITP_Qrad = Qrad;
}

/**
The probe runs inside the `tracer_diffusion` event of
`multicomponent-varprop.h`, where `TS` and `TG` hold the value of one phase.
This event only writes the result, so it can run at the end of the step.

`INT_TEMP_PROBE_EVERY` subsamples the output. Keep it at 1: the fast band of
[[tmax-is-a-function-of-dt]] is a single sample, and a subsample hides it. */

#ifndef INT_TEMP_PROBE_EVERY
# define INT_TEMP_PROBE_EVERY 1
#endif

event intres_output (i++, last) {
  if (i % INT_TEMP_PROBE_EVERY)
    return 0;

  if (pid() == 0) {
    static FILE * fr = fopen ("intres.dat", restarted ? "a" : "w");
    if (fr == NULL) {
      fprintf (stderr, "Error opening intres.dat\n");
      exit (1);
    }
    if (i == 0)
      fprintf (fr, "#t(1) dt(2) nint(3) res_max(4) res_l2(5) rel_max(6)"
                   " rel_l2(7) SG_max(8) SS_max(9) epsg_at_Smax(10)"
                   " epsg_at_resmax(11) epsg_min(12) nS1(13) nrel10(14)"
                   " TG_min(15) QS(16) QG(17) Qrad(18)\n");

    fprintf (fr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
             t, dt, ITP_nint, ITP_res_max, ITP_res_l2, ITP_rel_max,
             ITP_rel_l2, ITP_SG_max, ITP_SS_max, ITP_epsg_at_Smax,
             ITP_epsg_at_resmax, ITP_epsg_min, ITP_nS1, ITP_nrel10,
             ITP_TG_min, ITP_QS, ITP_QG, ITP_Qrad);
    fflush (fr);
  }
}
