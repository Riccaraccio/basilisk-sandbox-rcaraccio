#ifndef MULTICOMPONENT
  #define MULTICOMPONENT 1
#endif

#include "intgrad.h"

#ifdef EXPLICIT_DIFFUSION
  #include "diffusion-explicit.h"
#else
  #include "diffusion.h"
#endif

#include "common-phasechange.h"
#include "memoryallocation-varprop.h"
#include "int-temperature.h"
#include "int-concentration.h"
#include "multicomponent-properties.h"
#include "chemistry.h"

/**
The probe of the interface temperature balance. It changes no field. It is
compiled out unless the case sets `INT_TEMP_PROBE`. */

#if INT_TEMP_PROBE
# include "int-temperature-probe.h"
#endif

/**
The record of the Picard loop on the interface temperature. It changes no
field either. The loop itself is in the `tracer_diffusion` event below. */

#if INT_TEMP_PICARD
# include "int-temperature-picard.h"
#endif

/**
The record of the interface conductance, and the reason it must give back the
heat that it holds. It changes no field either. */

#if INT_TEMP_ROBIN
# include "int-temperature-robin.h"
#endif

/**
The tolerance of the two temperature solves, scaled to their own residual.
The inherited `TOLERANCE = 1e-5` asks the gas temperature for 7e-12 K, which
no solve can deliver, and the wasted iterations killed two runs. On by
default; set `INT_TEMP_TOL` to 0 to restore the inherited value. */

#ifndef INT_TEMP_TOL
# define INT_TEMP_TOL 1
#endif

#if INT_TEMP_TOL
# include "int-temperature-tol.h"
#endif

/**
## The gradient of the mass diffusion enthalpy term

`TS`, `TG`, `YGList_S`, `YGList_G` and the two mole fraction lists hold the
value of one phase each. Outside that phase they are 0. A centred stencil that
reaches into a cell without the phase therefore reads 0, not a temperature or
a mass fraction, and it gives a false gradient.

This function uses a neighbour only where the phase exists. With both
neighbours it gives the centred difference, which is what the previous version
gave in a full cell. With one neighbour it gives the one sided difference.
With no neighbour it gives 0.

Give `fS` for the solid side and `fG` for the gas side. */

/**
## The off switches of the three transport fixes

Each fix has a switch, so that a case can measure the fix against the previous
code. All three are on by default, and each one is independent of the others.

  `CORRECTIVE_CFL`      The Courant limit of the corrective flux. 0 removes
                        the limit and gives back the previous timestep.
                        Default 0.5.
  `CORRECTIVE_LIMITER`  The limited slope and the true Courant number of the
                        corrective flux. 0 gives back the unlimited centred
                        slope, and the mass flux in place of a velocity.
                        Default 1.
  `MDE_INTERFACE`       The mass diffusion enthalpy term in the interface
                        cells, with the phase weight and the phase aware
                        gradient. 0 gives back the term in the full cells
                        only, with the plain centred stencil. Default 1.

Set all three off to get the previous code exactly. */

#ifndef CORRECTIVE_LIMITER
# define CORRECTIVE_LIMITER 1
#endif

#ifndef MDE_INTERFACE
# define MDE_INTERFACE 1
#endif

/**
The term feeds the two temperature solves, and the cut cell branch reads
`TInt`. Without `SOLVE_TEMPERATURE` neither of them exists, so keep the
previous gate in that build. */

#if MDE_INTERFACE && !defined SOLVE_TEMPERATURE
# undef MDE_INTERFACE
# define MDE_INTERFACE 0
#endif

#ifdef MASS_DIFFUSION_ENTHALPY
foreach_dimension()
static double mde_gradient_x (Point point, scalar a, scalar ff)
{
#if MDE_INTERFACE
  bool vp = (ff[1] > F_ERR), vm = (ff[-1] > F_ERR);
  if (vp && vm) return (a[1] - a[-1])/(2.*Delta);
  if (vp)       return (a[1] - a[])/Delta;
  if (vm)       return (a[] - a[-1])/Delta;
  return 0.;
#else
  (void) ff;                          // the previous plain centred stencil
  return (a[1] - a[-1])/(2.*Delta);
#endif
}
#endif

/**
## The timestep limit of the corrective flux

`FICK_CORRECTED` moves each species with an extra velocity `u_c = phic/rho`,
and `MOLAR_DIFFUSION` adds `-D grad(MW_mix)/MW_mix` to it. That velocity is
not always small. With a heavy product in a light carrier, `MW_mix` changes by
a factor near 2 across the reaction front. Then `D d(ln MW_mix)/dx` reaches
0.2 m/s, which is more than the inflow velocity of the cases in `run/`.

The transport of that flux is explicit, and no event limited the timestep by
it. The scheme therefore ran at a Courant number above 1 at the front.

The `tracer_diffusion` event records `max(|u_c|/Delta)` in `corrective_uodx`.
This event turns that record into a limit on `dtmax`. The record is one step
old. That is safe, because the corrective velocity changes slowly.

Set `CORRECTIVE_CFL` to 0 to remove the limit and get the previous timestep. */

#ifdef FICK_CORRECTED
# ifndef CORRECTIVE_CFL
#  define CORRECTIVE_CFL 0.5
# endif

double corrective_uodx = 0.;    // max |u_c|/Delta of the last step
double corrective_dtmax = HUGE; // the limit that it gives

event stability (i++) {
  if (CORRECTIVE_CFL > 0. && corrective_uodx > 0.) {
    corrective_dtmax = CORRECTIVE_CFL/corrective_uodx;
    if (corrective_dtmax < dtmax)
      dtmax = corrective_dtmax;
  }
}
#endif

/**
## The probe of the gas temperature

`TG_PROBE` reports every gas cell whose temperature falls below
`TG_PROBE_TMIN`. It writes `tgprobe-<pid>.dat`.

The first question the file answers is WHERE the fall happens. `TGadv` is the
value that enters the diffusion solve, thus it carries the advection of the
step. `TGpost` is the value that leaves it. If `Tpre_phys` is already bad, the
advection makes it. If only `Tpost_phys` is bad, the source terms make it.

The next columns separate the sources. `qint` is the interface heat flux,
`qmde` is the mass diffusion enthalpy and `qrob` is the Robin correction.
`theta2` is the heat capacity of the cell, and it carries the factor `fG`,
which the interface flux does NOT carry. A small `fG` with a large `qint` is
the runaway. */

#ifdef TG_PROBE
scalar qint_dbg[], qmde_dbg[], qrob_dbg[], TGadv_dbg[];
# ifndef TG_PROBE_TMIN
#  define TG_PROBE_TMIN 100.
# endif
# ifndef TG_PROBE_MAX
#  define TG_PROBE_MAX 400
# endif
#endif

event reset_sources (i++) {
#ifdef SOLVE_TEMPERATURE
  foreach() {
    sST[] = 0.;
    sGT[] = 0.;
#ifdef TG_PROBE
    qint_dbg[] = 0.;
    qmde_dbg[] = 0.;
    qrob_dbg[] = 0.;
#endif
#if INT_TEMP_ROBIN
    betaST[] = 0.;
    betaGT[] = 0.;
#endif
  }
#endif

  reset (sGexpList, 0.);
  reset (sSexpList, 0.);
}


void update_mole_fields() {
  #ifdef MOLAR_DIFFUSION
  foreach() {
    double xG[NGS], yG[NGS];
    double MWmix;
    if (f[] > F_ERR) { // Internal gas phase
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_S[jj];
        yG[jj] = YG[];
      }
      mole_from_mass (xG, &MWmix, yG, NGS);
      // MWmixG_S[] = MWmix; // Already done in update_properties()
      for (int jj = 0; jj < NGS; jj++) {
        scalar XG = XGList_S[jj];
        XG[] = xG[jj];
      }
    }

    if (f[] < 1. - F_ERR) { // External gas phase
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        yG[jj] = YG[];
      }
      mole_from_mass (xG, &MWmix, yG, NGS);
      // MWmixG_G[] = MWmix; // Already done in update_properties()
      for (int jj = 0; jj < NGS; jj++) {
        scalar XG = XGList_G[jj];
        XG[] = xG[jj];
      }
    }
  }
  boundary (XGList_S);
  boundary (XGList_G); // Ensure boundary conditions are applied
  #endif
}

#ifdef SOLVE_TEMPERATURE
/**
## The interface heat source

This assembles the two interface heat sources `sST` and `sGT`, and, under
`INT_TEMP_ROBIN`, the two conductances `betaST` and `betaGT`. It was part of
the species source loop. It is a function of its own so that the Picard loop
of `INT_TEMP_PICARD` can call it again with an updated `TInt`.

The function adds to the four fields. The caller must set them to a known
value first. `event reset_sources` does that once per step; the Picard loop
does it again on each pass.

Caution: `TS` and `TG` must hold the value of one phase here, not the tracer
form. The event divides them at its start and multiplies them back at its
end. */

static void interface_temperature_sources (void)
{
  bool success = false;

  foreach() {
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      coord n = facet_normal (point, fS, fsS), p;
      double alpha = plane_alpha (fS[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      double bc = TInt[];
      double Strgrad = ebmgrad (point, TS, fS, fG, fsS, fsG, false, bc, &success);
      double Gtrgrad = ebmgrad (point, TG, fS, fG, fsS, fsG, true , bc, &success);

      n.x = fabs(n.x); n.y = fabs(n.y);

      double lambda1vh = n.x/(n.x+n.y)*lambda1v.x[] + n.y/(n.x+n.y)*lambda1v.y[];
      double lambda2vh = n.x/(n.x+n.y)*lambda2v.x[] + n.y/(n.x+n.y)*lambda2v.y[];

      double Sheatflux = lambda1vh*Strgrad;
      double Gheatflux = lambda2vh*Gtrgrad;

# ifdef AXI
      double aov = area*(y + p.y*Delta)/(Delta*y)*cm[];
# else
      double aov = area/Delta*cm[];
# endif

      sST[] += Sheatflux*aov;
      sGT[] += Gheatflux*aov;
#ifdef TG_PROBE
      qint_dbg[] = Gheatflux*aov;
#endif

/**
## The interface conductance on the diagonal

`ebmgrad` builds the gradient from the NEIGHBOURS of this cell along the
normal. The cell value `TG[]` is not in it. So the source that heats a cut
cell does not answer to the temperature of that cell: the only restoring
term is the internal diffusion, and the face fractions of a sliver make it
weak. The heat that one step delivers is then large against the heat
capacity `theta2 = cm*fG*rhoG*cpG`, which goes to zero with the gas
fraction. The probe of `int-temperature-probe.h` measures the ratio

    S = dt*lambda*h*aov/theta

and a measured run stopped at S of order 100: the gas temperature of one cut
cell went through zero and grew by a factor 2.8 per step for five steps.

Add a conductance `K` to the diagonal and give the same `K*T^n` back to the
source. The step then reads

    (theta/dt + K)(T^{n+1} - T^n) = div(D grad T^{n+1}) + src

so the change per step falls by `1/(1 + K*dt/theta)`, and the effective
exchange number becomes `dt*A/(theta + K*dt)`. This is the deferred
correction of a Dirichlet cut-cell condition: it keeps the accurate gradient
stencil and damps the path to it.

Caution: the two added terms cancel only when the field stops changing. On
its own this scheme therefore DESTROYS heat while the field moves. When `K`
fires the diagonal is exactly `A/SMAX`, so the cell keeps `SMAX/S` of the
heat that the interface gave it — one per cent at `SMAX = 1` and `S = 100`.
A measured run lost 20 to 25 K of surface temperature and 25 per cent of the
mass loss rate. `INT_TEMP_ROBIN_DEBT`, which is on by default, carries the
withheld heat to the next step and removes that loss. Read
`int-temperature-robin.h` before you change any of this.

Choose `K` to bring the effective exchange number down to
`INT_TEMP_ROBIN_SMAX` and no further:

    K = max (0, A/SMAX - theta/dt)

A cell already under the limit gets `K = 0` and is untouched, bit for bit.
Only the cells that the probe reports as unstable change.

**The value of `SMAX` decides whether that last sentence is true.** The first
version used `SMAX = 1`, which fired on 63 interface cells of about 90 and
rewrote the answer: `Tsurf` fell 2.1 K by t = 0.3 s and 25 K by t = 8 s, the
radial velocity at 1.5 mm changed sign, and the `full` ladder died at
t = 8.48. The default is now 20. It fires on 2 cells, it clamps only the
outliers that the probe reports, and it reproduces the run without the flag
to every printed digit over the first 0.3 s.

Caution: 20 is five times below the `S` of 100 that crashed a measured run,
but the crash case must be repeated at this value. Run `test-robinm` from the
t = 5 s dump and check that it passes t = 5.94 s before you trust the number.

`ebmgrad` is affine in the interface value, so `A` is exact. Do not
finite-difference it. */

# if INT_TEMP_ROBIN
#  ifndef INT_TEMP_ROBIN_SMAX
#   define INT_TEMP_ROBIN_SMAX 20.
#  endif
      double hS = ebmgrad (point, TS, fS, fG, fsS, fsG, false, 1., &success)
                - ebmgrad (point, TS, fS, fG, fsS, fsG, false, 0., &success);
      double hG = ebmgrad (point, TG, fS, fG, fsS, fsG, true,  1., &success)
                - ebmgrad (point, TG, fS, fG, fsS, fsG, true,  0., &success);

      /**
      The heat capacity of each phase, exactly as the `diffusion()` call
      builds `theta1` and `theta2`. Keep the two sites identical. The
      `VARCOEFF` block later divides `betaST` and `sST` by the same heat
      capacity that it divides `theta1` by, so `A/theta` is unchanged. */

      double theta1vh, theta2vh;
#  ifdef VARPROP
      theta1vh = fS[] > F_ERR ?
        porosity[]/fS[]*rhoGv_S[]*cpGv_S[] + (1. - porosity[]/fS[])*rhoSv[]*cpSv[] : 0.;
      theta2vh = rhoGv_G[]*cpGv_G[];
#  else
      theta1vh = fS[] > F_ERR ?
        porosity[]/fS[]*rhoG*cpG + (1. - porosity[]/fS[])*rhoS*cpS : 0.;
      theta2vh = rhoG*cpG;
#  endif
      double th1 = cm[]*max(fS[]*theta1vh, F_ERR);
      double th2 = cm[]*max(fG[]*theta2vh, F_ERR);

      double smax = INT_TEMP_ROBIN_SMAX;
      double KS = max (0., fabs(lambda1vh*hS)*aov/smax - th1/dt);
      double KG = max (0., fabs(lambda2vh*hG)*aov/smax - th2/dt);

      /**
      Keep the conductance. The debt update after the solve needs it to
      measure the heat that this step withheld. */

      KSf[] = KS;
      KGf[] = KG;

      betaST[] -= KS;
      betaGT[] -= KG;

      /**
      Give back the heat that the last step withheld. Without this term the
      conductance is a heat sink: the cell keeps only `SMAX/S` of what the
      interface gave it, and the rest is destroyed. With it, the steady
      heating rate is exact, because `(theta/dt + K) dT = q + K dT` reduces
      to `(theta/dt) dT = q`. See `int-temperature-robin.h`.

      `debtST` and `debtGT` are constant over a step, so a `INT_TEMP_PICARD`
      pass that rebuilds the source adds them again, which is correct. */

# ifndef INT_TEMP_ROBIN_DEBT
#  define INT_TEMP_ROBIN_DEBT 1
# endif
# if INT_TEMP_ROBIN_DEBT
      sST[] += KS*TS[] + debtST[];
      sGT[] += KG*TG[] + debtGT[];
#  ifdef TG_PROBE
      qrob_dbg[] = KG*TG[] + debtGT[];
#  endif
# else
      sST[] += KS*TS[];
      sGT[] += KG*TG[];
#  ifdef TG_PROBE
      qrob_dbg[] = KG*TG[];
#  endif
# endif
# endif
    }
  }
}

#endif

event tracer_diffusion (i++) {

  //Check the mass fractions Can be removed for performance
  check_and_correct_fractions (YGList_S, NGS, false);
  check_and_correct_fractions (YGList_G, NGS, true);
  check_and_correct_fractions (YSList,   NSS, false);

  foreach() {
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]/f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]/(1. - f[]) : 0.;
#endif

    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList_S[jj];
      YG[] = (f[] > F_ERR) ? YG[]/f[] : 0.;
    }
    
    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList_G[jj];
      YG[] = ((1. - f[]) > F_ERR) ? YG[]/(1. - f[]) : 0.;
    }
  }

#ifdef MOLAR_DIFFUSION
  update_mole_fields();
#endif


#ifdef SOLVE_TEMPERATURE
  //interface temperature first guess
  foreach() {
    TInt[] = 0.;
    if (f[] > F_ERR && f[] < 1. - F_ERR)
      TInt[] = (TS[] + TG[])/2;
  }

  #ifdef FIXED_INT_TEMP //Force interface temperature = TG0
  foreach()
    if (f[] > F_ERR && f[] < 1. - F_ERR)
      TInt[] = TG0;

  #elif TEMPERATURE_PROFILE
  double tv = TemperatureProfile_GetT(t);
  foreach() {
    if (f[] > F_ERR && f[] < 1. - F_ERR)
      TInt[] = tv;
    
    if (f[] < 1. - F_ERR)
      TG[] = tv;
  }

  #else //default: solve for interface temperature
  ijc_CoupledTemperature();
  #endif
#endif

  // first guess for species interface concentration
  foreach() {
    for (int jj=0; jj<NGS; jj++) {
      scalar YGInt = YGList_Int[jj];
      YGInt[] = 0.;
      if (f[] > F_ERR && f[] < 1. - F_ERR) {
        scalar YG_S = YGList_S[jj];
        scalar YG_G = YGList_G[jj];
        YGInt[] = (YG_G[] + YG_S[])/2;
        YGInt[] = clamp (YGInt[], 0., 1.);
      }
    }
  }

  //find the interface concentration for each species
  intConcentration();

#ifdef MOLAR_DIFFUSION // calculate the mole fractions at the interface
  foreach()
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      double xG[NGS], yG[NGS], MWmixInt;
      for (int jj=0; jj<NGS; jj++) {
        scalar YGInt = YGList_Int[jj];
        yG[jj] = YGInt[];
      }
      mole_from_mass (xG, &MWmixInt, yG, NGS);
      for (int jj=0; jj<NGS; jj++) {
        scalar XGInt = XGList_Int[jj];
        XGInt[] = xG[jj];
      }
    }
#endif

  //Calculate the source therm
  foreach() {
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      coord n = facet_normal (point, fS, fsS), p;
      double alpha = plane_alpha (fS[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      //Solid side
      double jS[NGS];
      for (int jj=0; jj<NGS; jj++) {
        scalar DmixG = DmixGList_S[jj];

        double rhoGvh_S;
        #ifdef VARPROP
        rhoGvh_S = rhoGv_S[];
        #else
        rhoGvh_S = rhoG;
        #endif
        
        #ifdef MOLAR_DIFFUSION
        scalar XG = XGList_S[jj];
        scalar XGInt = XGList_Int[jj];
        double Strgrad = ebmgrad (point, XG, fS, fG, fsS, fsG, false, XGInt[], &success);
        jS[jj] = (MWmixG_S[] > 0.) ?
          rhoGvh_S*DmixG[]*Strgrad*gas_MWs[jj]/MWmixG_S[] : 0.; // MW==0 in unfilled cells -> 0/0
        #else
        scalar YG    = YGList_S[jj];
        scalar YGInt = YGList_Int[jj];
        double Strgrad = ebmgrad (point, YG, fS, fG, fsS, fsG, false, YGInt[], &success);
        jS[jj] = rhoGvh_S*DmixG[]*Strgrad; 
        #endif
      }

      double jStot = 0.;
#ifdef FICK_CORRECTED
      for (int jj=0; jj<NGS; jj++)
        jStot += jS[jj];
#endif

      for (int jj=0; jj<NGS; jj++) {
        scalar sSexp = sSexpList[jj];
        scalar YGInt = YGList_Int[jj];
        jS[jj] -= jStot*YGInt[];
#ifdef AXI
        sSexp[] += jS[jj]*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
        sSexp[] += jS[jj]*area/Delta*cm[];
#endif
      }

      //Gas side
      double jG[NGS];
      for (int jj=0; jj<NGS; jj++) {
        scalar DmixG = DmixGList_G[jj];

        double rhoGvh_G;
#ifdef VARPROP
        rhoGvh_G = rhoGv_G[];
#else
        rhoGvh_G = rhoG;
#endif

#ifdef MOLAR_DIFFUSION
        scalar XG = XGList_G[jj];
        scalar XGInt = XGList_Int[jj];
        double Gtrgrad = ebmgrad (point, XG, fS, fG, fsS, fsG, true, XGInt[], &success);
        jG[jj] = (MWmixG_G[] > 0.) ?
          rhoGvh_G*DmixG[]*Gtrgrad*gas_MWs[jj]/MWmixG_G[] : 0.; // MW==0 in unfilled cells -> 0/0
#else
        scalar YG    = YGList_G[jj];
        scalar YGInt = YGList_Int[jj];
        double Gtrgrad = ebmgrad (point, YG, fS, fG, fsS, fsG, true, YGInt[], &success);
        jG[jj] = rhoGvh_G*DmixG[]*Gtrgrad;
#endif
      }

      double jGtot = 0.;
#ifdef FICK_CORRECTED
      for (int jj=0; jj<NGS; jj++)
        jGtot += jG[jj];
#endif

      for (int jj=0; jj<NGS; jj++) {
        scalar sGexp = sGexpList[jj];
        scalar YGInt = YGList_Int[jj];
        jG[jj] -= jGtot*YGInt[];
#ifdef AXI
        sGexp[] += jG[jj]*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
        sGexp[] += jG[jj]*area/Delta*cm[];
#endif
      }

    }
  }

#ifdef MASS_DIFFUSION_ENTHALPY

  /**
  ## The mass diffusion enthalpy source

  The term is `- sum_j cp_j (J_j - Y_j sum_k J_k) . grad T`. It is an explicit
  source of the two temperature solves. It is 0 when every `cp_j` is equal.

  Two properties of the previous version made it noisy at the front.

  1. The code applied it in full cells only (`f > 1 - F_ERR` on the solid side
     and `f < F_ERR` on the gas side). An interface cell got nothing. The term
     is largest in the cells next to the front, so each cell switched a large
     explicit source on, then off, then on again as the front passed it. The
     source now carries `fS[]` on the solid side and `fG[]` on the gas side.
     That is the same weight that `theta1` and `theta2` carry, so it is the
     consistent volume weight. A full cell keeps the previous value. An
     interface cell gets a fraction of the term instead of a step to 0.

  2. The code used the plain centred stencil `(a[1] - a[-1])/(2*Delta)`. But
     `TS`, `TG` and the two species lists hold the value of one phase only,
     and they are 0 outside that phase. A stencil that reaches into a cell
     without the phase reads 0 and gives a false gradient.
     `mde_gradient_x()` uses a neighbour only where the phase exists.

  A centred stencil is not correct in an interface cell either, and the phase
  test alone does not repair it. In an interface cell the code holds a phase
  AVERAGE, and that average belongs to the centroid of the phase, not to the
  centre of the cell. A difference between an interface cell and a full cell
  is therefore a difference of two values at two unknown positions. The error
  is of the order of the offset, which is of the order of `Delta`, so the
  gradient loses its order in the cells that this fix adds.

  An interface cell therefore uses `ebmgrad()` and the interface value, in the
  direction of the normal. That is the same method that the species source of
  this event uses, and that `int-temperature.h` uses for `TS` and `TG`. It
  keeps the geometry of the cut cell, and it needs no value from a cell that
  holds another phase. The term becomes the normal part of the product, which
  is the part that survives at a front. `ebmgrad()` gives both gradients along
  the same normal, so the sign of the product does not depend on the
  orientation. */

  foreach() {

    /**
    The weight of the two sides. `MDE_INTERFACE` on gives the phase fraction,
    which is the weight that `theta1` and `theta2` carry. `MDE_INTERFACE` off
    gives the previous gate: 1 in a full cell, and 0 in every other cell. */

#if MDE_INTERFACE
    double wS = fS[], wG = fG[];
    bool interfacial = (f[] > F_ERR && f[] < 1. - F_ERR);
#else
    double wS = (fS[] > 1. - F_ERR) ? 1. : 0.;
    double wG = (fG[] > 1. - F_ERR) ? 1. : 0.;
#endif

    if (wS > F_ERR) { //Internal gas phase
      double mdeGS = 0.;

#if MDE_INTERFACE
      if (interfacial) {

        /**
        The cut cell. `ebmgrad()` reads the gradient of the phase toward the
        interface, from the interface value. It never reads a cell of the
        other phase, and it does not assume that the stored value belongs to
        the centre of the cell. */

        bool ok = false;
        double gTn = ebmgrad (point, TS, fS, fG, fsS, fsG, false, TInt[], &ok);
        double jn[NGS], jntot = 0.;

        for (int jj=0; jj<NGS; jj++) {
          scalar Dmixv = DmixGList_S[jj];
# ifdef MOLAR_DIFFUSION
          scalar XG = XGList_S[jj];
          scalar XGInt = XGList_Int[jj];
          jn[jj] = (MWmixG_S[] > 0.) ?
            -rhoGv_S[]*Dmixv[]*gas_MWs[jj]/MWmixG_S[]*
             ebmgrad (point, XG, fS, fG, fsS, fsG, false, XGInt[], &ok) : 0.;
# else
          scalar YG = YGList_S[jj];
          scalar YGInt = YGList_Int[jj];
          jn[jj] = -rhoGv_S[]*Dmixv[]*
             ebmgrad (point, YG, fS, fG, fsS, fsG, false, YGInt[], &ok);
# endif
          jntot += jn[jj];
        }

        for (int jj=0; jj<NGS; jj++) {
          scalar YG = YGList_S[jj];
          scalar cpGv = cpGList_S[jj];
          mdeGS += cpGv[]*(jn[jj] - YG[]*jntot)*gTn;
        }
      }
      else
#endif
      {
      coord gTS = {0., 0., 0.};
      coord gYGj_S = {0., 0., 0.};
      coord gYGsum_S = {0., 0., 0.};

      foreach_dimension()
        gTS.x = mde_gradient_x (point, TS, fS);

      foreach_dimension() {
        for (int jj=0; jj<NGS; jj++) {
          scalar Dmixv = DmixGList_S[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList_S[jj];
          gYGsum_S.x -= (MWmixG_S[] > 0.) ?
            rhoGv_S[]*Dmixv[]*gas_MWs[jj]/MWmixG_S[]*mde_gradient_x (point, XG, fS) : 0.;
  # else
          scalar YG = YGList_S[jj];
          gYGsum_S.x -= rhoGv_S[]*Dmixv[]*mde_gradient_x (point, YG, fS);
  # endif
        }

        for (int jj=0; jj<NGS; jj++) {
          scalar YG = YGList_S[jj];
          scalar cpGv = cpGList_S[jj];
          scalar Dmixv = DmixGList_S[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList_S[jj];
          gYGj_S.x = (MWmixG_S[] > 0.) ?
            -rhoGv_S[]*Dmixv[]*gas_MWs[jj]/MWmixG_S[]*mde_gradient_x (point, XG, fS) : 0.;
  # else
          gYGj_S.x = -rhoGv_S[]*Dmixv[]*mde_gradient_x (point, YG, fS);
  # endif
          mdeGS += cpGv[]*(gYGj_S.x - YG[]*gYGsum_S.x)*gTS.x;
        }
      }
      }
      sST[] -= mdeGS*cm[]*wS;
    }

    if (wG > F_ERR) { //External gas phase
      double mdeGG = 0.;

#if MDE_INTERFACE
      if (interfacial) {

        /**
        The cut cell. `ebmgrad()` reads the gradient of the phase toward the
        interface, from the interface value. It never reads a cell of the
        other phase, and it does not assume that the stored value belongs to
        the centre of the cell. */

        bool ok = false;
        double gTn = ebmgrad (point, TG, fS, fG, fsS, fsG, true, TInt[], &ok);
        double jn[NGS], jntot = 0.;

        for (int jj=0; jj<NGS; jj++) {
          scalar Dmixv = DmixGList_G[jj];
# ifdef MOLAR_DIFFUSION
          scalar XG = XGList_G[jj];
          scalar XGInt = XGList_Int[jj];
          jn[jj] = (MWmixG_G[] > 0.) ?
            -rhoGv_G[]*Dmixv[]*gas_MWs[jj]/MWmixG_G[]*
             ebmgrad (point, XG, fS, fG, fsS, fsG, true, XGInt[], &ok) : 0.;
# else
          scalar YG = YGList_G[jj];
          scalar YGInt = YGList_Int[jj];
          jn[jj] = -rhoGv_G[]*Dmixv[]*
             ebmgrad (point, YG, fS, fG, fsS, fsG, true, YGInt[], &ok);
# endif
          jntot += jn[jj];
        }

        for (int jj=0; jj<NGS; jj++) {
          scalar YG = YGList_G[jj];
          scalar cpGv = cpGList_G[jj];
          mdeGG += cpGv[]*(jn[jj] - YG[]*jntot)*gTn;
        }
      }
      else
#endif
      {
      coord gTG = {0., 0., 0.};
      coord gYGj_G = {0., 0., 0.};
      coord gYGsum_G = {0., 0., 0.};

      foreach_dimension()
        gTG.x = mde_gradient_x (point, TG, fG);

      foreach_dimension() {
        for (int jj=0; jj<NGS; jj++) {
          scalar Dmixv = DmixGList_G[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList_G[jj];
          gYGsum_G.x -= (MWmixG_G[] > 0.) ?
            rhoGv_G[]*Dmixv[]*gas_MWs[jj]/MWmixG_G[]*mde_gradient_x (point, XG, fG) : 0.;
  # else
          scalar YG = YGList_G[jj];
          gYGsum_G.x -= rhoGv_G[]*Dmixv[]*mde_gradient_x (point, YG, fG);
  # endif
        }

        for (int jj=0; jj<NGS; jj++) {
          scalar YG = YGList_G[jj];
          scalar cpGv = cpGList_G[jj];
          scalar Dmixv = DmixGList_G[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList_G[jj];
          gYGj_G.x = (MWmixG_G[] > 0.) ?
            -rhoGv_G[]*Dmixv[]*gas_MWs[jj]/MWmixG_G[]*mde_gradient_x (point, XG, fG) : 0.;
  # else
          gYGj_G.x = -rhoGv_G[]*Dmixv[]*mde_gradient_x (point, YG, fG);
  # endif
          mdeGG += cpGv[]*(gYGj_G.x - YG[]*gYGsum_G.x)*gTG.x;
        }
      }
      }
      sGT[] -= mdeGG*cm[]*wG;
#ifdef TG_PROBE
      qmde_dbg[] = -mdeGG*cm[]*wG;
#endif
    }
  }
#endif //MASS_DIFFUSION_ENTHALPY

#ifdef SOLVE_TEMPERATURE

  /**
  Assemble the interface heat source. This used to sit inside the species
  loop above. The move is exact: `sST` and `sGT` are accumulators, the
  enthalpy block writes only bulk cells, and this writes only interface
  cells.

  `INT_TEMP_PICARD` keeps a copy of everything that does not depend on
  `TInt` — the spark of `spark.h` and the enthalpy of mass diffusion — so
  that each pass of the loop can rebuild the interface part alone. */

# if INT_TEMP_PICARD
  foreach() {
    sST_base[] = sST[];
    sGT_base[] = sGT[];
  }
# endif

  interface_temperature_sources();
#endif

#if defined VARPROP && !defined NO_EXPANSION
  update_divergence();
  // update_divergence_density();
#endif

#ifdef FICK_CORRECTED

  /**
  The largest `|u_c|/Delta` of this step, where `u_c = phic/rho` is the
  corrective velocity. The `stability` event of the next step reads it. */

  double uodx = 0.;
  face vector phicGtot[];
  foreach_face() {
    phicGtot.x[] = 0.;
    for (int jj=0; jj<NGS; jj++) {
      scalar DmixG = DmixGList_G[jj];
      double DmixGf = face_value(DmixG, 0);
      double rhoGf;
# ifdef VARPROP
      rhoGf = face_value(rhoGv_G, 0);
# else
      rhoGf = rhoG;
# endif

# ifdef MOLAR_DIFFUSION
      scalar XG = XGList_G[jj];
      double MWmixf = face_value(MWmixG_G, 0);
      phicGtot.x[] += (MWmixf > 0.) ? rhoGf*DmixGf*face_gradient_x (XG, 0)*gas_MWs[jj]/MWmixf : 0.;
# else
      scalar YG = YGList_G[jj];
      phicGtot.x[] += rhoGf*DmixGf*face_gradient_x (YG, 0);
# endif
    }
    phicGtot.x[] *= fsG.x[]*fm.x[];
  }

  face vector phicStot[];
  foreach_face() {
    phicStot.x[] = 0.;
    for (int jj=0; jj<NGS; jj++) {
      scalar DmixG = DmixGList_S[jj];
      double DmixGf = face_value(DmixG, 0);
      double rhoGf;
# ifdef VARPROP
      rhoGf = face_value(rhoGv_S, 0);
# else
      rhoGf = rhoG;
# endif

# ifdef MOLAR_DIFFUSION
      scalar XG = XGList_S[jj];
      double MWmixf = face_value(MWmixG_S, 0);
      phicStot.x[] += (MWmixf > 0.) ? rhoGf*DmixGf*face_gradient_x (XG, 0)*gas_MWs[jj]/MWmixf : 0.;
# else
      scalar YG = YGList_S[jj];
      phicStot.x[] += rhoGf*DmixGf*face_gradient_x (YG, 0);
# endif
    }
    phicStot.x[] *= fsS.x[]*fm.x[];
  }

  //Apply the Fick's law correction
  for (int jj=0; jj<NGS; jj++) {
    face vector phicjj[];
    foreach_face() {
      phicjj.x[] = phicGtot.x[];
#ifdef MOLAR_DIFFUSION
      scalar DmixG = DmixGList_G[jj];
      double DmixGf = face_value(DmixG, 0);
      double MWmixf = face_value(MWmixG_G, 0);

      double rhoGf;
# ifdef VARPROP
      rhoGf = face_value(rhoGv_G, 0);
# else
      rhoGf = rhoG;
# endif
      phicjj.x[] -= (MWmixf > 0.) ? rhoGf*DmixGf/MWmixf*face_gradient_x (MWmixG_G, 0)*fsG.x[]*fm.x[] : 0.;
#endif
    }

    /**
    Convert the corrective mass flux into a velocity before the transport.

    `tracer_fluxes()` reads its second argument as a velocity. It builds the
    Courant number `un = dt*uf/(fm*Delta)` from it, and it removes the slope
    with the factor `(1 - s*un)`. But `phicjj` is a mass flux in kg/m2/s, so
    `un` was too small by the density and the slope kept its full size at any
    Courant number.

    Divide by the face density here, and multiply the flux back after the
    call. The flux, and so the mass balance, is the same expression as
    before. Only `un` and the slope change.

    Caution: `gradients()` reads a `NULL` gradient as the unlimited centred
    slope, not as no slope. The previous `YG.gradient = NULL` therefore
    selected the least stable reconstruction, which is the opposite of what
    its comment says. `minmod2` is the limited one. */

    scalar YG = YGList_G[jj];

#if CORRECTIVE_LIMITER
    face vector rhocjj[];
#endif
    foreach_face (reduction(max:uodx)) {
      double rhoGf;
#ifdef VARPROP
      rhoGf = face_value (rhoGv_G, 0);
#else
      rhoGf = rhoG;
#endif
#if CORRECTIVE_LIMITER
      rhocjj.x[] = rhoGf;
      phicjj.x[] = (rhoGf > 0.) ? phicjj.x[]/rhoGf : 0.;
      if (fm.x[] > 0.)
        uodx = max (uodx, fabs (phicjj.x[])/(fm.x[]*Delta));
#else
      if (fm.x[] > 0. && rhoGf > 0.)   // record it, but do not change the flux
        uodx = max (uodx, fabs (phicjj.x[])/(rhoGf*fm.x[]*Delta));
#endif
    }

    double (* gradient_backup)(double, double, double) = YG.gradient; // we need to backup the gradient function
#if CORRECTIVE_LIMITER
    YG.gradient = minmod2; // NULL means the unlimited centred slope, not no slope
#else
    YG.gradient = NULL;    // the previous choice: the unlimited centred slope
#endif
    face vector flux[];
    tracer_fluxes (YG, phicjj, flux, dt, zeroc); //calculate the fluxes using the corrective velocity
    YG.gradient = gradient_backup; // restore the gradient function

#if CORRECTIVE_LIMITER
    // back to a mass flux, so that the balance is untouched
    foreach_face()
      flux.x[] *= rhocjj.x[];
#endif

    // apply the corrective fluxes
    foreach()
      foreach_dimension()
        YG[] += (rhoGv_G[] > 0.) ? dt/(rhoGv_G[])*(flux.x[] - flux.x[1])/(Delta*cm[]) : 0.;
  }

  for (int jj=0; jj<NGS; jj++) {
    face vector phicjj[];
    foreach_face() {
      phicjj.x[] = phicStot.x[];
#ifdef MOLAR_DIFFUSION
      scalar DmixG = DmixGList_S[jj];
      double DmixGf = face_value(DmixG, 0);
      double MWmixf = face_value(MWmixG_S, 0);

      double rhoGf;
# ifdef VARPROP
      rhoGf = face_value(rhoGv_S, 0);
# else
      rhoGf = rhoG;
# endif
      phicjj.x[] -= (MWmixf > 0.) ? rhoGf*DmixGf/MWmixf*face_gradient_x (MWmixG_S, 0)*fsS.x[]*fm.x[] : 0.;
#endif
  }

    /**
    Convert the corrective mass flux into a velocity before the transport.

    `tracer_fluxes()` reads its second argument as a velocity. It builds the
    Courant number `un = dt*uf/(fm*Delta)` from it, and it removes the slope
    with the factor `(1 - s*un)`. But `phicjj` is a mass flux in kg/m2/s, so
    `un` was too small by the density and the slope kept its full size at any
    Courant number.

    Divide by the face density here, and multiply the flux back after the
    call. The flux, and so the mass balance, is the same expression as
    before. Only `un` and the slope change.

    Caution: `gradients()` reads a `NULL` gradient as the unlimited centred
    slope, not as no slope. The previous `YG.gradient = NULL` therefore
    selected the least stable reconstruction, which is the opposite of what
    its comment says. `minmod2` is the limited one. */

    scalar YG = YGList_S[jj];

#if CORRECTIVE_LIMITER
    face vector rhocjj[];
#endif
    foreach_face (reduction(max:uodx)) {
      double rhoGf;
#ifdef VARPROP
      rhoGf = face_value (rhoGv_S, 0);
#else
      rhoGf = rhoG;
#endif
#if CORRECTIVE_LIMITER
      rhocjj.x[] = rhoGf;
      phicjj.x[] = (rhoGf > 0.) ? phicjj.x[]/rhoGf : 0.;
      if (fm.x[] > 0.)
        uodx = max (uodx, fabs (phicjj.x[])/(fm.x[]*Delta));
#else
      if (fm.x[] > 0. && rhoGf > 0.)   // record it, but do not change the flux
        uodx = max (uodx, fabs (phicjj.x[])/(rhoGf*fm.x[]*Delta));
#endif
    }

    double (* gradient_backup)(double, double, double) = YG.gradient; // we need to backup the gradient function
#if CORRECTIVE_LIMITER
    YG.gradient = minmod2; // NULL means the unlimited centred slope, not no slope
#else
    YG.gradient = NULL;    // the previous choice: the unlimited centred slope
#endif
    face vector flux[];
    tracer_fluxes (YG, phicjj, flux, dt, zeroc); //calculate the fluxes using the corrective velocity
    YG.gradient = gradient_backup; // restore the gradient function

#if CORRECTIVE_LIMITER
    // back to a mass flux, so that the balance is untouched
    foreach_face()
      flux.x[] *= rhocjj.x[];
#endif

    // apply the corrective fluxes
    foreach()
      foreach_dimension()
        YG[] += (rhoGv_S[] > 0.) ? dt/(rhoGv_S[])*(flux.x[] - flux.x[1])/(Delta*cm[]) : 0.;
  }

  corrective_uodx = uodx;
  #endif //FICK_CORRECTED

  scalar theta1[], theta2[];

#if TREE
  theta1.refine = fraction_refine;
  set_prolongation (theta1, fraction_refine);
  theta2.refine = fraction_refine;
  set_prolongation (theta2, fraction_refine);
#endif

  // Internal gas diffusion
  for (int jj=0; jj<NGS; jj++) {
    face vector DmixGf[];
    scalar DmixG = DmixGList_S[jj];
    foreach_face() {
      double rhoGvh_S;
#ifdef VARPROP
      rhoGvh_S = face_value(rhoGv_S, 0);
#else
      rhoGvh_S = rhoG;
#endif
      DmixGf.x[] = face_value(DmixG, 0)*rhoGvh_S*fsS.x[]*fm.x[];
    }

    foreach() {
#ifdef VARPROP
      theta1[] = cm[]*max(rhoGv_S[]*porosity[], F_ERR); // porosity is already multiplied by fS
#else
      theta1[] = cm[]*max(rhoG*porosity[], F_ERR); // porosity is already multiplied by fS
#endif
    }

    scalar YG = YGList_S[jj];
    scalar sSexp = sSexpList[jj];

#ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (YG, dt, D=DmixGf, theta=theta1);
#else
    diffusion (YG, dt, D=DmixGf, r=sSexp, theta=theta1);
#endif
  }

  //external diffusion
  for (int jj=0; jj<NGS; jj++) {
    face vector DmixGf[];
    scalar DmixG = DmixGList_G[jj];
    foreach_face() {
      double rhoGvh_G;
#ifdef VARPROP
      rhoGvh_G = face_value(rhoGv_G, 0);
#else
      rhoGvh_G = rhoG;
#endif
      DmixGf.x[] = face_value(DmixG, 0)*rhoGvh_G*fsG.x[]*fm.x[];
    }
    foreach() {
#ifdef VARPROP
      theta2[] = cm[]*max(fG[]*rhoGv_G[], F_ERR);
#else
      theta2[] = cm[]*max(fG[]*rhoG, F_ERR);
#endif
    }

    scalar YG = YGList_G[jj];
    scalar sGexp = sGexpList[jj];

#ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (YG, dt, D=DmixGf, theta=theta2);
#else
    diffusion (YG, dt, D=DmixGf, r=sGexp, theta=theta2);
#endif
  }

#ifdef SOLVE_TEMPERATURE

/**
## The Picard loop on the interface temperature

The scheme above is partitioned and does no iteration: it builds `TInt` from
the fields of step `n`, freezes the two interface fluxes into `sST` and `sGT`,
then solves each phase alone. The interface condition is therefore explicit
while the interior is implicit, and the lag of `TInt` is first order in `dt`.

`INT_TEMP_PICARD` repeats the pair of solves. Each pass restarts from the
fields of step `n`, rebuilds the interface source with the newest `TInt`, and
solves again. At the fixed point the flux balance holds with the new fields on
both sides, which is the fully implicit interface condition.

Each pass must rebuild the source, the conductance and the heat capacity,
because `diffusion()` destroys all three of its `r`, `beta` and `theta`
arguments (`$BASILISK/diffusion.h`). That is why the whole block is inside the
loop and not only the two solves.

`INT_TEMP_PICARD_MAXITER = 0` gives the present code exactly. Use it as the
inertness control.

Caution: `ijc_CoupledTemperature()` skips a cell whose `TS` or `TG` is not
positive, and a skipped cell keeps its old `TInt`. Such a cell adds nothing to
`dTInt_max` and so it looks converged when it is not. `picard.dat` reports the
count. Do not trust `dTInt_max` on a step whose count is not zero. */

/**
Keep the temperature before the solve. The debt update below needs it to
measure `T^{n+1} - T^n`. Nothing writes `TS` or `TG` between the call to
`interface_temperature_sources()` above and this point, so the snapshot
matches the fields that built the source. */

# if INT_TEMP_ROBIN
  foreach() {
    TS_rn[] = TS[];
    TG_rn[] = TG[];
  }
# endif

# if INT_TEMP_PICARD
#  ifndef INT_TEMP_PICARD_MAXITER
#   define INT_TEMP_PICARD_MAXITER 5
#  endif
#  ifndef INT_TEMP_PICARD_TOL
#   define INT_TEMP_PICARD_TOL 1e-2
#  endif
#  ifndef INT_TEMP_PICARD_OMEGA
#   define INT_TEMP_PICARD_OMEGA 1.
#  endif
#  if defined FIXED_INT_TEMP || defined TEMPERATURE_PROFILE
#   error "INT_TEMP_PICARD needs a solved TInt. FIXED_INT_TEMP and TEMPERATURE_PROFILE force it."
#  endif

  foreach() {
    TS_n[] = TS[];
    TG_n[] = TG[];
  }

  ITP_niter = 0.;
  ITP_dTInt = 0.;
  ITP_dTInt0 = 0.;
  ITP_nskip = 0.;

  for (int picard_m = 0; picard_m <= INT_TEMP_PICARD_MAXITER; picard_m++) {

  /**
  Undo the previous pass. `diffusion()` consumed `sST`, `sGT`, `betaST` and
  `betaGT`, so restore the part that does not depend on `TInt` and add the
  interface part again with the new `TInt`. */

    if (picard_m > 0) {
      foreach() {
        TS[] = TS_n[];
        TG[] = TG_n[];
        sST[] = sST_base[];
        sGT[] = sGT_base[];
#  if INT_TEMP_ROBIN
        betaST[] = 0.;
        betaGT[] = 0.;
#  endif
      }
      interface_temperature_sources();
    }
# endif // INT_TEMP_PICARD

  foreach_face() {
    lambda1f.x[] = face_value(lambda1v.x, 0)*fsS.x[]*fm.x[];
    lambda2f.x[] = face_value(lambda2v.x, 0)*fsG.x[]*fm.x[];
  }

  foreach() {
    double theta1vh, theta2vh;
# ifdef VARPROP
    theta1vh = fS[] > F_ERR ? porosity[]/fS[]*rhoGv_S[]*cpGv_S[] + (1. - porosity[]/fS[])*rhoSv[]*cpSv[] : 0.;
    theta2vh = rhoGv_G[]*cpGv_G[];
# else
    theta1vh = fS[] > F_ERR ? porosity[]/fS[]*rhoG*cpG + (1. - porosity[]/fS[])*rhoS*cpS : 0.;
    theta2vh = rhoG*cpG;
# endif

    theta1[] = cm[]*max(fS[]*theta1vh, F_ERR);
    theta2[] = cm[]*max(fG[]*theta2vh, F_ERR);
#ifdef TG_PROBE
    TGadv_dbg[] = TG[];   // the value that the advection of this step left
#endif
  }

#ifdef VARCOEFF
  foreach()
    porosity[] = (f[] > F_ERR) ? porosity[]/f[] : 0;

  foreach_face() {
    double ef = face_value(porosity, 0);
    lambda1f.x[] = (ef > F_ERR) ? lambda1f.x[] / (rhoG*cpG*ef + rhoS*cpS*(1. - ef)) : 0.;
    lambda2f.x[] = lambda2f.x[] / (rhoG*cpG);
  }

  foreach() {
    theta1[] = cm[] * max(fS[], F_ERR);
    theta2[] = cm[] * max(fG[], F_ERR);
  }

  foreach() {
    sST[] = (f[] > F_ERR) ? sST[] / (rhoG*cpG*porosity[] + rhoS*cpS*(1. - porosity[])) : 0.;
    sGT[] = sGT[] / (rhoG*cpG);
#if INT_TEMP_ROBIN
    betaST[] = (f[] > F_ERR) ? betaST[] / (rhoG*cpG*porosity[] + rhoS*cpS*(1. - porosity[])) : 0.;
    betaGT[] = betaGT[] / (rhoG*cpG);
#endif
  }

  foreach()
    porosity[] *= f[];
#endif

/**
## The tolerance of the two temperature solves

`INT_TEMP_TOL`, on by default, scales `TOLERANCE` to the residual of each
solve. The full derivation, the measured numbers and the columns of
`tsolve.dat` are in `int-temperature-tol.h`. The short version: the inherited
`TOLERANCE = 1e-5` asks the gas temperature for 7e-12 K, `poisson.h` then
raises `nrelax` for a target it can never meet, and two runs died of the
wasted iterations.

An earlier version of this used `rhoS*cpS` and covered the solid solve alone.
That was measured as a no-op at level 10 before the flame, because
`NITERMIN = 2` binds before the tolerance does. It was still the wrong scale,
and it left the gas solve, which is the one that fails, untouched.

Caution: measure any change here with several alternating runs, normalised by
CPU time. A single pair on a loaded machine once gave 43 per cent, which was
pure scatter; eleven proper runs gave 0.6 per cent.

Caution: with `INT_TEMP_PICARD` the outer loop cannot converge below what the
linear solve delivers. Keep `INT_TEMP_TOL_K` well under
`INT_TEMP_PICARD_TOL`, and re-measure `rel_max` after any change. */

# ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (TS, dt, D=lambda1f, r=sST, theta=theta1);
    diffusion_explicit (TG, dt, D=lambda2f, r=sGT, theta=theta2);
# else

  /**
  Read the two heat capacities BEFORE either solve. `diffusion()` overwrites
  `theta` in place. See `int-temperature-tol.h` for why the inherited
  `TOLERANCE` is the wrong number here. */

#  if INT_TEMP_TOL
    double th1max = 0., th2max = 0.;
    foreach (reduction(max:th1max) reduction(max:th2max)) {
      th1max = max (th1max, theta1[]);
      th2max = max (th2max, theta2[]);
    }
    double tol_save = TOLERANCE;
    double tolS = max (tol_save, th1max/dt*INT_TEMP_TOL_K);
    double tolG = max (tol_save, th2max/dt*INT_TEMP_TOL_K);
    ITT_tolS = tolS;
    ITT_tolG = tolG;
    mgstats mgS, mgG;
    mgG.i = 0; mgG.nrelax = 0; mgG.resa = 0.;
#  endif

#  if INT_TEMP_ROBIN
#   if INT_TEMP_TOL
    TOLERANCE = tolS;
    mgS = diffusion (TS, dt, D=lambda1f, r=sST, beta=betaST, theta=theta1);
#   else
    diffusion (TS, dt, D=lambda1f, r=sST, beta=betaST, theta=theta1);
#   endif
#   ifndef TEMPERATURE_PROFILE
#    if INT_TEMP_TOL
    TOLERANCE = tolG;
    mgG = diffusion (TG, dt, D=lambda2f, r=sGT, beta=betaGT, theta=theta2);
#    else
    diffusion (TG, dt, D=lambda2f, r=sGT, beta=betaGT, theta=theta2);
#    endif
#   endif
#  else
#   if INT_TEMP_TOL
    TOLERANCE = tolS;
    mgS = diffusion (TS, dt, D=lambda1f, r=sST, theta=theta1);
#   else
    diffusion (TS, dt, D=lambda1f, r=sST, theta=theta1);
#   endif
#   ifndef TEMPERATURE_PROFILE
#    if INT_TEMP_TOL
    TOLERANCE = tolG;
    mgG = diffusion (TG, dt, D=lambda2f, r=sGT, theta=theta2);
#    else
    diffusion (TG, dt, D=lambda2f, r=sGT, theta=theta2);
#    endif
#   endif
#  endif

#  if INT_TEMP_TOL
    TOLERANCE = tol_save;
    ITT_iS = mgS.i; ITT_nrelaxS = mgS.nrelax; ITT_resaS = mgS.resa;
    ITT_iG = mgG.i; ITT_nrelaxG = mgG.nrelax; ITT_resaG = mgG.resa;
#  endif

#ifdef TG_PROBE
  {
    static FILE * fpt = NULL;
    static int nt = 0;
    foreach (serial) {
      double gfr = 1. - f[];
      if (gfr > F_ERR && nt < TG_PROBE_MAX) {
        double Tpre  = TGadv_dbg[]/gfr;
        double Tpost = TG[]/gfr;
        if (Tpre < TG_PROBE_TMIN || Tpost < TG_PROBE_TMIN) {
          if (!fpt) {
            char nm[80];
            snprintf (nm, sizeof(nm), "tgprobe-%d.dat", pid());
            fpt = fopen (nm, "w");
            fprintf (fpt, "#t i x y level dt f fG theta2 rhoGv_G cpGv_G"
                          " TGadv TGpost Tpre_phys Tpost_phys sGT betaGT"
                          " qint qmde qrob lam2L lam2R\n");
          }
          fprintf (fpt, "%g %d %g %g %d %g %.17g %.17g %.17g %.17g %.17g"
                        " %.17g %.17g %.17g %.17g %.17g %.17g"
                        " %.17g %.17g %.17g %.17g %.17g\n",
                   t, i, x, y, level, dt, f[], fG[], theta2[],
                   rhoGv_G[], cpGv_G[],
                   TGadv_dbg[], TG[], Tpre, Tpost, sGT[], betaGT[],
                   qint_dbg[], qmde_dbg[], qrob_dbg[],
                   lambda2f.x[], lambda2f.x[1]);
          fflush (fpt);
          nt++;
        }
      }
    }
  }
#endif
# endif

# if INT_TEMP_PICARD

  /**
  The last pass is the answer. Do not rebuild `TInt` after it: nothing would
  use the new value, and `update_divergence()` already ran with the first one. */

    if (picard_m == INT_TEMP_PICARD_MAXITER)
      break;

    foreach()
      TInt_prev[] = TInt[];

    ijc_CoupledTemperature();

  /**
  Under-relaxation. The default weight is 1, which is no relaxation at all,
  so the branch costs one comparison per step. Lower the weight if
  `picard.dat` shows that `dTInt` does not fall from pass to pass. */

    if (INT_TEMP_PICARD_OMEGA != 1.) {
      double w = INT_TEMP_PICARD_OMEGA;
      foreach()
        if (f[] > F_ERR && f[] < 1. - F_ERR)
          TInt[] = w*TInt[] + (1. - w)*TInt_prev[];
    }

  /**
  The change of this pass, and the cells that `ijc_CoupledTemperature()` did
  not solve. Its guard is `f[] > F_ERR && f[] < 1.-F_ERR && TS[] > 0. &&
  TG[] > 0.`, so a cell that fails the last two keeps its old `TInt` and does
  not appear in `dTInt_max`. Count it, or the loop reports convergence in a
  cell it never touched. */

    double dTInt_max = 0.;
    double nint = 0., nskip = 0.;
    foreach (reduction(max:dTInt_max) reduction(+:nint) reduction(+:nskip))
      if (f[] > F_ERR && f[] < 1. - F_ERR) {
        nint += 1.;
        if (TS[] > 0. && TG[] > 0.)
          dTInt_max = max (dTInt_max, fabs (TInt[] - TInt_prev[]));
        else
          nskip += 1.;
      }

  /**
  Keep the change of the first pass as well. The contraction rate of the map
  is the ratio of the last change to the first, over the passes between them.
  A ratio taken between two steps measures nothing. */

    if (picard_m == 0)
      ITP_dTInt0 = dTInt_max;

    ITP_niter = picard_m + 1.;
    ITP_dTInt = dTInt_max;
    ITP_nint  = nint;
    ITP_nskip = nskip;

    if (dTInt_max < INT_TEMP_PICARD_TOL)
      break;
  }
# endif // INT_TEMP_PICARD

/**
The debt of this step: the heat that the conductance withheld, as a rate. The
next step adds it to the source. `KSf` and `KGf` hold the conductance of the
LAST pass, which is the pass that produced these temperatures.

The guard is the same one that `interface_temperature_sources()` uses, so a
cell that is no longer an interface cell gets a debt of zero. That costs one
step of heat, which is the same bound the carry itself has. */

# if INT_TEMP_ROBIN
  ITR_debt_max = 0.;
  ITR_debt_l1 = 0.;
  ITR_dTcap_max = 0.;
  ITR_ncap = 0.;

  foreach (reduction(max:ITR_debt_max) reduction(+:ITR_debt_l1)
           reduction(max:ITR_dTcap_max) reduction(+:ITR_ncap)) {
    debtST[] = 0.;
    debtGT[] = 0.;
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      double dTS = TS[] - TS_rn[], dTG = TG[] - TG_rn[];
      debtST[] = KSf[]*dTS;
      debtGT[] = KGf[]*dTG;

      if (KSf[] > 0. || KGf[] > 0.) {
        ITR_ncap += 1.;
        ITR_debt_max = max (ITR_debt_max,
                            max (fabs (debtST[]), fabs (debtGT[])));
        ITR_debt_l1 += (fabs (debtST[]) + fabs (debtGT[]))*dv();
        ITR_dTcap_max = max (ITR_dTcap_max, max (fabs (dTS), fabs (dTG)));
      }
    }
  }
# endif

/**
Measure the interface balance again, now with the new fields. `TS` and `TG`
still hold the value of one phase here; the block below puts them back into
tracer form. */

# if INT_TEMP_PROBE
  int_temperature_probe();
# endif
#endif

  //recover tracer form
  foreach() {
    for (scalar YG in YGList_S)
      YG[] = (f[] > F_ERR) ? YG[]*f[] : 0.;

    for (scalar YG in YGList_G)
      YG[] = ((1. - f[]) > F_ERR) ? YG[]*(1. - f[]) : 0.;
    
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]*f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]*(1. - f[]) : 0.;
    T[] = TS[] + TG[];
#endif
  }

  check_and_correct_fractions (YGList_S, NGS, false);
  check_and_correct_fractions (YGList_G, NGS, true);
}

/* 
This is actually a tracer_advection step.
We put it here so that it is executed after the default
tracer advection step.
*/
extern face vector ufsave;
face vector u_prime[];
event tracer_diffusion (i++,last) {

foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    f[] = (f[] < 1.-F_ERR) ? f[] : 1.;
    fS[] = f[]; fG[] = 1. - f[];
  }

  //Compute face gradients
  face_fraction (fS, fsS);
  face_fraction (fG, fsG);

  check_and_correct_fractions (YGList_S, NGS, false);
  check_and_correct_fractions (YGList_G, NGS, true);
  check_and_correct_fractions (YSList,   NSS, false);

#ifdef VARPROP
  update_properties();
#else
  update_properties_constant();
#endif

  // lose tracer form and extrapolate fields
  foreach() {
    porosity[] = (f[] > F_ERR) ? porosity[]/f[] : 0.;
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]/f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]/(1. - f[]) : 0.;

    TS[] = (f[] > F_ERR) ? TS[] : TG[];
    TG[] = (f[] < 1. - F_ERR) ? TG[] : TS[];
#endif

    for (int jj=0; jj<NGS; jj++) { 
      scalar YG_S = YGList_S[jj];
      scalar YG_G = YGList_G[jj];

      YG_S[] = (f[] > F_ERR) ? YG_S[]/f[] : 0.;
      YG_G[] = (f[] < 1. - F_ERR) ? YG_G[]/(1. - f[]) : 0.;

      YG_S[] = (f[] > F_ERR) ? YG_S[] : YG_G[];
      YG_G[] = (f[] < 1. - F_ERR) ? YG_G[] : YG_S[];
    }
  }

  advection_div(YGList_S, ufsave, dt);
  advection_div(YGList_G, ufsave, dt);

#ifdef SOLVE_TEMPERATURE
  foreach_face() {
    double ef = clamp(face_value(porosity, 0), 0., 1.);

    double rhoGvh_S, rhoSvh;
    double cpGvh_S, cpSvh;

    #ifdef VARPROP
    rhoGvh_S = face_value(rhoGv_S, 0); rhoSvh = face_value(rhoSv, 0);
    cpGvh_S = face_value(cpGv_S, 0); cpSvh = face_value(cpSv, 0);
    #else
    rhoGvh_S = rhoG; rhoSvh = rhoS;
    cpGvh_S = cpG; cpSvh = cpS;
    #endif

    double denom = rhoGvh_S*cpGvh_S*ef + rhoSvh*cpSvh*(1. - ef);
    u_prime.x[] = (fsS.x[] > F_ERR && denom > 0.) ? // denom==0 in unfilled solid faces -> 0/0
                  fsS.x[]*ufsave.x[]*(rhoGvh_S*cpGvh_S)/denom
                  : 0.;
  }

  advection_div({TS}, u_prime, dt);
# ifndef TEMPERATURE_PROFILE
  advection_div({TG}, ufsave, dt);
# endif
#endif

  // recover tracer form
  foreach() {
    porosity[] = (f[] > F_ERR) ? porosity[]*f[] : 0.;
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]*f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]*(1. - f[]) : 0.;
#endif

    for (int jj=0; jj<NGS; jj++) { 
      scalar YG_S = YGList_S[jj];
      scalar YG_G = YGList_G[jj];

      YG_S[] = (f[] > F_ERR) ? YG_S[]*f[] : 0.;
      YG_G[] = (f[] < 1. - F_ERR) ? YG_G[]*(1. - f[]) : 0.;
    }
  }
}
