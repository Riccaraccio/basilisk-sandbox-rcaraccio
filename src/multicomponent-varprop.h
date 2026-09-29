#ifndef MULTICOMPONENT
  #define MULTICOMPONENT 1
#endif

#include "intgrad.h"

#include "diffusion.h"

/**
`EXPLICIT_DIFFUSION` and `VARCOEFF` were removed on the `interface-vofbc`
branch. The interface heat flux now writes a term on the diagonal through the
`beta` argument of `diffusion()`. `diffusion_explicit()` applied `beta`
explicitly, which needs `dt < theta/|beta|` — the very constraint that this
work removes — and `VARCOEFF` rescaled each row by a heat capacity built from
constants while `theta` used the variable properties. Neither path had a live
user. */

#ifdef EXPLICIT_DIFFUSION
# error "EXPLICIT_DIFFUSION was removed. It cannot carry the interface diagonal."
#endif
#ifdef VARCOEFF
# error "VARCOEFF was removed. See the comment in multicomponent-varprop.h."
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
The record of the interface conductance, and the reason it must give back the
heat that it holds. It changes no field either. */


/**
## The transport of heat by the pore gas

The `tracer_diffusion` event advects `TS` with `u_prime`. `u_prime` is the velocity at which the gas that
flows through the pores carries the heat of the pseudo-phase:

  u_prime = fsS*uf*rhoG*cpG/(rhoG*cpG*eps + rhoS*cpS*(1 - eps))

The gas leaves the reaction front cold and flows out through the hot char,
so this term cools the char layer when the release rate rises.

The advection of `TS` with the solid velocity `ubf` is separate: `TS` is a
tracer of `f`, and `vof.h` moves it. */

/**
## The gas side does not move with the solid velocity

`TG` and `YGList_G` are tracers of `f`. The `vof` event of `shrinking.h` sets
`uf = ubf`, and `vof.h` then moves the tracers of `f` with `ubf`. `ubf` is
not zero in the gas, because the solver solves `psi` over the full domain.
The `tracer_diffusion` event then moves the gas side again with `ufsave`.

`ufsave` is the total volume flux `(1 - eps) ubf + eps v_g`. The projection
source does not contain `prod` or `zeta`. So `ufsave` already carries the
gas that fills the volume that the shrinkage frees. In the pure gas, a
second transport with `ubf` counts that flow two times.

`shrinking.h` therefore gives `ubf` only to the faces next to a cell with
solid. A face between two pure gas cells gets 0.

Caution: in a cut cell, the `ubf` transport of the gas side is not a double
count. It fills the gas volume that the interface frees with upwind gas. The
advection with `ufsave` moves the intrinsic value and does not see the change
of `1 - f`. Without the `ubf` transport, the freed volume of a sliver cell
takes the cold value of the cell itself, and the minimum `TG` fell by 40 to
70 K. So keep `ubf` on the faces of the cut cells. */

/**
The tolerance of the two temperature solves, scaled to their own residual.
The inherited `TOLERANCE = 1e-5` asks the gas temperature for 7e-12 K, which
no solve can deliver. See `int-temperature-tol.h`. */

# include "int-temperature-tol.h"

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
## The corrective flux and the mass diffusion enthalpy term

`CORRECTIVE_CFL` is the Courant limit of the corrective flux, 0.5 by
default. 0 removes the limit. The corrective flux uses a limited slope and
the true Courant number. The mass diffusion enthalpy term also acts in the
interface cells, with the phase weight and the phase aware gradient below.
Without `SOLVE_TEMPERATURE` the term acts in the full cells only, with the
plain centred stencil. */




/**
The term feeds the two temperature solves, and the cut cell branch reads
`TInt`. Without `SOLVE_TEMPERATURE` neither of them exists, so keep the
previous gate in that build. */


#ifdef MASS_DIFFUSION_ENTHALPY
foreach_dimension()
static double mde_gradient_x (Point point, scalar a, scalar ff)
{
#ifdef SOLVE_TEMPERATURE
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

#ifdef TG_PROBE
scalar qint_dbg[], qmde_dbg[], qrob_dbg[], TGadv_dbg[];
# ifndef TG_PROBE_TMIN
#  define TG_PROBE_TMIN 100.
# endif
# ifndef TG_PROBE_MAX
#  define TG_PROBE_MAX 400
# endif
#endif

#ifdef TG_PROBE

/**
`tg_stage_check()` reports every cell whose gas temperature is negative. It
tests the sign only, so it is correct in tracer form and in physical form.
Call it at each stage of the step. The tag of the first row names the stage
that makes the value negative. */

FILE * tgs_fp = NULL;
int tgs_n = 0;
FILE * tgq_fp = NULL;
int tgq_n = 0;

#define tg_stage_check(TAG) do {                                        \
  if (tgs_n < TG_PROBE_MAX) {                                           \
    foreach (serial) {                                                  \
      if (TG[] < 0. && tgs_n < TG_PROBE_MAX) {                          \
        if (!tgs_fp) {                                                  \
          char nm[80];                                                  \
          sprintf (nm, "tgstage-%d.dat", pid());                        \
          tgs_fp = fopen (nm, "w");                                     \
          fprintf (tgs_fp, "#tag t i x y level f TG TS T\n");           \
        }                                                               \
        fprintf (tgs_fp, "%s %g %d %g %g %d %.17g %.17g %.17g %.17g\n", \
                 TAG, t, i, x, y, level, f[], TG[], TS[], T[]);          \
        fflush (tgs_fp);                                                \
        tgs_n++;                                                        \
      }                                                                 \
    }                                                                   \
  }                                                                     \
} while (0)

/**
This event runs before every other `vof` event, because the same-name events
run in reverse order of the declaration. It therefore reports the state that
`adapt` left, before the advection of this step. */

event vof (i++) {
  tg_stage_check ("0-prevof");
}
#endif

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

This assembles the two interface heat sources `sST` and `sGT`.

The function adds to the two fields. The caller must set them to a known
value first. `event reset_sources` does that once per step.

Caution: `TS` and `TG` must hold the value of one phase here, not the tracer
form. The event divides them at its start and multiplies them back at its
end. */

static void interface_temperature_sources (void)
{
  bool success = false;
  foreach() {
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      coord n = interface_source_normal (point, fS, fsS), p;
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

      /**
      The interface heat flux. `update_divergence()` also reads `sST` and
      `sGT`, as the interface heat of the thermal expansion.

      Caution: `ebmgrad` builds the gradient from the NEIGHBOURS of the cell
      along the normal, so the source that heats a cut cell does not answer
      to the temperature of that cell. In a gas sliver the heat capacity
      `theta2 = cm*fG*rhoG*cpG` goes to zero with the gas fraction, and the
      exchange number `S = dt*lambda*h*aov/theta2` has no bound. A measured
      run reached `S` of order 100, and `TG` of one cell went through zero
      and then grew by a factor 2.8 per step. Nothing in this scheme bounds
      `TG` there, and the solid side has the same weakness when `fS` is
      small. The production cases write `TGmin_gas` and `nTGneg` as the
      early warning. The tag `oscillation-campaign-2026-09` keeps two
      remedies: a conductance on the diagonal (`INT_TEMP_ROBIN`) and a
      Dirichlet condition inside the operator (`INT_TEMP_VOFBC`). */

      sST[] += Sheatflux*aov;
      sGT[] += Gheatflux*aov;
#ifdef TG_PROBE
      qint_dbg[] = Gheatflux*aov;
#endif

    }
  }

}

#endif

/**
`DRHODT_BUDGET` measures the time level of the transport part of `drhodt`:
the explicit fluxes of `update_divergence()` against the implicit solves of
the event below. It does not change the run. See `drhodt-budget.h`. */

/**
## The transport part of `drhodt` from the implicit solves

When `DRI_ON` is 1 (see `multicomponent-properties.h`),
`update_divergence()` skips the diffusion fluxes and the interface sources.
The event below adds them after the solves, as the change that each solve
really made:

    theta*(X^{n+1} - X**)/dt = div (D grad X^{n+1}) + r + beta*X^{n+1}

`theta` is the capacity of the solve and carries `cm`, so the term is per
unit volume of the cell, as the flux form of `update_divergence()` is. The
weights and the divisors are the ones of `update_divergence()`, and the
divisors use the state before the solves, as there.

For the species, the expansion reads only `sum_j dY_j/MW_j`. The species
part therefore needs one sum for each phase, not one field for each species.

The corrected `drhodt` reaches the projection of the same step:
`project_sf()` reads `drhodt` in `advection_term` and in `projection`, and
both events run after `tracer_diffusion`. */

#if DRI_ON

/**
`dri_YS` and `dri_YG` hold `sum_j theta*(Y_j^{n+1} - Y_j**)/MW_j`. `dri_thY`
keeps the capacity of one species solve, because `diffusion()` overwrites
`theta` in place. `dri_TS`, `dri_TG`, `dri_th1` and `dri_th2` keep the
temperatures and the capacities before the temperature solves. `dri_cT` is
the temperature part of the change of `drhodt`. */

scalar dri_YS[], dri_YG[], dri_thY[];
# ifdef SOLVE_TEMPERATURE
scalar dri_TS[], dri_TG[], dri_th1[], dri_th2[], dri_cT[];
# endif

event defaults (i = 0) {
  dri_YS.nodump = dri_YG.nodump = dri_thY.nodump = true;
# ifdef SOLVE_TEMPERATURE
  dri_TS.nodump = dri_TG.nodump = true;
  dri_th1.nodump = dri_th2.nodump = dri_cT.nodump = true;
# endif
}

/**
The weights of the two phases in `drhodt`. Keep them identical to the last
loop of `update_divergence()`. */

static inline void dri_weights (double ff, double * wS, double * wG)
{
  *wS = ff > F_ERR ? 1. : 0.;
  *wG = ff < 1. - F_ERR ? 1. : 0.;
}

#endif // DRI_ON

#if DRHODT_BUDGET
# include "drhodt-budget.h"
#endif

/**
`SPECIES_CLAMP_PROBE` measures item TL-5: how far the species solves of the
event below push `Y` out of `[0,1]`, and the mass that the clamp at the end
of the event adds or removes. It does not change the run. See
`species-clamp-probe.h`, which sets the default of the flag to 0. */

#include "species-clamp-probe.h"

/**
## The properties after the solves

`PROPS_AFTER_SOLVES` at 1 runs the `properties` events once more at the end
of the event below, after the temperature and species solves. Then the
momentum, the Darcy drag and the projection of this step use the properties
of the new state, P^{n+1}, and the volume fraction after the VOF advection.
At 0 (the default) the code is the same as before this flag.

Caution: a call of `update_properties()` alone does not change `alphav`,
`rhov` or `mu`. The `properties` events of all the headers form one chain,
and `events.h` tests only the condition of the head of the chain. The head
is `event properties (i = 0)` of `multicomponent-properties.h`, because it is
declared last. So after step 0 the chain does not run in the step. It runs
only through `event ("properties")` in the `adapt` event of `centered.h`,
which ignores the conditions. `qcc -events` shows this. Thus at 0:

- `alphav`, `rhov` and `mu` come from the `adapt` event of the previous step.
  That is the state at the start of the step, P^n, with `f^n`.
- The Darcy drag reads `rhoGv_S` and `muGv_S` directly. These come from the
  `update_properties()` call of the `tracer_diffusion` event above, thus the
  state after the chemistry, P*.

This flag calls `event ("properties")`, which runs `update_properties()`
first and then the events that fill `alphav`, `rhov` and `mu`. The fields
are in tracer form at that point, which is the form that
`update_properties()` expects. The `adapt` event computes the properties
again after the grid changes. The cost is one more evaluation of the
properties in each step. */

#ifndef PROPS_AFTER_SOLVES
# define PROPS_AFTER_SOLVES 0
#endif

#if PROPS_AFTER_SOLVES && !defined(VARPROP)
# error "PROPS_AFTER_SOLVES needs VARPROP"
#endif

event tracer_diffusion (i++) {

#ifdef TG_PROBE
  tg_stage_check ("A-postvof");
#endif

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
      coord n = interface_source_normal (point, fS, fsS), p;
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
    The weight of the two sides is the phase fraction, which is the weight
    that `theta1` and `theta2` carry. Without `SOLVE_TEMPERATURE` the weight
    is 1 in a full cell and 0 in every other cell. */

#ifdef SOLVE_TEMPERATURE
    double wS = fS[], wG = fG[];
    bool interfacial = (f[] > F_ERR && f[] < 1. - F_ERR);
#else
    double wS = (fS[] > 1. - F_ERR) ? 1. : 0.;
    double wG = (fG[] > 1. - F_ERR) ? 1. : 0.;
#endif

    if (wS > F_ERR) { //Internal gas phase
      double mdeGS = 0.;

#ifdef SOLVE_TEMPERATURE
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

#ifdef SOLVE_TEMPERATURE
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
  cells. */


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

    face vector rhocjj[];
    foreach_face (reduction(max:uodx)) {
      double rhoGf;
#ifdef VARPROP
      rhoGf = face_value (rhoGv_G, 0);
#else
      rhoGf = rhoG;
#endif
      rhocjj.x[] = rhoGf;
      phicjj.x[] = (rhoGf > 0.) ? phicjj.x[]/rhoGf : 0.;
      if (fm.x[] > 0.)
        uodx = max (uodx, fabs (phicjj.x[])/(fm.x[]*Delta));
    }

    double (* gradient_backup)(double, double, double) = YG.gradient; // we need to backup the gradient function
    YG.gradient = minmod2; // NULL means the unlimited centred slope, not no slope
    face vector flux[];
    tracer_fluxes (YG, phicjj, flux, dt, zeroc); //calculate the fluxes using the corrective velocity
    YG.gradient = gradient_backup; // restore the gradient function

    // back to a mass flux, so that the balance is untouched
    foreach_face()
      flux.x[] *= rhocjj.x[];

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

    face vector rhocjj[];
    foreach_face (reduction(max:uodx)) {
      double rhoGf;
#ifdef VARPROP
      rhoGf = face_value (rhoGv_S, 0);
#else
      rhoGf = rhoG;
#endif
      rhocjj.x[] = rhoGf;
      phicjj.x[] = (rhoGf > 0.) ? phicjj.x[]/rhoGf : 0.;
      if (fm.x[] > 0.)
        uodx = max (uodx, fabs (phicjj.x[])/(fm.x[]*Delta));
    }

    double (* gradient_backup)(double, double, double) = YG.gradient; // we need to backup the gradient function
    YG.gradient = minmod2; // NULL means the unlimited centred slope, not no slope
    face vector flux[];
    tracer_fluxes (YG, phicjj, flux, dt, zeroc); //calculate the fluxes using the corrective velocity
    YG.gradient = gradient_backup; // restore the gradient function

    // back to a mass flux, so that the balance is untouched
    foreach_face()
      flux.x[] *= rhocjj.x[];

    // apply the corrective fluxes
    foreach()
      foreach_dimension()
        YG[] += (rhoGv_S[] > 0.) ? dt/(rhoGv_S[])*(flux.x[] - flux.x[1])/(Delta*cm[]) : 0.;
  }

  corrective_uodx = uodx;
  #endif //FICK_CORRECTED

  scalar theta1[], theta2[];

#ifdef TG_PROBE

/**
`tg_source_check()` reports the full energy balance of every cell whose gas
temperature is negative. Use it only inside this event, because `theta2` lives
here. `TGadv` is the value before the solve. */

#define tg_source_check(TAG) do {                                       \
  if (tgq_n < TG_PROBE_MAX) {                                           \
    foreach (serial) {                                                  \
      if (TG[] < 0. && tgq_n < TG_PROBE_MAX) {                          \
        if (!tgq_fp) {                                                  \
          char nm[80];                                                  \
          sprintf (nm, "tgsource-%d.dat", pid());                       \
          tgq_fp = fopen (nm, "w");                                     \
          fprintf (tgq_fp, "#tag t i x y level dt f fG theta2 sGT"      \
                           " betaGT qint qmde qrob TGadv TG TS\n");     \
        }                                                               \
        fprintf (tgq_fp, "%s %g %d %g %g %d %g %.17g %.17g %.17g %.17g" \
                         " %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",\
                 TAG, t, i, x, y, level, dt, f[], fG[], theta2[],       \
                 sGT[], betaGT[], qint_dbg[], qmde_dbg[], qrob_dbg[],   \
                 TGadv_dbg[], TG[], TS[]);                              \
        fflush (tgq_fp);                                                \
        tgq_n++;                                                        \
      }                                                                 \
    }                                                                   \
  }                                                                     \
} while (0)
#endif

#if TREE
  theta1.refine = fraction_refine;
  set_prolongation (theta1, fraction_refine);
  theta2.refine = fraction_refine;
  set_prolongation (theta2, fraction_refine);
#endif

#if DRI_ON
  foreach() {
    dri_YS[] = 0.;
    dri_YG[] = 0.;
  }
#endif
#if SPECIES_CLAMP_PROBE
  species_clamp_presolve();
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

#if DRI_ON
    foreach() {
      dri_thY[] = theta1[];
      dri_YS[] -= theta1[]*YG[]/gas_MWs[jj];
    }
#endif

    diffusion (YG, dt, D=DmixGf, r=sSexp, theta=theta1);

#if DRI_ON
    foreach()
      dri_YS[] += dri_thY[]*YG[]/gas_MWs[jj];
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

#if DRI_ON
    foreach() {
      dri_thY[] = theta2[];
      dri_YG[] -= theta2[]*YG[]/gas_MWs[jj];
    }
#endif

    diffusion (YG, dt, D=DmixGf, r=sGexp, theta=theta2);

#if DRI_ON
    foreach()
      dri_YG[] += dri_thY[]*YG[]/gas_MWs[jj];
#endif
  }

#if DRI_ON

  /**
  Add the species part of the transport to `drhodt`. The factor
  `MWmixG/rhoGv` is the one of `update_divergence()`, and the sign follows
  its last line. */

  if (dt > 0.)
    foreach() {
      double wS, wG;
      dri_weights (f[], &wS, &wG);
      double cS = (rhoGv_S[] > 0.) ? MWmixG_S[]/rhoGv_S[] : 0.;
      double cG = (rhoGv_G[] > 0.) ? MWmixG_G[]/rhoGv_G[] : 0.;
      drhodt[] -= (wS*cS*dri_YS[] + wG*cG*dri_YG[])/dt;
    }
#endif

#ifdef SOLVE_TEMPERATURE

/**
## The time level of the interface condition

The scheme is partitioned and does no iteration: it builds `TInt` from the
fields of step `n`, freezes the two interface fluxes into `sST` and `sGT`,
then solves each phase alone. The interface condition is therefore explicit
while the interior is implicit, and the lag of `TInt` is first order in `dt`.
A Picard loop over the pair of solves removed the lag to round-off, but the
median lag was 0.35 K and the loop cost 12.7 per cent more, so the code does
not iterate. */

/**
Keep the temperature before the solve. The debt update below needs it to
measure `T^{n+1} - T^n`. Nothing writes `TS` or `TG` between the call to
`interface_temperature_sources()` above and this point, so the snapshot
matches the fields that built the source. */




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

#if DRI_ON

  /**
  Keep the state and the capacities before the solves. */

  foreach() {
    dri_TS[] = TS[];
    dri_TG[] = TG[];
    dri_th1[] = theta1[];
    dri_th2[] = theta2[];
  }
#endif

#if DRHODT_BUDGET
  drhodt_budget_presolve (theta1, theta2);
#endif


/**
## The tolerance of the two temperature solves

`int-temperature-tol.h` scales `TOLERANCE` to the residual of each
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
pure scatter; eleven proper runs gave 0.6 per cent. */

#ifdef TG_PROBE
  tg_stage_check ("B-presolve");
#endif

  /**
  Read the two heat capacities BEFORE either solve. `diffusion()` overwrites
  `theta` in place. See `int-temperature-tol.h` for why the inherited
  `TOLERANCE` is the wrong number here. */

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

    TOLERANCE = tolS;
    mgS = diffusion (TS, dt, D=lambda1f, r=sST, theta=theta1);
#   ifndef TEMPERATURE_PROFILE
    TOLERANCE = tolG;
    mgG = diffusion (TG, dt, D=lambda2f, r=sGT, theta=theta2);
#   endif

    TOLERANCE = tol_save;
    ITT_iS = mgS.i; ITT_nrelaxS = mgS.nrelax; ITT_resaS = mgS.resa;
    ITT_iG = mgG.i; ITT_nrelaxG = mgG.nrelax; ITT_resaG = mgG.resa;

#if DRI_ON

  /**
  The temperature part of the transport, from the change of each solve. The
  weights and the divisors are the ones of `update_divergence()`, with the
  state before the solves. Keep it in `dri_cT` and add it to `drhodt`
  below. */

  if (dt > 0.)
    foreach() {
      double wS, wG;
      dri_weights (f[], &wS, &wG);
      double eps = f[] > F_ERR ? porosity[]/f[] : 0.;
      double cS = (dri_TS[]*rhoGv_S[]*cpGv_S[] > 0.) ?
        eps/(dri_TS[]*(rhoGv_S[]*cpGv_S[]*eps + rhoSv[]*cpSv[]*(1. - eps))) : 0.;
      double cG = (dri_TG[]*rhoGv_G[]*cpGv_G[] > 0.) ?
        1./(dri_TG[]*rhoGv_G[]*cpGv_G[]) : 0.;
      dri_cT[] = - (wS*cS*dri_th1[]*(TS[] - dri_TS[])
                    + wG*cG*dri_th2[]*(TG[] - dri_TG[]))/dt;
    }
  else
    foreach()
      dri_cT[] = 0.;
#endif

#if DRHODT_BUDGET
  drhodt_budget_postsolve();
#endif

#ifdef TG_PROBE
  tg_source_check ("S-postsolve");
#endif

#ifdef TG_PROBE
  {
    static FILE * fpt = NULL;
    static int nt = 0;
    foreach (serial) {
      double gfr = 1. - f[];
      if (gfr > F_ERR && nt < TG_PROBE_MAX) {
        double Tpre  = TGadv_dbg[];   // TG is physical inside this event
        double Tpost = TG[];
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


/**
Add the temperature part of the transport to `drhodt`. Nothing reads
`drhodt` between here and the projection except the probe. */

# if DRI_ON
  foreach()
    drhodt[] += dri_cT[];
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

#if SPECIES_CLAMP_PROBE
  species_clamp_postsolve();
#endif

  check_and_correct_fractions (YGList_S, NGS, false);
  check_and_correct_fractions (YGList_G, NGS, true);

#ifdef TG_PROBE
  tg_stage_check ("C-end433");
#endif

#if PROPS_AFTER_SOLVES
  event ("properties");
#endif
}

/* 
This is actually a tracer_advection step.
We put it here so that it is executed after the default
tracer advection step.
*/
extern face vector ufsave;
face vector u_prime[];

/**
## The velocity of the pore species

`uf` is the superficial velocity: the projection makes its divergence equal
to the gas source per unit volume of the cell. The diffusion solve of `YG_S`
uses the capacity `rho*eps*f` (`theta1`). So the equation of the pore species
in non-conservative form is

  eps*rho*dY/dt + rho*u.grad(Y) = div(rho*D*grad(Y)) + ...

and the consistent velocity is the interstitial velocity `u/eps`.
`advection_div` with `NO_ADVECTION_DIV` gives `-dt*u.grad(Y)` for any face
velocity, thus the velocity that it receives must be `u/eps`.

The event below therefore moves `YGList_S` with

  u_pore = ufsave/max (face_value (e1), PORE_EPS_MIN),  e1 = f*eps + 1 - f

`e1` is the one-fluid porosity that `centered-phasechange.h` also uses. It is
`eps` in a full solid cell and 1 in the gas, so the gas side does not change.
`PORE_EPS_MIN` stops a division by a small porosity.

Caution: `u_pore` is `1/eps` times larger than `uf` inside the particle. The
`stability` event below therefore adds a CFL limit on `u_pore`. It computes
the limit from `uf` and from the `f` and `porosity` of the start of the step,
because `uf` becomes `ufsave` in the `vof` event of `shrinking.h`. It only
lowers `dtmax`, and the `stability` events of `shrinking.h` and `centered.h`
run after it and use that value. `pore_dtmax` keeps the limit for output. */


#ifndef PORE_EPS_MIN
# define PORE_EPS_MIN 0.05
#endif

double pore_dtmax = HUGE; // the CFL limit of u_pore, for output only

event stability (i++) {

  /**
  `porosity` is in tracer form here, so `porosity + 1 - f` is `f*eps + 1 - f`.
  Only the faces with `e1 < 1` can give a limit below that of `uf`. The gas
  faces are left to `centered.h`. */

  scalar e1[];
  foreach()
    e1[] = porosity[] + 1. - f[];

  double dtp = HUGE;
  foreach_face (reduction(min:dtp)) {
    double ef = face_value (e1, 0);
    if (uf.x[] != 0. && ef < 1.) {
      double dtf = max (ef, PORE_EPS_MIN)*Delta*fm.x[]/fabs (uf.x[]);
      if (dtf < dtp)
        dtp = dtf;
    }
  }

  pore_dtmax = CFL*dtp;
  if (pore_dtmax < dtmax)
    dtmax = pore_dtmax;
}

event tracer_diffusion (i++,last) {

foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    f[] = (f[] < 1.-F_ERR) ? f[] : 1.;
    fS[] = f[]; fG[] = 1. - f[];
  }

#ifdef TG_PROBE
  tg_stage_check ("D-pre1557");
#endif

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

  {
    // porosity is intrinsic here, and 0 where f <= F_ERR
    scalar e1[];
    foreach()
      e1[] = f[]*porosity[] + 1. - f[];

    face vector u_pore[];
    foreach_face()
      u_pore.x[] = ufsave.x[]/max (face_value (e1, 0), PORE_EPS_MIN);

    advection_div(YGList_S, u_pore, dt);
  }
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

#ifdef TG_PROBE
  tg_stage_check ("E-end1557");
#endif
}
