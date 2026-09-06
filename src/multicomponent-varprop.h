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

event reset_sources (i++) {
#ifdef SOLVE_TEMPERATURE
  foreach() {
    sST[] = 0.;
    sGT[] = 0.;
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
# else
      sST[] += KS*TS[];
      sGT[] += KG*TG[];
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
  foreach() {
    if (f[] > 1. - F_ERR) { //Internal gas phase
      double mdeGS = 0.;
      coord gTS = {0., 0., 0.};
      coord gYGj_S = {0., 0., 0.};
      coord gYGsum_S = {0., 0., 0.};

      foreach_dimension()
        gTS.x = (TS[1] - TS[-1])/(2.*Delta);
      
      foreach_dimension() {
        for (int jj=0; jj<NGS; jj++) {
          scalar Dmixv = DmixGList_S[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList_S[jj];
          gYGsum_S.x -= (MWmixG_S[] > 0.) ?
            rhoGv_S[]*Dmixv[]*gas_MWs[jj]/MWmixG_S[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
  # else
          scalar YG = YGList_S[jj];
          gYGsum_S.x -= rhoGv_S[]*Dmixv[]*(YG[1] - YG[-1])/(2.*Delta);
  # endif
        }

        for (int jj=0; jj<NGS; jj++) {
          scalar YG = YGList_S[jj];
          scalar cpGv = cpGList_S[jj];
          scalar Dmixv = DmixGList_S[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList_S[jj];
          gYGj_S.x = (MWmixG_S[] > 0.) ?
            -rhoGv_S[]*Dmixv[]*gas_MWs[jj]/MWmixG_S[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
  # else
          gYGj_S.x = -rhoGv_S[]*Dmixv[]*(YG[1] - YG[-1])/(2.*Delta);
  # endif
          mdeGS += cpGv[]*(gYGj_S.x - YG[]*gYGsum_S.x)*gTS.x;
        }
      }
    sST[] -= mdeGS*cm[];
    }
   
    if (f[] < F_ERR) { //Internal gas phase
      double mdeGG = 0.;
      coord gTG = {0., 0., 0.};
      coord gYGj_G = {0., 0., 0.};
      coord gYGsum_G = {0., 0., 0.};

      foreach_dimension()
        gTG.x = (TG[1] - TG[-1])/(2.*Delta);
      
      foreach_dimension() {
        for (int jj=0; jj<NGS; jj++) {
          scalar Dmixv = DmixGList_G[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList_G[jj];
          gYGsum_G.x -= (MWmixG_G[] > 0.) ?
            rhoGv_G[]*Dmixv[]*gas_MWs[jj]/MWmixG_G[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
  # else
          scalar YG = YGList_G[jj];
          gYGsum_G.x -= rhoGv_G[]*Dmixv[]*(YG[1] - YG[-1])/(2.*Delta);
  # endif
        }

        for (int jj=0; jj<NGS; jj++) {
          scalar YG = YGList_G[jj];
          scalar cpGv = cpGList_G[jj];
          scalar Dmixv = DmixGList_G[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList_G[jj];
          gYGj_G.x = (MWmixG_G[] > 0.) ?
            -rhoGv_G[]*Dmixv[]*gas_MWs[jj]/MWmixG_G[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
  # else
          gYGj_G.x = -rhoGv_G[]*Dmixv[]*(YG[1] - YG[-1])/(2.*Delta);
  # endif
          mdeGG += cpGv[]*(gYGj_G.x - YG[]*gYGsum_G.x)*gTG.x;
        }
      }
    sGT[] -= mdeGG*cm[];
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

    scalar YG = YGList_G[jj];
    double (* gradient_backup)(double, double, double) = YG.gradient; // we need to backup the gradient function
    YG.gradient = NULL; //reset the gradient
    face vector flux[];
    tracer_fluxes (YG, phicjj, flux, dt, zeroc); //calculate the fluxes using the corrective velocity
    YG.gradient = gradient_backup; // restore the gradient function
    
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

    scalar YG = YGList_S[jj];
    double (* gradient_backup)(double, double, double) = YG.gradient; // we need to backup the gradient function
    YG.gradient = NULL; //reset the gradient
    face vector flux[];
    tracer_fluxes (YG, phicjj, flux, dt, zeroc); //calculate the fluxes using the corrective velocity
    YG.gradient = gradient_backup; // restore the gradient function
    
    // apply the corrective fluxes
    foreach()
      foreach_dimension()
        YG[] += (rhoGv_S[] > 0.) ? dt/(rhoGv_S[])*(flux.x[] - flux.x[1])/(Delta*cm[]) : 0.; 
  }
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
## The tolerance of the solid temperature solve

`TOLERANCE` is one number for every solver of the run, but each solver
compares it against a residual in its own units. The residual of
`diffusion()` is that of

    theta*(T^{n+1} - T^n)/dt = div(D grad T^{n+1}) + ...

so it carries `rho*cp*K/s`, that is W/m3. The case sets `TOLERANCE = 1e-5`
for the pressure solve, where the residual is a velocity divergence. Applied
to the solid temperature that value asks for about ten digits, because
`rhoS*cpS` is 2.79e6 here.

Scale it to the physics instead. Ask for a residual no larger than the one
that a temperature error of `INT_TEMP_TOL_K` over one step would make:

    TOLERANCE = rhoS*cpS*INT_TEMP_TOL_K/dt

Restore the old value straight after, because the next solver needs it.

The gas solve keeps the plain `TOLERANCE`. Do NOT copy this line to it with
`rhoG*cpG`: `run/test.c` never sets `cpG`, so it holds the default of 1
(`memoryallocation-varprop.h:39`) and the product is meaningless. Under
`VARPROP` the true gas heat capacity is the field `cpGv_G`, not the scalar.
A gas version needs a scale of its own, and it is not written yet.

This is OFF unless the case defines `INT_TEMP_TOL_K`. It is off because the
measurement below shows it does nothing here, not because it is wrong.

Measured 2026-09-04, level 10, the conductance build with the probe: eleven
runs of 200 s, one at a time, in alternating order, normalised by CPU time.

| `INT_TEMP_TOL_K` | steps per CPU second | answer |
|---|---|---|
| off (plain `TOLERANCE = 1e-5`) | 2.200 (sd 0.063) | reference |
| 1e-6 | 2.214 (sd 0.102) | **byte-identical** |

The two builds give **byte-identical** output over 390 steps: `max|dTsurf|`
and `max|d res_max|` are exactly zero, not merely small. The speed difference
of 0.6 % is far inside the run-to-run scatter, which reaches 14 % on repeats
of the same binary.

So the tolerance is not what stops that solve. `run/test.c` sets
`NITERMIN = 2`, `mgp_i` is exactly 2 in every run, and two multigrid cycles
already converge this operator. Nine orders of magnitude of tolerance change
nothing because the iteration floor binds first, not the tolerance.

Keep the flag: a case whose solid solve really does stall can use it, and it
costs nothing when it is off. Do not expect it to buy speed.

Caution: measure this kind of change with several alternating runs. A single
pair on a loaded machine gave 43 % here, which was pure scatter.

Caution: with `INT_TEMP_PICARD` the outer loop cannot converge below what the
linear solve delivers. Keep `INT_TEMP_TOL_K` well under
`INT_TEMP_PICARD_TOL`, and re-measure `rel_max` after any change. */

# ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (TS, dt, D=lambda1f, r=sST, theta=theta1);
    diffusion_explicit (TG, dt, D=lambda2f, r=sGT, theta=theta2);
# else
# ifdef INT_TEMP_TOL_K
    double tol_save = TOLERANCE;
    TOLERANCE = rhoS*cpS*INT_TEMP_TOL_K/dt;
# endif
#  if INT_TEMP_ROBIN
    diffusion (TS, dt, D=lambda1f, r=sST, beta=betaST, theta=theta1);
# ifdef INT_TEMP_TOL_K
    TOLERANCE = tol_save;
# endif
#   ifndef TEMPERATURE_PROFILE
    diffusion (TG, dt, D=lambda2f, r=sGT, beta=betaGT, theta=theta2);
#   endif
#  else
    diffusion (TS, dt, D=lambda1f, r=sST, theta=theta1);
# ifdef INT_TEMP_TOL_K
    TOLERANCE = tol_save;
# endif
#   ifndef TEMPERATURE_PROFILE
    diffusion (TG, dt, D=lambda2f, r=sGT, theta=theta2);
#   endif
#  endif
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
