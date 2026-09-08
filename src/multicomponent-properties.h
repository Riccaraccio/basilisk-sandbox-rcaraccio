/**
## Update Properties

Update the thermodynamic properties for the multicomponent phase
change model, and compute the lagrangian derivative of the density,
which is used as a sorce term for the velocity divergence, to
describe low Mach compressibility effects. */

#ifdef VARPROP
#include "solid-thermal-conductivity.h"

scalar rhoGv_G0[], rhoGv_S0[];
extern scalar porosity;
scalar DTDtS[], DTDtG[];
scalar * DYDtG_G = NULL;

/**
With `GAS_SOURCE_EXACT`, `chemistry.h` writes the expansion rate of the
gas-phase reactions here, as `cm[]*ln(rho_start/rho_end)/dt` per unit volume
of gas. `update_divergence()` adds it to `divu2`. The reaction part then does
not pass through `DYDtG_G` and `DTDtG`, and those two fields carry only the
diffusion and the interface terms. */

#if GAS_SOURCE_EXACT
scalar drhodt_chem[];
#endif
scalar * DYDtG_S = NULL;

/**
## The state of a newly uncovered gas cell

`GAS_STATE_FALLBACK` at 1 gives a cell that changes from solid to gas a
usable gas state, so that `update_properties()` can fill its properties. See
the block in the external gas branch below for the mechanism and for the
reason why the repair touches no conserved field.

`GAS_STATE_FALLBACK_FMIN` is the least gas fraction of a donor cell. A donor
below it is too close to the interface to carry a clean gas state.

Set `GAS_STATE_FALLBACK` to 0 to get the previous behaviour, bit for bit.
Define `GAS_STATE_FALLBACK_DEBUG` to write `gasfallback-<pid>.dat`, which
gives every cell that the repair examined and whether the repair succeeded. */

#ifndef GAS_STATE_FALLBACK
# define GAS_STATE_FALLBACK 1
#endif

#ifndef GAS_STATE_FALLBACK_FMIN
# define GAS_STATE_FALLBACK_FMIN 0.5
#endif

trace
void update_properties (void) {

  foreach()
    rhoGv_S0[] = rhoGv_S[]*f[] + (1.-f[])*rhoGv_G[]; //field looks nicer done in one field

  // Reset all the properties fields
  reset ({rhoGv_S, rhoGv_G, rhoSv,
          muGv_S, muGv_G,
          lambdaGv_S, lambdaGv_G, lambdaSv,
          cpGv_S, cpGv_G, cpSv}, 0.);
  reset (DmixGList_S, 0.);
  reset (DmixGList_G, 0.);
  reset ({MWmixG_S, MWmixG_G}, 0.);

  foreach() {
    ThermoState tsGh, tsSh;
    double Diff_coeff[NGS];
    if (f[] > F_ERR && TS[] > 0.) {
      // porosity fraction eps in [0,1]; the raw ratio can undershoot negative
      // in f~F_ERR sliver cells, which makes pow(.,4./3.) a NaN (FE_INVALID).
      double eps_frac = clamp (porosity[]/f[], 0., 1.);
      double xG[NGS], yG[NGS];
      double MWmixG;
      // Update internal gas properties
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList_S[jj];
        yG[jj] = YG[]/f[];
      }
      // empty pore gas: skip the fill; the gas fields stay at their reset-0 value,
      // which downstream consumers guard (MW>0 / denom>0) — do NOT fabricate a rho
      // from zero mole fractions.
      if (mole_from_mass (xG, &MWmixG, yG, NGS)) {
      MWmixG_S[] = MWmixG;

      tsGh.T = TS[]/f[];
      tsGh.P = Pref+p[];
      tsGh.x = xG;

      rhoGv_S[] = tpG.rhov (&tsGh);
      cpGv_S[] = tpG.cpv (&tsGh);
      lambdaGv_S[] = tpG.lambdav (&tsGh);
      muGv_S[] = tpG.muv (&tsGh);
      tpG.diff (&tsGh, Diff_coeff);
#ifdef MASS_DIFFUSION_ENTHALPY
      double cpG[NGS];
      tpG.cpvs (&tsGh, cpG);
      for(int jj=0; jj<NGS; jj++) {
        scalar cpGv = cpGList_S[jj];
        cpGv[] = cpG[jj];
      }
#endif // MASS_DIFFUSION_ENTHALPY

      for(int jj=0; jj<NGS; jj++) {
        scalar DmixGv = DmixGList_S[jj];
        #ifdef CONST_DIFF
        DmixGv[] = CONST_DIFF*pow(eps_frac, 4./3.);
        #else
        DmixGv[] = Diff_coeff[jj]*pow(eps_frac, 4./3.);
        #endif
      }
      } // internal gas filled

      // Update internal solid properties
      double xS[NSS], yS[NSS];
      double MWmixS;
      for (int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        yS[jj] = YS[]/f[];
      }
      // fully-consumed solid: skip; rhoSv/cpSv/lambda1v stay at reset 0
      // (guarded by downstream denom>0 checks).
      if (solid_mole_from_mass (xS, &MWmixS, yS, NSS)) {
      tsSh.T = TS[]/f[];
      tsSh.P = Pref+p[];
      tsSh.x = xS;

      rhoSv[] = tpS.rhov (&tsSh);
      cpSv[] = tpS.cpv (&tsSh);

      coord lambda_pf = pseudo_phase_thermal_conductivity(point, lambdaGv_S[],
                                                          eps_frac,
                                                          tsSh.T,
                                                          f);
      foreach_dimension()
        lambda1v.x[] = lambda_pf.x;
      } // internal solid filled
    }

#ifdef GASGATE_DEBUG
    if (f[] < 1. - F_ERR && !(TG[] > 0.)) {
      static FILE * fpg = NULL;
      static int ng = 0;
      if (ng < 200) {
        if (!fpg) {
          fpg = fopen ("gasgate-dbg.dat", "w");
          fprintf (fpg, "#t x y level f TG TS porosity rhoGv_G rhoGv_S\n");
        }
        fprintf (fpg, "%g %g %g %d %.17g %.17g %.17g %.17g %.17g %.17g\n",
                 t, x, y, level, f[], TG[], TS[], porosity[],
                 rhoGv_G[], rhoGv_S[]);
        fflush (fpg);
        ng++;
      }
    }
#endif

    if (f[] < 1. - F_ERR) {
      // Update external gas properties
      double xG[NGS], yG[NGS];
      double MWmixG;
      double gf = 1. - f[];
      double ytot = 0.;
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList_G[jj];
        yG[jj] = YG[]/gf;
        ytot += yG[jj];
      }
      double TGh = TG[]/gf;

#if GAS_STATE_FALLBACK
      /**
      ## The state of a cell that changes from solid to gas

      `TG` and `YGList_G` are in tracer form, thus `TG[]` is `(1-f)*TG` and
      `YG[]` is `(1-f)*Y`. Both stay near 0 while the cell is solid, which is
      correct. The body shrinks, `f` reaches 0, and the cell becomes gas. The
      recovered state is then `0/1`, which is 0.

      Without a repair the two gates below refuse the fill, every external gas
      property keeps the reset value 0, and `rhomix` in
      `variable-properties.h` becomes 0. The run stops on `1./rhomix`. At
      level 12 the cells are 4 times thinner than at level 10, thus many more
      of them change phase in each second, and the first stop comes at
      t = 0.02 instead of t = 5.94.

      Repair the LOCAL STATE only. Do not write `TG` or `YGList_G`. Those two
      are conserved fields. A write here adds enthalpy and species mass and
      breaks the balance. The properties are derived quantities, thus a
      fallback state changes no balance.

      The order of the fallback is: the solid side of the same cell for the
      temperature, then the neighbour that holds the most gas for whatever is
      still missing.

      Set `GAS_STATE_FALLBACK` to 0 to get the previous behaviour. With the
      flag at 0 this block is empty and the test below is the same test as
      before, because `gf` is larger than 0 here. */

      if (!(TGh > 0.) || !(ytot > 0.)) {
        double wbest = 0., Tdon = 0., ydon[NGS];
        for (int jj=0; jj<NGS; jj++)
          ydon[jj] = 0.;

        foreach_neighbor(1) {
          double gfn = 1. - f[];
          if (gfn > GAS_STATE_FALLBACK_FMIN && gfn > wbest && TG[] > 0.) {
            double ytn = 0.;
            for (int jj=0; jj<NGS; jj++) {
              scalar YG = YGList_G[jj];
              ytn += YG[];
            }
            if (ytn > 0.) {
              wbest = gfn;
              Tdon = TG[]/gfn;
              for (int jj=0; jj<NGS; jj++) {
                scalar YG = YGList_G[jj];
                ydon[jj] = YG[]/gfn;
              }
            }
          }
        }

        if (!(TGh > 0.)) {
          if (f[] > F_ERR && TS[] > 0.)
            TGh = TS[]/f[];
          else if (wbest > 0.)
            TGh = Tdon;
        }

        if (!(ytot > 0.) && wbest > 0.) {
          ytot = 0.;
          for (int jj=0; jj<NGS; jj++) {
            yG[jj] = ydon[jj];
            ytot += yG[jj];
          }
        }

#ifdef GAS_STATE_FALLBACK_DEBUG
        {
          static FILE * fpf = NULL;
          static int nf = 0;
          if (nf < 200) {
            if (!fpf) {
              char nm[80];
              snprintf (nm, sizeof(nm), "gasfallback-%d.dat", pid());
              fpf = fopen (nm, "w");
              fprintf (fpf, "#t x y level f TG TS ytot TGh wbest repaired\n");
            }
            fprintf (fpf, "%g %g %g %d %.17g %.17g %.17g %.17g %.17g %.17g %d\n",
                     t, x, y, level, f[], TG[], TS[], ytot, TGh, wbest,
                     (TGh > 0. && ytot > 0.) ? 1 : 0);
            fflush (fpf);
            nf++;
          }
        }
#endif
      }
#endif // GAS_STATE_FALLBACK

      // empty external gas: skip the fill (fields stay at reset 0, guarded downstream).
      if (TGh > 0. && mole_from_mass (xG, &MWmixG, yG, NGS)) {
      MWmixG_G[] = MWmixG;

      tsGh.T = TGh;
      tsGh.P = Pref+p[];
      tsGh.x = xG;

      rhoGv_G[] = tpG.rhov (&tsGh);
      muGv_G[] = tpG.muv (&tsGh);
      cpGv_G[] = tpG.cpv (&tsGh);
      lambdaGv_G[] = tpG.lambdav (&tsGh);
      tpG.diff (&tsGh, Diff_coeff);

#ifdef MASS_DIFFUSION_ENTHALPY
      double cpG[NGS];
      tpG.cpvs (&tsGh, cpG);
      for(int jj=0; jj<NGS; jj++) {
        scalar cpGv = cpGList_G[jj];
        cpGv[] = cpG[jj];
      }
#endif // MASS_DIFFUSION_ENTHALPY

      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix2v = DmixGList_G[jj];
# ifdef CONST_DIFF
        Dmix2v[] = CONST_DIFF;
# else
        Dmix2v[] = Diff_coeff[jj];
# endif
      }

      foreach_dimension()
        lambda2v.x[] = lambdaGv_G[];
      } // external gas filled
    }
  }
}

//Should  done in the default event but is executed before OS++ initialization otherwise
event init (i = 0) {
  DYDtG_G = NULL;
  DYDtG_S = NULL;

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[64];
    snprintf (name, sizeof(name), "DYDtG_%s_G", OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    a.nodump = true;
    DYDtG_G = list_append (DYDtG_G, a);
  }
  reset (DYDtG_G, 0.);
  
  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[64];
    snprintf (name, sizeof(name), "DYDtG_%s_S", OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    a.nodump = true;
    DYDtG_S = list_append (DYDtG_S, a);
  }
  reset (DYDtG_S, 0.);

#if TREE
  for (scalar s in {drhodt}) {
#if EMBED
    s.refine = refine_embed_linear;
    set_prolongation (s, refine_embed_linear);
#else
    s.refine  = refine_linear;
#endif
    set_restriction (s, restriction_volume_average);
  }
#endif
}

event properties (i = 0) {
  update_properties();
}

event cleanup (t = end)
{
  delete (DYDtG_G), free (DYDtG_G), DYDtG_G = NULL;
  delete (DYDtG_S), free (DYDtG_S), DYDtG_S = NULL;
}

event reset_sources (i++) {
  foreach() {
    DTDtG[] = 0.;
    DTDtS[] = 0.;
#if GAS_SOURCE_EXACT
    drhodt_chem[] = 0.;
#endif
  }

  reset (DYDtG_G, 0.);
  reset (DYDtG_S, 0.);
}

void update_divergence (void) {

//   // ENSURE THAT THE TRACER FORM IS LOST
//   /**
//   We define the variables used to compute the lagrangian derivative
//   on each level. */

  restriction ({T,TS,TG});
  restriction (YSList);
  restriction (YGList_G);
  restriction (YGList_S);
#ifdef MOLAR_DIFFUSION
  restriction (XGList_G);
  restriction (XGList_S);
#endif

//   /**
//   We calculate the Lagrangian derivative of the temperature fields. */

  face vector lambdagradTS[], lambdagradTG[];
  foreach_face() {
    lambdagradTS.x[] = face_value(lambda1v.x, 0)*face_gradient_x (TS, 0)*fm.x[]*fsS.x[];
    lambdagradTG.x[] = face_value(lambda2v.x, 0)*face_gradient_x (TG, 0)*fm.x[]*fsG.x[];
  }

  /**
  The interface heat source. Under `INT_TEMP_VOFBC` the source is split
  between `sST` and the diagonal `betaST`, so the whole term is
  `sST + betaST*TS`. Read only `sST` and the expansion loses the diagonal
  half, which for a sliver is the whole of it, and then the velocity and the
  mass loss rate are wrong with no message.

  `INT_TEMP_ROBIN` escapes this by accident: it adds `+KS*TS[]` to `sST` and
  `-KS` to `betaST`, and the two cancel at `TS = TS^n`. The exact split of
  `INT_TEMP_VOFBC` has no such cancellation, so the term must be written out
  here. */

  foreach() {
    foreach_dimension()
      DTDtS[] += (lambdagradTS.x[1] - lambdagradTS.x[])/Delta;
    DTDtS[] += sST[];

    foreach_dimension()
      DTDtG[] += (lambdagradTG.x[1] - lambdagradTG.x[])/Delta;
    DTDtG[] += sGT[];

#if INT_TEMP_VOFBC
    DTDtS[] += betaST[]*TS[];
    DTDtG[] += betaGT[]*TG[];
#endif
  }

  // EXTERNAL GAS PHASE
  /**
  We calculate the Lagrangian derivative for the chemical species mass
  fractions. */ 

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    scalar DmixGv = DmixGList_G[jj];
    scalar DYDtGjj = DYDtG_G[jj];

    face vector rhoDmixYGjj[];
    foreach_face() {
      double rhoGf = face_value (rhoGv_G, 0);
      double DmixGf = face_value (DmixGv, 0);
      rhoDmixYGjj.x[] = rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
    }

    scalar sgexp = sGexpList[jj];

    foreach() {
      foreach_dimension()
        DYDtGjj[] += (rhoDmixYGjj.x[1] - rhoDmixYGjj.x[])/Delta;
      DYDtGjj[] += sgexp[];
    }
  }

  /**
  We add diffusion correction contributions to the chemical species
  mass fraction derivatives. */

  face vector phicGtot[];
  foreach_face() {
    phicGtot.x[] = 0.;
#ifdef FICK_CORRECTED
    for (int jj=0; jj<NGS; jj++) {
      scalar DmixGv = DmixGList_G[jj];

      double rhoGf = face_value (rhoGv_G, 0);
      double DmixGf = face_value (DmixGv, 0);
# ifdef MOLAR_DIFFUSION
      double MWmixGf = face_value (MWmixG_G, 0);

      scalar XG = XGList_G[jj];
      phicGtot.x[] += (MWmixGf > 0.) ?
        rhoGf*DmixGf*gas_MWs[jj]/MWmixGf*face_gradient_x (XG, 0)*fm.x[]*fsG.x[] : 0.;
# else
      scalar YG = YGList_G[jj];
      phicGtot.x[] += rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
# endif // MOLAR_DIFFUSION
    }
#endif  // FICK_CORRECTED
  }

  for (int jj=0; jj<NGS; jj++) {
    face vector phicGjj[];
    foreach_face() {
      phicGjj.x[] = phicGtot.x[];
#ifdef MOLAR_DIFFUSION
      scalar DmixGv = DmixGList_G[jj];

      double rhoGf = face_value (rhoGv_G, 0);
      double DmixGf = face_value (DmixGv, 0);
      double MWmixGf = face_value (MWmixG_G, 0);

      phicGjj.x[] -= (MWmixGf > 0.) ?
        rhoGf*DmixGf/MWmixGf*face_gradient_x (MWmixG_G, 0)*fm.x[]*fsG.x[] : 0.;
#endif

      scalar YG = YGList_G[jj];
      phicGjj.x[] *= face_value (YG, 0);
    }

    scalar DYDtGjj = DYDtG_G[jj];

    foreach()
      foreach_dimension()
        DYDtGjj[] -= (phicGjj.x[1] - phicGjj.x[])/Delta;
  }

  // INTERNAL GAS PHASE
  /**
  We calculate the Lagrangian derivative for the chemical species mass
  fractions. */ 

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_S[jj];
    scalar DmixGv = DmixGList_S[jj];
    scalar DYDtGjj = DYDtG_S[jj];

    face vector rhoDmixYGjj[];
    foreach_face() {
      double rhoGf = face_value (rhoGv_S, 0);
      double DmixGf = face_value (DmixGv, 0);
      rhoDmixYGjj.x[] = rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsS.x[];
    }

    scalar ssexp = sSexpList[jj];

    foreach() {
      foreach_dimension()
        DYDtGjj[] += (rhoDmixYGjj.x[1] - rhoDmixYGjj.x[])/Delta;
      DYDtGjj[] += ssexp[];
    }
  }

  face vector phicStot[];
  foreach_face() {
    phicStot.x[] = 0.;
#ifdef FICK_CORRECTED
    for (int jj=0; jj<NGS; jj++) {
      scalar DmixGv = DmixGList_S[jj];

      double rhoGf = face_value (rhoGv_S, 0);
      double DmixGf = face_value (DmixGv, 0);
# ifdef MOLAR_DIFFUSION
      double MWmixGf = face_value (MWmixG_S, 0);

      scalar XG = XGList_S[jj];
      phicStot.x[] += (MWmixGf > 0.) ?
        rhoGf*DmixGf*gas_MWs[jj]/MWmixGf*face_gradient_x (XG, 0)*fm.x[]*fsS.x[] : 0.;
# else
      scalar YG = YGList_S[jj];
      phicStot.x[] += rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsS.x[];
# endif // MOLAR_DIFFUSION
    }
#endif  // FICK_CORRECTED
  }

  for (int jj=0; jj<NGS; jj++) {
    face vector phicSjj[];
    foreach_face() {
      phicSjj.x[] = phicStot.x[];
#ifdef MOLAR_DIFFUSION
      scalar DmixGv = DmixGList_S[jj];

      double rhoGf = face_value (rhoGv_S, 0);
      double DmixGf = face_value (DmixGv, 0);
      double MWmixGf = face_value (MWmixG_S, 0);

      phicSjj.x[] -= (MWmixGf > 0.) ?
        rhoGf*DmixGf/MWmixGf*face_gradient_x (MWmixG_S, 0)*fm.x[]*fsS.x[] : 0.;
#endif

      scalar YG = YGList_S[jj];
      phicSjj.x[] *= face_value (YG, 0);
    }

    scalar DYDtGjj = DYDtG_S[jj];

    foreach()
      foreach_dimension()
        DYDtGjj[] -= (phicSjj.x[1] - phicSjj.x[])/Delta;
  }

  // We calculate the one-field divergence by volume-averaging the liquid and the
  // gas-phase contributions.

  foreach() {
    double divu1 = 0., divu2 = 0.;

    // Add internal gas temperature contribution
    divu1 += (TS[]*rhoGv_S[]*cpGv_S[] > 0.) ?
      1./(TS[]*(rhoGv_S[]*cpGv_S[]*porosity[]/f[] + rhoSv[]*cpSv[]*(1-porosity[]/f[])))*DTDtS[] : 0.;

    // Add external gas temperature contribution
    divu2 += (TG[]*rhoGv_G[]*cpGv_G[] > 0.) ?
      1./(TG[]*rhoGv_G[]*cpGv_G[])*DTDtG[] : 0.;

    // Add internal gas chemical species contribution
    double divu1species = 0.;
    for (int jj=0; jj<NGS; jj++) {
      scalar DYDtGjj = DYDtG_S[jj];
      divu1species += 1./gas_MWs[jj]*DYDtGjj[];
    }
    divu1 += (rhoGv_S[] > 0.) ? MWmixG_S[]/rhoGv_S[]*divu1species : 0.;

    // Add external gas chemical species contribution
    double divu2species = 0.;
    for (int jj=0; jj<NGS; jj++) {
      scalar DYDtGjj = DYDtG_G[jj];
      divu2species += 1./gas_MWs[jj]*DYDtGjj[];
    }
    divu2 += (rhoGv_G[] > 0.) ? MWmixG_G[]/rhoGv_G[]*divu2species : 0.;

#if GAS_SOURCE_EXACT
    // Exact step mean of the expansion from the gas-phase reactions
    divu2 += drhodt_chem[];
#endif

    // Volume averaged contributions
    drhodt[] = divu1*f[] + divu2*(1. - f[]);

    // Adjust sign for internal convention
    drhodt[] *= -1.;
  }
}

void update_divergence_density_u (void) {

  scalar rhot[];
  foreach()
    rhot[] = rhoGv_S[]*f[] + rhoGv_G[]*(1. - f[]);
  
  scalar eps[];
  foreach()
    eps[] = porosity[] + (1. - f[]);
  
  scalar DrhoDt[];
  foreach()
    DrhoDt[] = (rhot[] - rhoGv_S0[])*eps[]/dt;
  
  vector grho[];
  gradients ({rhot}, {grho});

  foreach()
    foreach_dimension()
      DrhoDt[] += u.x[]*grho.x[];


  foreach(){
    DrhoDt[] = DrhoDt[]*cm[];

    double one_over_rho = (rhot[] > 0.) ? 1./rhot[] : 0.;

    if (iter > 1) {
      drhodt[] = DrhoDt[]*one_over_rho;
    }
  }
}

void update_divergence_density_uf (void) {
  
  scalar rhot[];
  foreach()
    rhot[] = rhoGv_S[]*f[] + rhoGv_G[]*(1. - f[]);
  
  scalar eps[];
  foreach()
    eps[] = porosity[] + (1. - f[]);
  scalar DrhoDt[];

  foreach()
    DrhoDt[] = (rhot[] - rhoGv_S0[])*eps[]/dt;
  
  face vector rhoGflux[];
  tracer_fluxes (rhot, uf, rhoGflux, dt, zeroc);

  foreach() 
    foreach_dimension()
      DrhoDt[] += (rhoGflux.x[1] - rhoGflux.x[] - rhot[]*(uf.x[1] - uf.x[]))/Delta;

  foreach(){
    DrhoDt[] = DrhoDt[]*cm[];

    double one_over_rho = (rhot[] > 0.) ? 1./rhot[] : 0.;

    if (iter > 1) {
      drhodt[] = DrhoDt[]*one_over_rho;
    }
  }
}

void update_divergence_density (void) {
  // Choose one of the two methods to compute drhodt
  update_divergence_density_u();
  // update_divergence_density_uf();
}
#endif
