/**
# Single-cell check of the two gas-expansion source forms

One cell, one step, no flow. The test does by hand what `chemistry.h` does per
cell: it integrates the gas batch reactor over `dt`, then builds the expansion
source in both forms and prints the two side by side.

    instantaneous   sources[j]  = MW_j*rgas_j            at the end state
                    sources[NGS] = QRgas + divq_rad      at the end state

    averaged        rho*(Y_end,j - Y_start,j)/dt
                    rho*cp*(T_end - T_start)/dt

Both feed `divu2` in `multicomponent-properties.h`:

    divu2 = DTDtG/(TG*rho*cp) + (MWmix/rho)*sum_j DYDtG_j/MW_j
    drhodt = -divu2      (for a pure-gas cell, f = 0)

so `divu2 > 0` is expansion and `drhodt < 0` is expansion.

The two forms must agree in sign. Analytically they are the same quantity:

    sum_j (MW_j*rgas_j)/MW_j = sum_j rgas_j = rho*d(1/MWmix)/dt

so the averaged species term is the exact time average of the instantaneous
one when rho does not change over the step. The test measures whether that
holds for real states and real timesteps.

The averaged form needs a weight for the increment. The exact quantity is

    (1/dt) * integral of rho(tau)*dY_j/dtau dtau

so the weight must represent `rho` over the step. The test measures three
weights: `rho_0` at the step start, which is what `chemistry.h` uses today,
the mean `rho_bar` of the start and end values, and `rho_end`.

## The density form

`update_divergence_density()` in `multicomponent-properties.h` does not build
the source from the reaction rates. It reads the change of the gas density
over the step and divides by the current density:

    density        drhodt = (rho_end - rho_0)/(rho_end*dt)

In one cell with no flow the advection term is zero, so this is the whole
form. The two densities come from the same ideal gas law as above. The column
`dens` reports this form as `update_divergence_density()` would produce it
when both densities belong to the same step.

## The exact expansion

At constant pressure the divergence of a gas parcel is `-d(ln rho)/dt`. Its
exact mean over the step needs no integration of the reactor:

    exact          divu2 = ln(rho_0/rho_end)/dt

The `reference` column is not this quantity. It is the exact mean of the
numerator only, divided by the same denominators that `divu2` in
`multicomponent-properties.h` uses in the default path. The ratio of the two,
in the column `ref/exact`, measures the error of these denominators. The three
density ratios `dens`, `avg` and `exact` differ only in the density that
divides the increment, and the difference is second order in
`1 - rho_end/rho_0`. In a burning cell that number is 0.3, so the choice is
not small.

## The denominators of the default path

The numerator and the denominators of `divu2` come from two time levels:

- `chemistry.h` weights the increment with `rhoGv_G` and `cpGv_G`. The last
  `update_properties()` before the chemistry is the one of the `adapt` event
  of the previous step, so these are the start values `rho_0` and `cp_0`.
- `update_divergence()` runs after the second `update_properties()` of the
  step, in the `tracer_diffusion` event of `multicomponent-varprop.h`. That
  call reads the state after the chemistry. So `TG`, `rhoGv_G`, `cpGv_G` and
  `MWmixG_G` hold the end values `T_end`, `rho_end`, `cp_end` and `MW_end`.
  In the real case `TG` also holds the advection of the step. One cell with
  no flow has no advection, so `T_end` is correct here.

An earlier version of this test divided by the start values `rho_0`, `cp_0`
and `MW_0`. It then showed an expansion under the exact value. The code does
not do that, and with the real denominators the default path gives an
expansion over the exact value. The columns `code/exact`, `heat` and `spec`
show this. `heat` and `spec` compare each part of the default path with the
same part of the exact form:

    exact heat part      ln(T_end/T_0)/dt
    exact species part   ln(MW_0/MW_end)/dt

The sum of the two is `ln(rho_0/rho_end)/dt`, because the pressure does not
change. The columns `1+x/2` and `1+x` give the first-order estimates of the
time-level review, with `x = (T_end - T_0)/T_0`. `GAS_SOURCE_EXACT` in
`chemistry.h` removes both errors, because it uses the exact form. */

#define F_ERR 1e-10
#ifndef NREF
# define NREF 400
#endif
#define SOLVE_TEMPERATURE 1
#define GAS_PHASE_REACTIONS 1

scalar f[];

#include "run.h"

attribute {
  scalar * tracers, c;
  bool inverse;
}

scalar p[];
scalar porosity[];
scalar zeta[];

double rhoG = 0.3;
double rhoS = 1500.;

/**
`variable-properties.h` declares these fields, but it also installs a
`properties` event that overwrites `mu`, `alphav` and `rhov` of the
Navier-Stokes solver. This test has no flow solver, so declare the fields here
and leave the event out. `VARPROP` matters to `reactors.h`: it makes the
reactor recompute `rhog` and `cpg` from the ideal gas law at every evaluation
of the right-hand side, which is what the real case does. */

#define VARPROP
scalar rhoGv_G[], rhoGv_S[], rhoSv[];
scalar muGv_G[], muGv_S[];
scalar lambdaGv_G[], lambdaGv_S[], lambdaSv[];
scalar cpGv_G[], cpGv_S[], cpSv[];

#include "memoryallocation-varprop.h"

event timestep (i++) {
  dtnext (DT);
}

#include "reactors.h"

/**
Build one state, integrate it, and report both source forms. */

/**
Reference: the true time average of the source over the step. Split `dt` into
`n` substeps, integrate each one, and evaluate the reactor at the end of each.
As `n` grows this converges to `(1/dt) * integral of the source dt`, which is
what both forms try to approximate. `sp` receives the species part
`sum_j sources_j/MW_j` and `hr` the heat-release part. */

static void reference_average (const double * ystart, double dtstep, int n,
                               double rho, double cp, double * sp, double * hr)
{
  double y[NGS + 1];
  for (int jj = 0; jj < NGS + 1; jj++)
    y[jj] = ystart[jj];

  double h = dtstep/n, accs = 0., acch = 0.;
  for (int k = 0; k < n; k++) {
    UserDataODE d;
    d.P = Pref; d.T = y[NGS]; d.sources = NULL; d.rhog = rho; d.cpg = cp;
    OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
                         NGS + 1, h, y, &d);

    double src[NGS + 1], dytmp[NGS + 1];
    UserDataODE e;
    e.P = Pref; e.T = y[NGS]; e.sources = src;
    gas_batch_nonisothermal_constantpressure (y, h, dytmp, &e);

    for (int jj = 0; jj < NGS; jj++)
      accs += src[jj]/gas_MWs[jj]*h;
    acch += src[NGS]*h;
  }
  *sp = accs/dtstep;
  *hr = acch/dtstep;
}

static void one_cell (double T0, double yfuel, double yO2, double dtstep,
                      const char * fuel)
{
  int ifuel = OpenSMOKE_IndexOfSpecies (fuel);
  int iO2   = OpenSMOKE_IndexOfSpecies ("O2");
  int iN2   = OpenSMOKE_IndexOfSpecies ("N2");

  double ystart[NGS + 1], yend[NGS + 1];
  for (int jj = 0; jj < NGS; jj++)
    ystart[jj] = 0.;
  ystart[ifuel] = yfuel;
  ystart[iO2]   = yO2;
  ystart[iN2]   = 1. - yfuel - yO2;
  ystart[NGS]   = T0;

  /**
  The gas density and heat capacity of the cell. `chemistry.h` reads these
  from `rhoGv_G`/`cpGv_G`, which the property update sets from the same
  ideal-gas relation, so build them the same way. */

  double x[NGS], MWmix;
  OpenSMOKE_MoleFractions_From_MassFractions (x, &MWmix, ystart);
  OpenSMOKE_GasProp_SetTemperature (T0);
  OpenSMOKE_GasProp_SetPressure (Pref);
  double rho = Pref*MWmix/(R_GAS*1000.*T0);
  double cp  = OpenSMOKE_GasProp_HeatCapacity (x);

  for (int jj = 0; jj < NGS + 1; jj++)
    yend[jj] = ystart[jj];

  UserDataODE data;
  data.P = Pref;
  data.T = ystart[NGS];
  data.sources = NULL;
  data.rhog = rho;
  data.cpg = cp;
  OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
                       NGS + 1, dtstep, yend, &data);

  /**
  The end-state density and heat capacity.

  Caution: do not read these from `data` after the solve. `reactors.h:256`
  starts with `UserDataODE data = *(UserDataODE *)args`, so the reactor works
  on a local copy. The `VARPROP` branch writes `data.rhog` and `data.cpg` in
  that copy only, and the caller never sees the new values. The same applies
  to `bin->rho = data.rhog` in the binning path of `chemistry.h`. So build the
  end state here from the ideal gas law, the same way as the start state. */

  double xe[NGS], MWmix_e;
  OpenSMOKE_MoleFractions_From_MassFractions (xe, &MWmix_e, yend);
  double TGe = yend[NGS];
  OpenSMOKE_GasProp_SetTemperature (TGe);
  OpenSMOKE_GasProp_SetPressure (Pref);
  double rho_end = Pref*MWmix_e/(R_GAS*1000.*TGe);

  /**
  `chemistry.h` builds the same value without a call to the property library,
  as `1/MW = sum_j Y_j/MW_j`. Check that the two agree. */

  double invMW = 0.;
  for (int jj = 0; jj < NGS; jj++)
    invMW += (yend[jj] > 0. ? yend[jj] : 0.)/gas_MWs[jj];
  double rho_end_arith = Pref/(R_GAS*1000.*TGe*invMW);
  if (fabs (rho_end_arith/rho_end - 1.) > 1e-6)
    fprintf (stderr, "# MISMATCH: arithmetic rho_end off by %.3e\n",
             rho_end_arith/rho_end - 1.);
  double cp_end  = OpenSMOKE_GasProp_HeatCapacity (xe);

  double rho_bar = 0.5*(rho + rho_end);

  /**
  Instantaneous form: one more call to the reactor at the converged end
  state, exactly as `gas_sources_instantaneous()` does. The call resets the
  temperature and the pressure of the property library, so read `cp_end`
  before it, as above. */

  double src[NGS + 1], dytmp[NGS + 1];
  UserDataODE d2;
  d2.P = Pref;
  d2.T = yend[NGS];
  d2.sources = src;
  gas_batch_nonisothermal_constantpressure (yend, dtstep, dytmp, &d2);

  double sp_inst = 0., sp_avg = 0., sp_bar = 0., sp_end = 0.;
  for (int jj = 0; jj < NGS; jj++) {
    double dY = (yend[jj] - ystart[jj])/dtstep/gas_MWs[jj];
    sp_inst += src[jj]/gas_MWs[jj];
    sp_avg  += rho*dY;                    // weight frozen at the step start
    sp_bar  += rho_bar*dY;                // trapezoidal weight
    sp_end  += rho_end*dY;                // weight at the step end
  }
  double dTdt_src = src[NGS];                              // W/m3
  double dTdt_avg = rho*cp*(yend[NGS] - ystart[NGS])/dtstep;
  double dTdt_bar = 0.5*(rho*cp + rho_end*cp_end)
                       *(yend[NGS] - ystart[NGS])/dtstep;
  double dTdt_bar0 = rho_bar*cp*(yend[NGS] - ystart[NGS])/dtstep;
  double dTdt_end = rho_end*cp_end*(yend[NGS] - ystart[NGS])/dtstep;

  /**
  The reference average, and a coarser one, so the convergence is visible. */

  double sp_ref, hr_ref, sp_ref2, hr_ref2;
  reference_average (ystart, dtstep, NREF, rho, cp, &sp_ref, &hr_ref);
  reference_average (ystart, dtstep, NREF/4, rho, cp, &sp_ref2, &hr_ref2);

  /**
  Every form shares the same denominators, so the ratio to the reference
  isolates the weight in the numerator, which is the only thing that changes
  between the forms. The denominators mirror `update_divergence()` in the
  default path: `TG`, `rhoGv_G`, `cpGv_G` and `MWmixG_G` hold the values of
  the state after the chemistry, thus `TGe`, `rho_end`, `cp_end` and
  `MWmix_e`. See "The denominators of the default path" above. */

  double TD = TGe, rhoD = rho_end, cpD = cp_end, MWD = MWmix_e;
  double divu2_inst = dTdt_src/(TD*rhoD*cpD) + MWD/rhoD*sp_inst;
  double divu2_avg  = dTdt_avg/(TD*rhoD*cpD) + MWD/rhoD*sp_avg;
  double divu2_bar  = dTdt_bar/(TD*rhoD*cpD) + MWD/rhoD*sp_bar;
  double divu2_bar0 = dTdt_bar0/(TD*rhoD*cpD) + MWD/rhoD*sp_bar;
  double divu2_end  = dTdt_end/(TD*rhoD*cpD) + MWD/rhoD*sp_end;
  double divu2_ref  = hr_ref/(TD*rhoD*cpD)   + MWD/rhoD*sp_ref;
  double divu2_ref2 = hr_ref2/(TD*rhoD*cpD)  + MWD/rhoD*sp_ref2;

  /**
  The two parts of the default path (`avg`), each against the same part of
  the exact form. */

  double heat_avg = dTdt_avg/(TD*rhoD*cpD);
  double spec_avg = MWD/rhoD*sp_avg;
  double heat_exact = log (TGe/T0)/dtstep;
  double spec_exact = log (MWmix/MWmix_e)/dtstep;
  double xT = (TGe - T0)/T0;

  /**
  The density form and the exact expansion. Both use only the two densities.
  The density form is `-DrhoDt/rhot` with the same-step fields, and `rhot`
  is the density at the end of the step. */

  double divu2_dens  = (rho/rho_end - 1.)/dtstep;
  double divu2_exact = log (rho/rho_end)/dtstep;

  fprintf (stderr,
           "%6.0f %6.3f %9.1e | %11.3e | %8.3f %7.4f %7.4f %7.4f %7.4f |"
           " %7.4f %7.4f %9.4f | %6.1f %7.4f %8.4f |"
           " %7.4f %7.4f %7.4f %7.4f %7.4f\n",
           T0, yfuel, dtstep,
           -divu2_ref,
           divu2_inst/divu2_ref, divu2_avg/divu2_ref,
           divu2_bar/divu2_ref, divu2_bar0/divu2_ref, divu2_end/divu2_ref,
           divu2_dens/divu2_ref, divu2_exact/divu2_ref,
           divu2_ref/divu2_exact,
           yend[NGS] - ystart[NGS], rho_end/rho,
           divu2_ref2/divu2_ref,
           divu2_avg/divu2_exact,
           heat_exact != 0. ? heat_avg/heat_exact : 0., 1. + xT/2.,
           spec_exact != 0. ? spec_avg/spec_exact : 0., 1. + xT);
}

int main() {
  TG0 = 1123.; TS0 = 300.;
  cpS = 1500.; cpG = 1100.;
  kinfolder = "biomass/dummy-solid-gas";
  DT = 1e-4;
  init_grid (1 << 2);
  run();
}

event init (i = 0) {
  restarted = true;

  /**
  `chemistry.h` does this in its own `init` event. This test does not include
  `chemistry.h`, so it must do it here, or `OpenSMOKE_ODESolver` runs with no
  workspace and segfaults. */

  OpenSMOKE_InitODESolver ();
  foreach() {
    f[] = 0.; porosity[] = 0.; zeta[] = 0.; p[] = 0.;
    for (int jj = 0; jj < NGS; jj++) {
      scalar YG = YGList_G[jj]; YG[] = (jj == OpenSMOKE_IndexOfSpecies ("N2"));
      scalar YGs = YGList_S[jj]; YGs[] = 0.;
      scalar YGi = YGList_Int[jj]; YGi[] = 0.;
      scalar sSexp = sSexpList[jj]; sSexp[] = 0.;
      scalar sGexp = sGexpList[jj]; sGexp[] = 0.;
    }
    for (int jj = 0; jj < NSS; jj++) { scalar YS = YSList[jj]; YS[] = 0.; }
    TG[] = TG0; TS[] = 0.; T[] = TG0; TInt[] = 0.;
  }

  fprintf (stderr, "# gas species in %s:", kinfolder);
  for (int jj = 0; jj < NGS; jj++)
    fprintf (stderr, " %s", OpenSMOKE_NamesOfSpecies (jj));
  fprintf (stderr, "\n\n");

  fprintf (stderr,
           "#   weight:  (none)  rho_0  mean(rho.cp)  rho_bar+cp_0"
           "  rho_end | density form, exact log form\n");
  fprintf (stderr,
           "#    T0  yfuel        dt |  drhodt_ref |      ratio to the"
           " reference       |  ratio to ref     ref/exact |"
           "     dT rhoE/rho ref conv | code/exact    heat   1+x/2"
           "    spec     1+x\n");
  fprintf (stderr,
           "#                        |    = -divu2 |     inst     avg     bar"
           "    bar0     end |    dens   exact           |"
           "                   |  default path against the exact form\n");

  const char * fuel = "CO";
  for (int k = 0; k < NGS; k++)
    if (!strcmp (OpenSMOKE_NamesOfSpecies (k), "CH4")) fuel = "CH4";

  double dts[] = {2e-4, 5e-5, 1e-5, 2e-6};
  for (int a = 0; a < 4; a++)
    one_cell (1900., 0.08, 0.10, dts[a], fuel);
  fprintf (stderr, "\n");
  double Ts[] = {1200., 1500., 1800., 1950., 2100.};
  for (int a = 0; a < 5; a++)
    one_cell (Ts[a], 0.08, 0.10, 2e-4, fuel);
  fprintf (stderr, "\n");
  double yf[] = {0.02, 0.05, 0.10, 0.20};
  for (int a = 0; a < 4; a++)
    one_cell (1900., yf[a], 0.15, 2e-4, fuel);

  exit (0);
}
