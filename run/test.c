/**
# The amplitude ladder: from the reduced case to `fatehi-combustion`

The slow oscillation of the release rate of the particle has the same
frequency in the reduced case and in the full case, but not the same
amplitude. At matched remaining mass the full case gives 3.5 to 4.1 per
cent and the reduced case gives 1.6 to 2.3 per cent. The probe temperature
differs by more: 23 K against 4.4 K. So the particle oscillates twice as
much, and the flame turns that oscillation into 2.5 times more temperature.

This case measures which ingredient carries the difference. The build with
no flags reproduces the reduced case `~/temp/expansion/test.c`. Each flag
adds one ingredient of `run/fatehi-combustion.c`.

| flag | reduced value | full value |
|---|---|---|
| `DA_VALUE` | 1e-12 | 1e-10 |
| `MOISTURE` | 0, pure biomass | 1, 6.1 % moisture and 0.4 % ash |
| `GRAVITY` | 0 | 1 |
| `SHAPE` | 0, sphere | 1, superquadric |
| `EMISSIVITY_DIBLASI` | 0, constant | 1, Di Blasi |
| the three transport flags | off | on |

`MAXLEVEL_VALUE` and `DT_VALUE` are separate: they test the grid and the
step, not an ingredient.

Compare the runs at matched remaining mass, never at matched time. The
cases burn at different rates, so the same instant is a different state.

Caution: `MOLAR_DIFFUSION`, `FICK_CORRECTED` and `MASS_DIFFUSION_ENTHALPY`
are existence-tested. Every site that this case reaches is `#ifdef`, 45 of
them in `src/`. A `#define MOLAR_DIFFUSION 0` therefore turns the flag ON.
The block below undefines them instead, so `-DMOLAR_DIFFUSION=0` means off,
and both `-DMOLAR_DIFFUSION` and `-DMOLAR_DIFFUSION=1` mean on. Do not
replace it with a `#ifndef ... 0` default. */

#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1

#if defined(MOLAR_DIFFUSION) && !MOLAR_DIFFUSION
# undef MOLAR_DIFFUSION
#endif

#if defined(FICK_CORRECTED) && !FICK_CORRECTED
# undef FICK_CORRECTED
#endif

#if defined(MASS_DIFFUSION_ENTHALPY) && !MASS_DIFFUSION_ENTHALPY
# undef MASS_DIFFUSION_ENTHALPY
#endif

/**
`TURN_OFF_HEAT_OF_REACTION` carries the same trap: `reactors.h` reads it with
`#ifdef` at lines 237 and 303, so `-DTURN_OFF_HEAT_OF_REACTION=0` would turn
it ON. Undefine it when the value is 0.

The flag keeps the kinetics of the species and zeroes the temperature source.
It therefore removes the endothermic arm of the loop and nothing else, which
is what separates the endothermic closure from the blowing closure. */

#if defined(TURN_OFF_HEAT_OF_REACTION) && !TURN_OFF_HEAT_OF_REACTION
# undef TURN_OFF_HEAT_OF_REACTION
#endif

/**
## Pure heating

`TURN_OFF_REACTIONS` removes every reaction and leaves the particle to heat
up. `chemistry.h` reads it at line 29 with `#ifndef`, so the flag compiles out
the whole `chemistry` event: no solid kinetics, no gas kinetics, and no stiff
solve. The run keeps the flow, the transport of the species, the two
temperature equations, the interface coupling and `darcy.h`.

The particle therefore keeps its mass and its composition. `omega` stays at
zero, so `gas_source` stays at zero, `zeta` stays at zero and the particle
neither shrinks nor blows. Columns 3, 6, 7 and 17 of `expansion.dat` must hold
0, and column 2 of `OutputData` must stay constant. Check those columns on the
first run: they are the proof that the reactions are off.

Caution: column 2 of `OutputData` does not hold exactly 1. The first
adaptation changes the discrete volume of the solid by about 0.6 per cent, so
the value falls to 0.9942 at the first output and holds there. Read the
constancy, not the value.

Use this build as the thermal reference of the ladder. It gives the transient
of pure heat conduction, so it separates the response of the temperature from
the loop of the reactions.

Caution: the flag carries the trap of `TURN_OFF_HEAT_OF_REACTION` in the
opposite direction. The test is `#ifndef`, so `-DTURN_OFF_REACTIONS=0` would
also turn the reactions OFF. The block below undefines it when the value is 0,
so `-DTURN_OFF_REACTIONS=0` means on, and both `-DTURN_OFF_REACTIONS` and
`-DTURN_OFF_REACTIONS=1` mean off.

Caution: no reaction makes H2O, so `T_H2O_weigthed_average` finds an empty
weight and returns `TG0` at each of its points. Columns 4, 5 and 6 of
`OutputData` and every column of `TemperatureProfile.dat` are constant in this
build. Read columns 10 to 12 of `OutputData`, which weight the same three
points uniformly and need no absorber, and `Tcore`, `Tbulk` and `Tsurf`,
columns 14 to 16 of `expansion.dat`. */

#if defined(TURN_OFF_REACTIONS) && !TURN_OFF_REACTIONS
# undef TURN_OFF_REACTIONS
#endif

/**
## No expansion of the gas

`NO_EXPANSION` removes the expansion of the gas from the divergence. Two
sites read it, and both use `#ifndef`:

`centered-phasechange.h:120` builds the source of the projection as
`div_source = gas_source + drhodt`. With the flag it keeps `gas_source`
alone. `multicomponent-varprop.h:1071` then does not call
`update_divergence()`, so `drhodt` is never written and it holds the 0 that
the `reset_sources` event gives it.

The particle still releases its mass, because `gas_source` carries the phase
change and it does not pass through `drhodt`. What goes away is the volume
which the gas gains when it heats up and when its composition changes. The
build therefore separates the blowing of the phase change from the expansion
of the gas around it.

`test-constprop` also has no expansion, but for another reason: it has no
`VARPROP`, so nothing writes `drhodt` at all. That build changes every
property as well. This one keeps the variable properties everywhere else, so
it is the clean test of the expansion alone.

Caution: the test is `#ifndef`, so `-DNO_EXPANSION=0` would also turn the
expansion OFF. The block below undefines the flag when the value is 0, so
`-DNO_EXPANSION=0` means the expansion is on, and both `-DNO_EXPANSION` and
`-DNO_EXPANSION=1` mean it is off.

Caution: the flag must reach `centered-phasechange.h`, so it stays above the
includes. */

#if defined(NO_EXPANSION) && !NO_EXPANSION
# undef NO_EXPANSION
#endif

/**
The log needs a value, and the six flags above carry none once they are
undefined. */

#ifdef NO_EXPANSION
# define NOEXP_ON 1
#else
# define NOEXP_ON 0
#endif

/**
The coupling of the interface. `INT_TEMP_VOFBC` and `INT_TEMP_PICARD` are
value-tested in `multicomponent-varprop.h`, so an undefined name reads as 0
and the two blocks below need no `defined()`.

The log must carry them. A build takes them from the Makefile rung only, and
`test-%.c: test.c` makes the link for ANY name, so a rung that the local
`Makefile` does not know still compiles and runs, with the flags of the
pattern rule alone. The header line is what tells the two apart. */

#if INT_TEMP_VOFBC
# define VOFBC_ON 1
#else
# define VOFBC_ON 0
#endif

#if INT_TEMP_PICARD
# define PICARD_ON 1
#else
# define PICARD_ON 0
#endif

#ifdef TURN_OFF_REACTIONS
# define NOREACT_ON 1
#else
# define NOREACT_ON 0
#endif

#ifdef TURN_OFF_HEAT_OF_REACTION
# define NOHEAT_ON 1
#else
# define NOHEAT_ON 0
#endif

#ifdef MOLAR_DIFFUSION
# define MOLAR_ON 1
#else
# define MOLAR_ON 0
#endif

#ifdef FICK_CORRECTED
# define FICK_ON 1
#else
# define FICK_ON 0
#endif

#ifdef MASS_DIFFUSION_ENTHALPY
# define MDE_ON 1
#else
# define MDE_ON 0
#endif

/**
The permeability of the reduced case. The full case leaves the default of
`darcy.h`, which is 1e-10, so the full particle is 100 times more
permeable. */

#ifndef DA_VALUE
# define DA_VALUE 1e-12
#endif

#ifndef MOISTURE
# define MOISTURE 0
#endif

#ifndef EMISSIVITY_DIBLASI
# define EMISSIVITY_DIBLASI 0
#endif

#ifndef GRAVITY
# define GRAVITY 0
#endif

#ifndef SHAPE
# define SHAPE 0 // 0 sphere, 1 superquadric
#endif

#ifndef MAXLEVEL_VALUE
# define MAXLEVEL_VALUE 10
#endif

/**
The two wavelet thresholds of the `adapt` event. They are separate from
`MAXLEVEL_VALUE` on purpose.

`maxlevel` sets how fine the mesh may become. These set how much of the
domain reaches it. A sweep of the thresholds at a FIXED `maxlevel` therefore
separates two errors that the level ladder confuses: the error of the
discretisation at the finest cell, and the error of the criterion that
decides which cells get one.

The 0.15 to 1.5 Hz band of `mdot` falls by about 4 per level, and it does not
answer to the timestep at all. If it also falls with these thresholds at
fixed `maxlevel`, the criterion carries it. If it does not, the finest cell
does, and only `maxlevel` will help. */

#ifndef ADAPT_T_TOL
# define ADAPT_T_TOL 5e0
#endif

#ifndef ADAPT_O_TOL
# define ADAPT_O_TOL 1.e-2
#endif

/**
The reduced reference uses 5e-4. Keep this value, or the unflagged build
does not reproduce `~/temp/expansion/avg`. */

#ifndef DT_VALUE
# define DT_VALUE 5e-4
#endif

/**
## The timestep pair

`Tmax` is a single-valued function of `dt`: 16 to 29 K for each halving, with
R^2 = 0.89 in `test-shape`. `dt` itself is quantised, because `dtnext()` snaps
the step so that the run lands on the four `t += 0.01` events, so `dt` can
only take the values 0.01/n.

`CFLNUM` makes a constant step possible. With `CFLNUM = 2` the CFL never
binds, so `DT_VALUE` sets every step, and it still caps a runaway.

Caution: `DT_VALUE` must divide 0.01 exactly, or `dtnext()` subdivides the
step and `dt` is not constant. Use 2e-4 = 0.01/50, not 1.75e-4.

Caution: `CFL` is assigned in the `defaults` event of
`navier-stokes/centered.h`, which runs after `main()`. So this value is set in
`event init` below, never in `main()`. `TOLERANCE` has no such event and stays
in `main()`.

Caution: at ignition the peak velocity reaches about 2 m/s and the CFL binds
even at `CFLNUM = 2`. Branch the fixed-step runs from a plateau snapshot, and
check that column 2 of `expansion.dat` is constant before you quote them. */

#ifndef CFLNUM
# define CFLNUM 0.8
#endif

/**
## The blowing sweep

In the pyrolysis case the Stefan flow is not a perturbation. The blowing
ratio against the free stream is about 1.06, the Peclet number `v_w R/alpha`
is 2.5 to 4.4, and the blockage of the conductive flux moves by a factor of 4
over one cycle. `UIN_VALUE` changes the ratio, so it tests whether the loop
closes through the blowing.

`v_w` follows the rate of pyrolysis, not `Uin`, so `Uin` changes the ratio and
not the Peclet number. A larger `Uin` thins the layer and weakens the
blockage; a smaller `Uin` strengthens it.

Caution: `Uin` also changes the supply of heat, so it changes the rate of
burn. Compare at matched remaining mass, never at matched time. */

#ifndef UIN_VALUE
# define UIN_VALUE 0.13
#endif

/**
## No inflow

`NO_INFLOW` removes the inlet. The particle then sits in a quiescent gas at
`TG0`: the domain is `10*D0` with its origin at the centre of the particle,
the gas temperature and composition are fixed at the right and top
boundaries, and gravity is off.

The switch is an integer on purpose. The preprocessor rejects a floating
constant in `#if`, so a test such as `#if !NO_INFLOW` stops the build of
every rung. Caution: do not test `UIN_VALUE` in `#if`. Test `NO_INFLOW`.

A build with the switch on sets `UIN_VALUE` to 0, so the log line reports the
velocity that the case really has. */

#ifndef NO_INFLOW
# define NO_INFLOW 0
#endif

#if NO_INFLOW
# undef UIN_VALUE
# define UIN_VALUE 0.
#endif

/**
## Pyrolysis only

`PYROLYSIS_ONLY` removes the oxidiser. It does not remove the reactions of
the gas: the Makefile adds `TURN_OFF_GAS_REACTIONS=1` to every
pyrolysis-only rung for that. Before that flag existed, the rung reproduced
`~/temp/10-fatehi/test.c`. The release still oscillates, 7 to 18
per cent peak-to-peak, so the flame does not close the loop. The case is much
cheaper and gives 6 cycles in 20 s, against 1.2 to 4.3 cycles in the runs
with a flame, so it is the correct case for the questions about the loop.

Both branches use `biomass/dummy-solid-gas`; see the note at `kinfolder`
below. The mechanism therefore holds O2 in either branch, but the
`PYROLYSIS_ONLY` build sets its mass fraction to zero everywhere. Keep the
lookups of O2 inside the `#else` branch anyway: they record which build wants
an oxidiser. */

#ifndef PYROLYSIS_ONLY
# define PYROLYSIS_ONLY 0
#endif

/**
## The shrinking of the particle

`ZETA_POLICY` selects the policy of `shrinking.h`, which splits the volume of
the phase change between the shrinkage of the solid and the release of gas.
The ladder uses `ZETA_REACTION`, which ties the shrinkage to the local
reaction rate. `ZETA_SWELLING` sets `zeta` to zero everywhere, so the
interface stays where it is: the solid keeps its volume, and the decomposition
raises the porosity instead.

`zeta` enters through one term only, the source of the Poisson equation of the
velocity potential, `prod = omega*f*zeta*cm/rhoS`
(`velocity-potential.h:47`). With `zeta = 0` that source is zero, the solid
velocity `ubf` is zero, and nothing advects the interface. `gas_source` does
not carry `zeta`, so the particle releases the same mass into the gas in
either policy. The two builds therefore separate the motion of the interface
from the release of mass.

A build with no reaction gives the same run under either policy. `omega` is
zero everywhere, so `ZETA_REACTION` finds no maximum and sets `zeta` to zero,
which is what `ZETA_SWELLING` sets. Use the pair as a control of the ladder:
the two runs must agree in every column.

`shrinking.h` declares the names, so this macro expands only in `main()`,
after the includes. */

#ifndef ZETA_POLICY
# define ZETA_POLICY ZETA_REACTION
#endif

#define ZETA_STR_(x) #x
#define ZETA_STR(x) ZETA_STR_(x)

/**
## A steady reaction rate

`OMEGA_CONST` overwrites the field `omega` after the `chemistry` event and
before `shrinking.h` reads it. The release of the particle then has no
oscillation, and the run answers one question: does the temperature still
oscillate when the release does not?

Read `src/CLAUDE.md` first. `omega` is the rate per cubic metre of solid
material, and it is the ONLY path by which the decomposition of the solid
reaches the divergence of the velocity. Three consumers read it, all of them
in the `phasechange` event of `shrinking.h`:

  `gas_source = -omega*(f - porosity)*(1/rhoG - 1/rhoS)`   the blowing
  `prod = omega*f*zeta*cm/rhoS`   the source of `psi`, hence `ubf`
  `set_zeta()`   under `ZETA_REACTION`, `zeta = omega/max(omega*f)`

The flag takes two modes.

`OMEGA_CONST 1` gives every cell of the solid the same value,
`OMEGA_CONST_VALUE`. The release is then uniform over the particle.

`OMEGA_CONST 2` multiplies the whole field by one number, so that the
integral `mdot` takes the value `OMEGA_CONST_MDOT` at every step. Prefer
this mode. It keeps the shape of the reaction front, and `zeta` does not
change at all, because the factor cancels in the ratio `omega/max(omega*f)`.
Mode 1 makes `omega` uniform, so `ZETA_REACTION` returns 1 over the whole
particle, and the solid shrinks as fast as it can. That is a second change
of the model, and it confuses the answer.

Caution: this build does not conserve mass, and it is a diagnostic build
only. The ODE of the chemistry has already run when this event fires, so the
solid loses its mass at the true rate while the gas receives the overwritten
rate. Mode 2 breaks the balance by the few per cent of the oscillation.
Mode 1 breaks it by much more. Never quote column 2 of `OutputData` from
this build.

Caution: the flag does not stop the release of the species. The pore gas
takes its products from `dy[]` of the reactor (`reactors.h:248`), not from
`omega`, so `YGList_S` still carries the true kinetics. The build therefore
separates the blowing from the composition. It does not remove every path
by which the reactions reach the gas.

`OMEGA_CONST_T0` holds the override back until the particle ignites. The
case ignites at about t = 6 s, so the default starts the override on the
plateau. The step of `mdot` at `T0` makes a transient of about 2 s. Discard
it before you read an amplitude.

## Where the two numbers come from

Both defaults are means of `~/test/new/test-base`. `OMEGA_CONST_MDOT` is the
mean of column 17 over 10 to 30 s. `OMEGA_CONST_VALUE` is the mass-weighted
mean of `omega` over the same window, which is `mdot*rhoS/Ms`, because
`omega` is the rate per cubic metre of solid material.

A window mean is only as good as the window, so the two quantities were
tested for a trend. The result decides which mode to run:

| window | mdot trend | mdot oscillation | omega trend | omega oscillation |
|---|---|---|---|---|
| 12-28 s | -4.9 % | 21.6 % | +76.5 % | 11.3 % |
| 15-25 s | -0.4 % | 13.6 % | +51.7 % | 8.6 % |
| 18-22 s | -1.0 % | 3.9 % | +20.5 % | 4.3 % |

`mdot` holds still on the plateau. The rate per unit mass does not: it rises
by 76 per cent over 12 to 28 s, because the solid which remains gets hotter
and each kilogram reacts faster. Mode 2 therefore fixes the quantity which
the case already holds nearly fixed, and it removes the oscillation and
little else. Mode 1 fixes a quantity which doubles over the run, so it is
about 40 per cent too high at the start of the window and 40 per cent too
low at the end.

Caution: read mode 2 over 12 to 28 s. The default `T0` gives a transient of
about 2 s, and the drift over that window is 5 per cent. Read mode 1 over 18
to 22 s only.

The units of `Ms0` were verified against a run: mode 1 at 83 gave
`mdot = 2.83047e-6` against the prediction `83*Ms0/rhoS = 2.83307e-6`. */

#ifndef OMEGA_CONST
# define OMEGA_CONST 0         // 0 off, 1 uniform value, 2 fixed mdot
#endif

#ifndef OMEGA_CONST_VALUE
# define OMEGA_CONST_VALUE 83. // mode 1 [kg/m3/s]
#endif

#ifndef OMEGA_CONST_MDOT
# define OMEGA_CONST_MDOT 1.53e-6 // mode 2, the target of column 17
#endif

#ifndef OMEGA_CONST_T0
# define OMEGA_CONST_T0 10.
#endif

/**
## Constant properties

`CONST_PROPERTIES` replaces `opensmoke-properties.h` with
`constant-properties.h`. The gas then carries one density, one viscosity, one
conductivity and one heat capacity over the whole domain, and the solid
carries one conductivity and one heat capacity. `run/POM.c` builds the same
stack this way.

The swap removes `VARPROP`, because `variable-properties.h` is what defines
it. That is not a small change of the model. Three consequences:

The gas no longer expands with the temperature. The density is a constant, so
its material derivative is zero, and the divergence carries only the source of
the phase change. The pair of builds therefore measures what the thermal
expansion of the gas contributes.

`multicomponent-properties.h` lies inside `#ifdef VARPROP`, so
`update_properties()` and `solid-thermal-conductivity.h` are both absent.
`constant-properties.h` writes `lambda1v` and `lambda2v` itself, from
`lambdaS`, `lambdaG` and the porosity. `lambdaSmodel` does not exist in this
build, so `main()` sets it in the other branch only.

Caution: `constant-properties.h` fixes the diffusion coefficient of every
species at 2.05e-5 m^2/s, the value for CO in N2 at 500 K, and the case cannot
change it. At 1123 K the varprop build gives about ten times more. Do not read
a difference in the transport of the species as a result of the constant
properties alone.

The six values below are the ones of the commented block in `main()`, for air
at about 1100 K. Each one takes a `-D` of its own. */

#ifndef CONST_PROPERTIES
# define CONST_PROPERTIES 0
#endif

#ifndef RHOG_CONST
# define RHOG_CONST 0.31     // air at 1100 K, 1 atm [kg/m3]
#endif

#ifndef MUG_CONST
# define MUG_CONST 4.5e-5    // [Pa s]
#endif

#ifndef LAMBDAG_CONST
# define LAMBDAG_CONST 0.08  // [W/m/K]
#endif

#ifndef CPG_CONST
# define CPG_CONST 1200.     // [J/kg/K]
#endif

#ifndef LAMBDAS_CONST
# define LAMBDAS_CONST 0.2   // [W/m/K]
#endif

#ifndef CPS_CONST
# define CPS_CONST 1500.     // [J/kg/K]
#endif

/**
## The four probes of the time level

Each probe answers one question of
`~/discretization-report/time-level-review.md` and of
`~/discretization-report/formulation-coherence.md`. Each one is diagnostic
only: with its flag at 0 the run is the same, bit for bit, as the run of the
code without it. No probe writes a field, a boundary condition or `dt`, and
no probe changes the order of an event that writes a field.

| flag | header | file | question |
|---|---|---|---|
| `CHEM_SPLIT_PROBE` | `chem-split-probe.h` | `chemsplit.dat` | TL-1, the split of the gas chemistry from the transport |
| `DRHODT_BUDGET` | `drhodt-budget.h` | `drhodtbudget.dat` | TL-2, the explicit fluxes of `drhodt` |
| `SHRINK_BUDGET` | `shrink-budget.h` | `shrinkbudget.dat` | NEW-1, NEW-4 and 9d, the shrinkage that the VOF sweep removes |
| `SPECIES_CLAMP_PROBE` | `species-clamp-probe.h` | `speciesclamp.dat` | TL-5, the clamp of the species after the implicit solves |

The four flags are ON in this case, because the next base run must carry
them. Build with `-DCHEM_SPLIT_PROBE=0`, `-DDRHODT_BUDGET=0`,
`-DSHRINK_BUDGET=0` or `-DSPECIES_CLAMP_PROBE=0` to turn one off. Each header holds the column layout and
the way to read it.

Cost, measured at level 8 over 0.5 s (1012 steps): the first three probes use
0.14 s of CPU together, 0.2 per cent of the run. `CHEM_SPLIT_PROBE` and
`DRHODT_BUDGET` work only on the step that starts at an output time
(`t += 0.01`). `SHRINK_BUDGET` works at every step, because its running
integrals need every step, and writes at the output times. Each probe
prints its own CPU time at the end of the log.

Cost of `SPECIES_CLAMP_PROBE`, measured on the same case: 0.018 s of CPU in
50 steps, 0.02 per cent of the run. It works only on the step that starts at
an output time, and it adds no field.

Caution: `DRHODT_BUDGET` and `SPECIES_CLAMP_PROBE` must be set before
`multicomponent-varprop.h`, which includes their headers and calls their
functions. So the flags stay in this block, above the includes. */

#ifndef CHEM_SPLIT_PROBE
# define CHEM_SPLIT_PROBE 1
#endif

#ifndef DRHODT_BUDGET
# define DRHODT_BUDGET 1
#endif

#ifndef SHRINK_BUDGET
# define SHRINK_BUDGET 1
#endif

#ifndef SPECIES_CLAMP_PROBE
# define SPECIES_CLAMP_PROBE 1
#endif

/**
## The end of the run

`TEND_VALUE` sets the end time. The ladder uses the default. A smoke test
uses a small value, for example `-DTEND_VALUE=0.5`.

Caution: keep the end time above the last time of every other event that
carries an upper limit, or read the note of `event stop` below. */

#ifndef TEND_VALUE
# define TEND_VALUE 40.
#endif

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"

/**
`multicomponent-varprop.h` defines `MULTICOMPONENT`, but it comes after
`constant-properties.h`, which reads the flag to declare the fields of the
diffusion coefficient. So this build defines it here, as `run/POM.c` does.

`opensmoke.h` comes with `opensmoke-properties.h` in the other branch. Without
it `common-phasechange.h` calls `OpenSMOKE_IndexOfSolidSpecies` before it sees
the declaration, and the build gives two implicit-declaration warnings. Include
it here. `qcc` includes each header once, so the later include in
`memoryallocation-varprop.h` does nothing. */

#if CONST_PROPERTIES
# define MULTICOMPONENT
# include "opensmoke.h"
# include "constant-properties.h"
#else
# include "opensmoke-properties.h"
#endif

#include "two-phase.h"

/**
`fatehi-combustion.c` puts `gravity.h` here, between `two-phase.h` and
`shrinking.h`. Same-name events run in reverse declaration order, so the
position of a module in this list changes when its events run. Keep the
order of the full case. */

#if GRAVITY && !NO_INFLOW
# include "gravity.h"
#endif

#include "superquadric.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
#include "divergence-budget.h"
#include "probe-cache.h"
#include "chem-split-probe.h"
#include "shrink-budget.h"

#include "view.h"

#if !NO_INFLOW

const double Uin = UIN_VALUE; //inlet velocity
/**
Caution: give `pf` the same conditions as `p`. Basilisk does not copy them,
so without these lines every boundary of `pf` is Neumann, and the projection
of `uf` has no solution. See the comment in `fatehi-combustion.c`. */

u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
pf[left]     = neumann (0.);
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);
psi[right]    = neumann (0.);

#else // NO_INFLOW

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);
psi[right]    = dirichlet (0.);

u.n[top]    = neumann (0.);
u.t[top]    = neumann (0.);
p[top]      = dirichlet (0.);
pf[top]     = dirichlet (0.);
psi[top]    = dirichlet (0.);

#endif

double tend = TEND_VALUE;
int maxlevel = MAXLEVEL_VALUE, minlevel = 2;
double solid_mass0 = 0.;
double D0 = 8e-3, H0 = 8e-3;

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

int main() {

  /**
  Caution: under MPI every rank shares this stderr. Guard the message with
  `pid() == 0`, or the log carries one copy per rank. */

  if (pid() == 0)
    fprintf (stderr, "# ladder: MOLAR=%d FICK=%d MDE=%d MOISTURE=%d GRAVITY=%d"
                     " SHAPE=%d DIBLASI=%d Da=%g DT=%g maxlevel=%d"
                     " CFL=%g Uin=%g PYRO=%d NOHEAT=%d NOREACT=%d"
                     " ZETA=%s CONSTP=%d NOEXP=%d VOFBC=%d PICARD=%d"
                     " TSADV=%d NOGASR=%d NOINFLOW=%d"
                     " OMEGACONST=%d OCVAL=%g OCMDOT=%g OCT0=%g nranks=%d\n",
             MOLAR_ON, FICK_ON, MDE_ON, MOISTURE, GRAVITY, SHAPE,
             EMISSIVITY_DIBLASI, (double) DA_VALUE, (double) DT_VALUE,
             MAXLEVEL_VALUE, (double) CFLNUM, (double) UIN_VALUE,
             PYROLYSIS_ONLY, NOHEAT_ON, NOREACT_ON,
             ZETA_STR(ZETA_POLICY), CONST_PROPERTIES, NOEXP_ON,
             VOFBC_ON, PICARD_ON, TS_PORE_ADVECTION,
             TURN_OFF_GAS_REACTIONS, NO_INFLOW,
             OMEGA_CONST, (double) OMEGA_CONST_VALUE,
             (double) OMEGA_CONST_MDOT, (double) OMEGA_CONST_T0, npe());

  /**
  `lambdaSmodel` comes with `solid-thermal-conductivity.h`, which
  `multicomponent-properties.h` includes inside `#ifdef VARPROP`. The constant
  build has neither, and `constant-properties.h` writes the conductivity of
  the two pseudo-phases itself. */

#if CONST_PROPERTIES
  rhoG    = RHOG_CONST;
  muG     = MUG_CONST;
  lambdaG = LAMBDAG_CONST;  cpG = CPG_CONST;
  lambdaS = LAMBDAS_CONST;
#else
  lambdaSmodel = L_TENWOLDE;
#endif

  TS0 = 300.; TG0 = 1123.;
  rhoS = 1550; cpS = 1800;
  eps0 = 0.2;

#if CONST_PROPERTIES
  cpS = CPS_CONST;
#endif

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_POLICY;

  DT = DT_VALUE;

  /**
  One mechanism for both branches.

  `biomass/dummy-solid` declares only `BIOMASS CHAR` as solid species. It has
  no `MOIST` and no `ASH`. `lambda_tenwolde` looks up both with the
  hard-error form (`solid-thermal-conductivity.h:171,173`), so a
  `PYROLYSIS_ONLY` build with that mechanism stopped at startup with "The
  requested species MOIST is not available".

  The oxidiser turns the gas phase off, not the mechanism. `PYROLYSIS_ONLY`
  sets `gas_start` to pure N2 and gives the boundaries the same value, so no
  reaction of the gas has an oxidiser. That does not stop the gas kinetics.
  No source file reads `GAS_PHASE_REACTIONS` since `3cbc1e9`, and three
  reactions of the scheme need no O2: TAR cracking, CH4 reforming and the
  water-gas shift. The Makefile therefore gives every pyrolysis-only rung
  `TURN_OFF_GAS_REACTIONS=1`, which removes them in the gas and in the
  pores.

  Two consequences. The mechanism carries 8 gas species against 3, so the
  species transport costs about 2.7 times more and a `PYROLYSIS_ONLY` run is
  slower than it was. And the solid mechanism adds the pair
  `MOIST => H2O` / `H2O => MOIST`, which stays active whatever `MOISTURE` is,
  so the pore water of pyrolysis can condense into `MOIST` even in a build
  that starts dry. Every flame case of the ladder already ran that way, so
  this makes the two branches consistent; but it means these runs are no
  longer comparable with `~/temp/{9,10,11}-fatehi`, which used
  `dummy-solid`. Compare new against new. */

  kinfolder = "biomass/dummy-solid-gas";
  shift_prod = true;

#if !NO_INFLOW
  L0 = 20*D0;
  origin (-L0/2, 0);
#else
  L0 = 10*D0;
#endif

#if EMISSIVITY_DIBLASI
  emissivity = emissivity_diblasi;
#else
  emissivity = emissivity_constant;
#endif

  Da = (coord){DA_VALUE, DA_VALUE};

#if GRAVITY && !NO_INFLOW
  /**
  `gravity.h` declares `coord G = {0.,0.,0.}`. Without this line the header
  is present, the acceleration event runs, and the gravity is zero. */

  G.x = -9.81;
#endif

  /**
  Start coarse. `event init` refines near the particle. */

  init_grid(1 << min (maxlevel, 8));

  TOLERANCE = 1e-5;
  NITERMIN = 2;

  run();
}

double r0;
event init (i = 0) {

  /**
  Caution: `navier-stokes/centered.h` assigns `CFL = 0.8` in its `defaults`
  event, which runs after `main()`. So the value belongs here. */

  CFL = CFLNUM;

  /**
  Refine near the particle BEFORE `fraction()`. `fraction()` computes the
  volume fraction on the grid that exists at this point.

  Caution: do not move this `refine()` to `main()`. `run()` calls
  `init_grid (N)` again (`$BASILISK/run.h:17`), and `init_grid` of the tree
  frees the grid. A `refine()` in `main()` does nothing, so the particle
  started on a uniform grid at level 8. The first adapt then built the finest
  cells from coarse PLIC lines: the solid lost 0.3 % at the first step, and the
  corner of the pellet smeared. `run/shrink-corner.c` measures this.

  The disc holds the particle and no more: the corner of a square pellet is at
  0.71 of its size. One cell holds 1104 fields, so a disc of 4 sizes at level
  11 would hold 2.3 GB, and the chemistry event of `i = 0` would run on all of
  it. The adapt of the first steps refines the gas near the particle.

  Caution: runs before this change took `solid_mass0` from the level 8 `f0`.
  Their normalized mass reads about 0.3 % lower. Compare new runs with new
  runs. */

  refine (circle (x, y, 0.75*max (D0, H0)) > 0. && level < maxlevel);

  scalar f0[];

#if SHAPE
  fraction (f0, superquadric (x, y, 20, 0.5*H0, 0.5*D0));
#else
  fraction (f0, circle(x, y, 0.5 * D0));
#endif

#if PYROLYSIS_ONLY
  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
#else
  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;
#endif

#if MOISTURE
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 0.935; // 93.5% biomass
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("MOIST")]   = 0.061; // 6.1% moisture
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("ASH")]     = 0.004; // 0.4% ash
#else
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;
#endif

  foreach()
    porosity[] = eps0*f0[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += f0[]*(1. - eps0)*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

#if !NO_INFLOW
  TG[left] = dirichlet (TG0);
#else
  TG[right] = dirichlet (TG0);
#endif

  TG[top] = dirichlet (TG0);

  /**
  The mechanism carries O2 in either branch. The `PYROLYSIS_ONLY` build gives
  it a mass fraction of zero at the boundaries, so the particle sees none. */

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
#if PYROLYSIS_ONLY
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      #if !NO_INFLOW
      YG[left] = dirichlet (1.);
      #else
      YG[right] = dirichlet (1.);
      #endif
      YG[top] = dirichlet (1.);
    }
#else
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      #if !NO_INFLOW
      YG[left] = dirichlet (0.765);
      #else
      YG[right] = dirichlet (0.765);
      #endif
      YG[top] = dirichlet (0.765);
    } else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
      #if !NO_INFLOW
      YG[left] = dirichlet (0.235);
      #else
      YG[right] = dirichlet (0.235);
      #endif
      YG[top] = dirichlet (0.235);
    }
#endif
    else {
      #if !NO_INFLOW
      YG[left] = dirichlet (0.);
      #else
      YG[right] = dirichlet (0.);
      #endif
      YG[top] = dirichlet (0.);
    }
  }

  if (restore (file = "last-snapshot", list = all)) {
    fprintf (stderr, "Restart file found!\n");
    restarted = true;
  } else {
    fprintf (stderr, "No restart file found, starting from scratch!\n");

    foreach() {
      f[] = f0[];
      porosity[] = eps0*f[];
    }
  }
}

/**
## The override of the reaction rate

Same-name events run in reverse declaration order, so this instance of
`phasechange` runs BEFORE the one of `shrinking.h`, which is declared
earlier. The `chemistry` event of `chemistry.h` attaches to the slot which
`shrinking.h` declares before its own `phasechange`, so it runs before this
one. The order of the step is therefore:

  reset_sources -> chemistry (writes omega) -> this event -> phasechange

Do not move this event, and do not rename it. Under another name it runs
after `shrinking.h` has already built `gas_source`, and the override does
nothing.

Caution: `omega` can carry a negative value, because the pair
`H2O => MOIST` deposits water back into the solid. Mode 2 therefore tests
the integral before it divides. */

#if OMEGA_CONST

/**
The rate which the kinetics asked for, before the override replaces it.
Column 17 of `expansion.dat` reads `omega` after this event, so in mode 2 it
holds `OMEGA_CONST_MDOT` by construction and it measures nothing. This
global carries the true rate to column 24, and it is the column which says
whether the solid still oscillates under a steady blowing. */

double mdot_true = 0.;

event phasechange (i++) {

  if (t < OMEGA_CONST_T0)
    return 0;

  mdot_true = 0.;
  foreach (reduction(+:mdot_true))
    mdot_true += omega[]*(f[] - porosity[])*dv();

#if OMEGA_CONST == 1

  foreach()
    omega[] = (f[] > F_ERR) ? OMEGA_CONST_VALUE : 0.;

#else // OMEGA_CONST == 2

  /**
  A particle which no longer reacts gives no scale. Leave the field as it
  is, and say so once. */

  if (mdot_true <= 1e-30) {
    static bool warned = false;
    if (pid() == 0 && !warned) {
      fprintf (stderr, "# OMEGA_CONST: mdot = %g at t = %g, no override\n",
               mdot_true, t);
      warned = true;
    }
    return 0;
  }

  double scale = OMEGA_CONST_MDOT/mdot_true;
  foreach()
    omega[] *= scale;

#endif
}

#endif // OMEGA_CONST

/**
The H2O-weighted path-averaged temperature of the full case: the mole
fraction over a quarter of the domain.

Caution: the weight follows `MOLAR_DIFFUSION`, because `XGList_G` only
exists when that flag is on. The full case always uses the mole fraction.
So the change of weight is part of the transport rung, and only
`test-transport` and `test-full` weight the probes as the full case does. */

double T_H2O_weigthed_average (double x_interp, int n_samples = 1 << (maxlevel - 1),
                               const double length = L0/4.) {
#ifdef MOLAR_DIFFUSION
  scalar YH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
#else
  scalar YH2O = YGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
#endif

  double numerator = 0., denominator = 0.;
  coord pos, box[2] = {{x_interp, 0.}, {x_interp, length}}, nn = {1, n_samples};
  foreach_region (pos, box, nn, reduction(+:numerator) reduction(+:denominator)) {
    double yH2O_local = interpolate_linear (point, YH2O, pos.x, pos.y, pos.z);
    numerator   += yH2O_local;
    denominator += yH2O_local / interpolate_linear (point, T, pos.x, pos.y, pos.z);
  }

  if (denominator <= 0) // avoid division by 0
    return TG0;

  return numerator/denominator;
}

/**
The same path average with a uniform H2O concentration.

The weight of `T_H2O_weigthed_average` cancels when the concentration is
uniform, and the average becomes the harmonic mean of the temperature along
the path:

  Tuni = N / sum (1/T)

This is the analogue of the weighted probe for a medium which carries the
same H2O everywhere. It needs no reaction to make the absorber, so it holds
a value in a `TURN_OFF_REACTIONS` build, where the weighted probe returns
`TG0`. The two forms differ only through the weight, so a pair of columns
separates the transport of H2O from the field of the temperature.

Caution: every sample counts the same, and the path is a quarter of the
domain, so the free stream at `TG0` dilutes this average. The weighted form
does not carry that dilution, because the plume holds the H2O. Compare the two
columns of one point against each other; do not read one for the other.

`T` is positive everywhere the case runs, but a sample outside the domain
returns `nodata`. The test below drops such a sample from both sums. */

double T_uniform_average (double x_interp, int n_samples = 1 << (maxlevel - 1),
                          const double length = L0/4.) {
  double numerator = 0., denominator = 0.;
  coord pos, box[2] = {{x_interp, 0.}, {x_interp, length}}, nn = {1, n_samples};
  foreach_region (pos, box, nn, reduction(+:numerator) reduction(+:denominator)) {
    double T_local = interpolate_linear (point, T, pos.x, pos.y, pos.z);
    if (T_local > 0. && T_local < nodata) {
      numerator   += 1.;
      denominator += 1./T_local;
    }
  }

  if (denominator <= 0) // avoid division by 0
    return TG0;

  return numerator/denominator;
}

/**
## Higher-order interpolation

`interpolate_linear` is bilinear on a 2x2 stencil. The function below is
biquadratic on the 3x3 stencil of the leaf that holds the sample: a Lagrange
parabola through the centres of the cells at -1, 0 and +1 in each direction.
It is exact for a quadratic field, where the bilinear form is exact only for a
linear field.

The parabola overshoots at a jump. `T` has a jump at the interface, and the
H2O field is near 0 at the edge of the plume, so an overshoot can give a
negative weight. The function therefore clips the result to the minimum and
maximum of the 9 values of the stencil. Where the field is smooth, the clip
does nothing.

Caution: the function is for cell-centred fields only (`T`, `u.x`, `u.y`,
species). It ignores `v.d`, so do not use it for a face field.

Caution: the case is axisymmetric, so the function is 2D only. On the axis
the row at `j = -1` holds the ghost values of the symmetry condition. */

static double interpolate_biquadratic (Point point, scalar v,
                                       double xp = 0., double yp = 0.)
{
  double xi = (xp - x)/Delta, eta = (yp - y)/Delta;
  double wx[3] = {0.5*xi*(xi - 1.), 1. - sq(xi), 0.5*xi*(xi + 1.)};
  double wy[3] = {0.5*eta*(eta - 1.), 1. - sq(eta), 0.5*eta*(eta + 1.)};

  double val = 0., vmin = HUGE, vmax = -HUGE;
  for (int ii = -1; ii <= 1; ii++)
    for (int jj = -1; jj <= 1; jj++) {
      double vij = v[ii,jj];
      val += wx[ii + 1]*wy[jj + 1]*vij;
      vmin = min (vmin, vij);
      vmax = max (vmax, vij);
    }

  return clamp (val, vmin, vmax);
}

/**
## The weighted average with the depth of the mesh

This is `T_H2O_weigthed_average` with two changes. It interpolates with
`interpolate_biquadratic`, and it records the level of the leaf at each
sample. The weight and the path are the same, so the result compares directly
with `TemperatureProfile.dat`.

The line holds 2 samples per cell at `maxlevel`, so every leaf along the line
gets at least one sample. `line_depth` summarises the levels:

  lmin, lmax  the coarsest and the finest leaf on the line
  lmean       the mean level over the samples
  fmax        the fraction of the samples on a leaf at `maxlevel`
  ycoarse     the smallest `y` where the leaf is coarser than `maxlevel`,
              which is how far from the axis the finest mesh reaches.
              It holds `length` when the whole line is at `maxlevel`.

If `fmax` or `ycoarse` moves in time at a fixed point, the mesh adapts along
that line, and the probe changes its resolution during the run.

The function is collective. Call it on every rank. Each sample is on one rank
only, and the reductions combine them. */

typedef struct {
  double lmin, lmax, lmean, fmax, ycoarse;
} line_depth;

double T_H2O_weighted_average_hq (double x_interp, line_depth * depth = NULL,
                                  int n_samples = 1 << (maxlevel - 1),
                                  const double length = L0/4.) {
#ifdef MOLAR_DIFFUSION
  scalar YH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
#else
  scalar YH2O = YGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
#endif

  double numerator = 0., denominator = 0.;
  double lmin = HUGE, lmax = -HUGE, lsum = 0., nsamp = 0., nfine = 0.;
  double ycoarse = length;

  coord pos, box[2] = {{x_interp, 0.}, {x_interp, length}}, nn = {1, n_samples};
  foreach_region (pos, box, nn,
                  reduction(+:numerator) reduction(+:denominator)
                  reduction(min:lmin) reduction(max:lmax)
                  reduction(+:lsum) reduction(+:nsamp) reduction(+:nfine)
                  reduction(min:ycoarse)) {
    double lev = point.level;
    lmin = min (lmin, lev);
    lmax = max (lmax, lev);
    lsum  += lev;
    nsamp += 1.;
    if (point.level >= maxlevel)
      nfine += 1.;
    else
      ycoarse = min (ycoarse, pos.y);

    double yH2O_local = interpolate_biquadratic (point, YH2O, pos.x, pos.y);
    double T_local    = interpolate_biquadratic (point, T, pos.x, pos.y);
    if (T_local > 0. && T_local < nodata) {
      numerator   += yH2O_local;
      denominator += yH2O_local/T_local;
    }
  }

  if (depth) {
    depth->lmin    = nsamp > 0. ? lmin : -1.;
    depth->lmax    = nsamp > 0. ? lmax : -1.;
    depth->lmean   = nsamp > 0. ? lsum/nsamp : -1.;
    depth->fmax    = nsamp > 0. ? nfine/nsamp : 0.;
    depth->ycoarse = ycoarse;
  }

  if (denominator <= 0) // avoid division by 0
    return TG0;

  return numerator/denominator;
}

event output (t += 0.01) {

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, restarted ? "a" : "w");
  if (fp == NULL) {
    fprintf (stderr, "Error opening OutputData\n");
    exit(1);
  }

  /**
  These three points keep the layout of `~/temp/expansion/test.c`, so that
  `OutputData` stays comparable with the runs of the archive. The header of
  that case named them 1, 2 and 4 mm, but it sampled 2, 6 and 11 mm from the
  surface. The names below are the distances the case actually samples.

  The five points of the full case go to `TemperatureProfile.dat`, which the
  same event writes below. Do not add them to `OutputData`:
  `slow_flicker.py` reads the first six columns of this file by position.

  Columns 10 to 12 hold the same three points with a uniform H2O
  concentration, which is the analogue of the weighted average for a medium
  which carries the same absorber everywhere. They go on the end for the same
  reason: a reader which stops at column 9 keeps working. Each column pairs
  with the weighted column of the same point, 10 with 4, 11 with 5, 12 with
  6.

  Caution: columns 8 and 9, `mgp_i` and `mgp_resa`, hold the same two globals
  as columns 11 and 12 of `expansion.dat`, at the same instant. They are a
  duplicate. They stay because they cost nothing and because a removal would
  move columns 10 to 12, which `explore.py` reads by position when a file
  carries no header line. */

  /**
  The six points below are the union of the three points of this file and the
  five points of `TemperatureProfile.dat`. The two files share the point at
  2 mm and the point at 11 mm. An earlier version of the case had one event
  for each file, so it computed those two path averages two times. This event
  writes both files and computes each path average one time. One path average
  costs `2^(maxlevel-1)` interpolations, so the merge saves 2 times 512
  interpolations at each output at level 10.

  Do not change the order of the columns of either file.
  `slow_flicker.py` reads the first six columns of `OutputData` and the first
  six columns of `TemperatureProfile.dat` by position. */

  double Tw[6], Tuni[3];
  double sample_points[6] = {H0/2 + 2e-3, H0/2 + 4e-3, H0/2 + 6e-3,
                             H0/2 + 8e-3, H0/2 + 11e-3, H0/2 + 15e-3};
  double uni_points[3] = {H0/2 + 2e-3, H0/2 + 6e-3, H0/2 + 11e-3};

  for (int ii = 0; ii < 6; ii++)
    Tw[ii] = T_H2O_weigthed_average (sample_points[ii]);

  for (int ii = 0; ii < 3; ii++)
    Tuni[ii] = T_uniform_average (uni_points[ii]);

  /**
  `Tavg` holds the three points of `OutputData`: 2, 6 and 11 mm. */

  double Tavg[3] = {Tw[0], Tw[2], Tw[4]};

  if (i == 0)
    fprintf (fp, "#t(1) Ms/Ms0(2) Tmax(3) Tavg_2mm(4) Tavg_6mm(5) Tavg_11mm(6)"
                 " dt(7) mgp_i(8) mgp_resa(9)"
                 " Tuni_2mm(10) Tuni_6mm(11) Tuni_11mm(12)\n");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[] - porosity[])*rhoS*dv();

  fprintf (fp, "%g %g %g %g %g %g %g %d %g %g %g %g\n",
           t, solid_mass/solid_mass0, probe_stats_T().max,
           Tavg[0], Tavg[1], Tavg[2], dt, mgp.i, mgp.resa,
           Tuni[0], Tuni[1], Tuni[2]);

  fflush(fp);

  /**
  ## The probes of the full case

  `run/fatehi-combustion.c` writes `TemperatureProfile.dat` with five
  H2O-weighted path averages at 2, 4, 8, 11 and 15 mm from the surface. This
  block writes the same five points in the same order, so that the amplitudes
  compare directly with the runs under `~/temp/fatehi` and `slow_flicker.py`
  reads them with the branch it already has. The temperature gap between the
  two cases is a factor 5, which is larger than the gap of the release rate,
  so this file carries the quantity that this campaign must reduce.

  `T_H2O_weigthed_average` is collective: it reduces over `foreach_region`.
  The loop above calls it on every rank, and this block writes on rank 0
  only. */

  if (pid() == 0) {
    static FILE * fpT = NULL;
    if (!fpT) {
      fpT = fopen ("TemperatureProfile.dat", restarted ? "a" : "w");
      if (fpT == NULL) {
        fprintf (stderr, "Error opening TemperatureProfile.dat\n");
        exit (1);
      }
      fprintf (fpT, "#t(1) T2mm(2) T4mm(3) T8mm(4) T11mm(5) T15mm(6)\n");
    }
    fprintf (fpT, "%g %g %g %g %g %g\n",
             t, Tw[0], Tw[1], Tw[3], Tw[4], Tw[5]);
    fflush (fpT);
  }
}

/**
## The high-order profile and the depth of the mesh

The same five lines as `TemperatureProfile.dat`, with
`T_H2O_weighted_average_hq`. Each line gives 6 columns: the temperature, then
`lmin`, `lmax`, `lmean`, `fmax` and `ycoarse`. The file is separate, so the
readers of `TemperatureProfile.dat` keep working.

Compare column `T` of a line here with the same line in
`TemperatureProfile.dat`. A large difference says that the interpolation
error is not small against the oscillation. */

event temperature_profile_hq (t += 0.01) {

  double sample_points[5] = {H0/2 + 2e-3, H0/2 + 4e-3, H0/2 + 8e-3,
                             H0/2 + 11e-3, H0/2 + 15e-3};
  const char * names[5] = {"2mm", "4mm", "8mm", "11mm", "15mm"};
  double Tavg[5];
  line_depth dep[5];

  for (int ii = 0; ii < 5; ii++)
    Tavg[ii] = T_H2O_weighted_average_hq (sample_points[ii], &dep[ii]);

  if (pid() == 0) {
    static FILE * fpH = NULL;
    if (!fpH) {
      fpH = fopen ("TemperatureProfileHQ.dat", restarted ? "a" : "w");
      if (fpH == NULL) {
        fprintf (stderr, "Error opening TemperatureProfileHQ.dat\n");
        exit (1);
      }
      if (!restarted) {
        fprintf (fpH, "#t(1)");
        for (int ii = 0; ii < 5; ii++) {
          int c = 2 + 6*ii;
          fprintf (fpH, " T%s(%d) lmin%s(%d) lmax%s(%d) lmean%s(%d)"
                        " fmax%s(%d) ycoarse%s(%d)",
                   names[ii], c, names[ii], c + 1, names[ii], c + 2,
                   names[ii], c + 3, names[ii], c + 4, names[ii], c + 5);
        }
        fprintf (fpH, "\n");
      }
    }
    fprintf (fpH, "%g", t);
    for (int ii = 0; ii < 5; ii++)
      fprintf (fpH, " %g %g %g %g %g %g", Tavg[ii], dep[ii].lmin, dep[ii].lmax,
               dep[ii].lmean, dep[ii].fmax, dep[ii].ycoarse);
    fprintf (fpH, "\n");
    fflush (fpH);
  }
}

/**
## Point probes around the particle

Each probe gives the local `u.x`, `u.y`, `T`, H2O and the level of the leaf
that holds the point. The values come from `interpolate_biquadratic`. The
level is the local depth of the mesh, so a jump in a probe signal at a change
of level is a grid event, not a physics event.

The H2O column follows the weight of `T_H2O_weigthed_average`: the mole
fraction with `MOLAR_DIFFUSION`, the mass fraction without it. The header
says which one the build writes.

The points are in the frame of `angular_profile`: the flow comes from the
left, `theta = 0` is the downstream pole on `+x`, and `r` is from the centre
of the particle. `R = D0/2 = 4 mm`.

| k | position | x [mm] | y [mm] |
|---|---|---|---|
| 0 | upstream stagnation, R + 1 mm | -5 | 0 |
| 1 | upstream, R + 3 mm | -7 | 0 |
| 2 | 135 deg, R + 1 mm | -3.54 | 3.54 |
| 3 | equator, R + 1 mm | 0 | 5 |
| 4 | equator, R + 3 mm | 0 | 7 |
| 5 | 45 deg, R + 1 mm | 3.54 | 3.54 |
| 6 | downstream, R + 1 mm | 5 | 0 |
| 7 | downstream, R + 4 mm | 8 | 0 |
| 8 | wake, R + 8 mm | 12 | 0 |

Caution: with `NO_INFLOW` the domain starts at `x = 0`. Points 0, 1 and 2
are then outside the domain, and their columns hold `nodata` and level -1.

The file is in long form: one line for each probe at each output. Select a
probe with its column 2. */

#define NPROBE 9

event probe_points (t += 0.01) {

  const double R = 0.5*D0;
  const double c45 = cos (pi/4.);
  coord probes[NPROBE] = {
    {-(R + 1e-3), 0.},
    {-(R + 3e-3), 0.},
    {-(R + 1e-3)*c45, (R + 1e-3)*c45},
    {0., R + 1e-3},
    {0., R + 3e-3},
    {(R + 1e-3)*c45, (R + 1e-3)*c45},
    {R + 1e-3, 0.},
    {R + 4e-3, 0.},
    {R + 8e-3, 0.}
  };

#ifdef MOLAR_DIFFUSION
  scalar YH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
#else
  scalar YH2O = YGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
#endif

  double ux[NPROBE], uy[NPROBE], Tp[NPROBE], yw[NPROBE], lev[NPROBE];

  for (int k = 0; k < NPROBE; k++) {
    double xp = probes[k].x, yp = probes[k].y;
    double vx = nodata, vy = nodata, vT = nodata, vw = nodata, vl = -1.;
    foreach_point (xp, yp, 0., reduction(min:vx) reduction(min:vy)
                   reduction(min:vT) reduction(min:vw) reduction(max:vl)) {
      vx = interpolate_biquadratic (point, u.x, xp, yp);
      vy = interpolate_biquadratic (point, u.y, xp, yp);
      vT = interpolate_biquadratic (point, T, xp, yp);
      vw = interpolate_biquadratic (point, YH2O, xp, yp);
      vl = point.level;
    }
    ux[k] = vx; uy[k] = vy; Tp[k] = vT; yw[k] = vw; lev[k] = vl;
  }

  if (pid() == 0) {
    static FILE * fpp = NULL;
    if (!fpp) {
      fpp = fopen ("probes.dat", restarted ? "a" : "w");
      if (fpp == NULL) {
        fprintf (stderr, "Error opening probes.dat\n");
        exit (1);
      }
      if (!restarted)
        fprintf (fpp, "#t(1) k(2) x(3) y(4) ux(5) uy(6) T(7) %s(8) level(9)\n",
#ifdef MOLAR_DIFFUSION
                 "xH2O"
#else
                 "yH2O"
#endif
                 );
    }
    for (int k = 0; k < NPROBE; k++)
      fprintf (fpp, "%g %d %g %g %g %g %g %g %g\n",
               t, k, probes[k].x, probes[k].y,
               ux[k], uy[k], Tp[k], yw[k], lev[k]);
    fflush (fpp);
  }
}

/**
Diagnostics for the phase-change expansion field, sampled every timestep.

Columns 3 to 5 compare div(uf) with gas_source only:
  Qsrc = int gas_source dV   the source computed by the chemistry
  Qdiv = int div(uf) dV      the expansion the velocity field actually carries
and 'resmax' = max |div(uf) + gas_source|. Note gas_source already carries
cm[], as does the discrete div(uf) below, so both integrate with sq(Delta).

Caution: 'resmax' is NOT the residual of the solver. The projection enforces
div(uf) = -div_source = -(gas_source + drhodt), so 'resmax' is approximately
max |drhodt|. Columns 20 to 23 come from `divergence-budget.h`. They read
div_source at the end of the step, before `adapt`. 'resds' is the residual of
the solver and must be approximately dt*mgp_resa (columns 2 and 12).

'ur' probes the radial velocity on a 45 degree ray just outside the particle:
this is the quantity the velocity vectors show.

Columns 17 to 19 carry the instruments that the amplitude campaign needs.

  mdot    the mass release rate of the particle, integrated:
          `omega*(f - porosity)` is the rate per unit volume, because `omega`
          is the rate per cubic metre of solid material and `(f - porosity)`
          is `f(1 - eps)`, the solid volume fraction. This is the metric of
          the campaign. Log it, do not differentiate `Ms/Ms0`: a numerical
          derivative of a sampled signal changes the amplitude, and the
          campaign compares amplitudes across nine runs.

          Caution: `omega` is positive where the solid decomposes, so `mdot`
          must be positive and must agree with `-d(Ms)/dt`. Check the sign on
          the first run before you trust the column.

  ncells  the number of leaf cells. An adaptation event changes it, so this
          column separates a grid event from a physics event. The pyrolysis
          case loses its spectral peak at level 11, so the mesh is a suspect
          for the excitation of the mode, not only for its resolution.

  nsolid  the number of cells with `f > F_ERR`. It says whether the reacting
          volume switches cells. The gas-source campaign falsified that for
          the fast band. It is free here.
*/

event probe_expansion (t += 0.01) {
  double Qsrc = 0., Qdiv = 0., resmax = 0.;

  foreach (reduction(+:Qsrc) reduction(+:Qdiv) reduction(max:resmax)) {
    double d = 0.;
    foreach_dimension()
      d += uf.x[1] - uf.x[];
    d /= Delta;                        // = cm*div(u), same weighting as gas_source

    Qsrc += gas_source[]*sq(Delta);
    Qdiv += d*sq(Delta);
    resmax = max (resmax, fabs (d + gas_source[]));
  }

  /**
  `angular_profile` needs the same integral for its column `un_pred`. Give it
  to the cache of `probe-cache.h`, so that the grid is swept one time only.
  This event runs before `angular_profile`, because the case declares it
  first. The cache holds the step index, so the order is not a condition:
  a reader of a later step computes the value again. */

  probe_div_uf_set (Qdiv);

  stats so = statsf (omega);

  double ur[3];
  for (int k = 0; k < 3; k++) {
    double r = 0.5*D0 + (k + 1)*0.5e-3, c = cos(pi/4.), s = sin(pi/4.);
    ur[k] = interpolate (u.x, r*c, r*s)*c + interpolate (u.y, r*c, r*s)*s;
  }

  /**
  Particle-side temperature. statsf(T).max is useless here: it sits on TG0
  because the particle starts at TS0 < TG0 and is always the coldest region,
  so a maximum only ever reports the inlet boundary condition.

  We use T[] (assigned once per step as TS[]+TG[], both in tracer form, so it
  is the volume-weighted mixture temperature) together with f[]. Both are
  unambiguous wherever this event runs, unlike TS[]/TG[] whose normalisation
  depends on the position within the timestep.

    Tcore  coldest point, i.e. the centre of the particle
    Tbulk  f-weighted mean over the solid, the bulk particle temperature
    Tsurf  mean over interface cells: where the Arrhenius feedback bites first
  */

  double Tcore = probe_stats_T().min;
  double Tbulk = 0., fvol = 0., Tsurf = 0., nsurf = 0.;
  double mdot = 0., nsolid = 0.;

  foreach (reduction(+:Tbulk) reduction(+:fvol)
           reduction(+:Tsurf) reduction(+:nsurf)
           reduction(+:mdot) reduction(+:nsolid)) {

    /**
    `omega` is zero outside the solid, so this sum needs no test. */

    mdot += omega[]*(f[] - porosity[])*dv();

    if (f[] > F_ERR) {
      nsolid += 1.;
      Tbulk += T[]*f[]*dv();
      fvol  += f[]*dv();
    }
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      Tsurf += T[];
      nsurf += 1.;
    }
  }
  Tbulk = fvol   > 0. ? Tbulk/fvol  : 0.;
  Tsurf = nsurf  > 0. ? Tsurf/nsurf : 0.;

  /**
  Everything above is collective (reductions / interpolate), so every rank
  holds the same values. Only rank 0 writes: basilisk wraps popen for MPI but
  not fopen, so an unguarded fprintf would emit one copy of each line per rank.
  */

  if (pid() == 0) {
    static FILE * fe = fopen ("expansion.dat", restarted ? "a" : "w");
    if (fe == NULL) {
      fprintf (stderr, "Error opening expansion.dat\n");
      exit(1);
    }
    if (i == 0)
      fprintf (fe, "#t(1) dt(2) Qsrc(3) Qdiv(4) resmax(5) omega_min(6) omega_max(7)"
                   " ur_0.5mm(8) ur_1mm(9) ur_1.5mm(10) mgp_i(11) mgp_resa(12)"
                   " mgpsf_i(13) Tcore(14) Tbulk(15) Tsurf(16)"
                   " mdot(17) ncells(18) nsolid(19)"
                   " Qds(20) Qdivb(21) resds(22) dsmax(23)"
#if OMEGA_CONST
                   " mdot_true(24)"
#endif
                   "\n");

    /**
    The new columns go on the end. `slow_flicker.py` reads this file with
    `load(path, 16)`, so it truncates to column 16 and keeps working. */

    fprintf (fe, "%g %g %g %g %g %g %g %g %g %g %d %g %d %g %g %g %g %ld %g"
                 " %g %g %g %g",
             t, dt, Qsrc, Qdiv, resmax, so.min, so.max,
             ur[0], ur[1], ur[2], mgp.i, mgp.resa, mgpsf.i,
             Tcore, Tbulk, Tsurf,
             mdot, grid->tn, nsolid,
             divb_Qds, divb_Qdiv, divb_resmax, divb_dsmax);
#if OMEGA_CONST
    fprintf (fe, " %g", mdot_true);
#endif
    fprintf (fe, "\n");
    fflush (fe);
  }
}

/**
## Angular profile along the particle surface

For a sphere blowing uniformly the *normal* velocity at the surface is the same
at every angle,

  un_pred = Q/(4 pi R^2),   Q = 2 pi Qdiv the total volumetric production,

which reduces to Qdiv/(2 R^2). |u| is *not* uniform even then, because the free
stream adds a tangential component that varies with angle. So 'un' answers "is
the blowing uniform?" while 'umag' is what the |u| screenshots show; logging
both separates the anomaly from the free-stream confound.

theta runs from the downstream pole (0 deg) through the equator (90) to the
upstream pole (180). Along each ray we locate the interface, sample the gas
just outside it, and take the *peak* of omega inside it -- peak rather than a
fixed depth, so the measurement follows the reaction front instead of sliding
off it as the front recedes. 'r_front' vs theta is then the front shape, and
'T_front' the temperature driving the local Arrhenius rate.

'T_gas' (column 10) is the temperature at the same point as 'un', just
outside the interface. The candidate mechanism of the slow oscillation is a
shield: a higher release rate raises the blowing, the blowing holds the hot
gas away from the surface, the surface cools, and the release falls. That
loop needs the temperature difference which drives the surface, and nothing
else in this case measures it. Read 'T_gas - T_front' against 'un'.

Every interpolate() here is collective and called an identical number of times
on every rank (the branch conditions depend only on reduced values), so the
arrays agree across ranks and only rank 0 writes.
*/

#define NANG 24        // angular samples over [0,pi]
#define NRAD 32        // radial samples along each ray

double profile_offset = 0.15;  // gas-side sampling offset, in units of R

event angular_profile (t += 0.01) {
  const double R = 0.5*D0, dr = profile_offset*R;

  double th[NANG], ri[NANG], un[NANG], um[NANG], om[NANG], rf[NANG], Tf[NANG];
  double Tg_out[NANG];

  for (int k = 0; k < NANG; k++) {
    double theta = (k + 0.5)*pi/NANG, c = cos(theta), s = sin(theta);

    /**
    Locate the interface along the ray. Marching outward to the first f < 0.5
    keeps this correct if the interface ever starts moving again.
    */

    double rint = R;
    for (int j = 0; j < NRAD; j++) {
      double rr = (j + 0.5)*1.5*R/NRAD;
      double ff = interpolate (f, rr*c, rr*s);
      if (ff != nodata && ff < 0.5) { rint = rr; break; }
    }

    /**
    Peak reaction rate along the ray, and the temperature where it peaks.
    */

    double ommax = 0., rfront = 0.;
    for (int j = 0; j < NRAD; j++) {
      double rr = (j + 0.5)*rint/NRAD;
      double o = interpolate (omega, rr*c, rr*s);
      if (o != nodata && fabs(o) > fabs(ommax)) {
        ommax = o; rfront = rr;
      }
    }

    /**
    Gas side, just outside the interface.
    */

    double ux = interpolate (u.x, (rint + dr)*c, (rint + dr)*s);
    double uy = interpolate (u.y, (rint + dr)*c, (rint + dr)*s);
    double Tg = interpolate (T, (rint + dr)*c, (rint + dr)*s);

    th[k] = theta*180./pi;
    ri[k] = rint;
    un[k] = ux*c + uy*s;                 // outward normal component = blowing
    um[k] = sqrt (sq(ux) + sq(uy));
    om[k] = ommax;
    rf[k] = rfront;
    Tf[k] = interpolate (T, rfront*c, rfront*s);
    Tg_out[k] = Tg;
  }

  /**
  Total production, for the uniform-blowing reference. Same weighting as in
  probe_expansion: the discrete divergence carries cm[], so it integrates
  with sq(Delta). `probe_expansion` has already computed this integral in
  this step, so `probe_div_uf()` returns its value and sweeps no grid. */

  double Qdiv = probe_div_uf();
  double un_pred = Qdiv/(2.*sq(R));

  if (pid() == 0) {
    static FILE * fa = fopen ("angular.dat", restarted ? "a" : "w");
    if (fa == NULL) {
      fprintf (stderr, "Error opening angular.dat\n");
      exit(1);
    }
    if (i == 0)
      fprintf (fa, "#t(1) theta_deg(2) r_int(3) un(4) umag(5) omega_max(6)"
                   " r_front(7) T_front(8) un_pred(9) T_gas(10)\n");

    /**
    `T_gas` goes on the end, because `slow_flicker.py` reads this file with
    `load(path, 9)` and truncates to column 9. */

    for (int k = 0; k < NANG; k++)
      fprintf (fa, "%g %g %g %g %g %g %g %g %g %g\n",
               t, th[k], ri[k], un[k], um[k], om[k], rf[k], Tf[k], un_pred,
               Tg_out[k]);
    fflush (fa);
  }
}

#if TREE
event adapt (i++) {
  scalar oxidiser = YGList_G[OpenSMOKE_IndexOfSpecies ("O2")];

  adapt_wavelet_leave_interface ({T, oxidiser}, {f},
    (double[]){ADAPT_T_TOL, ADAPT_O_TOL}, maxlevel, minlevel, 2);
  

  // Unrefine for outflow condition
  unrefine (x > L0*0.4);
}
#endif

/**
The movie shows the isolines of the TAR mass fraction in the two halves, over
the temperature and over the oxidiser. `TAR_G + TAR_S` is the mixture value,
as `O2_G + O2_S` is for the oxidiser. The five levels go from 0.02 to 0.1 in
steps of 0.02. Change `TAR_MIN`, `TAR_MAX` and `TAR_NISO` if the plume holds
less or more TAR. */

#ifndef TAR_MIN
# define TAR_MIN 0.02
#endif

#ifndef TAR_MAX
# define TAR_MAX 0.1
#endif

#ifndef TAR_NISO
# define TAR_NISO 5
#endif

event movie (t += 1) {
  clear();
  view (theta=0, phi=0, psi=-pi/2., width = 1080, height = 1080);
  squares ("T", min = 300, max = 2000, spread = -1, linear = true);
  isoline ("TAR_G + TAR_S", n = TAR_NISO, min = TAR_MIN, max = TAR_MAX,
           lw = 1., lc = {1., 1., 1.});
  draw_vof ("f", lw = 1.5);
  mirror ({0, 1}) {
    squares ("O2_G + O2_S", min = 0., max = 0.235, spread = -1, linear = true);
    isoline ("TAR_G + TAR_S", n = TAR_NISO, min = TAR_MIN, max = TAR_MAX,
             lw = 1., lc = {1., 1., 1.});
    draw_vof ("f", lw = 1.5);
  }
  save ("movie.mp4");
}

event dump (t = 1; t += 1) {
  dump("last-snapshot");
}

/**
## The end of the run

The event returns 1, so `events()` stops the loop at once. A bare
`event stop (t = tend)` also ends this case, because no other event carries
an upper limit on `t`; but it first runs one more full step, whose result no
file holds. The `return` also protects the case against a new event with a
condition such as `t <= X`, which would hold the loop open past `tend`. */

event stop (t = tend) {
  return 1;
}

/**
~~~gnuplot
~~~
**/
