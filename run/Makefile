include $(BASILISK)/Makefile.defs

CFLAGS += -O2 -disable-dimensions -I../src  -I$(HOME)/basilisk/OpenSMOKEppInterface/src  

recirculation.tst: CFLAGS+= -I$(HOME)/basilisk/basilisk-sandbox-ecipriano/src
marangoni-deff.tst: CFLAGS+= -I$(HOME)/basilisk/basilisk-sandbox-ecipriano/src
marangoni-sweep.tst: CFLAGS+= -I$(HOME)/basilisk/basilisk-sandbox-ecipriano/src
marangoni-point.tst: CFLAGS+= -I$(HOME)/basilisk/basilisk-sandbox-ecipriano/src
marangoni-hill.tst: CFLAGS+= -I$(HOME)/basilisk/basilisk-sandbox-ecipriano/src




# Flame-jitter experiment on the fatehi configuration (see gas-source.c).
# One source; each case is a symbolic link that this Makefile creates on
# demand, so only gas-source.c lives in the repository. Makefile.defs already
# declares %.c precious, so the links are not deleted as intermediates.
#
# Run gas-source-ignite first, then copy its ignition-snapshot into each case
# directory, so every case branches from the same state.

# Ask make for a variable, so a script can reuse the exact include paths
# of this Makefile instead of repeating them:  make print-CFLAGS
.PHONY: print-%
print-%:
	@echo '$($*)'

GSAB_CASES = ignite inst avg zmix zeta dt exact exact0 dtfix dtvar
GSAB_LINKS = $(GSAB_CASES:%=gas-source-%.c)

gas-source-%.c: gas-source.c
	ln -sf $< $@

.PHONY: gas-source-links gas-source-unlink
gas-source-links: $(GSAB_LINKS)
gas-source-unlink:
	rm -f $(GSAB_LINKS)

GSAB_REF = -DGAS_SOURCE_AVERAGED=0
GSAB_AVG = -DGAS_SOURCE_AVERAGED=1

# Fixed-timestep pair. `gas-source-sdiag` measured a period-2 oscillation in
# `Qrho` with `corr(dQrho, d(dt)) = +0.84`, so `dt` is a suspect in the loop
# source -> divergence -> velocity -> CFL -> dt -> source. The pair separates
# "dt varies" from "dt is smaller": both run at the same mean step, and only
# `dtfix` holds it constant.
#
# GSAB_DTFIX must sit below the smallest CFL-limited step, otherwise the cap
# does not bind and the test is void. The measured window gave a mean step of
# 2.52e-4 s at CFL 0.5, with a spread of 8.4 percent. Check that column 2 of
# probe.dat is constant; if it is not, lower this value.
#
# GSAB_DTVAR_CFL scales the CFL-limited step to the same mean:
#   0.5 * 1.75e-4 / 2.52e-4 = 0.347
# Both are overridable:  make gas-source-dtfix.tst GSAB_DTFIX=1.5e-4
GSAB_DTFIX = 1.75e-4
GSAB_DTVAR_CFL = 0.347

gas-source-ignite.s:   CFLAGS += $(GSAB_REF) -DIGNITION=1 -DTEND=12
gas-source-ignite.tst: CFLAGS += $(GSAB_REF) -DIGNITION=1 -DTEND=12

gas-source-inst.s:     CFLAGS += $(GSAB_REF)
gas-source-inst.tst:   CFLAGS += $(GSAB_REF)

gas-source-avg.s:      CFLAGS += -DGAS_SOURCE_AVERAGED=1
gas-source-avg.tst:    CFLAGS += -DGAS_SOURCE_AVERAGED=1

gas-source-zmix.s:     CFLAGS += $(GSAB_REF) -DFLAME_PRINT_TIME=0.005
gas-source-zmix.tst:   CFLAGS += $(GSAB_REF) -DFLAME_PRINT_TIME=0.005

gas-source-zeta.s:     CFLAGS += $(GSAB_REF) -DZETA_POLICY=ZETA_CONST
gas-source-zeta.tst:   CFLAGS += $(GSAB_REF) -DZETA_POLICY=ZETA_CONST

gas-source-dt.s:       CFLAGS += $(GSAB_REF) -DCFLNUM=0.125
gas-source-dt.tst:     CFLAGS += $(GSAB_REF) -DCFLNUM=0.125

# Exact expansion source with the divergence-source filter, and without it.
gas-source-exact.s:    CFLAGS += -DGAS_SOURCE_EXACT=1
gas-source-exact.tst:  CFLAGS += -DGAS_SOURCE_EXACT=1

gas-source-exact0.s:   CFLAGS += -DGAS_SOURCE_EXACT=1 -DGAS_SOURCE_FILTER_PASSES=0
gas-source-exact0.tst: CFLAGS += -DGAS_SOURCE_EXACT=1 -DGAS_SOURCE_FILTER_PASSES=0

# Fixed step. CFLNUM=2 never binds at this DTMAX, so DTMAX sets every step.
# It still caps a runaway, which an unlimited CFL would not.
GSAB_DTFIX_FLAGS = $(GSAB_AVG) -DDTMAX=$(GSAB_DTFIX) -DCFLNUM=2
gas-source-dtfix.s:    CFLAGS += $(GSAB_DTFIX_FLAGS)
gas-source-dtfix.tst:  CFLAGS += $(GSAB_DTFIX_FLAGS)

# The control: the same mean step, still set by the CFL.
GSAB_DTVAR_FLAGS = $(GSAB_AVG) -DCFLNUM=$(GSAB_DTVAR_CFL)
gas-source-dtvar.s:    CFLAGS += $(GSAB_DTVAR_FLAGS)
gas-source-dtvar.tst:  CFLAGS += $(GSAB_DTVAR_FLAGS)


# ---------------------------------------------------------------------------
# Amplitude ladder on the fatehi configuration (see test.c).
#
# The reduced case and the full case oscillate at the same frequency and at
# different amplitudes: at matched remaining mass the release rate gives 3.5
# to 4.1 percent in the full case against 1.6 to 2.3 percent in the reduced
# case, and the probe temperature gives 23 K against 4.4 K. Each case below
# adds ONE ingredient of fatehi-combustion.c to the reduced case, so that the
# amplitude can be attributed to that ingredient.
#
# One source; each case is a symbolic link that this Makefile creates on
# demand, as for gas-source.c above. Makefile.defs declares %.c precious, so
# the links survive, and .gitignore keeps them out of the repository.
#
# Read the results at matched remaining mass, never at matched time. The
# cases burn at different rates, so the same instant is a different state.
#
#   make test-base.tst test-da.tst test-moist.tst ...
#   python3 ~/oscillation-report/analysis/slow_flicker.py test-da --t0 10 --t1 30
#
# Caution: verify the header line of the log before you quote a case. It
# prints what the build actually enabled, not what the case name promises.

TLAD_CASES = base da moist grav transport emis shape full l11 l11noshrink\
             dtfix dtvar pyro uin2 uin05 noheat heat \
             heatnoshrink pyronoshrink basenoshrink fullnoshrink \
             constprop constpropnoshrink noexp pyrotsoff pyronoshrinktsoff l9 \
             probe probefl probem \
             robin robinfl robinm \
             picard picardfl picardm picardnr \
             rdt1 rdt2 pdt1 pdt2 \
             fullprobe fullrobin fullpicard \
             robindb robinls fullrobinls fullrobin20 fullrobininf \
             vofbc vofbcfl vofbcbase vofbcm vofbcnp fullvofbc vofbctrans \
             atol2 atol4 toloff \
             molar fick mde \
             transportold transportcfl transportlim transportmde fullold \
             transrdt2 transoldrdt2 \
             frozen frozenloose frozenfull

TLAD_LINKS = $(TLAD_CASES:%=test-%.c)

test-%.c: test.c
	ln -sf $< $@

.PHONY: test-links test-unlink
test-links: $(TLAD_LINKS)
test-unlink:
	rm -f $(TLAD_LINKS)

# Every rung of this ladder is read against the runs in ~/test/solid, and those
# were made when PIN_SOLID_INTERIOR was 1 by default. The default became 0 on
# 2026-09-09, so the ladder sets it back to 1 here.
#
# Caution: do not remove this line while the archive is the control. Without
# the pin the wavelet criterion coarsens the body of the particle, and it
# coarsens it differently for each case, so a comparison between two rungs
# measures the grid as much as the physics.

TLAD_PIN = -DPIN_SOLID_INTERIOR=1

test-%.s:   CFLAGS += $(TLAD_PIN)
test-%.tst: CFLAGS += $(TLAD_PIN)

# The three transport flags travel together: the full case has all three.
TLAD_TRANSPORT = -DMOLAR_DIFFUSION -DFICK_CORRECTED -DMASS_DIFFUSION_ENTHALPY

# Every ingredient at the full-case value. This build must approach
# fatehi-combustion.c. If it does not, one ingredient is still missing.
TLAD_FULL = $(TLAD_TRANSPORT) -DDA_VALUE=1e-10 -DMOISTURE=1 -DGRAVITY=1 \
            -DSHAPE=1 -DEMISSIVITY_DIBLASI=1

# test-base takes no flags: the pattern rule creates the link, and the
# defaults of test.c reproduce ~/temp/expansion/test.c.

test-da.s:          CFLAGS += -DDA_VALUE=1e-10
test-da.tst:        CFLAGS += -DDA_VALUE=1e-10

test-moist.s:       CFLAGS += -DMOISTURE=1
test-moist.tst:     CFLAGS += -DMOISTURE=1

test-grav.s:        CFLAGS += -DGRAVITY=1
test-grav.tst:      CFLAGS += -DGRAVITY=1

test-transport.s:   CFLAGS += $(TLAD_TRANSPORT)
test-transport.tst: CFLAGS += $(TLAD_TRANSPORT)

test-emis.s:        CFLAGS += -DEMISSIVITY_DIBLASI=1
test-emis.tst:      CFLAGS += -DEMISSIVITY_DIBLASI=1

test-shape.s:       CFLAGS += -DSHAPE=1
test-shape.tst:     CFLAGS += -DSHAPE=1

test-full.s:        CFLAGS += $(TLAD_FULL)
test-full.tst:      CFLAGS += $(TLAD_FULL)

# The grid test. The pyrolysis-only case loses its spectral peak at level 11,
# so the mode there is not grid-converged. The flame case has never been run
# at level 11. Compare this one against test-base at matched remaining mass.
test-l11.s:         CFLAGS += -DMAXLEVEL_VALUE=11
test-l11.tst:       CFLAGS += -DMAXLEVEL_VALUE=11

# A variant of the grid test that doesn't shrink the grid.
test-l11noshrink.s: CFLAGS += -DMAXLEVEL_VALUE=11 -DNO_SHRINK=1
test-l11noshrink.tst: CFLAGS += -DMAXLEVEL_VALUE=11 -DNO_SHRINK=1

# ---------------------------------------------------------------------------
# Loop block, added 2026-09-04. See ~/oscillation-report/solid-ladder/.
#
# Two bands, and they need different cases.
#
#   (a) Is the slow band an artifact of the step?      dtfix, dtvar
#   (b) What closes the loop, the blowing or the heat
#       of reaction?                                   pyro, uin2, uin05, noheat
#   (c) Does the mode survive refinement with a flame? l9 (with base and l11)
#
# Question (b) does NOT need a flame. The release oscillates 7 to 18 per cent
# peak-to-peak in ~/temp/{9,10,11}-fatehi, which carry no oxidiser. Those runs
# are much cheaper and give 6 cycles in 20 s, against 1.2 to 4.3 cycles in the
# runs with a flame. So (b) runs with PYROLYSIS_ONLY=1.
#
# Caution: quote peak-to-peak, and say so. The rms is 3.5 to 5 times smaller.
# Caution: band-pass, never a Savitzky-Golay window near the period you want.
#          A window of 2 s deletes the 0.3 Hz band. Use a corner of 0.01 to
#          0.02 Hz for the release; 0.05 Hz costs about 30 per cent of it.

# --- (a) the timestep pair -------------------------------------------------
# The pair separates "dt varies" from "dt is smaller": both hold the same mean
# step and only dtfix keeps it constant.
#
# TLAD_DTFIX must divide 0.01 exactly, or dtnext() subdivides the step to land
# on the four `t += 0.01` events and dt is not constant. 2e-4 = 0.01/50.
# CFLNUM=2 stops the CFL from binding, so DT_VALUE sets every step; it still
# caps a runaway.
#
# Verify column 2 of expansion.dat is constant before you quote dtfix. If it
# is not, lower TLAD_DTFIX to the next divisor of 0.01.
#
# TLAD_DTVAR_CFL scales the CFL-limited step to the same mean. Measure the
# mean step of test-base over the same window and set
#   TLAD_DTVAR_CFL = 0.8 * TLAD_DTFIX / <mean dt of test-base>
# The value below assumes a mean of 4.55e-4:  0.8 * 2e-4 / 4.55e-4 = 0.35
# Both are overridable:  make test-dtfix.tst TLAD_DTFIX=1e-4
#
# Predictions. In dtfix the single-sample Tmax spikes must disappear, and the
# slow band must not move. If the spikes stay, Tmax = f(dt) is wrong. If the
# slow band moves, it carries a component of the step.
#
# Caution: a run from scratch is NOT constant. Measured on 2026-09-04, the
# first steps come out at 1.93e-4 and 1.96e-4, which are 0.01/52 and 0.01/51,
# so the CFL binds through the impulsive start even at CFLNUM=2.
#
# So BRANCH BOTH RUNS FROM A PLATEAU SNAPSHOT. `event dump (t = 1; t += 1)`
# writes `last-snapshot`, and `event init` restores it, so copy the file of a
# finished run into the case directory before you start.
#
# CORRECTED 2026-09-07. An earlier version of this comment said that both runs
# die because the impulsive start of a run from scratch puts the true CFL near
# 2.8. That is wrong, and it cost a repeat of the same failed pair.
#
# The evidence: both logs say `Restart file found!`, so both DID branch from a
# snapshot, and both die at t = 6.47 and t = 6.49, not at the start. The state
# before the crash is healthy - dt nominal, resmax 1 to 2.5, Tcore 302.8 K,
# Tsurf 729 K rising smoothly - and the one thing that moves is omega_max,
# which climbs from 3690 to 5693 in 0.11 s. Both die ON THE IGNITION RAMP.
#
# It is not a timestep instability either. test-dtvar runs at CFLNUM=0.35,
# which is MORE conservative than the 0.8 of test-base, with a smaller step
# (3.30e-4 against 4.09e-4) and a lower residual (2.45 against 2.73) at the
# same instant, and it dies while test-base lives.
#
# THE FIX IS INT_TEMP_ROBIN, and the proof already exists. test-rdt1 carries
# the SAME DT=2e-4 CFL=2 as test-dtfix plus TLAD_ROBIN, and it crosses t=6.45,
# 6.60, 7.0 and 8.0 and runs to t=20.4. So test-dtfix + Robin IS test-rdt1:
# do not re-run it. The rung that is missing is the varying-step partner at a
# matched mean step. Robin is near enough to neutral to allow the comparison:
# over the matched-mass window test-rdt1 and test-base agree to 0.25 per cent
# on mdot and 4.3 K on Tsurf.
#
# Caution: test-rdt1, test-rdt2, test-pdt1 and test-pdt2 were all BUILT
# 2026-09-04 18:23, which is BEFORE the heat-loss fix of d2d5f2c
# (2026-09-06 16:59). They ran with INT_TEMP_ROBIN_SMAX=1 and no
# INT_TEMP_ROBIN_DEBT. Re-run the block before you quote it against test-base.
# To check any run: `stat -c %y <case>/warn` is the BUILD time, because the
# .tst harness writes that file from the compiler. The mtime of the data files
# only says when the run last wrote, and a running binary never picks up a
# source change.
TLAD_DTFIX     = 2e-4
TLAD_DTVAR_CFL = 0.35

test-dtfix.s:       CFLAGS += -DDT_VALUE=$(TLAD_DTFIX) -DCFLNUM=2
test-dtfix.tst:     CFLAGS += -DDT_VALUE=$(TLAD_DTFIX) -DCFLNUM=2

test-dtvar.s:       CFLAGS += -DCFLNUM=$(TLAD_DTVAR_CFL)
test-dtvar.tst:     CFLAGS += -DCFLNUM=$(TLAD_DTVAR_CFL)

# --- (b) what closes the loop ----------------------------------------------
# test-pyro is the reference of this group and reproduces ~/temp/10-fatehi.
# The gate of the frozen gas cells (src/chemistry.h).
#
# `FROZEN_CELL_GATE` is OFF by default. In a cell that does not react, the
# stiff solve of the 88-species scheme still costs about 100 evaluations of
# the right-hand side. The gate spends one evaluation to find that out, takes
# the explicit update, which is exact at that size, and skips the solve.
#
# `test-frozen` is the rung to compare against `test-base`. Both must give the
# same answer: at the default tolerance of 1e-15 a 2-D A/B at level 8 was
# equal to the ungated run in every column of every output row.
#
#   make test-base.tst test-frozen.tst
#   diff test-base/OutputData-* test-frozen/OutputData-*
#
# `test-frozenloose` opens the tolerance to 1e-10. It is faster, 1.43x against
# 1.28x on a level-7 fatehi with the 88-species scheme, but it perturbs the
# trajectory: 21 of 62 output rows moved in the same A/B. Use it only when the
# run does not have to compare against another run.
#
# `test-frozenfull` is the gate on top of the full flag set, to check that the
# gate and the transport flags do not interact.
#
# Caution: compare two builds only inside one batch of runs. The absolute time
# of the case swings 15 per cent from one run to the next on the same binary,
# while the ratio inside a batch repeats to 1 per cent. Wall-clock claims made
# across batches are meaningless.

TLAD_FROZEN = -DFROZEN_CELL_GATE=1
TLAD_FROZEN_LOOSE = -DFROZEN_CELL_GATE=1 -DFROZEN_CELL_YTOL=1e-10 \
                    -DFROZEN_CELL_TTOL=1e-6

test-frozen.s:          CFLAGS += $(TLAD_FROZEN)
test-frozen.tst:        CFLAGS += $(TLAD_FROZEN)

test-frozenloose.s:     CFLAGS += $(TLAD_FROZEN_LOOSE)
test-frozenloose.tst:   CFLAGS += $(TLAD_FROZEN_LOOSE)

test-frozenfull.s:      CFLAGS += $(TLAD_FULL) $(TLAD_FROZEN)
test-frozenfull.tst:    CFLAGS += $(TLAD_FULL) $(TLAD_FROZEN)

TLAD_PYRO = -DPYROLYSIS_ONLY=1

test-pyro.s:        CFLAGS += $(TLAD_PYRO)
test-pyro.tst:      CFLAGS += $(TLAD_PYRO)

# The blowing sweep. v_w follows the rate of pyrolysis, not Uin, so Uin moves
# the blowing ratio (about 1.06 now) and not the Peclet number. A larger Uin
# thins the layer and weakens the blockage; a smaller Uin strengthens it.
#   uin2   ratio about 0.53, the amplitude must FALL
#   uin05  ratio about 2.1,  the amplitude must RISE
# If the amplitude does not move with Uin, the blowing does not close the loop.
#
# Caution: Uin also changes the supply of heat, so it changes the rate of
# burn. Compare at matched remaining mass, never at matched time.
test-uin2.s:        CFLAGS += $(TLAD_PYRO) -DUIN_VALUE=0.26
test-uin2.tst:      CFLAGS += $(TLAD_PYRO) -DUIN_VALUE=0.26

test-uin05.s:       CFLAGS += $(TLAD_PYRO) -DUIN_VALUE=0.065
test-uin05.tst:     CFLAGS += $(TLAD_PYRO) -DUIN_VALUE=0.065

# Removes the endothermic arm and nothing else. The reaction absorbs 165 to
# 203 kJ for each kg of biomass, so it is a real sink. If the oscillation
# stays, the loop is the blowing.
#
# Caution: the flag is existence-tested in reactors.h. test.c undefines it
# when the value is 0, so -DTURN_OFF_HEAT_OF_REACTION=0 means off, and both
# the bare form and =1 mean on. Do not pass =0 to a build that lacks that
# block.
test-noheat.s:      CFLAGS += $(TLAD_PYRO) -DTURN_OFF_HEAT_OF_REACTION
test-noheat.tst:    CFLAGS += $(TLAD_PYRO) -DTURN_OFF_HEAT_OF_REACTION

# Pure heating: no reaction at all. The particle keeps its mass and its
# composition, and the case gives the transient of pure heat conduction. It is
# the thermal reference of the ladder, and it is cheap because it runs no
# stiff solve.
#
# Check the columns that prove the reactions are off: columns 3, 6, 7 and 17
# of expansion.dat hold 0, and column 2 of OutputData stays constant.
#
# Caution: column 2 of OutputData does not hold exactly 1. The first adaptation
# changes the discrete volume of the solid by about 0.6 per cent, so the value
# falls to 0.9942 at the first output and holds there.
#
# Caution: chemistry.h tests the flag with #ifndef, so -DTURN_OFF_REACTIONS=0
# would also turn the reactions off. test.c undefines it when the value is 0.
# Do not pass =0 to a build that lacks that block.
#
# Caution: no reaction makes H2O, so the H2O-weighted probes return TG0. Read
# Tcore, Tbulk and Tsurf, columns 14 to 16 of expansion.dat.
test-heat.s:        CFLAGS += $(TLAD_PYRO) -DTURN_OFF_REACTIONS
test-heat.tst:      CFLAGS += $(TLAD_PYRO) -DTURN_OFF_REACTIONS

# The fixed-interface block. ZETA_SWELLING sets zeta to zero,
# so the source of the velocity potential is zero, the solid velocity ubf is
# zero, and nothing advects the interface. The particle keeps its volume and
# the decomposition raises the porosity instead.
#
# gas_source does not carry zeta, so the release of mass into the gas is the
# same in both policies. The pair therefore separates the motion of the
# interface from the release of mass.
#
# test-heatnoshrink is a control, not an experiment. With no reaction omega is
# zero, so ZETA_REACTION already gives zeta = 0, which is what ZETA_SWELLING
# gives. The run must equal test-heat in every column:
#
#   make test-heat.tst test-heatnoshrink.tst
#   diff test-heat/OutputData-10 test-heatnoshrink/OutputData-10
#
# The other three release mass, so the two policies give two different runs.
# Each one pairs with the rung of the same name:
#
#   test-pyronoshrink   against test-pyro    pyrolysis only, no flame
#   test-basenoshrink   against test-base    the reduced case with a flame
#   test-fullnoshrink   against test-full    every ingredient of the full case
#
# Caution: the particle keeps its volume, so it does not lose surface as it
# decomposes, and it burns at a different rate. Compare at matched remaining
# mass, never at matched time.
TLAD_NOSHRINK = -DZETA_POLICY=ZETA_SWELLING

test-heatnoshrink.s:    CFLAGS += $(TLAD_PYRO) -DTURN_OFF_REACTIONS $(TLAD_NOSHRINK)
test-heatnoshrink.tst:  CFLAGS += $(TLAD_PYRO) -DTURN_OFF_REACTIONS $(TLAD_NOSHRINK)

test-pyronoshrink.s:    CFLAGS += $(TLAD_PYRO) $(TLAD_NOSHRINK)
test-pyronoshrink.tst:  CFLAGS += $(TLAD_PYRO) $(TLAD_NOSHRINK)

# The transport of heat by the pore gas, off. TS_PORE_ADVECTION=0 removes the
# advection of TS with u_prime in multicomponent-varprop.h and nothing else.
# The pore gas leaves the reaction front cold and cools the char on its way
# out, so a rise of the release cools the surface at once. In
# test-pyronoshrink that is corr(mdot, Tsurf) = -0.92 at zero lag, with a slow
# band of 0.80 per cent of mdot.
#
#   test-pyronoshrinktsoff  against test-pyronoshrink. The cleanest pair: no
#                           flame, no motion of the interface. If the band and
#                           the -0.92 both go, the loop runs through this term.
#   test-pyrotsoff          against test-pyro. Tells whether the amplifier of
#                           the moving interface needs the term too.
#
# Caution: test-pyro is a 09-04 build. Rebuild it next to test-pyrotsoff, or
# the pair measures a week of source changes as well.
#
# Caution: a build without the term is a different physical model. It
# measures the mechanism. It is not a fix.
#
# Check the header line: TSADV=0.
TLAD_TSOFF = -DTS_PORE_ADVECTION=0

test-pyronoshrinktsoff.s:   CFLAGS += $(TLAD_PYRO) $(TLAD_NOSHRINK) $(TLAD_TSOFF)
test-pyronoshrinktsoff.tst: CFLAGS += $(TLAD_PYRO) $(TLAD_NOSHRINK) $(TLAD_TSOFF)

test-pyrotsoff.s:           CFLAGS += $(TLAD_PYRO) $(TLAD_TSOFF)
test-pyrotsoff.tst:         CFLAGS += $(TLAD_PYRO) $(TLAD_TSOFF)

test-basenoshrink.s:    CFLAGS += $(TLAD_NOSHRINK)
test-basenoshrink.tst:  CFLAGS += $(TLAD_NOSHRINK)

test-fullnoshrink.s:    CFLAGS += $(TLAD_FULL) $(TLAD_NOSHRINK)
test-fullnoshrink.tst:  CFLAGS += $(TLAD_FULL) $(TLAD_NOSHRINK)

# The constant-property pair. CONST_PROPERTIES replaces opensmoke-properties.h
# with constant-properties.h, so the gas carries one density, one viscosity,
# one conductivity and one heat capacity over the whole domain, and the solid
# carries one conductivity and one heat capacity. run/POM.c builds the same
# stack this way.
#
# The swap removes VARPROP, because variable-properties.h is what defines it.
# The density is then a constant, so the gas does not expand with the
# temperature and the divergence carries only the source of the phase change.
# Read the pair against test-base and test-basenoshrink.
#
# The six values are the ones of the commented block in main(), for air at
# about 1100 K. Each takes a -D of its own:
#   RHOG_CONST 0.31   MUG_CONST 4.5e-5   LAMBDAG_CONST 0.08
#   CPG_CONST 1200    LAMBDAS_CONST 0.2  CPS_CONST 1500
#
# Caution: constant-properties.h fixes the diffusion coefficient of every
# species at 2.05e-5 m^2/s, CO in N2 at 500 K, and the case cannot change it.
# At 1123 K the varprop build gives about ten times more. Do not read a
# difference in the transport of the species as a result of the constant
# properties alone.
TLAD_CONSTPROP = -DCONST_PROPERTIES=1

test-constprop.s:           CFLAGS += $(TLAD_CONSTPROP)
test-constprop.tst:         CFLAGS += $(TLAD_CONSTPROP)

test-constpropnoshrink.s:   CFLAGS += $(TLAD_CONSTPROP) $(TLAD_NOSHRINK)
test-constpropnoshrink.tst: CFLAGS += $(TLAD_CONSTPROP) $(TLAD_NOSHRINK)

# The base case with no expansion of the gas. NO_EXPANSION drops drhodt from
# the source of the projection (centered-phasechange.h:120) and stops the call
# to update_divergence() (multicomponent-varprop.h:1071), so drhodt keeps the
# 0 that reset_sources gives it.
#
# The particle still releases its mass: gas_source carries the phase change
# and it does not pass through drhodt. What goes away is the volume the gas
# gains as it heats up and as its composition changes. Read this one against
# test-base.
#
# test-constprop also has no expansion, because it has no VARPROP and nothing
# writes drhodt. That build changes every property as well. This one keeps the
# variable properties, so it is the clean test of the expansion alone.
#
# Caution: the two sites test the flag with #ifndef, so -DNO_EXPANSION=0 would
# also turn the expansion off. test.c undefines it when the value is 0. Do not
# pass =0 to a build that lacks that block.
test-noexp.s:       CFLAGS += -DNO_EXPANSION
test-noexp.tst:     CFLAGS += -DNO_EXPANSION

# --- (c) the grid sweep with a flame ---------------------------------------
# Without a flame the amplitude falls with refinement: 17.9, 8.0, 7.2 per cent
# peak-to-peak at levels 9, 10 and 11. That is a damped mode that the noise of
# the grid excites. With l9, test-base and test-l11 the same sweep exists for
# the flame case. If the amplitude holds there while it falls without a flame,
# the feedback of the flame is what sustains the mode.
test-l9.s:          CFLAGS += -DMAXLEVEL_VALUE=9
test-l9.tst:        CFLAGS += -DMAXLEVEL_VALUE=9

# ---------------------------------------------------------------------------
# Probe block, added 2026-09-04. See src/int-temperature-probe.h.
#
# The probe measures the error of the partitioned coupling between TS and TG.
# It changes no field, so a probe run and its plain twin must give the same
# trajectory. Verify that before you quote the probe.
#
# It writes intres.dat, one line per step. Two columns decide the question:
#
#   SG_max   the heat that the frozen interface source delivers in one step,
#            divided by the heat capacity of the gas in that cell. A value
#            above 1 means the source overshoots the cell.
#   rel_max  the interface energy balance, evaluated again after the two
#            diffusion solves, divided by the largest flux. It is zero before
#            the solves by construction.
#
# Read the header of src/int-temperature-probe.h for the decision table.
#
# Caution: the file holds one line per step, about 150 bytes. A 30 s run at
# dt = 2.5e-4 gives about 18 MB. Do not subsample it: the fast band is a
# single sample and INT_TEMP_PROBE_EVERY would hide it.
TLAD_PROBE = -DINT_TEMP_PROBE=1

# probe    pyrolysis only. The cheapest case that carries the slow band, and
#          it gives 6 cycles in 20 s. Run this one first.
# probefl  the reduced case with a flame. The reference of the campaign.
# probem   moisture. This build stopped at t = 5.94 s with the sliver SIGFPE
#          of the TG runaway. Restart it serially from a t = 5 dump and the
#          probe reports SG_max and TG_min through the failure, in about
#          3 minutes. This is the fastest decisive test of the three.
test-probe.s:       CFLAGS += $(TLAD_PYRO) $(TLAD_PROBE)
test-probe.tst:     CFLAGS += $(TLAD_PYRO) $(TLAD_PROBE)

test-probefl.s:     CFLAGS += $(TLAD_PROBE)
test-probefl.tst:   CFLAGS += $(TLAD_PROBE)

test-probem.s:      CFLAGS += -DMOISTURE=1 $(TLAD_PROBE)
test-probem.tst:    CFLAGS += -DMOISTURE=1 $(TLAD_PROBE)

# The interface conductance on the diagonal of the diffusion operator. It
# keeps the probe on, so the same intres.dat columns compare the two builds
# line by line. The pair to run first is probem against robinm: probem stops
# at t = 5.94 s, so robinm must pass that time.
#
# Caution: on its own the conductance DESTROYS heat while the field moves.
# The cell keeps only SMAX/S of the interface heat, which was one per cent in
# a measured run: Tsurf fell 20 to 25 K and mdot 25 per cent, on every ladder.
# INT_TEMP_ROBIN_DEBT carries the withheld heat to the next step and removes
# that loss. It is 1 by default. Read src/int-temperature-robin.h.
TLAD_ROBIN = -DINT_TEMP_ROBIN=1 $(TLAD_PROBE)

# The lossy scheme, for the A/B comparison only. Do not run physics with it.
TLAD_ROBINLOSSY = -DINT_TEMP_ROBIN_DEBT=0 $(TLAD_ROBIN)

test-robindb.s:     CFLAGS += $(TLAD_ROBIN)
test-robindb.tst:   CFLAGS += $(TLAD_ROBIN)

test-robinls.s:     CFLAGS += $(TLAD_ROBINLOSSY)
test-robinls.tst:   CFLAGS += $(TLAD_ROBINLOSSY)

test-fullrobinls.s:   CFLAGS += $(TLAD_FULL) $(TLAD_ROBINLOSSY)
test-fullrobinls.tst: CFLAGS += $(TLAD_FULL) $(TLAD_ROBINLOSSY)

# The adapt-threshold ladder, at a FIXED maxlevel of 10. The level ladder
# (l9/base/l11) changes the finest cell AND the share of the domain that
# reaches it. These targets change only the second. If the 0.15 to 1.5 Hz band
# of mdot falls here too, the adapt criterion carries the error and the fix is
# cheap. If it stays flat, only maxlevel helps.
#
# Read the amplitude the same way for all of them, over a window of at least
# 19 s. Shorter windows hold about 10 bins in the band and scatter by 4x.
#
# They carry the full ladder. Two reasons. The band is 7.43 per cent there
# against 2.85 in the bare build, so the signal is 2.6 times larger against
# the same run-to-run scatter. And test-fullprobe is already the matched
# reference: it is the same flags with the thresholds at their defaults, and
# it measures 7.43 per cent, the same as test-full, so the probe is inert.
#
# Caution: compare them against test-fullprobe, NOT against test-base. The
# full ladder changes Da, moisture, gravity, shape, emissivity and transport
# all at once.
# The A/B for the temperature-solve tolerance. One binary difference: a tiny
# INT_TEMP_TOL_K makes max(TOLERANCE, theta/dt*K) collapse back to the
# inherited TOLERANCE, so test-toloff reproduces the old behaviour AND still
# writes tsolve.dat. Compare the iteration columns, which are deterministic;
# do not compare wall clock without several alternating runs.
test-toloff.s:    CFLAGS += -DINT_TEMP_TOL_K=1e-30 $(TLAD_FULL) $(TLAD_PROBE)
test-toloff.tst:  CFLAGS += -DINT_TEMP_TOL_K=1e-30 $(TLAD_FULL) $(TLAD_PROBE)

TLAD_ATOL = $(TLAD_FULL) $(TLAD_PROBE)

test-atol2.s:    CFLAGS += -DADAPT_T_TOL=2.5 -DADAPT_O_TOL=5.e-3 $(TLAD_ATOL)
test-atol2.tst:  CFLAGS += -DADAPT_T_TOL=2.5 -DADAPT_O_TOL=5.e-3 $(TLAD_ATOL)

test-atol4.s:    CFLAGS += -DADAPT_T_TOL=1.25 -DADAPT_O_TOL=2.5e-3 $(TLAD_ATOL)
test-atol4.tst:  CFLAGS += -DADAPT_T_TOL=1.25 -DADAPT_O_TOL=2.5e-3 $(TLAD_ATOL)

# The SMAX ladder. The residual error of the conductance is a LAG of about
# S/SMAX steps in the cut cells, so it falls as SMAX rises. s20 caps the
# exchange number at 20, well under the 100 that crashed; sinf never fires
# and must reproduce the no-robin run exactly.
test-fullrobin20.s:   CFLAGS += $(TLAD_FULL) -DINT_TEMP_ROBIN_SMAX=20. $(TLAD_ROBIN)
test-fullrobin20.tst: CFLAGS += $(TLAD_FULL) -DINT_TEMP_ROBIN_SMAX=20. $(TLAD_ROBIN)

test-fullrobininf.s:   CFLAGS += $(TLAD_FULL) -DINT_TEMP_ROBIN_SMAX=1e30 $(TLAD_ROBIN)
test-fullrobininf.tst: CFLAGS += $(TLAD_FULL) -DINT_TEMP_ROBIN_SMAX=1e30 $(TLAD_ROBIN)

test-robin.s:       CFLAGS += $(TLAD_PYRO) $(TLAD_ROBIN)
test-robin.tst:     CFLAGS += $(TLAD_PYRO) $(TLAD_ROBIN)

test-robinfl.s:     CFLAGS += $(TLAD_ROBIN)
test-robinfl.tst:   CFLAGS += $(TLAD_ROBIN)

test-robinm.s:      CFLAGS += -DMOISTURE=1 $(TLAD_ROBIN)
test-robinm.tst:    CFLAGS += -DMOISTURE=1 $(TLAD_ROBIN)

# The Picard loop on the interface temperature. It repeats the pair of
# temperature solves until TInt stops moving, so the interface condition
# becomes implicit instead of lagged. It keeps the probe and the conductance
# on, so intres.dat compares line by line with a robin build.
#
# Read the header of src/int-temperature-picard.h for how to read picard.dat.
#
# The loop and the conductance fix different failures. The conductance bounds
# the response of a cut cell to the interface source, so it stops the sliver
# runaway. The loop removes the lag of TInt. Measured 2026-09-04, pre-flame:
# rel_max falls from 3.8e-3 to 7e-16, at 2 passes median and 12.7 % of wall
# clock. test-picardnr drops the conductance, to test whether the loop needs
# it; in that window it does not.
#
# INT_TEMP_PICARD_MAXITER=0 must reproduce the robin build bit for bit. That
# is the inertness control; run it before you quote any picard result.
# ---------------------------------------------------------------------------
# The interface value as a Dirichlet condition of the solve, added 2026-09-08.
# See src/plicbc.h and src/int-temperature-vofbc.h.
#
# This needs the two patches in src/basilisk-patches/ applied to $BASILISK.
# Check first, or the build fails on a missing `flux` argument:
#
#   grep -c 'flux = flux' $(BASILISK)/diffusion.h     # must print 1
#
# Run them in this order. Each answers one question and the later ones are
# pointless if an earlier one fails.
#
#   1. test-vofbcm   the crash. test-probem dies at t = 5.94 s. This build
#                    must pass it with no INT_TEMP_ROBIN and no SMAX. It needs
#                    about 7 s of the 40, so it is the cheapest decisive run.
#   2. test-vofbc    against test-probe, on the pyrolysis-only case. Does the
#                    answer move where there was no pathology? Read rel_max of
#                    intres.dat and Tsurf of expansion.dat.
#   3. test-vofbcnp  the same build with the Picard loop OFF. The hook fixes
#                    the lag of the NEIGHBOUR values; Picard fixes the lag of
#                    TInt itself. This rung measures whether the second still
#                    matters once the first is gone.
#   4. test-fullvofbc  only after 1 to 3. Compare at matched remaining mass,
#                    never at matched time.
#
# Caution: INT_TEMP_ROBIN is an #error alongside INT_TEMP_VOFBC. They both
# claim the interface term and only one can have it.
TLAD_VOFBC = -DINT_TEMP_VOFBC=1 -DINT_TEMP_PICARD=1 $(TLAD_PROBE)

test-vofbc.s:       CFLAGS += $(TLAD_PYRO) $(TLAD_VOFBC)
test-vofbc.tst:     CFLAGS += $(TLAD_PYRO) $(TLAD_VOFBC)

test-vofbcfl.s:     CFLAGS += $(TLAD_VOFBC)
test-vofbcfl.tst:   CFLAGS += $(TLAD_VOFBC)

test-vofbcbase.s:   CFLAGS += -DINT_TEMP_VOFBC=1 $(TLAD_PROBE)
test-vofbcbase.tst: CFLAGS += -DINT_TEMP_VOFBC=1 $(TLAD_PROBE)

test-vofbcm.s:      CFLAGS += -DMOISTURE=1 $(TLAD_VOFBC)
test-vofbcm.tst:    CFLAGS += -DMOISTURE=1 $(TLAD_VOFBC)

# The Picard loop off, everything else the same as test-vofbc.
test-vofbcnp.s:     CFLAGS += $(TLAD_PYRO) -DINT_TEMP_VOFBC=1 $(TLAD_PROBE)
test-vofbcnp.tst:   CFLAGS += $(TLAD_PYRO) -DINT_TEMP_VOFBC=1 $(TLAD_PROBE)

test-fullvofbc.s:   CFLAGS += $(TLAD_FULL) $(TLAD_VOFBC)
test-fullvofbc.tst: CFLAGS += $(TLAD_FULL) $(TLAD_VOFBC)

# The three transport flags with the coupling, and nothing else. test-vofbcm
# passes and test-fullvofbc stops at t = 10.34, and five ingredients separate
# them: Da 1e-12 -> 1e-10, gravity, the shape, the emissivity of Di Blasi, and
# these three flags. No other VOFBC run carries them, and they act on the same
# diffusion solves as INT_TEMP_VOFBC.
#
# This rung must pass t = 10.34. If it stops, the transport flags carry the
# failure of test-fullvofbc. If it runs, they do not, and the other four
# ingredients hold it.
#
# Caution: the transport flags do nothing before the flame lights. A paired
# run is bit for bit identical over 560 rows up to the ignition at t = 5.66.
# test-fullvofbc writes its first pf warning at line 5 of the log, which is
# the FIRST step, so a clean run here does not clear the transport flags of
# every charge; it clears them of the collapse. Read the first ten lines of
# the log of this run for pf warnings before you read the end of it.
test-vofbctrans.s:   CFLAGS += $(TLAD_TRANSPORT) $(TLAD_VOFBC)
test-vofbctrans.tst: CFLAGS += $(TLAD_TRANSPORT) $(TLAD_VOFBC)

TLAD_PICARD = -DINT_TEMP_PICARD=1 $(TLAD_ROBIN)

test-picard.s:      CFLAGS += $(TLAD_PYRO) $(TLAD_PICARD)
test-picard.tst:    CFLAGS += $(TLAD_PYRO) $(TLAD_PICARD)

test-picardfl.s:    CFLAGS += $(TLAD_PICARD)
test-picardfl.tst:  CFLAGS += $(TLAD_PICARD)

test-picardm.s:     CFLAGS += -DMOISTURE=1 $(TLAD_PICARD)
test-picardm.tst:   CFLAGS += -DMOISTURE=1 $(TLAD_PICARD)

# The loop without the conductance. Measured 2026-09-04 over 516 steps to
# t = 0.25: it converges just as well as test-picard, at 2 passes median.
# So the loop does NOT need the conductance in a window with no extreme
# sliver. Keep this case: the crash restart is what decides the question.
test-picardnr.s:    CFLAGS += -DINT_TEMP_PICARD=1 $(TLAD_PROBE)
test-picardnr.tst:  CFLAGS += -DINT_TEMP_PICARD=1 $(TLAD_PROBE)

# --- the dt-dependence pair, the payoff test of the Picard loop -------------
#
# The loop exists to remove the lag of TInt, which is first order in dt. If it
# works, the Tmax difference between a step and a half step must shrink.
# Measure the SAME pair with and without the loop and compare the two gaps:
#
#   Tmax(rdt1) - Tmax(rdt2)   the gap without the loop
#   Tmax(pdt1) - Tmax(pdt2)   the gap with it
#
# The recorded sensitivity is 16 to 29 K per halving. A loop that removes the
# lag must cut that gap. If the gap does not move, the lag of TInt is not what
# makes Tmax follow the step, and the next suspects are update_divergence(),
# which keeps the first-pass TInt, and the mass-diffusion enthalpy.
#
# Both values divide 0.01 exactly: 2e-4 = 0.01/50, 1e-4 = 0.01/100. CFLNUM=2
# stops the CFL from binding, so DT_VALUE sets every step.
#
# Caution: START ALL FOUR FROM THE SAME PLATEAU SNAPSHOT. A run from scratch
# is not constant in dt, because the CFL binds through the impulsive start and
# through ignition. Copy last-snapshot of a finished run into each of the four
# case directories before you submit. Verify column 2 of expansion.dat is
# constant before you quote any of them.
TLAD_DT1 = 2e-4
TLAD_DT2 = 1e-4

test-rdt1.s:        CFLAGS += $(TLAD_ROBIN) -DDT_VALUE=$(TLAD_DT1) -DCFLNUM=2
test-rdt1.tst:      CFLAGS += $(TLAD_ROBIN) -DDT_VALUE=$(TLAD_DT1) -DCFLNUM=2

test-rdt2.s:        CFLAGS += $(TLAD_ROBIN) -DDT_VALUE=$(TLAD_DT2) -DCFLNUM=2
test-rdt2.tst:      CFLAGS += $(TLAD_ROBIN) -DDT_VALUE=$(TLAD_DT2) -DCFLNUM=2

test-pdt1.s:        CFLAGS += $(TLAD_PICARD) -DDT_VALUE=$(TLAD_DT1) -DCFLNUM=2
test-pdt1.tst:      CFLAGS += $(TLAD_PICARD) -DDT_VALUE=$(TLAD_DT1) -DCFLNUM=2

test-pdt2.s:        CFLAGS += $(TLAD_PICARD) -DDT_VALUE=$(TLAD_DT2) -DCFLNUM=2
test-pdt2.tst:      CFLAGS += $(TLAD_PICARD) -DDT_VALUE=$(TLAD_DT2) -DCFLNUM=2

# ---------------------------------------------------------------------------
# The full case with the modifications, added 2026-09-04.
#
# WHY THE FULL CASE. Measured at matched remaining mass (Ms/Ms0 0.87 to 0.55,
# band 0.10 to 0.45 Hz, restart seams removed), the full case is about twice
# the reduced case on every instrument that the campaign quotes:
#
#   case            mdot pp   Tsurf pp   T2mm..T15mm pp
#   test-base         5.0 %     4.66 K   8.9 to 10.6 K
#   test-full         7.6 %    11.32 K   12.4 to 19.5 K
#   test-transport   10.6 %    12.29 K   27.1 to 40.2 K
#
# So the transport rung carries most of the gap, and the full case is where a
# fix has to prove itself. Caution: the H2O weight of the probes follows
# MOLAR_DIFFUSION, so only transport and full weight them as the full case
# does. Compare the probe columns of these five cases with test-full and
# test-transport only, never with test-base.
#
# THE REFERENCE. test-fullprobe is test-full plus the probe. The probe changes
# no field: test-probefl and test-base give 65.0 against 65.6 per cent on
# mdot over t = 2.0 to 6.4 s. So test-fullprobe reproduces test-full and adds
# intres.dat. Run it FIRST. Nothing below can be read without it, because no
# full-case run has ever written SG_max, and the cut-cell over-drive that
# test-fullrobin removes has only been measured in the reduced case
# (test-probe: SG_max mean 133, 57 of 95 interface cells above 1).
#
# THE TWO MODIFICATIONS, and what each one tests.
#
#   fullrobin   The conductance on the diagonal. It caps the effective
#               exchange number of a cut cell at INT_TEMP_ROBIN_SMAX = 1.
#               A cell already under the limit gets K = 0 and is untouched bit
#               for bit, and the two added terms cancel at steady state, so
#               the fixed point does not move. A change in the slow band is
#               therefore attributable to the cut-cell over-drive alone.
#               This is the run to make after fullprobe.
#
#   fullpicard  The conductance plus the Picard loop on TInt. Expect this to
#               be WORSE than fullrobin, not better. Measured on the moisture
#               pair over t = 5.05 to 10.85 s: test-robinm 10.4 per cent
#               against test-picardm 23.5 per cent. The loop drives rel_max to
#               7e-16 and still raises the slow band. Run it to confirm that
#               the lag of TInt is not the source, and do not expect a fix.
#
# REFUTED, do not add a zeta rung. The gas-source campaign already tested
# ZETA_CONST against ZETA_REACTION on the full configuration, branched from
# one ignition snapshot, both to t = 20. Window t = 13 to 20, band 0.10 to
# 0.45 Hz, from ~/temp/gas-source/gas-source-{inst,zeta}/probe.dat:
#
#   case   zeta        mdot pp   Tp2 pp   Tmax pp   r_flame pp   omega_max pp
#   inst   6 REACTION   18.7 %   35.4 K    17.5 K       219 um         18.5 %
#   zeta   2 CONST      14.3 %   50.8 K    19.9 K       242 um         34.7 %
#
# ZETA_CONST does not remove the slow band, and it makes the probe
# temperature and the peak reaction rate worse. The global maximum of
# set_zeta is not the source. gas-source.c keeps a ZETA_POLICY flag if the
# question ever needs re-opening; test.c does not need one.
#
# ORDER: fullprobe, then fullrobin, then fullpicard.
#
# Caution: the full case is the slowest of the ladder. test-full reached
# t = 26.1 s where test-base reached 36.5 s. The matched-mass window needs
# Ms/Ms0 down to 0.55, which test-full reaches at t = 24.9 s. Give these runs
# the wall clock to pass t = 25, or the window holds under two cycles and the
# band amplitude is not measurable.
test-fullprobe.s:   CFLAGS += $(TLAD_FULL) $(TLAD_PROBE)
test-fullprobe.tst: CFLAGS += $(TLAD_FULL) $(TLAD_PROBE)

test-fullrobin.s:   CFLAGS += $(TLAD_FULL) $(TLAD_ROBIN)
test-fullrobin.tst: CFLAGS += $(TLAD_FULL) $(TLAD_ROBIN)

test-fullpicard.s:  CFLAGS += $(TLAD_FULL) $(TLAD_PICARD)
test-fullpicard.tst:CFLAGS += $(TLAD_FULL) $(TLAD_PICARD)


# ---------------------------------------------------------------------------
# Transport fix block, added 2026-09-07. See src/multicomponent-varprop.h.
#
# WHAT THE LADDER MEASURED. The transport rung doubles both bands. On mdot the
# slow band goes up 2.1 times and the fast band 2.5 times, so no new mode
# appears and the whole front spectrum gets louder. The front geometry does
# not move: the reaction zone is 0.260 mm against 0.258 mm, and the Arrhenius
# gain is 4.49 against 4.59 per cent per K. But the mean step falls from
# 4.24e-4 to 3.68e-4 s, and its coefficient of variation doubles from 8.5 to
# 15.8 per cent. A doubled scatter of the step is the signature of an explicit
# term that fights a stability limit that the solver does not know about.
#
# WHAT THE CODE DOES. Three properties of the three flags produce it. Each one
# now has an off switch in src/multicomponent-varprop.h, and each switch is
# independent of the other two.
#
#   CORRECTIVE_CFL      MOLAR_DIFFUSION and FICK_CORRECTED move each species
#                       with an extra velocity u_c = phic/rho. MW_mix changes
#                       by a factor near 2 across the front, so
#                       D d(ln MW_mix)/dx reaches 0.2 m/s against an inflow of
#                       0.13 m/s. That transport is explicit, and no event
#                       limited the step by it, so the scheme ran above a
#                       Courant number of 1 at the front. 0 removes the limit.
#                       Default 0.5.
#
#   CORRECTIVE_LIMITER  tracer_fluxes() reads its second argument as a
#                       velocity. It builds un = dt*uf/(fm*Delta) and it
#                       removes the slope with the factor (1 - s*un). The code
#                       gave it a mass flux in kg/m2/s, so un was too small by
#                       the density and the slope kept its full size at any
#                       Courant number. The same call set YG.gradient = NULL,
#                       and gradients() reads NULL as the unlimited centred
#                       slope, not as no slope. 0 gives both back. Default 1.
#
#   MDE_INTERFACE       MASS_DIFFUSION_ENTHALPY was applied in full cells
#                       only (f > 1 - F_ERR and f < F_ERR). An interface cell
#                       got nothing, and the term is largest in the cells next
#                       to the front. So each cell switched a large explicit
#                       source on, then off, then on again as the front passed
#                       it. 0 gives that gate back. Default 1.
#
#                       Caution: an interface cell does NOT use the centred
#                       stencil. It holds a phase average, and that average
#                       belongs to the centroid of the phase, not to the
#                       centre of the cell, so a centred difference there is
#                       first order at best. The interface cells use
#                       ebmgrad() and the interface value, in the direction of
#                       the normal, as int-temperature.h does for TS and TG.
#                       The full cells keep the centred stencil.
#
# Caution: none of the three switches changes the mass balance. The corrective
# flux keeps the same expression, because the code divides by the face density
# before tracer_fluxes() and multiplies the flux back after it. The heat
# source keeps the volume weight that theta1 and theta2 carry, so a full cell
# holds the previous value. Only the reconstruction, the Courant number and
# the interface cells change.
#
# THE CONTROL. test-transportold is the previous code: all three switches off.
# Read every case of this block against it, never against test-base, because
# test-base carries no transport flag at all.
#
# ORDER: transportold, then transport, then the three that remove one fix
# each, then fullold. Give fullold the wall clock of a full case: test-full
# reached t = 26.1 s where test-base reached 36.5 s, and the matched-mass
# window needs t = 24.9 s.
#
#   make test-transportold.tst test-transport.tst
#   python3 ~/oscillation-report/analysis/slow_flicker.py test-transportold
#
# Caution: read the header line of each log before you quote a case. It prints
# what the build enabled, not what the case name promises.

TLAD_LEGACY = -DCORRECTIVE_CFL=0 -DCORRECTIVE_LIMITER=0 -DMDE_INTERFACE=0

# --- which of the three flags amplifies ------------------------------------
# TLAD_TRANSPORT carries all three flags together, so the ladder has never
# said which one carries the factor 2. These three cases say it. Each one adds
# ONE flag to test-base. Compare them with test-base and with test-transport.
#
# Caution: the H2O weight of the probe columns follows MOLAR_DIFFUSION. Only
# test-molar and test-transport weight them as test-transport does, so compare
# the probe columns of test-fick and test-mde with test-base alone.

test-molar.s:         CFLAGS += -DMOLAR_DIFFUSION
test-molar.tst:       CFLAGS += -DMOLAR_DIFFUSION

test-fick.s:          CFLAGS += -DFICK_CORRECTED
test-fick.tst:        CFLAGS += -DFICK_CORRECTED

test-mde.s:           CFLAGS += -DMASS_DIFFUSION_ENTHALPY
test-mde.tst:         CFLAGS += -DMASS_DIFFUSION_ENTHALPY

# --- what each fix is worth ------------------------------------------------
# All four carry the three transport flags. test-transport has every fix on.
# Each case below removes ONE fix, so the difference against test-transport is
# the value of that fix alone.

test-transportold.s:  CFLAGS += $(TLAD_TRANSPORT) $(TLAD_LEGACY)
test-transportold.tst:CFLAGS += $(TLAD_TRANSPORT) $(TLAD_LEGACY)

test-transportcfl.s:  CFLAGS += $(TLAD_TRANSPORT) -DCORRECTIVE_CFL=0
test-transportcfl.tst:CFLAGS += $(TLAD_TRANSPORT) -DCORRECTIVE_CFL=0

test-transportlim.s:  CFLAGS += $(TLAD_TRANSPORT) -DCORRECTIVE_LIMITER=0
test-transportlim.tst:CFLAGS += $(TLAD_TRANSPORT) -DCORRECTIVE_LIMITER=0

test-transportmde.s:  CFLAGS += $(TLAD_TRANSPORT) -DMDE_INTERFACE=0
test-transportmde.tst:CFLAGS += $(TLAD_TRANSPORT) -DMDE_INTERFACE=0

# --- the full case ---------------------------------------------------------
# The control for test-full and test-fullprobe. The full case is where a fix
# has to prove itself: it holds 7.6 per cent on mdot against 5.0 per cent for
# test-base, and the transport rung carries most of that gap.

test-fullold.s:       CFLAGS += $(TLAD_FULL) $(TLAD_LEGACY)
test-fullold.tst:     CFLAGS += $(TLAD_FULL) $(TLAD_LEGACY)

# The forced-step transport run, and its control.
#
# At DT=1e-4 the corrective Courant limit never fires: the measured corrective
# velocity is 0.333 m/s, which gives a limit of 2.35e-4 at CORRECTIVE_CFL=0.5,
# and DT is far below it. So CORRECTIVE_CFL is INERT in both rungs below, and
# the only fixes that act are CORRECTIVE_LIMITER and MDE_INTERFACE.
#
# Caution: do NOT read test-transrdt2 against the archived test-transport.
# That run was built 2026-09-04 00:53, before the transport fixes of 44ddaf6,
# so the difference would mix the fixes with the step. test-transoldrdt2 is
# the control: the same step, the same source, every fix off. The pair gives
# the step alone in the old code, which is the question the campaign asks.

test-transrdt2.s:     CFLAGS += $(TLAD_TRANSPORT) -DDT_VALUE=$(TLAD_DT2) -DCFLNUM=2
test-transrdt2.tst:   CFLAGS += $(TLAD_TRANSPORT) -DDT_VALUE=$(TLAD_DT2) -DCFLNUM=2

test-transoldrdt2.s:  CFLAGS += $(TLAD_TRANSPORT) $(TLAD_LEGACY) -DDT_VALUE=$(TLAD_DT2) -DCFLNUM=2
test-transoldrdt2.tst:CFLAGS += $(TLAD_TRANSPORT) $(TLAD_LEGACY) -DDT_VALUE=$(TLAD_DT2) -DCFLNUM=2

test-l11full.s:   CFLAGS += $(TLAD_FULL) $(TLAD_ROBIN) -DMAXLEVEL_VALUE=11
test-l11full.tst: CFLAGS += $(TLAD_FULL) $(TLAD_ROBIN) -DMAXLEVEL_VALUE=11

TLAD_TGDBG = -DTG_PROBE -DTG_PROBE_TMIN=250. -DTG_PROBE_MAX=2000 \
             -DPROPERTIES_TRAP=0 -DPHASE_AWARE_PROPERTIES=1 \
             -DGAS_STATE_FALLBACK_DEBUG

test-l12full.s:   CFLAGS += $(TLAD_FULL) $(TLAD_ROBIN) -DMAXLEVEL_VALUE=12 $(TLAD_TGDBG)
test-l12full.tst: CFLAGS += $(TLAD_FULL) $(TLAD_ROBIN) -DMAXLEVEL_VALUE=12 $(TLAD_TGDBG)
