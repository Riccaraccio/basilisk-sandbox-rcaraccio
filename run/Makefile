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

# Fixed step. The `stability` event of vof.h sets CFL to 0.5 at each step
# when it is larger, so CFLNUM=2 runs as CFL 0.5. DTMAX sets the step only
# where the step at CFL 0.5 is larger than DTMAX. The measured window above
# (2.52e-4 at CFL 0.5) is larger, but check column 2 of probe.dat.
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
             constprop constpropnoshrink noexp noexpnoshrink pyrotsoff pyronoshrinktsoff l9 \
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
             frozen frozenloose frozenfull \
             pyronouinnoshrink pyronouinnoshrinkl9 ytol \
             l12 l11full l12full \
             dt2 dt4 ttol10 ttol2p5 ttol1p25 otol5e3 \
             omegac omegau

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

# Moisture, compared against test-base. The build of 2026-09-04 stopped at
# t = 5.94 s with the TG runaway in a sliver cell. The target has no interface
# treatment on purpose, so that the pair with test-base stays clean. Since
# 34ee613 (2026-09-08) TG_FGMIN acts in every build. Re-run it as it is. If it
# stops at t = 5.94 s again, TG_FGMIN alone does not hold the sliver.
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

# The grid test, one more level. Without this target the pattern rule builds
# test-l12 at the default level 10, and the run is a copy of test-base.
# The statistics go to OutputData-12. Check the header line: maxlevel=12.
test-l12.s:         CFLAGS += -DMAXLEVEL_VALUE=12
test-l12.tst:       CFLAGS += -DMAXLEVEL_VALUE=12

# A variant of the grid test that doesn't shrink the grid.
# NO_SHRINK is not a flag of test.c. The fixed interface is TLAD_NOSHRINK.
test-l11noshrink.s:   CFLAGS += -DMAXLEVEL_VALUE=11 $(TLAD_NOSHRINK)
test-l11noshrink.tst: CFLAGS += -DMAXLEVEL_VALUE=11 $(TLAD_NOSHRINK)

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
# CORRECTED 2026-09-22. CFLNUM=2 does not stop the CFL from binding. The
# `stability` event of vof.h (line 140) sets CFL to 0.5 at each step when it
# is larger, so CFLNUM=2 runs as CFL 0.5. DT_VALUE sets the step only where
# the step at CFL 0.5 is larger: u < 0.5*Delta/dt, thus 0.39 m/s at
# dt = 2e-4 on the finest level.
#
# Verify column 2 of expansion.dat is constant before you quote dtfix. If it
# is not, lower TLAD_DTFIX to the next divisor of 0.01.
#
# TLAD_DTVAR_CFL scales the CFL-limited step to the same mean. Measure the
# mean step of test-base over the same window and set
#   TLAD_DTVAR_CFL = 0.5 * TLAD_DTFIX / <mean dt of test-base>
# Use 0.5 and not the 0.8 of CFLNUM in test-base: vof.h sets CFL to 0.5 at
# each step, so test-base runs at 0.5. With a mean of 4.55e-4 the formula
# gives 0.5 * 2e-4 / 4.55e-4 = 0.22. The value 0.35 below came from the old
# formula with 0.8. It is not changed here: the user must choose the value.
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
# which is MORE conservative than the 0.5 at which test-base really runs
# (vof.h sets CFL to 0.5 when CFLNUM is larger), with a smaller step
# (3.30e-4 against 4.09e-4) and a lower residual (2.45 against 2.73) at the
# same instant, and it dies while test-base lives.
#
# UPDATED 2026-09-15. Both stops are the TG runaway in a sliver cell, and the
# two builds (2026-09-07 11:46) came before TG_FGMIN (34ee613, 2026-09-08),
# which now acts in every build. The pair is read against test-base, so it
# keeps NO interface treatment. Re-run both as they are.
#
# If they stop at t = 6.5 s again, TG_FGMIN does not hold the sliver. Then the
# fallback is INT_TEMP_ROBIN: test-rdt1 carries the same DT=2e-4 CFL=2 as
# test-dtfix plus TLAD_ROBIN, so test-dtfix + Robin is test-rdt1. The pair is
# then test-rdt1 against a test-dtvar with TLAD_ROBIN, not against test-base.
#
# Caution: test-rdt1, test-rdt2, test-pdt1 and test-pdt2 on disk were BUILT
# 2026-09-04 18:23, which is BEFORE the heat-loss fix of d2d5f2c
# (2026-09-06 16:59). They ran with INT_TEMP_ROBIN_SMAX=1 and no
# INT_TEMP_ROBIN_DEBT. A rebuild now gets SMAX=20 and the debt term. Re-run
# the block before you quote it against test-base.
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

# --- the timestep ladder, added 2026-09-16 -----------------------------------
# The effect of the step on T at 2 mm and on the mass loss, against test-base.
# Each rung divides the cap DT_VALUE by 2. The CFL number does not follow:
# vof.h sets CFL to 0.5 at each step when it is larger, so test-base runs at
# 0.5, not 0.8. From test-base to test-dt2 the CFL goes from 0.5 to 0.4, thus
# a factor of 0.8, not 0.5. From test-dt2 to test-dt4 it goes from 0.4 to 0.2.
# test-base runs near its cap (mean step 4.4e-4 against DT_VALUE 5e-4), so a
# smaller CFL alone does not halve the step: test-dtvar gives only 1.4, which
# is 0.5/0.35.
#
#   test-base   DT_VALUE 5e-4     CFLNUM 0.8
#   test-dt2    DT_VALUE 2.5e-4   CFLNUM 0.4
#   test-dt4    DT_VALUE 1.25e-4  CFLNUM 0.2
#
# The step still adapts, so the runs start from scratch and go through the
# ignition like test-base. Do not use CFLNUM=2 here: test-dtfix dies on the
# ignition ramp at t = 6.47. Both values divide 0.01 exactly.
# Cost: about 2 and 4 times the steps of test-base.
test-dt2.s:         CFLAGS += -DDT_VALUE=2.5e-4 -DCFLNUM=0.4
test-dt2.tst:       CFLAGS += -DDT_VALUE=2.5e-4 -DCFLNUM=0.4

test-dt4.s:         CFLAGS += -DDT_VALUE=1.25e-4 -DCFLNUM=0.2
test-dt4.tst:       CFLAGS += -DDT_VALUE=1.25e-4 -DCFLNUM=0.2

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

# Pyrolysis only: no oxidiser AND no gas kinetics. PYROLYSIS_ONLY removes the
# O2. It does not remove the gas reactions: since 3cbc1e9 (2026-08-26) no
# source file reads GAS_PHASE_REACTIONS, and the gas kinetics always run.
# The scheme holds three gas reactions that need no O2: TAR cracking
# (TAR => 3H2 + 5CO + CH4, nine moles out for one in), CH4 reforming and the
# water-gas shift. TURN_OFF_GAS_REACTIONS removes them, in the gas cells and
# in the pores.
#
# Caution: every pyro run in ~/test/solid up to 2026-09-11 ran WITH the gas
# kinetics. Rebuild a rung before you compare it with a run of this
# definition. Check the header line: PYRO=1 NOGASR=1.
TLAD_PYRO = -DPYROLYSIS_ONLY=1 -DTURN_OFF_GAS_REACTIONS=1

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

# No inflow, fixed interface, pyrolysis only. NO_INFLOW is an integer switch:
# the preprocessor rejects a floating constant in #if, so test.c cannot test
# UIN_VALUE there. Check the header line: Uin=0 ZETA=ZETA_SWELLING NOINFLOW=1.
test-pyronouinnoshrink.s:   CFLAGS += $(TLAD_PYRO) -DNO_INFLOW=1 $(TLAD_NOSHRINK)
test-pyronouinnoshrink.tst: CFLAGS += $(TLAD_PYRO) -DNO_INFLOW=1 $(TLAD_NOSHRINK)

test-pyronouinnoshrinkl9.s:   CFLAGS += $(TLAD_PYRO) -DNO_INFLOW=1 -DMAXLEVEL_VALUE=9 $(TLAD_NOSHRINK)
test-pyronouinnoshrinkl9.tst: CFLAGS += $(TLAD_PYRO) -DNO_INFLOW=1 -DMAXLEVEL_VALUE=9 $(TLAD_NOSHRINK)

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
# Caution: build the two rungs of each pair in the same batch, or the pair
# measures the source changes between the builds as well. The runs of
# 2026-09-11 were built before TURN_OFF_GAS_REACTIONS joined TLAD_PYRO, so
# they carry the gas kinetics. A rebuild now does not.
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

test-noexpnoshrink.s:       CFLAGS += -DNO_EXPANSION $(TLAD_NOSHRINK)
test-noexpnoshrink.tst:     CFLAGS += -DNO_EXPANSION $(TLAD_NOSHRINK)

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

# The tolerance of the gas species solves. YTOL_PROBE does not change the run:
# it solves each species again every YTOL_EVERY steps, against a reference and
# a set of trial tolerances, and writes ysolve.dat and yerr.dat. Read the
# header of src/species-tol-probe.h for the columns.
test-ytol.s:    CFLAGS += -DYTOL_PROBE=1 $(TLAD_FULL)
test-ytol.tst:  CFLAGS += -DYTOL_PROBE=1 $(TLAD_FULL)

# --- the adaptation tolerance, bare flags, added 2026-09-16 ------------------
# The effect of the refinement criterion on T at 2 mm, against test-base, at
# the same maxlevel. The defaults are ADAPT_T_TOL = 5 K and ADAPT_O_TOL = 1e-2
# (the O2 mass fraction). The T rungs change only the T threshold; otol5e3
# changes only the O2 threshold, so the two criteria stay separate.
#
#   test-ttol10     T 10 K      coarser than base
#   test-base       T 5 K
#   test-ttol2p5    T 2.5 K
#   test-ttol1p25   T 1.25 K
#   test-otol5e3    O2 5e-3     T at the default
#
# Caution: these are NOT test-atol2 and test-atol4. Those carry the full set
# and change both thresholds; their reference is test-fullprobe.
test-ttol10.s:      CFLAGS += -DADAPT_T_TOL=10.
test-ttol10.tst:    CFLAGS += -DADAPT_T_TOL=10.

test-ttol2p5.s:     CFLAGS += -DADAPT_T_TOL=2.5
test-ttol2p5.tst:   CFLAGS += -DADAPT_T_TOL=2.5

test-ttol1p25.s:    CFLAGS += -DADAPT_T_TOL=1.25
test-ttol1p25.tst:  CFLAGS += -DADAPT_T_TOL=1.25

test-otol5e3.s:     CFLAGS += -DADAPT_O_TOL=5.e-3
test-otol5e3.tst:   CFLAGS += -DADAPT_O_TOL=5.e-3

# The steady-release pair. Both rungs overwrite `omega` after the chemistry
# event, so the particle blows at a rate which does not oscillate. Read them
# against test-base over 12 to 30 s. If the slow band of the probe
# temperature survives, the loop does not close through the release rate.
#
# test-omegac rescales the field to a fixed mdot. It keeps the shape of the
# front and it leaves `zeta` unchanged, so it is the clean rung. Prefer it.
#
# test-omegau makes `omega` uniform. Under ZETA_REACTION that also sets
# `zeta` to 1 over the whole particle, so it changes the shrinkage as well.
# It is the control of the shape, not of the amplitude.
#
# Caution: neither rung conserves the mass of the solid. The ODE has already
# run when the override fires. Do not read column 2 of OutputData.
#
# The values are the mean of ~/test/new/test-base over 10 to 30 s. Override
# them on the command line:  make test-omegac.tst TLAD_OMEGA_MDOT=1.59e-6
#
# Caution: mdot holds still on the plateau, the rate per unit mass does not.
# Over 12 to 28 s the trend of mdot is -4.9 per cent and the trend of the
# mass-weighted omega is +76.5 per cent. So read test-omegac over 12 to 28 s,
# and read test-omegau over 18 to 22 s only. The table is in test.c.

TLAD_OMEGA_MDOT = 1.53e-6
TLAD_OMEGA_VAL  = 83.
TLAD_OMEGA_T0   = 10.

TLAD_OMEGAC = -DOMEGA_CONST=2 -DOMEGA_CONST_MDOT=$(TLAD_OMEGA_MDOT) \
              -DOMEGA_CONST_T0=$(TLAD_OMEGA_T0)
TLAD_OMEGAU = -DOMEGA_CONST=1 -DOMEGA_CONST_VALUE=$(TLAD_OMEGA_VAL) \
              -DOMEGA_CONST_T0=$(TLAD_OMEGA_T0)

test-omegac.s:      CFLAGS += $(TLAD_OMEGAC)
test-omegac.tst:    CFLAGS += $(TLAD_OMEGAC)

test-omegau.s:      CFLAGS += $(TLAD_OMEGAU)
test-omegau.tst:    CFLAGS += $(TLAD_OMEGAU)

test-omegaunoshrink.s:      CFLAGS += $(TLAD_OMEGAU) $(TLAD_NOSHRINK)
test-omegaunoshrink.tst:    CFLAGS += $(TLAD_OMEGAU) $(TLAD_NOSHRINK)

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
# does not stop the CFL from binding: vof.h sets CFL to 0.5 at each step, so
# DT_VALUE sets the step only where u < 0.5*Delta/dt.
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
# Caution: the table above comes from the code before the transport fixes of
# 44ddaf6 (2026-09-07). With the fixes the ranking changes. Measured
# 2026-09-14 at matched mass 0.9 to 0.5, band 0.10 to 0.45 Hz: test-transport
# gives 4.68 per cent on mdot, below test-base at 5.94 per cent, and
# test-transportold (fixes off) gives 8.67 per cent. The transport rung is no
# longer the loudest ingredient.
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
#               exchange number of a cut cell at INT_TEMP_ROBIN_SMAX, which is
#               20 since d2d5f2c (2026-09-06). The same commit added
#               INT_TEMP_ROBIN_DEBT, which carries the withheld heat to the
#               next step. The first runs used SMAX = 1 and no debt term, and
#               they lost heat. A cell already under the limit gets K = 0 and
#               is untouched bit for bit. A change in the slow band is
#               therefore attributable to the cut-cell over-drive alone.
#               This is the run to make after fullprobe.
#
#               It is also the L10 point of the grid ladder of the full case:
#               test-l11full and test-l12full carry the same flags.
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
# Caution: do NOT read test-transrdt2 against test-transportbkp. That run was
# built 2026-09-04 00:53, before the transport fixes of 44ddaf6, so the
# difference would mix the fixes with the step. The folder test-transport now
# holds the rebuild of 2026-09-07 14:56, with the fixes. test-transoldrdt2 is
# the control: the same step, the same source, every fix off. The pair gives
# the step alone in the old code, which is the question the campaign asks.

test-transrdt2.s:     CFLAGS += $(TLAD_TRANSPORT) -DDT_VALUE=$(TLAD_DT2) -DCFLNUM=2
test-transrdt2.tst:   CFLAGS += $(TLAD_TRANSPORT) -DDT_VALUE=$(TLAD_DT2) -DCFLNUM=2

test-transoldrdt2.s:  CFLAGS += $(TLAD_TRANSPORT) $(TLAD_LEGACY) -DDT_VALUE=$(TLAD_DT2) -DCFLNUM=2
test-transoldrdt2.tst:CFLAGS += $(TLAD_TRANSPORT) $(TLAD_LEGACY) -DDT_VALUE=$(TLAD_DT2) -DCFLNUM=2

# The grid ladder of the full case: test-fullrobin (L10), test-l11full (L11)
# and test-l12full (L12). All three carry TLAD_FULL and TLAD_ROBIN, so the
# level is the only difference. Do not use test-full as the L10 point: it has
# no conductance.
#
# Caution: the three must share every flag. If test-l12full stops at
# t = 0.02 s on 1./rhomix, add -DPHASE_AWARE_PROPERTIES=1 to all three
# targets, not to test-l12full alone.
test-l11full.s:   CFLAGS += $(TLAD_FULL) $(TLAD_ROBIN) -DMAXLEVEL_VALUE=11
test-l11full.tst: CFLAGS += $(TLAD_FULL) $(TLAD_ROBIN) -DMAXLEVEL_VALUE=11

test-l12full.s:   CFLAGS += $(TLAD_FULL) $(TLAD_ROBIN) -DMAXLEVEL_VALUE=12
test-l12full.tst: CFLAGS += $(TLAD_FULL) $(TLAD_ROBIN) -DMAXLEVEL_VALUE=12

# The debug switches of the gas state. For a diagnostic build only, append
# $(TLAD_TGDBG) to the two lines of the target, then remove it again. Do not
# pass CFLAGS on the command line: that replaces the flags of the target.
#
# Caution: PROPERTIES_TRAP=0 turns off the FPE trap over the properties loop,
# so a bad cell carries inf into the projection and the run continues. Never
# quote a result of such a build. The run of 2026-09-08 in ~/test/solid
# was built this way.
TLAD_TGDBG = -DTG_PROBE -DTG_PROBE_TMIN=250. -DTG_PROBE_MAX=2000 \
             -DPROPERTIES_TRAP=0 -DPHASE_AWARE_PROPERTIES=1 \
             -DGAS_STATE_FALLBACK_DEBUG


# --- the restart pairs of the time-level fixes, added 2026-09-22 ------------
# Each rung restores the snapshot of test-base at t = 15 s and runs to
# t = 21 s with ONE flag changed. The dt policy is the same as test-base
# (DT_VALUE 5e-4, CFLNUM 0.8, which vof.h clamps to 0.5).
#
#   test-sref     no flag: the reference. It has the same restart as the
#                 other rungs, so the transient of the restart cancels.
#   test-sdri     DRHODT_IMPLICIT=1     TL-2: drhodt from the solves
#   test-sgse     GAS_SOURCE_EXACT=1    TL-3: the log form of the expansion
#   test-sprops   PROPS_AFTER_SOLVES=1  TL-4: P^{n+1} in the momentum
#   test-subf2    GAS_UBF_ADVECTION=2   item 7: the ubf mask
#
# What to read, against test-sref:
#   - chemsplit.dat: Tmax_end, Tmid and the gap Tmax_chem - Tmax_pre.
#   - drhodtbudget.dat: in test-sdri, L1d/L1r must fall to about 0 in
#     region H. With the flag on, the probe compares the new drhodt against
#     the operator of the solve, so the ratio is the residual of the fix.
#   - shrinkbudget.dat: Cdiff (the drift of NEW-1) must not change.
#   - speciesclamp.dat: the clamp mass on the plateau (TL-5).
#   - OutputData and TemperatureProfile.dat: the mean Tmax and the slow band.
#     Caution: 6 s holds only 2 periods of the slow band at 0.35 Hz. Quote
#     only the ratios of the band against test-sref, not the absolute
#     values. For a firm result on the band, run to t = 27 s or more.
#
# Procedure:
#   1. Build every rung from the source of f50229d or later.
#   2. For each rung:
#        mkdir -p <rung> && cp <snapshot of t = 15 s> <rung>/last-snapshot
#   3. Run with the rank count of test-base (8).
#   4. Check that the log says "Restart file found!".
#
# restore() matches the fields by name. DRHODT_IMPLICIT adds 4 fields that
# the snapshot does not hold. They start at 0, and the code sets them again
# in each step, so the restart is valid.
#
# Caution: each rung overwrites its last-snapshot at t = 16 s. Keep the copy
# of t = 15 s somewhere else.
#
# Cost: each rung costs about the same as test-base over 6 s. test-sprops
# costs about 5 % more.
TLAD_SRST_TEND = 21.
TLAD_SRST      = -DTEND_VALUE=$(TLAD_SRST_TEND)

test-sref.s:        CFLAGS += $(TLAD_SRST)
test-sref.tst:      CFLAGS += $(TLAD_SRST)

test-sdri.s:        CFLAGS += $(TLAD_SRST) -DDRHODT_IMPLICIT=1
test-sdri.tst:      CFLAGS += $(TLAD_SRST) -DDRHODT_IMPLICIT=1

test-sgse.s:        CFLAGS += $(TLAD_SRST) -DGAS_SOURCE_EXACT=1
test-sgse.tst:      CFLAGS += $(TLAD_SRST) -DGAS_SOURCE_EXACT=1

test-sprops.s:      CFLAGS += $(TLAD_SRST) -DPROPS_AFTER_SOLVES=1
test-sprops.tst:    CFLAGS += $(TLAD_SRST) -DPROPS_AFTER_SOLVES=1

test-subf2.s:       CFLAGS += $(TLAD_SRST) -DGAS_UBF_ADVECTION=2
test-subf2.tst:     CFLAGS += $(TLAD_SRST) -DGAS_UBF_ADVECTION=2

# --- the Strang split of the gas chemistry (TL-1), added 2026-09-22 --------
# GAS_CHEMISTRY_STRANG=1 integrates the gas reactor over dt/2 before the
# transport and over dt/2 after the solves. The solid and the pore gas stay
# on the full dt. See the doc block of the flag in src/chemistry.h.
#
# Each rung restores the snapshot of test-base at t = 15 s and runs to
# t = 17 s, as the constant-dt ladder (test-sdt200/100/050).
#
#   test-sstr      the CFL dt (about 4e-4)    compare with test-sref
#   test-sstr200   DT_VALUE 2e-4              compare with test-sdt200
#   test-sstr100   DT_VALUE 1e-4              compare with test-sdt100
#   test-sstrg     as test-sstr, with FROZEN_CELL_GATE=1. Compare its cost
#                  with test-sstr. The gate tests each half with its own dt/2.
#
# CFLNUM 0.5 is the value that vof.h applies to test-sref too (it clamps 0.8
# to 0.5), so the CFL dt of test-sstr is comparable with test-sref.
#
# What to read, against the Lie rung of the same dt:
#   - chemsplit.dat: Tmax_end (col 6) is the state after the second half.
#     Tmax_tr (col 14) is the state after the transport, Q2 (col 15) the
#     heat of the second half. Qchem (col 13) holds only the first half.
#     Tmid (col 7) is not a predictor under the split.
#   - OutputData: Tmax, the gas temperatures, the solid mass.
#   - The log: the cost of the chemistry event. Expect about 1.6 to 1.9
#     times the gas chemistry of the Lie rung.
# The agent estimated Tmax_end at about 2305, 2320 and 2330 K for the three
# dt. The split stays first order in the full step, because the diffusion
# solves are backward Euler. The rungs measure the real gain.
#
# Procedure:
#   1. Build every rung from feccbb4 or later.
#   2. For each rung: mkdir -p <rung> && cp <snapshot of t = 15 s> <rung>/last-snapshot
#   3. Regenerate Makefile.tests on the login node before qsub.
#   4. Run on 8 ranks. Check that the log says "Restart file found!" and
#      that the "# ladder:" line says STRANG=1.
#
# Caution: each rung overwrites its last-snapshot at t = 16 s. Keep the copy
# of t = 15 s somewhere else.
TLAD_STRANG = -DGAS_CHEMISTRY_STRANG=1 -DCFLNUM=0.5 -DTEND_VALUE=17.

test-sstr.s:        CFLAGS += $(TLAD_STRANG)
test-sstr.tst:      CFLAGS += $(TLAD_STRANG)

test-sstr200.s:     CFLAGS += $(TLAD_STRANG) -DDT_VALUE=2e-4
test-sstr200.tst:   CFLAGS += $(TLAD_STRANG) -DDT_VALUE=2e-4

test-sstr100.s:     CFLAGS += $(TLAD_STRANG) -DDT_VALUE=1e-4
test-sstr100.tst:   CFLAGS += $(TLAD_STRANG) -DDT_VALUE=1e-4

test-sstrg.s:       CFLAGS += $(TLAD_STRANG) -DFROZEN_CELL_GATE=1
test-sstrg.tst:     CFLAGS += $(TLAD_STRANG) -DFROZEN_CELL_GATE=1

# --- the full runs of the base case and of the best flag set, 2026-09-22 ---
# Both runs start at t = 0 and run to TEND_VALUE (40 s by default), with the
# dt policy of test-base (DT_VALUE 5e-4, the CFL limit). They cover the
# ignition and the sliver phase, which the restart rungs from t = 15 s do
# not cover.
#
#   test-fbest       the best flag set of the time-level campaign, on the
#                    base case. Compare with test-base.
#   test-fbestfull   the same flag set on the full case (TLAD_FULL: the full
#                    transport, Da 1e-10, moisture, gravity, the shape and
#                    the emissivity of Di Blasi). Compare with test-full.
#
# The best flag set (TLAD_BEST):
#     GAS_SOURCE_EXACT=1      TL-3: the exact log form of the gas expansion.
#                             The default form expands 16 to 39 % too much.
#     DRHODT_IMPLICIT=1       TL-2: drhodt from the solves. At dt 4.5e-4 it
#                             gives the far field and mdot of dt 1e-4.
#     GAS_UBF_ADVECTION=2     item 7: the ubf mask. Mode 0 makes TG of the
#                             sliver cells 40 to 50 K deeper.
#     GAS_CHEMISTRY_STRANG=1  TL-1: the Strang split of the gas chemistry.
#                             Its restart ladder (test-sstr*) has not
#                             finished. Expected gain about +28 K of Tmax.
#     FROZEN_CELL_GATE=1      the gate lowers the cost of the two calls of
#                             the solver for each step. At its tolerance
#                             (1e-15) it does not change the result.
#
# Not in TLAD_BEST: PROPS_AFTER_SOLVES (no effect, +4.6 % cost) and
# FICK_CORRECTED (not tested alone). test-fbestfull gets FICK_CORRECTED from
# TLAD_TRANSPORT. The fixes that are on by default stay on in both runs, and
# shift_prod stays true.
#
# Caution: test-base and test-full ran on older code. They are not exact
# references for the fixes that are on by default since then.
#
# What to compare: Tmax and the gap of chemsplit.dat, the solid mass and the
# drift (shrinkbudget.dat), mdot and Tsurf (expansion.dat), the clamp
# (speciesclamp.dat), the slow band (OutputData, TemperatureProfile.dat), dt
# and the total cost. Under Strang, read chemsplit.dat as the doc block of
# the Strang rungs above says.
#
# Caution: the Strang split costs about 1.6 to 1.9 times the gas chemistry
# without the gate. The full case is the slowest case of the ladder, so
# give test-fbestfull the longest wall time. Do not restart these runs from
# the t = 15 s snapshot: remove any last-snapshot from the folder before
# the first submission.
#
# Regenerate Makefile.tests on the login node before qsub. Run on 8 ranks,
# as test-base.
TLAD_BEST = -DGAS_SOURCE_EXACT=1 -DDRHODT_IMPLICIT=1 -DGAS_UBF_ADVECTION=2 \
            -DGAS_CHEMISTRY_STRANG=1 -DFROZEN_CELL_GATE=1


test-fbest.s:       CFLAGS += $(TLAD_BEST)
test-fbest.tst:     CFLAGS += $(TLAD_BEST)

test-fbestfull.s:   CFLAGS += $(TLAD_FULL) $(TLAD_BEST)
test-fbestfull.tst: CFLAGS += $(TLAD_FULL) $(TLAD_BEST)

# --- the projection spikes, added 2026-09-23 --------------------------------
# The point probes show one-sample spikes of the velocity, of up to 8 % of
# the local speed. At a spike the residual of the projection (mgp.resa) is
# 2 to 20 times its median, and the grid often changes by more cells. The
# residual stays under TOLERANCE/dt^2, so the solver accepts the step. See
# the section "The projection spikes" in test.c.
#
# Each run restores the snapshot of test-base at t = 15 s and runs to
# t = 21 s, as the restart rungs above. Each run writes projstep.dat, one
# line for each step (PROJECTION_PROBE=1).
#
#   test-sprobe   the reference: test-sref with the probe
#   test-sptol    TOLERANCE 1e-7 instead of 1e-5. Tests the residual.
#   test-snit4    NITERMIN 4 instead of 2. A cheaper candidate fix.
#   test-fixgrid  run/test-fixgrid.c: no adapt. It refines a box around the
#                 particle (x -12 to +24 mm, y 0 to 12 mm) to level 10 one
#                 time, after the restore. Tests the adapt.
#
# What to read, against test-sprobe:
#   - probes.dat: the count of the one-sample spikes at the probes 0 and 2.
#   - projstep.dat: mgp_tol (col 6) against PROJ_TOLERANCE, dumax (col 14)
#     and its position, and nf + nc (col 12 + col 13) of the adapt.
#     A line normally holds the adapt of step i - 1 (adapt_i = i - 1).
#   - The cost: the cycles mgp_i (col 4) and the wall time.
#   - OutputData and expansion.dat: the slow band must not change. If it
#     changes, the spikes feed the slow band too.
# Decision:
#   - The spikes stop in test-sptol or test-snit4: the residual of the
#     projection is the cause. The cheaper of the two is the fix.
#   - The spikes stop in test-fixgrid only: the adapt is the cause, through the
#     prolongation of u and not through the residual.
#   - The spikes stay in all runs: neither is the cause.
#
# Procedure:
#   1. Build every run from the commit of this section or later.
#   2. For each run:
#        mkdir -p <run> && cp <snapshot of t = 15 s> <run>/last-snapshot
#   3. Regenerate Makefile.tests on the login node before qsub.
#   4. Run on 8 ranks, as test-sref. Check that the log says "Restart file
#      found!" and that the "# ladder:" line gives the PROJTOL, NITERMIN,
#      PROJPROBE and FIXEDGRID values of the run. The log of test-fixgrid
#      must also say "# fixed grid".
#
# Add projstep.dat to the list of files of sync-solid.sh.
#
# Cost: test-sprobe about test-sref. test-sptol and test-snit4 add projection
# cycles: expect a few per cent more. test-fixgrid about 2 to 2.5
# times test-sref.
# projstep.dat grows by about 300 bytes for each step (about 5 MB for 6 s).
TLAD_SPIKE = $(TLAD_SRST) -DPROJECTION_PROBE=1

test-sprobe.s:      CFLAGS += $(TLAD_SPIKE)
test-sprobe.tst:    CFLAGS += $(TLAD_SPIKE)

test-sptol.s:       CFLAGS += $(TLAD_SPIKE) -DPROJ_TOLERANCE=1e-7
test-sptol.tst:     CFLAGS += $(TLAD_SPIKE) -DPROJ_TOLERANCE=1e-7

test-snit4.s:       CFLAGS += $(TLAD_SPIKE) -DPROJ_NITERMIN=4
test-snit4.tst:     CFLAGS += $(TLAD_SPIKE) -DPROJ_NITERMIN=4

# test-fixgrid.c is a real source, not a link to test.c. This explicit rule
# with an empty recipe stops the pattern rule "test-%.c: test.c" from
# replacing the file with a link. Do not remove it, and do not add
# test-fixgrid to TLAD_CASES: "make test-unlink" removes those files.
# .gitignore has an exception for this file. The pattern rules give the pin.
test-fixgrid.c: ;

test-fixgrid.s:     CFLAGS += $(TLAD_SPIKE)
test-fixgrid.tst:   CFLAGS += $(TLAD_SPIKE)

# --- the pulses next to the interface, added 2026-09-24 ---------------------
# Result of the runs above: the residual does not cause the spikes.
# test-sptol and test-snit4 keep them. The spikes are pulses of about 4 ms,
# and each one starts in a gas cell next to the interface. The first
# test-fixgrid found no snapshot and ran from t = 0, so it did not test the
# adapt on 15 to 21 s.
#
# projstep.dat now has columns 27 to 49 (see the doc of projection_probe in
# test.c):
#   27-39  the cell of the largest |du| in the zone next to the interface,
#          with f, the change of f, gas_source, drhodt, the change of
#          div_source, porosity and the f of its 3x3 block
#   40-44  the cell of the largest change of div_source in that zone
#   45-49  the outlet cell (75, 15) mm: u and its change, with the sign
# Columns 1 to 26 do not change.
#
#   test-sprobe2         as test-sprobe, with the new columns. Restart at
#                        t = 15 s. Does a cut cell that toggles start the
#                        pulses? Is the outlet cell a mode of the solver?
#   test-sfixgrid        test-fixgrid from the snapshot of t = 15 s. Its
#                        source is a link to test-fixgrid.c. It stops at once
#                        if it finds no last-snapshot
#                        (FIXED_GRID_NEED_RESTART=1). Compare with
#                        test-sprobe2. Does the adapt cause the pulses?
#   test-omegaunoshrinkp test-omegaunoshrink with the probe, from t = 0 to
#                        t = 21 s. No snapshot of t = 18 s exists, because
#                        the dump overwrites last-snapshot every second. It
#                        looks for the source of the burst at 18.6 to 19.6 s.
#                        Caution: identical builds do not give identical
#                        trajectories, so the burst can come at another time
#                        or not come at all.
#
# What to read:
#   - A pulse: du_ni (col 27) large for about 10 steps. On its first line,
#     is |df_ni| (col 32) or |dds_ni| (col 35) large? Is the cell a cut cell
#     or the neighbour of one (ncut_ni, col 39)?
#   - The outlet: the sign of dux_out (col 47) from step to step.
#   - test-sfixgrid against test-sprobe2: the count of pulses on 15 to 21 s.
#
# Procedure:
#   1. Build from the commit of this section or later.
#   2. For test-sprobe2 and test-sfixgrid:
#        mkdir -p <run> && cp <snapshot of t = 15 s of test-base> <run>/last-snapshot
#      Do not reuse the folder test-fixgrid: its files hold t = 0 to 11 and
#      the new run would add to them.
#   3. For test-omegaunoshrinkp: an empty folder, no last-snapshot.
#   4. Regenerate Makefile.tests on the login node before qsub. Run on 8
#      ranks. Check "Restart file found!" in the log of the two restarts.
#
# Cost: test-sprobe2 about test-sref. test-sfixgrid about 2 times
# test-sprobe (3.1 s per step in test-fixgrid). test-omegaunoshrinkp about
# test-base to t = 21 s. projstep.dat grows by about 550 bytes for each step.

test-sfixgrid.c: test-fixgrid.c
	ln -sf $< $@

test-sprobe2.s:          CFLAGS += $(TLAD_SPIKE)
test-sprobe2.tst:        CFLAGS += $(TLAD_SPIKE)

test-sfixgrid.s:         CFLAGS += $(TLAD_SPIKE) -DFIXED_GRID_NEED_RESTART=1
test-sfixgrid.tst:       CFLAGS += $(TLAD_SPIKE) -DFIXED_GRID_NEED_RESTART=1

test-omegaunoshrinkp.s:   CFLAGS += $(TLAD_OMEGAU) $(TLAD_NOSHRINK) \
                                    -DPROJECTION_PROBE=1 -DTEND_VALUE=21.
test-omegaunoshrinkp.tst: CFLAGS += $(TLAD_OMEGAU) $(TLAD_NOSHRINK) \
                                    -DPROJECTION_PROBE=1 -DTEND_VALUE=21.

# ---------------------------------------------------------------------------
# Smear of the pellet corner under shrinkage (see shrink-corner.c).
#
# Level 11 of fatehi-combustion.c puts solid in gas cells at the corner by
# t = 4 s. Only the VOF advection with ubf and the adaptation change f, so
# each case below switches ONE of the suspects of that defect. Run them serial:
# the diagnostics are not written for MPI.
#
#   make shrink-corner-l11.tst shrink-corner-l10.tst ...
#   make shrink-corner-all
#
# Compare column 13 (err, the L1 distance to the exact shape) of debris.dat and
# the maps corner-<t>.txt. static is the control: no shrinkage, so any defect in
# that case comes from the adaptation alone.
#
# INIT_FIX: fatehi-combustion.c calls refine() in main(), but run() calls
# init_grid (N) again and frees that grid, so production starts at level 8.
# Every case keeps that start (INIT_FIX=0) except init, initstatic and
# noadapt, which refine before fraction(). Without adaptation a level 8 start
# never refines, so noadapt needs INIT_FIX=1.

SHCO_CASES = l11 l10 noshift const uniform noadapt static init initstatic
SHCO_LINKS = $(SHCO_CASES:%=shrink-corner-%.c)

shrink-corner-%.c: shrink-corner.c
	ln -sf $< $@

.PHONY: shrink-corner-links shrink-corner-all
shrink-corner-links: $(SHCO_LINKS)
shrink-corner-all: $(SHCO_CASES:%=shrink-corner-%.tst)

shrink-corner-l10.s:       CFLAGS += -DMAXLEVEL=10
shrink-corner-l10.tst:     CFLAGS += -DMAXLEVEL=10

shrink-corner-noshift.s:   CFLAGS += -DSHIFT_PROD=0
shrink-corner-noshift.tst: CFLAGS += -DSHIFT_PROD=0

shrink-corner-const.s:     CFLAGS += -DZETA_POLICY=ZETA_CONST
shrink-corner-const.tst:   CFLAGS += -DZETA_POLICY=ZETA_CONST

shrink-corner-uniform.s:   CFLAGS += -DOMEGA_SHAPE=0
shrink-corner-uniform.tst: CFLAGS += -DOMEGA_SHAPE=0

shrink-corner-noadapt.s:   CFLAGS += -DADAPT=0 -DINIT_FIX=1
shrink-corner-noadapt.tst: CFLAGS += -DADAPT=0 -DINIT_FIX=1

shrink-corner-static.s:    CFLAGS += -DMDOT=0
shrink-corner-static.tst:  CFLAGS += -DMDOT=0

shrink-corner-init.s:      CFLAGS += -DINIT_FIX=1
shrink-corner-init.tst:    CFLAGS += -DINIT_FIX=1

shrink-corner-initstatic.s:   CFLAGS += -DINIT_FIX=1 -DMDOT=0
shrink-corner-initstatic.tst: CFLAGS += -DINIT_FIX=1 -DMDOT=0
