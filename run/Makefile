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

TLAD_CASES = base da moist grav transport emis shape full l11 \
             dtfix dtvar pyro uin2 uin05 noheat l9 \
             probe probefl probem \
             robin robinfl robinm \
             picard picardfl picardm picardnr \
             rdt1 rdt2 pdt1 pdt2 \
             fullprobe fullrobin fullpicard \
             robindb robinls fullrobinls fullrobin20 fullrobininf \
             atol2 atol4

TLAD_LINKS = $(TLAD_CASES:%=test-%.c)

test-%.c: test.c
	ln -sf $< $@

.PHONY: test-links test-unlink
test-links: $(TLAD_LINKS)
test-unlink:
	rm -f $(TLAD_LINKS)

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
# so the CFL binds through the impulsive start even at CFLNUM=2. Ignition does
# the same: the peak velocity reaches about 2 m/s there, and a constant step of
# 2e-4 would put the true CFL near 2.8, which is unstable.
#
# So BRANCH BOTH RUNS FROM A PLATEAU SNAPSHOT. On the plateau the peak velocity
# is about 0.28 m/s, the CFL-limited step is about 1.1e-3 at CFLNUM=2, and
# DT_VALUE binds with room to spare. `event dump (t = 1; t += 1)` already
# writes `last-snapshot`, and `event init` restores it, so copy the file of a
# finished run into the case directory before you start.
TLAD_DTFIX     = 2e-4
TLAD_DTVAR_CFL = 0.35

test-dtfix.s:       CFLAGS += -DDT_VALUE=$(TLAD_DTFIX) -DCFLNUM=2
test-dtfix.tst:     CFLAGS += -DDT_VALUE=$(TLAD_DTFIX) -DCFLNUM=2

test-dtvar.s:       CFLAGS += -DCFLNUM=$(TLAD_DTVAR_CFL)
test-dtvar.tst:     CFLAGS += -DCFLNUM=$(TLAD_DTVAR_CFL)

# --- (b) what closes the loop ----------------------------------------------
# test-pyro is the reference of this group and reproduces ~/temp/10-fatehi.
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
