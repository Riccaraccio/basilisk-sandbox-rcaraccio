include $(BASILISK)/Makefile.defs

CFLAGS += -O2 -disable-dimensions -I../src  -I$(HOME)/basilisk/OpenSMOKEppInterface/src  

marangoni-deff.tst: CFLAGS+= -I$(HOME)/basilisk/basilisk-sandbox-ecipriano/src
marangoni-sweep.tst: CFLAGS+= -I$(HOME)/basilisk/basilisk-sandbox-ecipriano/src
marangoni-point.tst: CFLAGS+= -I$(HOME)/basilisk/basilisk-sandbox-ecipriano/src

# Ask make for a variable, so a script can reuse the exact include paths
# of this Makefile instead of repeating them:  make print-CFLAGS
.PHONY: print-%
print-%:
	@echo '$($*)'
