include $(BASILISK)/Makefile.defs

CFLAGS += -O2 -disable-dimensions -I../src  -I$(HOME)/basilisk/OpenSMOKEppInterface/src  

# The tests that link OpenSMOKE need its shared library at run time.
export LD_LIBRARY_PATH := $(HOME)/basilisk/OpenSMOKEppInterface/build:$(LD_LIBRARY_PATH)
