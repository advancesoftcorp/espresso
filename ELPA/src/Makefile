# Makefile for Modules

include ../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG).

ELPAMODS = \
elpa1.o \
elpa2.o \
elpa2_kernels.o \
blockedQR.o

all : libelpa.a

libelpa.a: $(ELPAMODS)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@
	mv $@ ../    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *.F90 *.mod *.L

include make.depend
