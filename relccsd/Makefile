MODS= -DMOD_DFT
MODS_IBM=,,-DMOD_DFT,,,
#
# VERSION : $Revision: 1.9.12.2 $
# DATE    : $Date: 2004/09/01 10:19:13 $
# FILE    : relccsd/Makefile.in
#

include ../Makefile.config

F_OBJS =ccdriv.o ccints.o cceqns.o  ccintm.o ccgetv.o ccsort.o \
      ccgenp.o ccrpa.o  ccxblas.o ccpam.o ccgrad.o ccmain.o ccfspc.o waio.o

F_CPPS = $(F_OBJS:.o=.i)

pre:	$(F_CPPS)

ftnchek:pre
	ftnchek $(CHEKFLAGS) $(F_CPPS) > ftnchek.out

all:	$(F_OBJS) $(C_OBJS)

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i

dirac.x:
	cd .. ; $(MAKE) dirac.x

dirmpi.x:
	cd .. ; $(MAKE) dirmpi.x

depend:
	makedepend -- $(DEPEND_FLAGS) -- $(F_OBJS:.o=.F) $(C_OBJS:.o=.c)
# DO NOT DELETE THIS LINE -- make depend depends on it.

