MODS= -DMOD_DFT
MODS_IBM=,,-DMOD_DFT,,,
#
# VERSION : $Revision: 1.7.18.1 $
# DATE    : $Date: 2004/09/01 10:19:08 $
# FILE    : dirrci/Makefile.in
#

include ../Makefile.config

F_OBJS = cidirr.o cigenp.o cigosc.o cimain.o

F_CPPS = $(F_OBJS:.o=.i)

pre:	$(F_CPPS)

ftnchek:pre
	ftnchek $(CHEKFLAGS) $(F_CPPS) > ftnchek.out

all:	$(F_OBJS)

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i

dirac.x:
	cd .. ; $(MAKE) dirac.x

dirmpi.x:
	cd .. ; $(MAKE) dirmpi.x

depend:
	makedepend -- $(DEPEND_FLAGS) -- $(F_OBJS:.o=.F) 
# DO NOT DELETE THIS LINE -- make depend depends on it.

