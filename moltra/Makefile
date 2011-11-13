MODS= -DMOD_DFT
MODS_IBM=,,-DMOD_DFT,,,
#
# VERSION : $Revision: 1.5.6.1 $
# DATE    : $Date: 2004/09/01 10:19:11 $
# FILE    : moltra/Makefile.in
#

include ../Makefile.config

F_OBJS =  trainp.o trapam.o tradr1.o tradr2.o tradr3.o traone.o \
          traout.o traprp.o tradr4.o tradr5.o

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

