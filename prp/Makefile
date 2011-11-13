MODS= -DMOD_DFT
MODS_IBM=,,-DMOD_DFT,,,
#
# VERSION : $Revision: 1.12.2.2 $
# DATE    : $Date: 2004/09/01 10:19:12 $
# FILE    : prp/Makefile.in
#

include ../Makefile.config

F_OBJS =paminp.o pamprp.o pamexp.o pamgrd.o pamogrd.o pamfck.o pamset.o \
	pamxlr.o pamrvc.o pamrsp.o rspmpg.o mp2lag.o pamxqr.o

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

