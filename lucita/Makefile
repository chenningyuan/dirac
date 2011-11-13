MODS= -DMOD_DFT
MODS_IBM=,,-DMOD_DFT,,,
#
# VERSION : $Revision: 1.3.18.1 $
# DATE    : $Date: 2004/09/01 10:19:10 $
# FILE    : lucita/Makefile.in
#

include ../Makefile.config

F_OBJS=ciinfo.o  density.o  diagonal.o  diagonalize.o  doubledir.o  \
       fock.o  interface.o  io.o  lucia_cc.o  lucia_icci.o  memory.o \
       program.o  property.o  sigma.o  strings.o  transform.o  utils.o  \
       vecmat.o

F_CPPS = $(F_OBJS:.o=.i)

pre:	$(F_CPPS)

ftnchek:pre
	ftnchek $(CHEKFLAGS) $(F_CPPS) > ftnchek.out

dirac.x:
	cd .. ; $(MAKE) dirac.x

dirmpi.x:
	cd .. ; $(MAKE) dirmpi.x

all:	$(F_OBJS)

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i


depend:
	makedepend -- $(DEPEND_FLAGS) -- $(F_OBJS:.o=.F)
# DO NOT DELETE THIS LINE -- make depend depends on it.


