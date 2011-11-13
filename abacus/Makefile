MODS= -DMOD_DFT
MODS_IBM=,,-DMOD_DFT,,,
#
# VERSION : $Revision: 1.10.4.1 $
# DATE    : $Date: 2004/09/01 10:19:06 $
# FILE    : abacus/Makefile.in
#

include ../Makefile.config

F_OBJS = her1ave.o    her1car.o    her1drv.o    her1int.o\
         her1odc.o    her1out.o    her1pro.o    her1sym.o\
         her2ave.o    her2dir.o    her2drv.o    her2el1.o\
         her2el2.o    her2fck.o    her2her.o    her2odc.o\
         her2out.o    her2sym.o    herbas.o     herdirect.o\
         herdrv.o     hergam.o     hergeo.o     hergp.o\
         herrdn.o     hersol.o     hersym.o     hertst.o\
         her2gab.o    hermem.o     her1cat.o \
         abanuc.o     abatro.o     abaopt.o abaortvec.o \
         abarint.o abavrml.o abawalk.o abaop2.o

SLAVE_OBJS = herpar.o

F_CPPS = $(F_OBJS:.o=.i) $(SLAVE_OBJS:.o=.i)

pre:	$(F_CPPS)

ftnchek:pre
	ftnchek $(CHEKFLAGS) $(F_CPPS) > ftnchek.out

all:	$(F_OBJS)

slave: $(SLAVE_OBJS)

clean :	
	$(RM) -f *.o *~ *.f *.lst *.i

dirac.x:
	cd .. ; $(MAKE) dirac.x

dirmpi.x:
	cd .. ; $(MAKE) dirmpi.x

depend:
	makedepend -- $(DEPEND_FLAGS) -- $(F_OBJS:.o=.F) $(SLAVE_OBJS:.o=.F)
# DO NOT DELETE THIS LINE -- make depend depends on it.

