
include ../../Makefile.config

F_OBJS = ccfspc.o 

all:	$(F_OBJS) 

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i

dirac.x:
	cd .. ; $(MAKE) dirac.x

dirmpi.x:
	cd .. ; $(MAKE) dirmpi.x

dirpvm.x:
	cd .. ; $(MAKE) dirpvm.x

depend:
	makedepend -- $(DEPEND_FLAGS) -- $(F_OBJS:.o=.F)
# DO NOT DELETE THIS LINE -- make depend depends on it.

