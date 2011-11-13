MODS= -DMOD_DFT
MODS_IBM=,,-DMOD_DFT,,,
#
# VERSION : $Revision: 1.11.18.1 $
# DATE    : $Date: 2004/09/01 10:19:11 $
# FILE    : pdpack/Makefile.in
#

include ../Makefile.config

#
# Delete the entry libpdlib.a(gpblas.o) if you supply your own BLAS !
#

F_OBJS=libpdlib.a(arhpack.o) libpdlib.a(dge.o) libpdlib.a(dsp.o)\
       libpdlib.a(linextra.o) libpdlib.a(printpkg.o)\
       libpdlib.a(gpblas.o) libpdlib.a(eispack.o)

SRC=arhpack.F dge.F dsp.F linextra.F printpkg.F \
    gpblas.F  eispack.F

F_CPPS = $(F_OBJS:.o=.i)

all :	libpdlib.a

pre:	$(F_CPPS)

ftnchek:pre
	ftnchek $(CHEKFLAGS) $(F_CPPS) > ftnchek.out

libpdlib.a: $(F_OBJS)

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.a

dirac.x:
	cd .. ; $(MAKE) dirac.x

dirmpi.x:
	cd .. ; $(MAKE) dirmpi.x

depend:
	makedepend -- $(DEPEND_FLAGS) -- $(SRC:.o=.F) 
# DO NOT DELETE THIS LINE -- make depend depends on it.

