MODS= -DMOD_DFT
MODS_IBM=,,-DMOD_DFT,,,

include ../Makefile.config

F_OBJS = pamdft.o dftden.o dftdrv.o dftfun.o dftgrd.o grtest.o \
	Lebedev-Laikov.o

C_OBJS = functionals.o emb-functional.o fun-slater.o fun-becke.o fun-vwn.o fun-lyp.o \
	fun-pw86x.o fun-p86c.o fun-kin_tf.o fun-kin_vw.o

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

fun-tester: fun-tester.c fun-slater.c fun-vwn.c fun-becke.c fun-lyp.c functionals.c fun-lb94.c \
	emb-functional.c fun-p86c.c fun-pw86x.c fun-kin_tf.c fun-kin_vw.c
	$(CC) -g -O -I../include -DNO_BACKWARD_COMP $^ -o $@  -lm

depend:
	makedepend -- $(DEPEND_FLAGS) -- $(F_OBJS:.o=.F) $(C_OBJS:.o=.c)
# DO NOT DELETE THIS LINE -- make depend depends on it.

