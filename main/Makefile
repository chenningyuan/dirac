MODS= -DMOD_DFT
MODS_IBM=,,-DMOD_DFT,,,
include ../Makefile.config

C_OBJS = diracmalloc.o alloc_mem.o get_memsize.o print_logo.o ini_mpi.o\
         end_mpi.o csectid.o main.o exec_dirac.o mpi_errtrap.o pr_timing.o\
         get_time.o 

F_OBJS = dirac.o

F_CPPS = $(F_OBJS:.o=.i)

pre:	$(F_CPPS)

ftnchek:pre
	ftnchek $(CHEKFLAGS) $(F_CPPS) > ftnchek.out

all:	$(F_OBJS) $(C_OBJS)

clean:
	$(RM) -f *.o *~ *.f *.lst *.i

dirmpi.x:
	cd .. ; $(MAKE) dirmpi.x

dirac.x:
	cd .. ; $(MAKE) dirac.x

depend:
	makedepend -- $(DEPEND_FLAGS) -- $(F_OBJS:.o=.F) $(C_OBJS:.o=.c)
# DO NOT DELETE THIS LINE -- make depend depends on it.

