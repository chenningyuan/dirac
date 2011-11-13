MODS= -DMOD_DFT
MODS_IBM=,,-DMOD_DFT,,,
#
# VERSION : $Revision: 1.15.4.1 $
# DATE    : $Date: 2004/09/01 10:19:09 $
# FILE    : gp/Makefile.in
#

include ../Makefile.config

F_OBJS =  gphjj.o    gpjeppe.o     gptrygve.o     mempkg.o \
	  qpack.o    vcraypack.o   gpsaue.o       gpluuk.o \
	  gpjth.o    gplucia.o     gpkrmc.o

C_OBJS = gpc.o

MXM_OBJS = mxmpack.o

F_CPPS = $(F_OBJS:.o=.i)

pre:	$(F_CPPS)

ftnchek:pre
	ftnchek $(CHEKFLAGS) $(F_CPPS) > ftnchek.out

all:	setcomp $(F_OBJS) $(C_OBJS)

setcomp:
	( \
	   rm -f gpc.o; \
	   rm -f compinfo.h ; \
	   echo "#define WHO_COMPILED \"${LOGNAME}\"" > compinfo.h ; \
	   echo "#define WHERE_COMPILED \"${HOST}\"" >> compinfo.h ; \
	)

mxm:	$(MXM_OBJS)

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i


dirac.x:
	cd .. ; $(MAKE) dirac.x

dirmpi.x:
	cd .. ; $(MAKE) dirmpi.x

depend: setcomp
	makedepend -- $(DEPEND_FLAGS) -- $(F_OBJS:.o=.F) $(C_OBJS:.o=.c)
# DO NOT DELETE THIS LINE -- make depend depends on it.

