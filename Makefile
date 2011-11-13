MODS=  -DMOD_DFT   
MODS_IBM=  -DMOD_DFT   
include ./Makefile.config

#
# Note that only object files from fortran files (.F) should be included in OBJS,
# object files from C files (.c) should be in C_OBJS.
# 

OBJS = \
dirac/dircnv.o dirac/dirgp.o  dirac/dirone.o dirac/dirv1c.o \
dirac/dirout.o dirac/dirana.o dirac/dirrdn.o dirac/dirscf.o \
dirac/dirset.o dirac/dirtra.o dirac/dirtwo.o dirac/dirden.o \
dirac/dirgrd.o dirac/dirsol.o dirac/dirgeo.o dirac/dirspf.o \
dirac/diropen.o dirac/diraos.o \
main/dirac.o \
abacus/her1ave.o    abacus/her1car.o    abacus/her1cat.o    abacus/her1drv.o    abacus/her1int.o\
abacus/her1odc.o    abacus/her1out.o    abacus/her1pro.o    abacus/her1sym.o\
abacus/her2ave.o    abacus/her2dir.o    abacus/her2drv.o    abacus/her2el1.o\
abacus/her2el2.o    abacus/her2fck.o    abacus/her2her.o    abacus/her2odc.o\
abacus/her2out.o    abacus/her2sym.o    abacus/herbas.o\
abacus/herdirect.o  abacus/herdrv.o     abacus/hergam.o     abacus/hergeo.o\
abacus/hergp.o      abacus/herrdn.o     abacus/hersol.o\
abacus/hersym.o     abacus/hertst.o     abacus/her2gab.o    abacus/hermem.o\
abacus/abanuc.o     abacus/abatro.o     abacus/abaortvec.o abacus/abarint.o \
abacus/abaopt.o     abacus/abavrml.o abacus/abawalk.o abacus/abaop2.o \
gp/gphjj.o          gp/gpjeppe.o        gp/gptrygve.o       gp/mempkg.o\
gp/qpack.o          gp/vcraypack.o      gp/gpsaue.o         gp/gpluuk.o \
gp/gpjth.o          gp/gplucia.o        gp/gpkrmc.o \
prp/paminp.o        prp/pamexp.o        prp/pamprp.o        prp/pamset.o \
prp/pamgrd.o        prp/pamogrd.o       prp/pamfck.o        \
prp/pamxlr.o        prp/pamrsp.o        prp/pamrvc.o        prp/rspmpg.o \
prp/mp2lag.o        prp/pamxqr.o        \
moltra/tradr1.o moltra/tradr2.o moltra/tradr3.o moltra/tradr4.o moltra/tradr5.o \
moltra/trainp.o moltra/traone.o moltra/traout.o moltra/trapam.o \
moltra/traprp.o \
relccsd/ccpam.o  relccsd/ccdriv.o relccsd/cceqns.o relccsd/ccints.o \
relccsd/ccsort.o relccsd/ccgenp.o relccsd/ccgetv.o relccsd/ccintm.o \
relccsd/ccrpa.o  relccsd/ccxblas.o relccsd/ccgrad.o relccsd/ccfspc.o \
relccsd/waio.o  \
dirrci/cidirr.o dirrci/cigenp.o dirrci/cigosc.o \
lucita/ciinfo.o  lucita/density.o  lucita/diagonal.o  lucita/diagonalize.o  lucita/doubledir.o  \
lucita/fock.o  lucita/interface.o  lucita/io.o  lucita/lucia_cc.o  lucita/lucia_icci.o  lucita/memory.o \
lucita/program.o  lucita/property.o  lucita/sigma.o  lucita/strings.o  lucita/transform.o  lucita/utils.o  \
lucita/vecmat.o

#
# All .o files for C routines in .c files are listed separately from OBJS
# (otherwise you'll get a lot of irrelevant error messages from ftnchek)
#

C_OBJS = \
main/alloc_mem.o    main/csectid.o \
main/diracmalloc.o  main/end_mpi.o    main/exec_dirac.o \
main/get_memsize.o  main/ini_mpi.o    main/main.o \
main/mpi_errtrap.o  main/pr_timing.o  main/print_logo.o \
main/get_time.o \
gp/gpc.o

#

OBJSLAVE = abacus/herpar.o

PDLIB = -Lpdpack -lpdlib

#
# RELCCSD_OBJS and RELCCSD_MAIN are for stand-alone relccsd.x and mpiccsd.x
#

RELCCSD_OBJS = \
relccsd/ccdriv.o relccsd/cceqns.o relccsd/ccints.o \
relccsd/ccsort.o relccsd/ccgenp.o relccsd/ccgetv.o relccsd/ccintm.o \
relccsd/ccrpa.o  relccsd/ccxblas.o relccsd/ccgrad.o  \
relccsd/ccfspc.o  relccsd/waio.o 

RELCCSD_MAIN = relccsd/ccmain.o

#
# DIRRCI_OBJS and DIRRCI_MAIN are for stand-alone dirrci.x
#

DIRRCI_OBJS = \
dirrci/cidirr.o dirrci/cigenp.o dirrci/cigosc.o 

DIRRCI_MAIN = dirrci/cimain.o

#
# GP_OBJS are for stand-alone programs relccsd.x, mpiccsd.x, and dirrci.x
#         and for some 'make utils' programs.
#

GP_OBJS = \
gp/gphjj.o          gp/gpjeppe.o        gp/gptrygve.o       gp/mempkg.o\
gp/qpack.o          gp/vcraypack.o      gp/gpsaue.o         gp/gpluuk.o \
gp/gpjth.o          gp/gplucia.o        gp/gpkrmc.o         gp/gpc.o


# dft files written in fortran:

DFTOBJS = \
dft/pamdft.o dft/dftden.o      dft/dftdrv.o     dft/dftgrd.o dft/grtest.o dft/Lebedev-Laikov.o

# dft files written in C:

DFTOBJS_C = \
dft/functionals.o dft/emb-functional.o dft/fun-slater.o dft/fun-vwn.o dft/fun-becke.o dft/fun-lyp.o dft/fun-pw86x.o dft/fun-p86c.o dft/fun-kin_tf.o dft/fun-kin_vw.o

OBJSMXM = gp/mxmpack.o

TOT_FORTRANOBJS = $(OBJS) $(ESROBJS) $(CCFOCKSPOBJS) \
	$(DFTOBJS) $(DNFOBJS) $(ERIOBJS) $(XPPOBJS) $(KRMC_FOBJS)

F_CPPS = $(TOT_FORTRANOBJS:.o=.i)

dirac.x: MAIN_OBJ DIR_OBJ PRP_OBJ HER_OBJ GP_OBJ MXM_OBJ PD_OBJ \
	TRA_OBJ CCM_OBJ CI_OBJ LUCITA_OBJ   DFT_OBJ      
	$(F77) $(LDFLAGS) -o $(INSTALLDIR)/dirac.x \
	$(OBJS) $(C_OBJS) \
	   $(DFTOBJS) $(DFTOBJS_C)      \
	$(OBJSMXM) $(PDLIB) $(LIB_LIST) 
#
# Run ftnchek on Dirac
#
pre:	$(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) $(F_CPPS) > dirac.ftnchek
dirmpi.x: MAIN_OBJ DIR_OBJ PRP_OBJ HER_OBJ GP_OBJ MXM_OBJ PD_OBJ \
	TRA_OBJ CCM_OBJ LUCITA_OBJ CI_OBJ   DFT_OBJ       SLAVE_OBJ
	$(F77) $(LDFLAGS) -o $(INSTALLDIR)/dirac.x \
	$(OBJS) $(C_OBJS) $(OBJSLAVE) \
	  $(DFTOBJS) $(DFTOBJS_C)      \
	$(OBJSMXM) $(LIB_LIST) $(PDLIB) \
	$(MPI_LIB_PATH) $(MPI_LIB)

relccsd.x: GP_OBJ PD_OBJ \
	CCM_OBJ
	$(F77) $(LDFLAGS) -o $(INSTALLDIR)/relccsd.x \
	$(RELCCSD_MAIN) $(RELCCSD_OBJS) $(GP_OBJS) $(OBJSMXM) \
	$(LIB_LIST) $(PDLIB)

mpiccsd.x: GP_OBJ PD_OBJ \
	CCM_OBJ
	$(F77) $(LDFLAGS) -o $(INSTALLDIR)/relccsd.x \
	$(RELCCSD_MAIN) $(RELCCSD_OBJS) $(GP_OBJS) $(OBJSMXM) \
	$(LIB_LIST) $(PDLIB) \
	$(MPI_LIB_PATH) $(MPI_LIB)


dirrci.x: GP_OBJ PD_OBJ \
	CI_OBJ
	$(F77) $(LDFLAGS) -o $(INSTALLDIR)/dirrci.x \
	$(DIRRCI_MAIN) $(DIRRCI_OBJS) $(GP_OBJS) $(OBJSMXM) \
	$(LIB_LIST) $(PDLIB)


utils	: GP_OBJ PD_OBJ UTIL_OBJ 
	$(F77) $(LDFLAGS) -o $(INSTALLDIR)/twofit.x \
	utils/twofit.o abacus/hergp.o $(GP_OBJS) $(OBJSMXM) $(PDLIB) $(LIB_LIST) 
	$(F77) $(LDFLAGS) -o $(INSTALLDIR)/labread.x \
	utils/labread.o abacus/hergp.o $(GP_OBJS) $(OBJSMXM) $(PDLIB) $(LIB_LIST) 
	$(F77) $(LDFLAGS) -o $(INSTALLDIR)/rspread.x \
	utils/rspread.o 
	$(F77) $(LDFLAGS) -o $(INSTALLDIR)/polfit.x \
	utils/polfit.o abacus/hergp.o $(GP_OBJS) $(OBJSMXM) $(PDLIB) $(LIB_LIST) 


clean   :
	cd dirac    ; $(MAKE) clean
	cd abacus   ; $(MAKE) clean
	cd prp      ; $(MAKE) clean
	cd main ; $(MAKE) clean
	cd moltra   ; $(MAKE) clean
	cd relccsd  ; $(MAKE) clean
	cd dirrci  ; $(MAKE) clean
	cd lucita  ; $(MAKE) clean
	cd dft      ; $(MAKE) clean
	cd gp       ; $(MAKE) clean
	cd pdpack   ; $(MAKE) clean
	$(RM) -f *~
#
#  Export target
#
export:
	./export
#
#  Objects
#

MAIN_OBJ:
	cd main ; $(MAKE) all

DIR_OBJ :
	cd dirac ; $(MAKE) all

PRP_OBJ :
	cd prp; $(MAKE) all

TRA_OBJ :
	cd moltra; $(MAKE) all

CCM_OBJ :
	cd relccsd ; $(MAKE) all

SLAVE_OBJ :
	cd abacus  ; $(MAKE) slave

CI_OBJ:
	cd dirrci ; $(MAKE) all

LUCITA_OBJ:
	cd lucita ; $(MAKE) all


DFT_OBJ :
	cd dft; $(MAKE) all

HER_OBJ :
	cd abacus ; $(MAKE) all

GP_OBJ :
	cd gp ; $(MAKE) all

MXM_OBJ :
	cd gp ; $(MAKE) mxm


PD_OBJ:
	cd pdpack ; $(MAKE) all

UTIL_OBJ:
	cd utils ; $(MAKE) all


depend:
	cd dirac; $(MAKE) depend
	cd prp; $(MAKE) depend
	cd abacus; $(MAKE) depend
	cd gp; $(MAKE) depend
	cd pdpack; $(MAKE) depend
	cd main; $(MAKE) depend
	cd relccsd; $(MAKE) depend
	cd moltra; $(MAKE) depend
	cd dirrci; $(MAKE) depend
	cd lucita; $(MAKE) depend
	cd dft; $(MAKE) depend
