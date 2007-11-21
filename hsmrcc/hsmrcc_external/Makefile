#################################################
#
# This is Makefile for NEWHSMR program
# This can be compiled with GNU Make
# (not ordinary Make)
#
# Makefile generates standalone hsmrcc.exe program.
#
# For the DIRAC enviroment adapted by Miro Ilias, Tel Aviv, 2007
#
##################################################

# Use ksh or bash for shell
#SHELL = /bin/bash

#-----------------------------------------------
# control chars for bold display
# cs : startbold ce : endbold
#-----------------------------------------------

cs = [1m 
ce = [m

###########################################
#
#          FORTRAN90 compiler 
#
#    ... from DIRAC
#
###########################################

##### include variables from DIRAC Makefile.config !!!  ####
include ../../Makefile.config

#FC      = pgf90
#FC      = gfortran
FC      =  $(F77)

###########################################
#
#              C compiler 
#
#   ... from DIRAC
#
###########################################
CCo      = $(CC)
#CCo      = $(F77)

##############################################################################
###
###        set up libraries from DIRAC
###
###  NOTE: in Makefile.config $(LIB_LIST) the -llapack must be before -lblas
###
##############################################################################
FCLIBS  = $(LIB_LIST) 
#FCLIBS  = -lacml
#FCLIBS  = -lblas -llapack
#FCLIBS  = $(LIB_LIST)  -lblas
#FCLIBS  = $(LIB_LIST)  -lacml
#FCLIBS  = $(LIB_LIST) -lblas -llapack

### linking flags ###
#LKFLAGS = -Mnomain
#LKFLAGS = 
LKFLAGS = $(LDFLAGS) 

CCLIBS  = 

##########################################################
#
#      These are options for specific compilers
#
##########################################################

# guessed default, if no compiler in the list ....
        FFLAGS  =  -I$(moddir) -module $(moddir) $(OPTFLAGS) 

ifeq ($(FC),pgf90)
### Portland group
        FFLAGS  =  -I$(moddir) -module $(moddir) $(OPTFLAGS) -Mnomain
       # FFLAGS  =  -I$(moddir) -module $(moddir) $(OPTFLAGS) -pg -Mnomain -Mbounds 
  # ... -pg flags causing pgdbg error stop !!!
       # FFLAGS  =  -I$(moddir) -module $(moddir) $(OPTFLAGS) -pg -Mnomain
endif

ifeq ($(FC),gfortran)
## gfortran/gcc
        FFLAGS  = -I$(moddir)  -J$(moddir) $(OPTFLAGS)
endif

ifeq ($(FC),ifort)
### Intel Fortran 
        FFLAGS  =  -I$(moddir) -module $(moddir) $(OPTFLAGS) 
endif

# ... others combinations...
#FFLAGS  = -Mdclchk -i8 -i8storage -I$(moddir) -module $(moddir) -fast
#FFLAGS  = -I$(moddir)  -fdefault-integer-8 -J$(moddir) $(LDFLAGS)
#FFLAGS  = -I$(moddir) -J$(moddir) $(LDFLAGS)
#FFLAGS  = -I$(moddir) -module $(moddir) $(LDFLAGS)  -Mdclchk 
#FFLAGS  = -I$(moddir)  -fdefault-integer-8 -J$(moddir)
#CFLAGS  = -I$(moddir) -i8 -i8storage -fast
#CFLAGS  = -I$(moddir) -fdefault-integer-8 

CFLAGS  = -I$(moddir) $(C_OPTFLAGS)

# These are options for generating a debuggable code
#FFLAGS  = -Mdclchk -i8 -i8storage -I$(moddir) -module $(moddir) -g -C
#CFLAGS  = -I$(moddir) -i8 -i8storage -g


##########################################################################
#
#
#  Directories for modules, Fortran90/C source files and object files.
#
#
###########################################################################

# ... directory with module source files and module files  ...
moddir = modules
# ... directory with source files  ...
srcdir = source
# ... directory with objects ...
objdir = objects

# specifying location of
# Fortran Source Files

vpath
vpath %.o   $(objdir)
#vpath %.f90 $(moddir) $(srcdir)
vpath %.F90 $(moddir) $(srcdir)
vpath %.c   $(srcdir)
vpath %.h   $(moddir)

#################################
####     Object Files        ####
#################################

modobjs := parameters.o   options.o      orbmod.o    iomodule.o   \
           tindexmod.o    timingmod.o    tamps.o     vacmod.o     \
           intgenmod.o    convmod.o      shift.o     heffmod.o    \
           directdiags.o  renormdiags.o  vinfmod.o   gmsmod.o     \
           syminf.o       diagonalize.o  hmatrix.o   symblock.o   \
           s2matele.o     prntmat.o      n_choose_r.o

f90objs := initguess.o   inihsbase.o   endhsmr.o   makezero.o   \
           sethsmropts.o printinfo.o   endenom.o   diagrams.o   \
           hsmrcc.o      inihsiter.o   runinfo.o   maproots.o   \
           startiter.o   compensate.o  genhsmrpars.o

c89objs := getopt.o hsmrmain.o

#######################
####               ####
#### DIRAC imports ####
####               ####
#######################
# ... no dirac objects/libraries ! 
#dirac_objs :=
#dirac_objs := ../../gp/gptrygve.o  ../../gp/qpack.o ../../gp/gphjj.o  
dirac_objs :=  ../../gp/qpack.o ../../gp/gptrygve.o  ../../gp/gphjj.o 
#dirac_libs := ../../pdpack/libpdlib.a
dirac_libs := 

objects := $(modobjs) $(f90objs) $(c89objs)

###############################################
###      name of the executable file        ###
###############################################

### to fix main dependence form the ../../Makefile.config file
dirac.x : $(executable)

executable = hsmrcc.exe

### default : executable file ###
default : $(executable) 
	@echo "Independent executable file " $(cs) $(executable) $(ce) " is ready." 

##################
###  linking   ###
##################

$(executable) : $(objects)
	@echo -n "Linking all objects files in" $(cs)$(objdir)$(ce) "into executable ..."
	@echo  
	@echo -n "...   where used system libraries are " $(cs) $(FCLIBS) $(CCLIBS) $(ce) 
	@echo  
	@echo -n "...   used compiler linking flags are " $(cs) $(LKFLAGS) $(FFLAGS) $(ce) 
	@echo  
	@echo -n "...      and used DIRAC resources are " $(cs) $(dirac_objs) $(dirac_libs) $(ce) 
	@echo  
	@$(FC) -o $@ $(LKFLAGS) $(FFLAGS) $(addprefix $(objdir)/,$(objects)) \
                     $(dirac_objs) $(dirac_libs) $(FCLIBS) $(CCLIBS)
	@echo "  is done."
	@echo "The independedr executable file generated is ==> " $(cs)$(@F)$(ce) " <=="

##################
###  cleaning  ###
##################
clean :
	@echo -n "Cleaning up objects files in" $(cs)$(objdir)$(ce) " directory ..."
	@rm -f $(addprefix $(objdir)/,$(objects))
	@echo "  Done."
	@echo -n "Cleaning up module files in" $(cs)$(moddir)$(ce) " directory ..."
	@rm -f $(moddir)/*.mod
	@echo "  Done."
	@echo -n "Cleaning up executable " $(cs)$(executable)$(ce) " in main hsmrcc_external directory ..."
	@rm -f $(executable) 
	@echo "  Done."

### dependencies other than implicit dependencies

###################################################################
#              Intermodule Dependencies
###################################################################

$(f90objs) : $(modobjs)
tindexmod.o   : parameters.o
timingmod.o   : iomodule.o
iomodule.o    : parameters.o options.o
syminf.o      : iomodule.o
orbmod.o      : parameters.o iomodule.o
vacmod.o      : parameters.o iomodule.o maxmin.h gmsmod.o symblock.o s2matele.o
tamps.o       : parameters.o iomodule.o
gmsmod.o      : parameters.o iomodule.o options.o syminf.o orbmod.o
intgenmod.o   : parameters.o iomodule.o vacmod.o timingmod.o maxmin.h gmsmod.o
heffmod.o     : parameters.o orbmod.o options.o iomodule.o vacmod.o tamps.o tindexmod.o timingmod.o gmsmod.o
convmod.o     : parameters.o options.o iomodule.o heffmod.o
shift.o       : parameters.o options.o orbmod.o iomodule.o convmod.o gmsmod.o
directdiags.o : parameters.o iomodule.o tamps.o timingmod.o tind.h maxmin.h
renormdiags.o : parameters.o options.o iomodule.o vacmod.o tamps.o tindexmod.o heffmod.o timingmod.o
vinfmod.o     : parameters.o options.o iomodule.o gmsmod.o shift.o intgenmod.o convmod.o
s2matele      : parameters.o iomodule.o gmsmod.o
diagonalize.o : iomodule.o
hmatrix.o     : iomodule.o

########################################################################
# Subroutines Depending on Modules
########################################################################

diagrams.o    : parameters.o iomodule.o tamps.o heffmod.o directdiags.o convmod.o renormdiags.o timingmod.o
hsmrcc.o      : convmod.o intgenmod.o options.o iomodule.o shift.o heffmod.o vinfmod.o
endenom.o     : parameters.o orbmod.o vacmod.o maxmin.h tind.h
makezero.o    : parameters.o vacmod.o tind.h
compensate.o  : parameters.o orbmod.o vacmod.o shift.o tind.h maxmin.h
inihsbase.o   : parameters.o options.o iomodule.o intgenmod.o convmod.o shift.o vacmod.o timingmod.o vinfmod.o gmsmod.o
inihsiter.o   : iomodule.o tamps.o directdiags.o convmod.o heffmod.o shift.o
printinfo.o   : options.o iomodule.o timingmod.o heffmod.o convmod.o
sethsmropts.o : parameters.o options.o iomodule.o shift.o convmod.o intgenmod.o 
endhsmr.o     : parameters.o options.o iomodule.o heffmod.o gmsmod.o directdiags.o convmod.o tamps.o shift.o vinfmod.o vacmod.o
initguess.o   : parameters.o options.o iomodule.o tamps.o heffmod.o tind.h
runinfo.o     : parameters.o options.o iomodule.o intgenmod.o gmsmod.o orbmod.o shift.o convmod.o
getopt.o      : getopt.h
hsmrmain.o    : sethsmropts.o getopt.o getopt.h
maproots.o    : parameters.o heffmod.o iomodule.o convmod.o
startiter.o   : parameters.o heffmod.o tamps.o timingmod.o
genhsmrpars.o : parameters.o iomodule.o orbmod.o

# The following FLAGS are for debug purpose
# these can be uncommented if any group of
# files require debug information.

#directdiags.o : FFLAGS = -u -I$(moddir) -g2 -ladebug -C

# pattern rules for compilation

%.o : %.F90
	@printf "Compiling F90 File ==> %s %s %s \t %s" $(cs) $@ $(ce) "....."
	@echo -n "using" $(cs) $(FC) $(FFLAGS) "-c -o" $(ce)
	@$(FC) $(FFLAGS) -c -o $(objdir)/$@ $<
	@printf "%s\n" " ...Done."

%.o : %.c
	@printf "Compiling C89 File ==> %s %s %s \t %s" $(cs) $@ $(ce) "....."
	@echo -n "using" $(cs) $(CCo) $(CFLAGS) "-c -o" $(ce)
	@$(CCo) $(CFLAGS) -c -o $(objdir)/$@ $<
	@printf "%s\n" " ...Done."

