#!/bin/sh
#######################################################################################################
#
#
#           Script for copying and launching the DIRAC suite on grid CE
#
#
#######################################################################################################

# check the parameter - VO name

if [[ $1 != "voce" && $1 != "compchem" && $1 != "isragrid" && $1 != "osg" ]]; then
 echo -e "\n wrong parameter - VO name : $1 "
 exit -12
else
 VO=$1
 echo -e "OK, you specified properly the VO=$VO"
 #exit 0
fi

unset THIS_HOSTNAME
THIS_HOSTNAME=`hostname -f`
if [[ -z $THIS_HOSTNAME ]]; then
 THIS_HOSTNAME=`hostname`
fi

echo -e "\n --- Hostname: $THIS_HOSTNAME "; 
echo -e "Server specification: \c";  uname -a
######################################################
#     get number of CPUs on the given node
######################################################
nprocs_all=`cat /proc/cpuinfo | grep processor | wc -l`
nprocs_PBS=`cat $PBS_NODEFILE | wc -l`
nprocs=$nprocs_all
echo -e "\n==== CPU \c"; cat /proc/cpuinfo | grep "model name"  | tail -n 1
echo -e "==== Number of all CPUs     :  $nprocs   ";
echo -e "==== Number of allowed CPUs :  $nprocs_PBS  ";
if [[ -n $PBS_NODEFILE ]]; then
  echo -e "PBS_NODEFILE=$PBS_NODEFILE"
  echo -e "PBS_NUM_PPN=$PBS_NUM_PPN"
fi 

echo -e "\n ---- server's total memory :\c"; cat /proc/meminfo | grep MemTotal
# check space
echo -e "\n--- Server's filesystem (df): "; df -h .
echo -e "--- My current directory :  $PWD "
echo -e "--- Occupied space by me (du) : \c"; du -h -s

# get info about CE #
echo -e "\n querry attribute SMPSize:"
lcg-info --list-ce --query "CE=$CE_ID" --attrs SMPSize
echo -e "\n querry attribute PlatformArch:"
lcg-info --list-ce --query "CE=$CE_ID" --attrs PlatformArch
echo -e "querry attributes TotalCPUs, Memory, ClockSpeed, MaxCPUTime :"
lcg-info --list-ce --query "CE=$CE_ID" --attrs TotalCPUs,Memory,ClockSpeed,MaxCPUTime


echo -e "\n\n ****  server environment:  set command  **** \n" 
set

echo -e "\n\n Is there cmake ? which cmake="; which cmake
echo -e "Python present : \c"; which python; python -V

###########################
### Check files on SE: ####
###########################
echo -e  "\n Files on my vo=$VO space, lcg-ls -l lfn://grid/$VO/ilias:"
lcg-ls -l lfn://grid/$VO/ilias || { echo -e "...previous command failed, once again... "; lcg-ls -l lfn://grid/$VO/ilias; }
echo -e "\n  File attributes lcg-lr  lfn://grid/$VO/ilias/dirac_grid_suite.tgz:"
lcg-lr    lfn://grid/$VO/ilias/dirac_grid_suite.tgz || { echo -e "...previous command failed, once again... "; lcg-lr  lfn://grid/$VO/ilias/dirac_grid_suite.tgz;  }
echo " ------------------------------------------------------ "
#####################################
###  download file to the current CE 
#####################################
echo -e "\n Going to download from SE, lcg-cp lfn://grid/$VO/ilias/dirac_grid_suite.tgz  file://$PWD/dirac_grid_suite.tgz...."
lcg-cp -v lfn://grid/$VO/ilias/dirac_grid_suite.tgz  file://$PWD/dirac_grid_suite.tgz || { echo "Error exit,downloading of dirac_grid_suite.tg failed" ; exit -2; }
echo -e "\nDownloaded file attributes (ls -lt dirac_grid_suite.tgz) :"
ls -lt dirac_grid_suite.tgz || { echo -e "No dirac_grid_suite.tgz in local directory found, error exit"; exit -3;}
echo " ---------------------------------------------------------------------------------------------------------------------"

##########################################
##  make it executable and finally launch
##########################################
echo -e "\n Job started at : \c"; date

# directory with all static executables - dirac.x, cmake, ctest, mpirun - and related directories
 export BUILD=build_Intel_OpenMPI_MKL_ILP64_STATIC
 export BUILD_SERIAL=build_Intel_serial_MKL_ILP64_STATIC
 export BUILD_GNU_SERIAL=build_GNU_serial_ILP64_STATIC

 unset OPAL_PREFIX
 export OPAL_PREFIX=$PWD/$BUILD
 echo -e "Variable OPAL_PREFIX=$OPAL_PREFIX"
 export PATH=$PWD/$BUILD:$PATH
 #export PATH=$PWD/$BUILD_SERIAL:$PWD/$BUILD:$PATH
 echo -e "Modified PATH=$PATH"
 echo -e "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"


# here exit
#exit

echo -e "\n\n ---- Unpacking dowloaded DIRAC suite, tar : \c "; which tar;  tar --version
echo -e "performing tar xfz dirac_grid_suite.tgz ...\c"
tar xfz dirac_grid_suite.tgz
echo -e "done; see the content of $PWD :"; ls -lt

#echo -e "\n making executables:"
#ls -lt $BUILD/dirac.x  $BUILD/cmake  $BUILD/ctest $BUILD/mpirun  pam  runtest pamadm setup
#chmod  u+x  $BUILD/dirac.x  $BUILD/cmake  $BUILD/ctest $BUILD/mpirun  pam  runtest pamadm setup
#ls -lt $BUILD/dirac.x  $BUILD/cmake  $BUILD/ctest $BUILD/mpirun  pam  runtest pamadm setup

# create diracrc file for paralle run on the grid CE
cp maintenance/grid_computing/diracrc_grid  diracrc
echo -e "\n ---- diracrc created: "; cat diracrc

echo -e "\n --- Launching simple pam test  --- \n "; 
python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --outcmo --mpi=$nprocs
#python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --outcmo --mpi=8
#python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --outcmo 
#python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --scratch=$PWD --outcmo 

#####################################################################
#           specify tests for two runtest runs
#####################################################################

regular_subset='fscc fscc_highspin fscc_intensities krci_energy krci_properties_omega_tdm krci_properties_perm_dipmom cosci_energy acmoin gaunt basis_automatic_augmentation basis_input bss_energy density_at_nuclei dft_ac dft_ao_eval dft_response_symmetry fde_static-vemb_dipole lucita_large lucita_short  mcscf_energy  molecular_mean_field_hamiltonian_schemes mp2_energy mp2_properties  projection_analysis response_lao_shield  symmetry_recognition symmetry_scf_energy visual_lao_shielding ecp'

benchmark_subset='cc_benchmark cc_linear fscc_intensities gaunt mcscf_energy moltra_4or6'


#############################################################
#               check own static cmake stuff
#############################################################

  echo -e "\n Performing pamadm step:"
  python ./pamadm -a

  cp maintenance/grid_computing/diracrc_build_grid  $BUILD/diracrc
  cp maintenance/grid_computing/diracrc_build_grid  $BUILD_SERIAL/diracrc
  cp maintenance/grid_computing/diracrc_build_grid  $BUILD_GNU_SERIAL/diracrc
  echo -e "\n ---- own diracrc in build-directories created: "; cat $BUILD/diracrc

#######################################################################################
#unset DO_CMAKE
DO_CMAKE=1
if [[ -n $DO_CMAKE ]]; then
  echo -e "\n\n===== Performing false configuration setup in order to enable tests running with static mpirun,dirac.x ===== \n"
  cd $BUILD
  echo -e "\n I am in build-directory  $PWD \n"
  echo -e "--- Directory $BUILD content ---:"; ls -lt
# FIX THIS  ERROR:
# CMake Error: The current CMakeCache.txt directory /tmp/ilias_OUWzTUf41_ntFwYFY3Pk7Q/build/CMakeCache.txt is different than the directory /home/ilias/qch_work/software/dirac_git/trunk_current/build where CMakeCache.txt was created. This may result in binaries being created in the wrong place. If you are not sure, reedit the CMakeCache.txt
#CMake Error: The source "/tmp/ilias_OUWzTUf41_ntFwYFY3Pk7Q/CMakeLists.txt" does not match the source "/home/ilias/qch_work/software/dirac_git/trunk_current/CMakeLists.txt" used to generate cache.  Re-run cmake with a different source directory.
  mv CMakeCache.txt CMakeCache.txt_STORED
#note: make  -DENABLE_MPI=OFF in order to avoid cmake's MPI configuration on an unknown grid's CE
  ./cmake  -DVERBOSE_OUTPUT=ON -DBUILDNAME="$THIS_HOSTNAME.grid_ompi_intel_mkl_static" -DENABLE_MPI=OFF  -DRUNTEST_NR_PROC=$nprocs -DRUNTEST_TIMEOUT="00:18:00" -DDART_TESTING_TIMEOUT=99999 -DRUNTEST_BINARY=$PWD/dirac.x ..
# make ExperimentalConfigure
  make ExperimentalStart
  make ExperimentalTest
  make ExperimentalSubmit
## return back into main directory
  cd ..
  echo -e "\n I am in main directory $PWD \n"
  date
fi

  echo -e "\n\n --- Launching parallel runtest - OpenMPI+Intel+MKL - with selected tests --- \n "; date 
  python ./runtest  --verbose --backup --binary=$PWD/$BUILD/dirac.x --mpi=$nprocs --tests="$regular_subset"

  echo -e "\n\n --- Launching benchmark parallel runtest - OpenMPI+Intel+MKL --- \n "; date 
  python ./runtest --benchmark --verbose --backup --binary=$PWD/$BUILD/dirac.x --mpi=$nprocs  --tests="$benchmark_subset"


if [[ -n $DO_CMAKE ]]; then
  echo -e "\n\n===== Performing false configuration setup in order to enable tests running with static Intel, serial dirac.x ===== \n"
### Again, do the buildup using serial standalone executable ###
# http://software.intel.com/en-us/articles/recommended-settings-for-calling-intelr-mkl-routines-from-multi-threaded-applications/
# set MKL threads
  let nprocsm1=$nprocs-1
  export MKL_NUM_THREADS=$nprocsm1; echo -e "MKL threads, MKL_NUM_THREADS=$MKL_NUM_THREADS " 
  cd $BUILD_SERIAL
  echo -e "\n I am in build-directory $PWD \n"
  echo -e "--- Directory $BUILD_SERIAL content ---:"; ls -lt
  mv CMakeCache.txt CMakeCache.txt_STORED
  #note: make  -DENABLE_MPI=OFF in order to avoid cmake's MPI configuration on an unknown grid's CE
 ./cmake  -DVERBOSE_OUTPUT=ON -DBUILDNAME="$THIS_HOSTNAME.grid_serial_intel_mkl_static" -DENABLE_MPI=OFF -DRUNTEST_TIMEOUT="00:18:00" -DDART_TESTING_TIMEOUT=99999 -DRUNTEST_BINARY=$PWD/dirac.x ..
  make ExperimentalStart
  make ExperimentalTest
  make ExperimentalSubmit
## return back into main directory
  cd ..
  echo -e "\n I am again in the main directory $PWD \n"
  date
fi

  echo -e "\n\n --- Launching selected serial runtest - Intel+MKL --- \n "; date 
  python ./runtest --verbose --backup --binary=$PWD/$BUILD_SERIAL/dirac.x --tests="$regular_subset"

  echo -e "\n\n --- Launching benchmark serial runtest - Intel+MKL --- \n "; date 
  python ./runtest --all --benchmark --verbose --backup --binary=$PWD/$BUILD_SERIAL/dirac.x --tests="$benchmark_subset"

### Again, do the buildup using serial standalone executable ###
if [[ -n $DO_CMAKE ]]; then
  echo -e "\n\n===== Performing false configuration setup in order to enable tests running with static, GNU, serial dirac.x ===== \n"
  cd $BUILD_GNU_SERIAL
  echo -e "\n I am in build-directory $PWD \n"
  echo -e "--- Directory $BUILD_GNU_SERIAL content ---:"; ls -lt
  mv CMakeCache.txt CMakeCache.txt_STORED
  #note: make  -DENABLE_MPI=OFF in order to avoid cmake's MPI configuration on an unknown grid's CE
  ./cmake  -DVERBOSE_OUTPUT=ON -DBUILDNAME="$THIS_HOSTNAME.grid_serial_GNU_static" -DENABLE_MPI=OFF -DRUNTEST_TIMEOUT="00:18:00" -DDART_TESTING_TIMEOUT=99999 -DRUNTEST_BINARY=$PWD/dirac.x ..
  make ExperimentalStart
  make ExperimentalTest
  make ExperimentalSubmit
## return back into main directory
  cd ..
  echo -e "\n I am again in the main directory $PWD \n"
  date
fi

  echo -e "\n\n --- Launching serial runtest - GNU - with selected tests  --- \n "; date
  python ./runtest  --verbose --backup --binary=$PWD/$BUILD_GNU_SERIAL/dirac.x --tests="$regular_subset"

  echo -e "\n\n --- Launching benchmark serial runtest - GNU  --- \n "; date 
  python ./runtest --all --benchmark --verbose --backup --binary=$PWD/$BUILD_GNU_SERIAL/dirac.x --tests="$benchmark_subset"


# ----------------------------------------- #
# ---- try own buildup on the grid CE ----- #
# ----------------------------------------- #

 echo -e "\n\n\n ==== Try own buildup on the host server (in build directory) ==== \n "

 unset FC
 if [[ -n `which pgf90` ]]; then
   FC=`which pgf90`
   echo -e "\n pgf90 ? $FC ";
   pgf90 -V
 elif [[ -n `which ifort` ]]; then
   FC=`which ifort`
   echo -e "\n ifort ? $FC"; 
   ifort --version
 elif [[ -n `which lfc` ]]; then
   FC=`which lfc`
   lfc --version
   echo -e "\n lfc ? $FC "; 
 elif [[ -n  `which gfortran44` ]]; then
   FC=`which gfortran44`
   echo -e "\n gfortran44 ? $FC "; 
   gfortran44 --version
 elif [[ -n  `which gfortran` ]]; then
   FC=`which gfortran`
   echo -e "\n gfortran ? $FC "; 
   gfortran --version
 else
   echo -e "\n No suitable Fortran compiler found..."
 fi

 unset CC
 if [[ -n `which pgcc` ]]; then
   CC=`which pgcc`
   pgcc -V
   echo -e "\n pgcc ? $CC ";
 elif [[ -n `which icc` ]]; then
   CC=`which icc`
   echo -e "\n icc ? $CC"; 
   icc --version
 elif [[ -n `which fcc` ]]; then
   CC=`which fcc`
   echo -e "\n fcc ? $CC ";
   fcc --version 
 elif [[ -n  `which gcc44` ]]; then
   CC=`which gcc44`
   echo -e "\n gcc44 ? $CC "; 
   gcc44 --version
 elif [[ -n  `which gcc` ]]; then
   CC=`which gcc`
   echo -e "\n gcc ? $CC "; 
   gcc --version
 else
   echo -e "\n No suitable C compiler found..."
 fi

 unset CXX
 if [[ -n `which pgCC` ]]; then
   CXX=`which pgCC`
   pgCC -V
   echo -e "\n pgCC ? $CXX ";
 elif [[ -n `which icpc` ]]; then
   CXX=`which icpc`
   echo -e "\n iccpc ? $CXX"; 
   icpc --version
 elif [[ -n `which fcc` ]]; then
   CXX=`which fcc`
   echo -e "\n fcc ? $CXX ";
   fcc --version 
 elif [[ -n  `which g++44` ]]; then
   CXX=`which g++44`
   echo -e "\n g++44 ? $CXX "; 
   g++44 --version
 elif [[ -n  `which g++` ]]; then
   CXX=`which g++`
   echo -e "\n g++ ? $CXX "; 
   g++ --version
 else
   echo -e "\n No suitable C++ compiler found !"
 fi

 if [[ -n $MKLROOT ]]; then
    echo -e "\n ===  Found Intel MKL library, MKLROOT=$MKLROOT === \n"
    export MATH_ROOT=$MKLROOT ; echo "MATH_ROOT variable activated"
 fi

##########################################################
## check first envirinment variables related to OpenMPI ##
##########################################################
 if [[ -n $MPI_OPENMPI_PATH ]]; then
   echo -e "\n MPI_OPENMPI_PATH=$MPI_OPENMPI_PATH"
   ls -lt $MPI_OPENMPI_PATH
 fi


### check presence of parallel MPI compilers ###
 MPISWITCH='off'; MPISWITCH_F90='off'; MPISWITCH_CC='off'; MPISWITCH_CXX='off';
 if [[ -n `which mpif90` ]]; then
    echo -e "\n Is  mpif90 ? \c"; which mpif90
    mpif90 --version
    MPISWITCH_F90='on'
 fi
 if [[ -n `which mpicc` ]]; then
    echo -e "\n Is  mpicc ? \c";  which mpicc
    mpicc --version
    MPISWITCH_CC='on'
 fi
 if [[ -n `which mpiCC` ]]; then
    echo -e "\n Is  mpiCC ? \c";  which mpiCC
    mpiCC --version
    MPISWITCH_CXX='on'
 fi
 if [[ $MPISWITCH_F90 == 'on' && $MPISWITCH_CC == 'on' && $MPISWITCH_CXX == 'on' ]]; then
    if [[ -n `which mpirun` ]]; then
        echo -e "\n Is  mpirun ? \c";  which mpirun
        MPISWITCH='on'
    fi
 fi

 BUILD_NAME="$THIS_HOSTNAME.gridCE_ownbuild"
 echo "Assigned BUILD NAME=$BUILD_NAME"
 dirac_file=$PWD/build/dirac.x

 if [[ -n $FC  && -n $CC ]]; then
    echo "buildup with selected FC=$FC and CC=$CC"
    python ./setup --fc=$FC --cc=$CC -D VERBOSE_OUTPUT=ON --int64=on --mpi=off -D BUILDNAME=$BUILD_NAME -D RUNTEST_TIMEOUT="00:40:00" -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_BINARY=$dirac_file
 elif [[ -n $FC  && -z $CC ]]; then
    echo "buildup with selected FC=$FC"
    python ./setup --fc=$FC  -D VERBOSE_OUTPUT=ON --int64=on --mpi=off -D BUILDNAME=$BUILD_NAME -D RUNTEST_TIMEOUT="00:40:00" -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_BINARY=$dirac_file
 elif [[ -z $FC  && -n $CC ]]; then
    echo "buildup with selected CC=$CC"
    python ./setup --cc=$CC -D VERBOSE_OUTPUT=ON --int64=on --mpi=off -D BUILDNAME=$BUILD_NAME -D RUNTEST_TIMEOUT="00:40:00" -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_BINARY=$dirac_file
 elif [[ -z $FC  && -z $CC ]]; then
    echo "buildup with cmake automatic selection of compilers"
    python ./setup -D VERBOSE_OUTPUT=ON --int64=on --mpi=off -D BUILDNAME=$BUILD_NAME -D RUNTEST_TIMEOUT="00:40:00" -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_BINARY=$dirac_file
 fi

 if [ -d "build" ]; then
   cd build
   make Experimental
   #make ExperimentalStart
   #make ExperimentalConfigure
   #make ExperimentalBuild
   #make ExperimentalSubmit

       # return to the main directory
   cd ..
 fi

 unset BUILD_LOCAL
 if [ -f "$dirac_file" ]; then
   echo -e "\n\n --- Launching  runtest - selected tests - of own dirac.x buildup on CE --- \n "; date
   echo -e "dirac.x file: \c "; ls -lt $dirac_file
   python ./runtest --verbose --backup --binary=$dirac_file --tests="$regular_subset"
   BUILD_LOCAL=build
 else
   echo -e "\n own buildup of dirac.x failed, no tests run \n"
 fi 


##############################################################
#
# pack selected files to get them back from CE
#
##############################################################
echo -e "\n --------------------------------- \n "; 
# delete old tar-ball first
ls -lt dirac_grid_suite.tgz
echo -e "\n deleting the old dirac_grid_suite.tgz..."
rm dirac_grid_suite.tgz
echo "check files..."
ls -lt
echo -e "\n --- Packing all wanted stuff back from the grid CE  --- "; 
tar --version
echo -e "\n we have to pack (ls -lt) :"
ls -lt
echo " "
#tar cvzf dirac_suite_back.tgz test *
tar czf dirac_grid_suite_back.tgz test $BUILD_LOCAL
echo -e "\n selected files of the DIRAC suite packed back , ls -lt dirac_grid_suite_back.tgz:";ls -lt dirac_grid_suite_back.tgz

############################
#  print out some info .. 
############################
echo -e "\n -------------------------------------------------------"
echo -e "Job finished at: \c"; date
echo -e "I am here: \c"; pwd
echo -e "My space: "; df -h .
echo -e "\n All files:"; ls -lt
echo -e "\n Total occupied space: \c"; du -h -s .

exit 0
