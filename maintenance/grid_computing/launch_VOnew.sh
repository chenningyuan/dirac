#!/bin/sh
#######################################################################################################
#
#
#           Script for copying and launching the DIRAC suite on grid CE
#
#
#######################################################################################################

# check the parameter - VO name
if [[ $1 != "voce" && $1 != "compchem" && $1 != "isragrid" && $1 != "osg" && $1 != "sivvp.slovakgrid.sk" ]]; then
 echo -e "\n wrong parameter - VO name : $1 "
 exit 12
else
 VO=$1
 echo -e "OK, you specified properly the VO=$VO, continuing \n"
fi

# include all external functions
if [ -e "UtilsCE.sh" ]
then
  source ./UtilsCE.sh
  #declare -F # list function
else
  echo "Source file UtilsCE not found! Error exit !"
  exit 333
fi

print_CE_info
querry_CE_attributes $VO
check_file_on_SE $VO
# download & unpack tar-file onto CE - MUST be successfull or exit
download_from_SE $VO "DIRAC_grid_suite.tgz"

# get number of procs #
unset nprocs
get_nprocs_CE nprocs
#RETVAL=$?; [ $RETVAL -ne 0 ] && exit 5
echo -e "\n Number of #CPU obtained from the function: $nprocs \n"

##########################################
#      set build dirs and paths          #
##########################################
# directory with all static executables - dirac.x, cmake, ctest, mpirun - and related directories
export BUILD=build_intel_openmpi_static
export BUILD_SERIAL=build_serial_intel_static

unset OPAL_PREFIX
export OPAL_PREFIX=$PWD/$BUILD
echo -e "Variable OPAL_PREFIX=$OPAL_PREFIX"
export PATH=$PWD/$BUILD:$PATH
#export PATH=$PWD/$BUILD_SERIAL:$PWD/$BUILD:$PATH
echo -e "Modified PATH=$PATH"
echo -e "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

unpack_DIRAC
#RETVAL=$?; [ $RETVAL -ne 0 ] && exit 6

#############################################################
#           check own static cmake, mpirun stuff
#############################################################
echo -e "\n\n"
echo -e "own mpirun in PATH ?\c"; which mpirun; mpirun --version

echo -e "\n --- Launching simple parallel pam test  --- \n "; 
python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --outcmo --mpi=$nprocs --dirac=$BUILD/dirac.x

#####################################################################
#           specify tests for two runtest runs
#####################################################################

regular_subset='fscc fscc_highspin fscc_intensities krci_energy krci_properties_omega_tdm krci_properties_perm_dipmom cosci_energy acmoin gaunt basis_automatic_augmentation basis_input bss_energy density_at_nuclei dft_ac dft_ao_eval dft_response_symmetry fde_static-vemb_dipole lucita_large lucita_short  mcscf_energy  molecular_mean_field_hamiltonian_schemes mp2_energy mp2_properties  projection_analysis response_lao_shield  symmetry_recognition symmetry_scf_energy visual_lao_shielding ecp'

benchmark_subset='cc_benchmark cc_linear fscc_intensities gaunt mcscf_energy moltra_4or6'

  echo -e "\n\n --- Going to launch parallel runtest - OpenMPI+Intel+MKL - with quick tests suite  --- \n "; date 
 # python ./runtest  --verbose --backup --builddir=$PWD/$BUILD --mpi=$nprocs --tests="$regular_subset"
 # python ./runtest  --verbose --backup --builddir=$PWD/$BUILD --mpi=$nprocs --quick

  echo -e "\n\n --- Going to launching selected serial runtest - Intel+MKL - with quick test suite --- \n "; date 
  python ./runtest --verbose --backup --builddir=$PWD/$BUILD_SERIAL  --tests="cosci_energy"
  #python ./runtest --verbose --backup --builddir=$PWD/$BUILD_SERIAL  --quick

  #echo -e "\n\n --- Launching benchmark serial runtest - Intel+MKL --- \n "; date 
  #python ./runtest --all --benchmark --verbose --backup --binary=$PWD/$BUILD_SERIAL/dirac.x --tests="$benchmark_subset"

##############################################################
#                                                            #
#      pack selected files to get them back from CE          #
#                                                            #
##############################################################
echo -e "\n --------------------------------- \n "; 
# delete old tar-ball first
ls -lt DIRAC_grid_suite.tgz
echo -e "\n deleting the old DIRAC_grid_suite.tgz..."
rm dirac_grid_suite.tgz
echo "check files..."
ls -lt
echo -e "\n --- Packing all wanted stuff back from the grid CE  --- "; 
tar --version
echo -e "\n we have to pack (ls -lt) :"
ls -lt
echo " "

tar czf DIRAC_grid_suite_back.tgz test *.out
echo -e "\n selected directories/files of the DIRAC suite packed back, ls -lt DIRAC_grid_suite_back.tgz:";ls -lt DIRAC_grid_suite_back.tgz

# upload final tarball onto SE so that you can dowload it later
# upload_to_SE $VO

#############################################
#### flush out some good-bye message ... ####
#############################################
final_message

exit 0

