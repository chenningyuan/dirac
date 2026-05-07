:orphan:
 

Installation on Zeus HPC cluster in Lille (Last tested 05/05/2026)
=============================================================================================

Configuration and compilation
-----------------------------

Make sure to be in dirac folder before running any setup script. A GNU build (with openblas?) will work just fine, but we will demonstrate with intel compilers.
An example of script, with intel compilers with MKL, to set up DIRAC is:

.. code-block:: bash

   module purge
   module load cmake/3.29.2
   module load intel/2024/compilers
   
   # if compiling with cuda, also define e.g.
   
   export WRAP=NOWRAP
   export BUILD_TYPE=OPT
   export EXA_TALSH_ONLY=YES
   export CMAKE_Fortran_COMPILER=/share/compilers/intel/2024/compiler/2024.0/bin/ifx
   export CMAKE_C_COMPILER=/share/compilers/intel/2024/compiler/2024.0/bin/icx
   export CMAKE_CXX_COMPILER=/share/compilers/intel/2024/compiler/2024.0/bin/icpx
   export TOOLKIT=INTEL
   export EXA_OS=LINUX
   export BLASLIB=MKL
   export EXA_GPU=NONE
   export MPILIB=NONE
   export FOLDER=build_oneapi2024
   rm -rf $FOLDER
   
   ./setup --cmake-options="-DCMAKE_VERBOSE_MAKEFILE=ON -DENABLE_PCMSOLVER=OFF -DENABLE_PELIB=OFF" --fc=ifx --cc=icx --cxx=icpx $FOLDER
   cd $FOLDER
   make -j4

Testing the installation
-----------------------------

After successfully building DIRAC, you can test the installation by using ctest:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --nodes=1
   #SBATCH --ntasks-per-node=1
   #SBATCH --cpus-per-task=12
   #SBATCH --time=10:00:00
   #SBATCH --mem=100G
   
   ########################################################
   # modules to use

   module purge
   module load cmake/3.29.2
   module load intel/2024/compilers
   
   ########################################################
   # DIRAC setup
   #
   export DIRACBIN=<PATH-INSIDE-THE-BUILD-DIRECTORY>
   export DIRACPAM=${DIRACBIN}/pam
   
   export OMP_NUM_THREADS=12
   export OMP_DYNAMIC=FALSE
   export OMP_MAX_ACTIVE_LEVELS=3
   export OMP_THREAD_LIMIT=256
   export OMP_WAIT_POLICY=PASSIVE
   export OMP_PROC_BIND="spread,spread,spread"
   export MKL_NUM_THREADS=12
   export MKL_DYNAMIC=FALSE
   
   #
   ########################################################
   
   cd ${DIRACBIN}
   
   export DIRAC_MPI_COMMAND=""

   # For only executing tests related to exacorr only...
   ctest -VV -R exacorr

