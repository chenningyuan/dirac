:orphan:
 

Installation on Irene HPC cluster (TGCC) (Last tested 05/05/2026)
=============================================================================================

The documentation website can be found here: https://hpc.cea.fr/tgcc-public/en/html/tgcc-public.html.

Configuration and compilation
-----------------------------

Make sure to be in dirac folder before running any setup script.
An example of script, working on AMD Rome and V100  nodes, is as follow:

.. code-block:: bash

    module purge
    module load datadir/<DATADIR>
    module load feature/mkl/multi-threaded
    module load intel/20.0.0
    module load mkl/20.0.0
    module load git
    module load cmake/3.26.4
    module load flavor/cuda/nvhpc-227
    module load cuda/11.7
    module load python3/3.8.10
    module load hdf5/1.12.0
    
    export EXA_GPU=NVIDIA
    export WRAP=NOWRAP
    export TOOLKIT=INTEL
    export BUILD_TYPE=OPT
    export MPILIB=NONE
    export BLASLIB=MKL
    export PATH_BLAS_MKL=${MKLROOT}/lib/intel64
    export PATH_CUDA=/ccc/products/NVHPC-22.7/system/default/Linux_x86_64/22.7/cuda/11.7/
    export GPU_CUDA=CUDA
    export CUDA_SM_ARCH=70
    export EXA_OS=LINUX
    export CMAKE_Fortran_COMPILER=/ccc/products/ifort-20.0.0/system/default/20.0.0/bin/intel64/ifort
    export CMAKE_C_COMPILER=/ccc/products/icc-20.0.0/system/default/20.0.0/bin/intel64/icc
    export CMAKE_CXX_COMPILER=/ccc/products/icc-20.0.0/system/default/20.0.0/bin/intel64/icpc
    export PATH_CUDA_INC="/ccc/products/nvhpc-22.7/system/default/Linux_x86_64/22.7/cuda/11.7/include -I/ccc/products/nvhpc-22.7/system/default/Linux_x86_64/22.7/cuda/11.7/math_libs/include/"
    export PATH_CUDA_LIB="/ccc/products/nvhpc-22.7/system/default/Linux_x86_64/22.7/cuda/11.7/lib64 -L/ccc/products/nvhpc-22.7/system/default/Linux_x86_64/22.7/cuda/11.7/math_libs/lib64"
    export CMAKE_EXE_LINKER_FLAGS='-L/ccc/products/nvhpc-22.7/system/default/Linux_x86_64/22.7/math_libs/11.7/targets/x86_64-linux/lib -L/ccc/products/nvhpc-22.7/system/default/Linux_x86_64/22.7/cuda/11.7/lib64 -lcublas -lcudart -lnvToolsExt'
    
    FOLDER=build_mkl20_cuda117_V100_64_test
    rm -rf $FOLDER
    ./setup --cmake-options="-DCMAKE_EXE_LINKER_FLAGS='-L/ccc/products/nvhpc-22.7/system/default/Linux_x86_64/22.7/math_libs/11.7/targets/x86_64-linux/lib -L/ccc/products/nvhpc-22.7/system/default/Linux_x86_64/22.7/cuda/11.7/lib64 -lcublas -lcudart -lnvToolsExt' -DCMAKE_VERBOSE_MAKEFILE=ON -DENABLE_PCMSOLVER=OFF -DENABLE_PELIB=OFF -DENABLE_STIELTJES=OFF -DEXATENSOR_GIT_REPO_LOCATION=/ccc/work/cont003/<DATADIR>/gomes/computer_center_environment/relmbdev/ExaTENSOR" --fc=ifort --cc=icc --cxx=icpc $FOLDER
    cd $FOLDER
    make -j4

You should replace <DATADIR> with your local installation of EXATENSOR on Irene. (This tutorial has to be updated to show how to do so.)

Testing the installation
-----------------------------

After successfully building DIRAC, you can test the installation by using ctest:

.. code-block:: bash

    module purge
    module load datadir/<DATADIR>
    module load feature/mkl/multi-threaded
    module load intel/20.0.0
    module load mkl/20.0.0
    module load git
    module load cmake/3.26.4
    module load flavor/cuda/nvhpc-227
    module load cuda/11.7
    module load python3/3.8.10
    module load hdf5/1.12.0
    
    export EXA_GPU=NVIDIA
    
    # code for A100 gpus
    export CUDA_SM_ARCH=70
    
    ########################################################
    # DIRAC setup
    #
    export DIRACBIN=<ABSOLUTE-PATH>/build_mkl20_cuda117_V100_64_test
    export DIRACPAM=${DIRACBIN}/pam
    
    export OMP_NUM_THREADS=36
    export OMP_DYNAMIC=FALSE
    export OMP_MAX_ACTIVE_LEVELS=3
    export OMP_THREAD_LIMIT=256
    export OMP_WAIT_POLICY=PASSIVE
    export MKL_NUM_THREADS=36
    export MKL_DYNAMIC=FALSE
    export OMP_PROC_BIND="spread,spread,spread"
    
    # Define BASDIR
    export BASEDIR=/ccc/cont003/home/unilille/gomes/dirac
    export BASDIR_PATH=$BASEDIR/basis:$BASEDIR/basis_dalton:$BASEDIR/basis_ecp
    export DIRAC_TMPDIR=$<DATADIR>_CCCSCRATCHDIR/testing-dirac25_serial_intel20_cuda117_v100l
    mkdir -p $DIRAC_TMPDIR
    
    #
    ########################################################
    
    cd ${DIRACBIN}
    
    export DIRAC_MPI_COMMAND=""
    export TALSH_GPUS=1
    ctest -VV -R exacorr_talsh
