conda deactivate
module purge
module load arch/a100
module load cmake/3.31.4
module load python/3.10.4
module load gcc/11.4.1
#module load openblas/0.3.20 
module load cuda/12.0.0
module load hdf5/1.12.0

# if compiling with cuda, also define e.g.
export CMAKE_EXE_LINKER_FLAGS="-I/gpfslocalsys/cuda/12.0.0/include -L/gpfslocalsys/cuda/12.1.0/lib64 -lcudart -lcublas -lnvToolsExt -lcuda"
export EXA_GPU=NVIDIA
export GPU_CUDA=CUDA
export PATH_CUDA=/gpfslocalsys/cuda/12.0.0
export CUDA_SM_ARCH=80

# EXACC variables
export WRAP=NOWRAP
export BUILD_TYPE=OPT
export EXA_TALSH_ONLY=YES
export CMAKE_Fortran_COMPILER=/usr/bin/gfortran
export CMAKE_C_COMPILER=/usr/bin/gcc
export CMAKE_CXX_COMPILER=/usr/bin/g++
export TOOLKIT=GNU
export EXA_OS=LINUX
export BLASLIB=NONE #OPENBLAS
export PATH_BLAS_OPENBLAS=/gpfslocalsup/spack_soft/openblas/0.3.20/gcc-8.5.0-aegtps5t7pceoowdakecseayj2uxl4te/lib
unset PATH_BLAS_OPENBLAS
export FCFLAGS="-I/gpfslocalsys/cuda/12.0.0/include -L/gpfslocalsys/cuda/12.0.0/lib64 -lcudart -lcublas -lnvToolsExt -lcuda"
export FFLAGS="-I/gpfslocalsys/cuda/12.0.0/include -L/gpfslocalsys/cuda/12.0.0/lib64 -lcudart -lcublas -lnvToolsExt -lcuda"

# setup command
#export FOLDER=build_gcc11_cuda12_openblas_i4_a100
export FOLDER=build_gcc11_cuda12_i4_a100
rm -rf $FOLDER
./setup --cmake-options="-DCMAKE_VERBOSE_MAKEFILE=ON -DENABLE_PCMSOLVER=OFF -DENABLE_PELIB=OFF" --fc=/usr/bin/gfortran --cc=/usr/bin/gcc --cxx=/usr/bin/g++ $FOLDER
cd $FOLDER
make -j6
