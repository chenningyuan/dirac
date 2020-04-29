#ExaCorr setup on Summit:
export WRAP=NOWRAP
export TOOLKIT=GNU
export BUILD_TYPE=OPT
export MPILIB=OPENMPI
export BLASLIB=ESSL
export GPU_CUDA=CUDA
export GPU_SM_ARCH=70
export EXA_OS=LINUX
export WITH_CUTT=NO

export PATH_OPENMPI=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-8.1.1/spectrum-mpi-10.3.1.2-20200121-chae23sgwacfeot7vxkpfboz6wao2c33

export PATH_BLAS_ESSL=/sw/summit/essl/6.1.0-2/essl/6.1/lib64
export PATH_IBM_XL_CPP=/sw/summit/xl/16.1.1-5/xlC/16.1.1/lib
export PATH_IBM_XL_FOR=/sw/summit/xl/16.1.1-5/xlf/16.1.1/lib
export PATH_IBM_XL_SMP=/sw/summit/xl/16.1.1-5/xlsmp/5.1.1/lib

export PATH_CUDA=/sw/summit/cuda/10.1.243

export PATH_CUTT=/ccs/home/div/src/cutt

#Commands for compiling:
module load gcc/8.1.1
module load spectrum-mpi
module load essl
module load cuda
module load git
module load cmake
module load python
module load netlib-lapack

python setup --mpi --omp --type=release --cmake-options="-DENABLE_EXATENSOR=ON -DENABLE_LAPLACE=OFF -DENABLE_PCMSOLVER=OFF -DENABLE_STIELTJES=OFF -DENABLE_GEN1INT=OFF -DENABLE_PELIB=OFF" build_gcc811_essl_cuda_release
