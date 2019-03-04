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

export PATH_OPENMPI=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-8.1.0/spectrum-mpi-10.2.0.10-20181214-ywwu2vqdnrntrxmgek5vfbw42h26kmkg

export PATH_BLAS_ESSL=/sw/summit/essl/6.1.0-2/essl/6.1/lib64
export PATH_IBM_XL_CPP=/sw/summit/xl/16.1.1-1/xlC/16.1.1/lib
export PATH_IBM_XL_FOR=/sw/summit/xl/16.1.1-1/xlf/16.1.1/lib
export PATH_IBM_XL_SMP=/sw/summit/xl/16.1.1-1/xlsmp/5.1.1/lib

export PATH_CUDA=/sw/summit/cuda/9.2.148

export PATH_CUTT=/ccs/home/div/src/cutt

#python setup --mpi --omp --type=release --cmake-options="-DENABLE_EXATENSOR=ON -DENABLE_PCMSOLVER=OFF -DENABLE_STIELTJES=OFF -DENABLE_GEN1INT=OFF -DENABLE_PELIB=OFF" dirac_gcc8_essl_cuda_opt
