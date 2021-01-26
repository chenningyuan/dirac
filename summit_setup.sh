#!/bin/bash

source ~/.modules

python setup --mpi --omp --type=release --cmake-options="-DENABLE_EXATENSOR=ON -DENABLE_LAPLACE=OFF -DENABLE_PCMSOLVER=OFF -DENABLE_STIELTJES=OFF -DENABLE_GEN1INT=OFF -DENABLE_PELIB=OFF" build_gcc811_essl_cuda_release
