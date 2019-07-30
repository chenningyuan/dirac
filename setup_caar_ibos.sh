# 1) Be sure to have the correct environment in your bashrc:
# for caar branch:
#export WRAP=NOWRAP
#export TOOLKIT=GNU
#export BUILD_TYPE=OPT # DEV gives faster build, but much slower code, OPT for faster code but slower build.
#export MPILIB=MPICH # or OPENMPI
#export BLASLIB=ATLAS
#export PATH_BLAS_ATLAS=/usr/local/lib
#export EXA_OS=LINUX # macos : replace by NO_LINUX
#export PATH_MPICH=/usr/local/ # or with OPENMPI: PATH_OPENMPI=/usr/local/ (or /opt/local, it really depends on your installation)
#export GPU_CUDA=NOCUDA

# 1) copy the file dirac_laptop into ~/.diracrc

# 2) execute:
./setup --fc=mpifort --cc=mpicc --cxx=mpicxx --python=python3 --mpi --omp --type=release --cmake-options="-DENABLE_EXATENSOR=ON -DENABLE_GENIAO=ON -DENABLE_PCMSOLVER=OFF -DENABLE_STIELTJES=OFF -DENABLE_UTILITIES=OFF -DENABLE_GEN1INT=OFF -DENABLE_PELIB=OFF -DENABLE_RELCCSD_STANDALONE=OFF" build-caar-ibos
#serial
#./setup --cc=gcc-8 --cxx=g++-8 --python=python3 --omp --cmake-options="-DENABLE_EXATENSOR=ON -DENABLE_GENIAO=ON -DENABLE_PCMSOLVER=OFF -DENABLE_STIELTJES=OFF -DENABLE_UTILITIES=OFF -DENABLE_GEN1INT=OFF -DENABLE_PELIB=OFF -DENABLE_RELCCSD_STANDALONE=OFF" build-caar-ibos-serial


# 3) Execute link-librairy.sh in the build directory. Once the build is done, copy pam_from_Luuk_caar-ibos into pam.

#link-librairy.sh is a file containing:
#!/bin/bash  
#echo -e "$(cat ./CMakeFiles/dirac.x.dir/link.txt)-lstdc++"   > ./CMakeFiles/dirac.x.dir/link.txt
#echo -e "$(cat ./CMakeFiles/exacorr.x.dir/link.txt)-lstdc++" > ./CMakeFiles/exacorr.x.dir/link.txt

# 4) do: 
# $> cd build
# $> make

# 5) for the program to work, replace pam by pam_from_luuk_caar-ibos:
# $> cp ../pam_from_luuk_caar-ibos pam

# That's all.
