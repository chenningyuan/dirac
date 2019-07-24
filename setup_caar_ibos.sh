./setup --fc=mpifort --cc=mpicc --cxx=mpicxx --python=python3 --mpi --omp --type=release --cmake-options="-DENABLE_EXATENSOR=ON -DENABLE_GENIAO=ON -DENABLE_PCMSOLVER=OFF -DENABLE_STIELTJES=OFF -DENABLE_UTILITIES=OFF -DENABLE_GEN1INT=OFF -DENABLE_PELIB=OFF -DENABLE_RELCCSD_STANDALONE=OFF" build-caar-ibos
#serial
#./setup --cc=gcc-8 --cxx=g++-8 --python=python3 --omp --cmake-options="-DENABLE_EXATENSOR=ON -DENABLE_GENIAO=ON -DENABLE_PCMSOLVER=OFF -DENABLE_STIELTJES=OFF -DENABLE_UTILITIES=OFF -DENABLE_GEN1INT=OFF -DENABLE_PELIB=OFF -DENABLE_RELCCSD_STANDALONE=OFF" build-caar-ibos-serial


# Then, execute link-librairy.sh in the build directory. Once the build is done, copy pam_from_Luuk_caar-ibos into pam. Add -hosts localhost to os.environ['DIRAC_MPI_COMMAND'] = "mpiexec -n 4".
# Or maybe just set it in the bashrc ; export DIRAC_MPI_COMMAND="mpiexec -n 4 -hosts localhost"
