cd /home/USER/software/DIRAC
python ./pamadm -a --builddir=build
python ./setup -D BUILDNAME="open_mpi-gnu.i8-cygwin" --int64 --mpi --cc=mpicc --cxx=mpicxx --fc=mpifort build
cd build
ctest --track Miro -D ExperimentalConfigure
ctest -j $(cat /proc/cpuinfo | grep processor | wc -l) --track Miro -D ExperimentalBuild
ctest -j $(expr $(cat /proc/cpuinfo | grep processor | wc -l) / 2) --track Miro -D ExperimentalTest
ctest --track Miro -D ExperimentalSubmit