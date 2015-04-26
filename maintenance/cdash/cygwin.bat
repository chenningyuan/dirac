@echo off

C:\cygwin64\bin\sh --login -c "cd /home/milias/test/ && git clone git@gitlab.com:dirac/dirac.git trunk_cloned && cd /home/milias/test/trunk_cloned  && git submodule update --init && python setup --int64 --mpi --cc=mpicc --cxx=mpicxx --fc=mpifort build && cd build && make -j 4 && ctest -j $(expr $(cat /proc/cpuinfo | grep processor | wc -l) * 2) -L short && logout"