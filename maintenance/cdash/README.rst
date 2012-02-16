
Contributing to the DIRAC Dashboard
-----------------------------------

Before you start, set --scratch=/your/scratch/path in ~/.diracrc.

In the case of parallel DIRAC installation, you have to set up more
parameters in the diracrc (=~/.diracrc, or .diracrc or diracrc) file, like
--machfile=$PBS_NODEFILE
--scratch=/your/scratch/path
--global-scratch-disk

It is very easy to contribute to the DIRAC Dashboard:
http://repo.ctcc.no/CDash/index.php?project=DIRAC

You don't need to install anything extra except CMake.
It is as easy as typing::

  $ ./setup [--flags] -D BUILDNAME="my-build-name"
  $ cd build
  $ make Experimental

The result of "make Experimental" appears automatically on the Dashboard.

Make Experimental takes a long time (2 hours or more). This is because we have
many very slow tests. For experimenting with make experimental you can cheat
and remove many lines in cmake/Tests.cmake (and do not rerun pamadm).  Later
you can regenerate the full cmake/Tests.cmake.

Once "make Experimental" works you can create a script to run "make Nightly".
To get inspired have a look at the cdash example recipes in the repository.


  Setting effective automatic buildup and tests
  ----------------------------------------------

It is recommended to use "git clone" scheme for repetitive downloading, building up 
and testing the entire DIRAC suite.

We give here few examples of the crontab buildup prescriptions on various servers.

Note: you have to have ~/.diracrc file ready with flags

-- pd.uniza.sk --

05 00 * * 1,3,5,7  DIRAC=/home/ilias/QCH_Work/qch_softs/dirac_space/trunk_cloned; if [ -d "$DIRAC" ]; then echo "Deleting old directory $DIRAC"; /bin/rm -rf $DIRAC; fi; git clone git@repo.ctcc.no:dirac.git $DIRAC; cd $DIRAC; ./setup -DVERBOSE_OUTPUT=ON --int64 -D DART_TESTING_TIMEOUT=99999 -D RUNTEST_TIMEOUT="00:40:00" -D BUILDNAME="pd buildup"; cd build; make Experimental

05 00 * * 2,4,6 export PATH=/home/ilias/bin/openmpi_gfor_gcc_ILP64/bin:$PATH; export LD_LIBRARY_PATH=/home/ilias/bin/openmpi_gfor_gcc_ILP64/lib:$LD_LIBRARY_PATH; DIRAC=/home/ilias/QCH_Work/qch_softs/dirac_space/trunk_cloned; if [ -d "$DIRAC" ]; then echo "Deleting old directory $DIRAC"; /bin/rm -rf $DIRAC; fi; git clone git@repo.ctcc.no:dirac.git $DIRAC; cd $DIRAC; ./setup -DVERBOSE_OUTPUT=ON --int64 --fc=mpif90 --cc=mpicc --cxx=mpiCC -D BUILDNAME="pd_uniza.linux.openmpi_gfortran_gcc_atlas.ilp64.opt" --math="-L/home/ilias/bin/atlas-lapack-ilp64/ATLAS/atlas_lapack_ILP64/lib -llapack -lcblas -lf77blas -latlas" --int64 -D RUNTEST_NR_PROC=2 -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_TIMEOUT="00:40:00"; cd build; make Experimental

-- frpd2.utcpd.sk --

05 00 * * 1,2,3,4,5,6,7  export PATH=/home/ilias/bin/openmpi_i32lp64_static/bin:$PATH; export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LIBRARY_PATH; export LD_LIBRARY_PATH=/home/ilias/bin/openmpi_i32lp64_static/lib:$LD_LIBRARY_PATH; DIRAC=/home/ilias/Work/qch/qch_software/dirac/trunk_cloned; if [ -d "$DIRAC" ]; then echo "Deleting old directory $DIRAC"; /bin/rm -rf $DIRAC; fi; git clone git@repo.ctcc.no:dirac.git $DIRAC; cd $DIRAC; ./setup --int64 --fc=mpif90 --cc=mpicc --cxx=mpiCC -D BUILDNAME="linux.OpenMPI_gfortran.gcc.ownmath.ilp64_full_STATIC" -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON -D DART_TESTING_TIMEOUT=99999 --static --internal-math -D RUNTEST_NR_PROC=2; cd build; make Experimental

03 02 * * 1,2,3,4,5,6,7  DIRAC=/home/ilias/Work/qch/qch_software/dirac/trunk_cloned; if [ -d "$DIRAC" ]; then echo "Deleting old directory $DIRAC"; /bin/rm -rf $DIRAC; fi; git clone git@repo.ctcc.no:dirac.git $DIRAC; cd $DIRAC; ./setup --int64 -D BUILDNAME="linux.gfortran.gcc.usr64_libs.i32lp64" -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON -D DART_TESTING_TIMEOUT=99999; cd build; make Nightly


