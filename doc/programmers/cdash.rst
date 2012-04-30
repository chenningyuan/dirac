
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
To get inspired have a look at the cdash example recipes in the repository
under maintenance/cdash. Then put it under crontab. Here is an example::

  # m  h  dom mon dow command
    05 00 *   *   *   /home/bast/nightly/smokey_debian-6.0_gfortran-4.4.5
    00 02 *   *   *   killall -9 dirac.x
    05 02 *   *   *   /home/bast/nightly/smokey_debian-6.0_gfortran-4.4.5_i8
    00 04 *   *   *   killall -9 dirac.x
    05 04 *   *   *   /home/bast/nightly/smokey_debian-6.0_gfortran-4.4.5_i8_mpi
    00 06 *   *   *   killall -9 dirac.x

Examples of working buildup scripts for wide range of compilers, libraries etc 
are in the DIRAC subdirectory *maintenance/cdash/*.

Alternatively, instead of keeping separate buildup script(s)  
the user can put all necessary commands into compact crontab prescription::

 MAILTO=<your@email> # if server is set for email notification
 PATH=/home/ilias/bin/openmpi_gfor_gcc_ILP64/bin:/usr/local/bin:/usr/bin:/bin:/usr/games:.:/home/ilias/bin
 LD_LIBRARY_PATH=/home/ilias/bin/openmpi_gfor_gcc_ILP64/lib
 DIRAC=/home/ilias/QCH_Work/qch_softs/dirac_space/trunk_cloned
 LBS="-L/home/ilias/bin/atlas-lapack-ilp64/ATLAS/atlas_lapack_ILP64/lib -llapack -lcblas -lf77blas -latlas"
 SETUP_CMD='./setup -DVERBOSE_OUTPUT=ON --int64 --fc=mpif90 --cc=mpicc --cxx=mpiCC -D BUILDNAME="linux.openmpi_gfortran_gcc_atlas.ilp64.opt" --explicit-libs=$LBS --int64 -D RUNTEST_NR_PROC=2 -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_TIMEOUT="00:40:00"'
 BUILD_CMD="echo $DIRAC; if [ -d $DIRAC ]; then echo 'Deleting old directory $DIRAC'; /bin/rm -rf $DIRAC; fi; echo 'cloning fresh Dirac'; git clone git@repo.ctcc.no:dirac.git $DIRAC; cd $DIRAC; echo 'launching setup'; $SETUP_CMD; cd build; make Experimental"
 # m  h  dom mon dow command
 05 00 * * 1,2,3,4,5,6,7 /bin/bash -c "$BUILD_CMD"

