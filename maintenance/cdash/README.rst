
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

We give here few examples of the crontab prescrition of buildup on various servers.

pd.uniza.sk:
05 00 * * 1,2,3,4,5,6,7 DIRAC=/home/ilias/QCH_Work/qch_softs/dirac_space/trunk_cloned; if [ -d "$DIRAC" ]; then echo "Deleting old directory $DIRAC"; /bin/rm -rf $DIRAC; fi; git clone git@repo.ctcc.no:dirac.git $DIRAC; cd $DIRAC; ./setup -DVERBOSE_OUTPUT=ON --int64 --debug; cd build; make Experimental

