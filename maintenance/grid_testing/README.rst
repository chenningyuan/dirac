Launching  DIRAC in grid environment
======================================

Testing of the DIRAC program suite in the grid enviroment.

DIRAC is packed and the tar-ball is copied onto VO location,
from, where assigned CE unpacks it and performs calculations.

Python script is not yet used for launching jdl-initiated script due to possibility
that not all CE are equipped with appropriate python version. So we resort to the bash shell.

The aim of this folder is to provide user with some teplates for so-called 
grid-computing.

Making dirac.x static
---------------------

Don't forget to include "--static" to your setup flags.

Making Python static
--------------------

 http://stackoverflow.com/questions/1150373/compile-the-python-interpreter-statically

Download Python-2.7.2

  $./configure --prefix=/home/ilias/bin/python_static  LDFLAGS="-static -static-libgcc" CPPFLAGS="-static"

put "*static*" into Modules/Setup.local

Currently I (MI) wasn't able to get static version due to compilation problems...

------------------
Data management
------------------

Accesible nodes/storage space:
lcg-infosites -vo compchem se
lcg-infosites -vo compchem all

VOCE:
-----
What you have in your VOCE space:
lcg-ls -l  lfn://grid/voce/ilias/

COMPCHEM:
---------
What you have in your COMPCHEM space:
 list of all subdirectories: lcg-ls lfn://grid/compchem

 your own space:  lcg-ls lfn://grid/compchem/ilias
                  lcg-ls lfn://grid/compchem/ilias/dirac_space

 copy executable to your compchem space:
lcg-cr -d se.grid.unipg.it  -l  lfn://grid/compchem/ilias/dirac.x --vo compchem  dirac.x
guid:9198887f-a283-486e-af4f-fb0ed4649dab

even better:
lcg-cr -d se.grid.unipg.it  -l  lfn://grid/compchem/ilias/dirac.x --vo compchem  ../../build/dirac.x
guid:e3b25fad-ad2a-4a01-9b9b-dd2472e86f0c

Delete from the SE space:
 lcg-del -a   lfn:/grid/compchem/ilias/dirac.x

Working with the "gLite" middleware
-----------------------------------

Some "gLite" howtos :
 http://egee-uig.web.cern.ch/egee-uig/production_pages/SimpleJobCycle.html
 http://iag.iucc.ac.il/workshop/complex_jobs.htm

Submit job script: 
::
 glite-wms-job-submit -o <JOB_ID_file> -a submit.jdl

Get job status:
:: 
 glite-wms-job-status  -i <JOB_ID_file>

Get job files back:
::
 glite-wms-job-output -i <JOB_ID_file>

