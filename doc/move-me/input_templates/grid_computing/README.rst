 DIRAC in grid environment
============================

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

 copy big executable to your compchem space:
ilias@grid3.ui.savba.sk:~/qch_work/software/dirac_git/trunk/input_templates/grid_computing/.lcg-cr -d se.grid.unipg.it  -l  lfn://grid/compchem/ilias/dirac.x --vo compchem  dirac.x
guid:9198887f-a283-486e-af4f-fb0ed4649dab

even better:
ilias@grid3.ui.savba.sk:~/qch_work/software/dirac_git/trunk/input_templates/grid_computing/.lcg-cr -d se.grid.unipg.it  -l  lfn://grid/compchem/ilias/dirac.x --vo compchem  ../../build/dirac.x
guid:e3b25fad-ad2a-4a01-9b9b-dd2472e86f0c

copy pam:
ilias@grid3.ui.savba.sk:~/qch_work/software/dirac_git/trunk/input_templates/grid_computing/.lcg-cr -d se.grid.unipg.it  -l lfn://grid/compchem/ilias/pam  --vo compchem  ../../pam
guid:1f3f225c-27ef-4f5a-ae87-6d2706ba2198

Delete from the SE space:
 lcg-del -a   lfn:/grid/compchem/ilias/dirac.x

Working with the "gLite" middleware
---------------------------------------

Some "gLite" howtos 
 http://egee-uig.web.cern.ch/egee-uig/production_pages/SimpleJobCycle.html
 http://iag.iucc.ac.il/workshop/complex_jobs.htm

Submit job script:
glite-wms-job-submit -a submit.jdl

Get job status:
glite-wms-job-status  <JOB_ID>

Get job files:
glite-wms-job-output  <JOB_ID>

