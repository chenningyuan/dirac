 DIRAC in grid environment
============================

The aim of this folder is to provide user with some teplates for so-called 
grid-computing.

Making dirac.x static
---------------------

Don't forget to include "--static" to your setup flags.

Making Python static
--------------------

Download Python-2.7.2

  $./configure --prefix=/home/ilias/bin/python_static  LDFLAGS="-static -static-libgcc" CPPFLAGS="-static"

put "*static*" into Modules/Setup.local


Working with "gLite" middleware
--------------------------------

Some "gLite" howtos 
 http://egee-uig.web.cern.ch/egee-uig/production_pages/SimpleJobCycle.html
 http://iag.iucc.ac.il/workshop/complex_jobs.htm

What you have in your VO space:
lcg-ls -l  lfn://grid/voce/ilias/

Submit job script:
glite-wms-job-submit -a submit.jdl

Get job status:
glite-wms-job-status  <JOB_ID>

Get job files:
glite-wms-job-output  <JOB_ID>


