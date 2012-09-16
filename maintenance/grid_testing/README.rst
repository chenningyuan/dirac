Launching DIRAC in the grid environment
=======================================

Testing of the DIRAC program suite in the grid enviroment.

DIRAC is packed and the tar-ball is copied onto chosen VO location,
from where an assigned CE unpacks it and performs own calculations.

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
 http://bugs.python.org/msg161300

Download Python-2.7.3

::

  $./configure --prefix=/home/ilias/bin/python_static  LDFLAGS="-static -static-libgcc" CPPFLAGS="-static"

put "\*static\*" into Modules/Setup.dist and type ::

  $ make -i install
  $ make install

------------------------------
Working with your certificate
------------------------------

Initialize: 
::
  voms-proxy-init --voms voce -hours 24 -vomslife 24:00 --out ~/voce_cert; export X509_USER_PROXY=~/voce_cert; voms-proxy-info --all
  voms-proxy-init --voms compchem -hours 240 --out ~/compchem_cert; export X509_USER_PROXY=~/compchem_cert
  voms-proxy-init --voms osg  -hours 24 -vomslife 24:00 --out ~/osg_cert; export X509_USER_PROXY=~/osg_cert; voms-proxy-info --all

  myproxy-init  --pshost myproxy.cnaf.infn.it # compchem
  myproxy-init  --pshost px.ui.savba.sk       # voce
  myproxy-init -c 504 -k miro2 -l ilias -x -R "/C=SK/O=SlovakGrid/O=UMB/CN=Miroslav Ilias" # initialize for 504 hours

Repeat for the certificate regeneration:
::
  myproxy-get-delegation -v -l ilias -k miro2 -a /tmp/x509up_u501
  voms-proxy-init --noregen --voms voce


List of all valid proxy...
::
 voms-proxy-info --all

Destroy proxy: 
::
  voms-proxy-destroy
  myproxy-destroy -v -l ilias -k miro2


------------------
VO Data management
-------------------

Accesible nodes/storage space:

::

  lcg-infosites -vo compchem se
  lcg-infosites -vo compchem all

What you have in your VO space:

::

  lcg-ls -l  lfn://grid/voce/ilias/
  lcg-ls lfn://grid/compchem/ilias
  lcg-ls lfn://grid/compchem
 
 copy executable to your compchem space:
:: 

  lcg-cr -d se.grid.unipg.it  -l  lfn://grid/compchem/ilias/dirac.x --vo compchem  dirac.x
  guid:9198887f-a283-486e-af4f-fb0ed4649dab

even better:

::

  lcg-cr -d se.grid.unipg.it  -l  lfn://grid/compchem/ilias/dirac.x --vo compchem  ../../build/dirac.x
  guid:e3b25fad-ad2a-4a01-9b9b-dd2472e86f0c

Delete data from your SE space:

::

  lcg-del -a   lfn:/grid/compchem/ilias/dirac.x

Copy files from SE onto current directory:

:: 

 lcg-cp  lfn://grid/voce/ilias/dirac_current.tgz      file://$PWD/dirac_current.tgz
 lcg-cp  lfn://grid/compchem/ilias/dirac_current.tgz  file://$PWD/dirac_current.tgz


-----------------------------------
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

Get job files back

::

 glite-wms-job-output -i <JOB_ID_file>

Querry computing elements on specific software:

::

 lcg-info  --list-ce  --query 'LRMS=pbs' --vo voce
 lcg-info  --list-ce  --query 'LRMS=pbs' --vo compchem
 lcg-info  --list-ce  --query 'LRMS=pbs' --vo osg
 lcg-info  --list-ce --query 'TotalCPUs>=8' --vo voce
 lcg-info  --list-ce --query 'TotalCPUs>=24,FreeCPUs>=5' --vo compchem
 lcg-info  --list-ce --query 'TotalCPUs>=24,FreeCPUs>=5,FreeJobSlots>=2' --vo voce
 lcg-info --list-ce  --query 'CE=*' --attrs EstRespTime,TotalCPUs,Memory,ClockSpeed,Cluster --vo voce
 lcg-info --list-ce  --query 'CE=*' --attrs EstRespTime,MaxCPUTime,TotalCPUs,Memory,ClockSpeed,MaxTotalJobs,Cluster  --vo voce
 lcg-info --list-ce  --query 'CE=*' --attrs EstRespTime,MaxCPUTime,TotalCPUs,Memory,ClockSpeed,Cluster,VMemory   --vo compchem
 lcg-info --list-ce --attrs MaxWCTime --vo voce
 lcg-info --list-ce --attrs RunningJobs,FreeCPUs,MaxWCTime,MaxCPUTime --vo voce
 lcg-info --list-ce --attrs PlatformArch --vo voce
 lcg-info --list-ce --query 'PlatformArch=x86_64' --vo voce

Querry tag attributes :

::

 lcg-info --list-ce --query 'Tag=*MPICH*' --attrs 'CE' --vo voce
 lcg-info --list-ce --query 'Tag=*GCC*'   --attrs 'CE' --vo voce

