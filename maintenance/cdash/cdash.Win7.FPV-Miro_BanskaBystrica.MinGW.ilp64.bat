@ECHO OFF
cls
rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits
rem echo "Hello, building up DIRAC on Windows 7 PC Desktop"

set DIRAC=C:\Users\milias\Documents\Work\QCH\Dirac_git\trunk
set BUILD=build
set BUILD_STATIC=%DIRAC%\build_gnu_i8_static

rem enable keyring....done automatically, ie Windows 7 computer has the private key loaded
rem "C:\Program Files\PuTTY\pageant.exe"  C:\Users\milias\.ssh\id_rsa_milias_fpv.ppk

rem go to the local DIRAC directory
cd %DIRAC%

rem download fresh Dirac - MUST HAVE full path of git
C:\Git\bin\git.exe pull origin release-12

python pamadm -a 

python setup -D BUILDNAME="Miro-PC-Desktop-Win7" -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --int64 --force %BUILD%

cd %BUILD%
mingw32-make Experimental

rem return to the top of the DIRAC directory
cd %DIRAC%
