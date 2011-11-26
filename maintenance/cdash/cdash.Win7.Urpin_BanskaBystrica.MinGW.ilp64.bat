@ECHO OFF

cls

rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits

echo "Hello, building up DIRAC on Windows 7 server Urpin, KCH FPV UMB"

set DIRAC=C:\Users\milias\Documents\Work\QCH\dirac_git\trunk
set BUILD=%DIRAC%\build_gnu_i8_debug

echo %BUILD%

rem go to the local DIRAC directory
cd %DIRAC%

rem enable keyring....
C:\Users\milias\pageant_link.lnk C:\Users\milias\Documents\milias_private_key.ppk

rem download fresh Dirac - MUST HAVE full path of git
C:\Git\bin\git pull

python pamadm -a 

rem delete/recreate build directory

rmdir /S /Q %BUILD%
mkdir %BUILD%
cd %BUILD%

rem own cmake configure step
cmake -G "MinGW Makefiles" -DENABLE_64BIT_INTEGERS=1 -DCMAKE_BUILD_TYPE=Debug -DBUILDNAME="Win7.Urpin_BB.MinGW.ownmath.debug.ilp64" -DVERBOSE_OUTPUT=ON %DIRAC%  -Wno-dev

rem make clean
rem make dirac.x
mingw32-make Experimental
rem mingw32-make  clean
rem mingw32-make dirac_interfaces

rem for now do only the compilation and submit info about it
rem mingw32-make -j 2
rem mingw32-make ExperimentalSubmit

rem return back
rem cd %DIRAC%\maintenance\cdash

rem return to the top of the DIRAC directory
cd %DIRAC%