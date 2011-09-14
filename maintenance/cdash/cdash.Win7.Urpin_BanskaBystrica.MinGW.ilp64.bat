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
rem git pull origin master

python pamadm -a 

rem delete/recreate build directory

rmdir /S /Q %BUILD%
mkdir %BUILD%
cd %BUILD%

rem own cmake configure step
cmake -G "MinGW Makefiles"  -DENABLE_64BIT_INTEGERS=1 -DCMAKE_BUILD_TYPE=Debug -DBUILDNAME="Win7.Urpin_BB.MinGW.gfortran.gcc.ownmath.ilp64" -DVERBOSE_OUTPUT=ON %DIRAC%

rem make clean
rem make dirac.x
rem mingw32-make Experimental
rem mingw32-make  clean
rem mingw32-make dirac_interfaces

rem for now do only the compilation and submit info about it
mingw32-make -j 2
rem mingw32-make ExperimentalSubmit

rem return back
rem cd %DIRAC%\maintenance\cdash
