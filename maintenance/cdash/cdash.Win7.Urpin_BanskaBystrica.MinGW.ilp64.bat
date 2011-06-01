@ECHO OFF

cls

rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits

echo "Hello, building up DIRAC on Windows 7 server Urpin, KCH FPV UMB"

set DIRAC=C:\Users\milias\Documents\Work\QCH\svn-dirac\trunk

rem go to the local DIRAC directory
cd %DIRAC%

svn up
python pamadm -a 

rem delete/recreate directory
rmdir /S /Q %DIRAC%\build
mkdir %DIRAC%\build
cd %DIRAC%\build

rem own cmake configure step
cmake -G "MinGW Makefiles"  -DENABLE_64BIT_INTEGERS=1 -DCMAKE_BUILD_TYPE=Debug  %DIRAC%

rem make clean
rem make dirac.x
rem mingw32-make Experimental
rem mingw32-make  clean
mingw32-make -j 2
