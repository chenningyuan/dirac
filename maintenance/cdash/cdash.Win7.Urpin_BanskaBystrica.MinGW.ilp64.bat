@ECHO OFF

cls

rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits

echo "hello, building up DIRAC on Urpin, KCH FPV UMB"

set DIRAC=C:\Users\milias\Documents\Work\QCH\svn-dirac\trunk

rem go to the local DIRAC directory
cd %DIRAC%

rem radovan: bad idea, use svn co, not svn up
svn up

python pamadm -a

rem configure step
cmake -G "MinGW Makefiles"  -DENABLE_64BIT_INTEGERS=1 -DCMAKE_BUILD_TYPE=Debug  %DIRAC%

cd %DIRAC%\build
rem make clean
rem make dirac.x
rem mingw32-make Experimental
mingw32-make -j 2


