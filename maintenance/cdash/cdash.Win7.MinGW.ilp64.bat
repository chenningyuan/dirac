@ECHO OFF

cls

rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits

echo "hello, building up DIRAC"
set DIRAC=C:\Users\dv6-2145ec\Dirac-Windows\trunk

rem go to the local DIRAC directory
cd %DIRAC%

rem radovan: bad idea, use svn co, not svn up
svn up

python pamadm -a --xcfun 

rem configure step
cmake -G "MinGW Makefiles" "FC=gfortran"  -DENABLE_64BIT_INTEGERS=1 -DCMAKE_BUILD_TYPE=Debug  %DIRAC%

cd %DIRAC%\build
rem make clean
rem make dirac.x
mingw32-make dirac.x


