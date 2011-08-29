@ECHO OFF

cls

rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits

echo "hello, building up DIRAC"
rem set DIRAC=C:\Users\dv6-2145ec\Dirac-Windows\gfortran_gcc_ilp64_ownmath\trunk
rem set DIRAC=C:\Users\dv6-2145ec\Documents\MiroIlias_osobne\dirac_software\git_downloaded\trunk

rem go to the local DIRAC directory
cd %DIRAC%

git pull origin master

python pamadm -a 

rem delete/recreate directory
rmdir /S /Q %DIRAC%\build
mkdir %DIRAC%\build
cd %DIRAC%/build
rem configure step
cmake -G "MinGW Makefiles" "FC=gfortran"  -DENABLE_64BIT_INTEGERS=1 -DCMAKE_BUILD_TYPE=Debug -DVERBOSE_OUTPUT=ON %DIRAC%

rem make clean
rem make dirac.x
mingw32-make dirac.x


