@ECHO OFF

cls

rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits

echo "Hello, building up DIRAC on Windows 7 server Urpin, KCH FPV UMB"

set DIRAC=C:\Users\milias\Documents\Work\QCH\dirac_git\trunk
set BUILD=%DIRAC%\build_gnu_i8
set BUILD_STATIC=%DIRAC%\build_gnu_i8_static

echo %BUILD%

rem go to the local DIRAC directory
cd %DIRAC%

rem enable keyring....
rem "C:\Program Files (x86)\PuTTY\pageant.exe"  C:\Users\milias\Documents\milias_private_key.ppk

rem download fresh Dirac - MUST HAVE full path of git
C:\Git\bin\git.exe pull --rebase

python pamadm -a 

rem delete/recreate build directory

rmdir /S /Q %BUILD%
mkdir %BUILD%
cd %BUILD%

rem own cmake configure step - MinGW must be in path !
cmake -Wno-dev -G "MinGW Makefiles" "FC=gfortran" -DENABLE_64BIT_INTEGERS=1 -DCMAKE_BUILD_TYPE=Release -DBUILDNAME="Win7.Urpin_BB.MinGW.ownmath.debug.ilp64" -DVERBOSE_OUTPUT=ON %DIRAC%  

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
