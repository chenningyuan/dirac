@ECHO OFF
cls
rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits
echo "Hello, building up DIRAC on Windows 7 server Urpin FPV-UMB"

set DIRAC=C:\Users\milias\Documents\Work\QCH\dirac_git\trunk
set BUILD=%DIRAC%\build
set BUILD_STATIC=%DIRAC%\build_gnu_i8_static

echo %BUILD%
rem go to the local DIRAC directory
cd %DIRAC%

rem enable keyring....
rem "C:\Program Files (x86)\PuTTY\pageant.exe"  C:\Users\milias\Documents\milias_private_key.ppk

rem download fresh Dirac - MUST HAVE full path of git
rem C:\Git\bin\git.exe pull 

python pamadm -a 

rem delete the old build directory
rmdir /S /Q %BUILD%

python setup -D BUILDNAME="Urpin-Win7-MinGW-ILP64" -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --build=%BUILD% -G "MinGW Makefiles" --int64

cd %BUILD%
mingw32-make Experimental

rem return to the top of the DIRAC directory
cd %DIRAC%
