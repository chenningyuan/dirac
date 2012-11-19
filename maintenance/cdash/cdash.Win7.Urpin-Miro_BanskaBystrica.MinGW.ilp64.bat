@ECHO OFF
cls
rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits
echo "Hello, building up DIRAC on Windows 7 server Urpin FPV-UMB"

set DIRAC=C:\Users\milias\Documents\Work\QCH\dirac_git\trunk

set BUILD=build_MinGW64_ilp64_serial
set BUILD_DBG=build_MinGW64_ilp64_serial_debug
set BUILD_STATIC=build_MinGW64_ilp64_serial_debug_static

rem go to the local DIRAC directory
cd %DIRAC%

rem enable keyring....
rem "C:\Program Files (x86)\PuTTY\pageant.exe"  C:\Users\milias\Documents\milias_private_key.ppk

rem download fresh Dirac - MUST HAVE full path of git
C:\Git\bin\git.exe pull origin release-12

python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe"  -D BUILDNAME="Urpin_Win7-MinGW64-ILP64-release" -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON -D RUNTEST_BINARY=%DIRAC%\%BUILD%\dirac.x.exe  --generator="MinGW Makefiles" --int64=on --force %BUILD%
cd %BUILD%
mingw32-make Experimental

rem return to the top of the DIRAC directory
cd %DIRAC%
python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe"  -D BUILDNAME="Urpin_Win7-MinGW64-ILP64-debug-release" --debug -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON -D RUNTEST_BINARY=%DIRAC%\%BUILD_DBG%\dirac.x.exe  --generator="MinGW Makefiles" --int64=on --force %BUILD_DBG%
cd %BUILD_DBG%
mingw32-make Experimental

rem return to the top of the DIRAC directory
cd %DIRAC%
python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe"  -D BUILDNAME="Urpin_Win7-MinGW64-ILP64-debug-static-release" --debug --static -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON -D RUNTEST_BINARY=%DIRAC%\%BUILD_STATIC%\dirac.x.exe  --generator="MinGW Makefiles" --int64=on --force %BUILD_STATIC%
cd %BUILD_STATIC%
mingw32-make Experimental
