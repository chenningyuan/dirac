@ECHO OFF
cls
rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits
echo "Hello, building up DIRAC on Windows 7 server Urpin FPV-UMB"

set DIRAC=C:\Users\milias\Documents\Work\QCH\dirac_git\trunk
set DIRACx=C:/Users/milias/Documents/Work/QCH/dirac_git/trunk

set BUILD=build_MinGW64_ilp64_serial
set BUILD_DBG=build_MinGW64_ilp64_serial_debug
set BUILD_STATIC=build_MinGW64_ilp64_serial_debug_check_static
set BUILD_I32LP64=build_MinGW64_i32lp64_serial
set BUILD_I32LP64_CHECK=build_MinGW64_i32lp64_serial_check

rem go to the local DIRAC directory
cd %DIRAC%

rem enable keyring....
rem "C:\Program Files (x86)\PuTTY\pageant.exe"  C:\Users\milias\Documents\milias_private_key.ppk

rem download fresh Dirac - MUST HAVE full path of git
C:\Git\bin\git.exe pull origin release-12
python pamadm --release

python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe"  -D BUILDNAME="Urpin_Win7-MinGW64-ILP64-release"  -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON  --generator="MinGW Makefiles" --int64=on --force %BUILD%
cd %BUILD%
cmake ..
mingw32-make clean
mingw32-make Experimental

rem return to the top of the DIRAC directory
rem cd %DIRAC%
rem python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe"  -D BUILDNAME="Urpin_Win7-MinGW64-ILP64-debug-release" --debug  -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --int64=on --force %BUILD_DBG%
rem cd %BUILD_DBG%
rem mingw32-make Experimental

C:\Git\bin\git.exe pull origin release-12
rem return to the top of the DIRAC directory
cd %DIRAC%
python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe"  -D BUILDNAME="Urpin_Win7-MinGW64-ILP64-debug-check-static-release" --debug --static --check -D DART_TESTING_TIMEOUT=99999 -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --int64=on --force %BUILD_STATIC%
cd %BUILD_STATIC%
cmake ..
mingw32-make clean
mingw32-make Experimental

C:\Git\bin\git.exe pull origin release-12
rem return to the top of the DIRAC directory
cd %DIRAC%
python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe"  -D BUILDNAME="Urpin_Win7-MinGW64-I32LP64-release12" -D DART_TESTING_TIMEOUT=99999 -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON  --generator="MinGW Makefiles" --int64=off --force %BUILD_I32LP64%
cd %BUILD_I32LP64%
mingw32-make Experimental

C:\Git\bin\git.exe pull origin release-12
rem return to the top of the DIRAC directory
cd %DIRAC%
python setup --check -D GIT_EXECUTABLE="C:\Git\bin\git.exe" -D BUILDNAME="Urpin_Win7-MinGW64-I32LP64-check_release12" -D DART_TESTING_TIMEOUT=99999 -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON  --generator="MinGW Makefiles" --int64=off --force %BUILD_I32LP64_CHECK%
cd %BUILD_I32LP64_CHECK%
cmake ..
mingw32-make clean
mingw32-make Experimental

