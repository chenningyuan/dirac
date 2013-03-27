@ECHO OFF

rem
rem True DOS written script for the DIRAC buildup !
rem

cls

rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits
rem echo "Hello, building up DIRAC on Windows 7 server Urpin FPV-UMB"

set DIRAC=C:\Users\milias\Documents\Work\QCH\dirac_git\trunk

set BUILD=build_MinGW64_ilp64_serial
set BUILD_DBG=build_MinGW64_ilp64_serial_debug
set BUILD_STATIC=build_MinGW64_ilp64_serial_debug_check_static
set BUILD_I32LP64=build_MinGW64_i32lp64_serial
set BUILD_I32LP64_CHECK=build_MinGW64_i32lp64_serial_check

rem go to the local DIRAC directory
cd %DIRAC%

rem initialize GIT_SSH to make PuTTy keyring work - 
rem set GIT_SSH=C:\Program Files (x86)\PuTTY\plink.exe
rem echo "I have GIT_SSH=%GIT_SSH%"
rem echo "I have SVN_SSH=%SVN_SSH%"
rem echo "I have PLINK_PROTOCOL=%PLINK_PROTOCOL%"

rem enable keyring....
rem "C:\Program Files (x86)\PuTTY\pageant.exe"  C:\Users\milias\.ssh\milias_private_key.ppk

rem check the server connectivity:
rem plink -agent git@repo.ctcc.no

rem download fresh Dirac - MUST HAVE full path of git
rem C:\Git\bin\git.exe pull origin release-12
C:\Git\bin\git.exe pull
rem python pamadm -a --builddir=%BUILD%

python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe"  -D BUILDNAME="Urpin_Win7-MinGW64-ILP64-master"  -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON  --generator="MinGW Makefiles" --int64=on --force %BUILD%
cd %BUILD%
cmake ..
mingw32-make clean
mingw32-make Experimental

rem return to the top of the DIRAC directory
rem cd %DIRAC%
rem python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe"  -D BUILDNAME="Urpin_Win7-MinGW64-ILP64-debug-release" --debug  -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --int64=on --force %BUILD_DBG%
rem cd %BUILD_DBG%
rem mingw32-make Experimental

rem C:\Git\bin\git.exe pull 
rem return to the top of the DIRAC directory
rem cd %DIRAC%
rem python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe"  -D BUILDNAME="Urpin_Win7-MinGW64-ILP64-debug-check-static-master" --debug --static --check -D DART_TESTING_TIMEOUT=99999 -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --int64=on --force %BUILD_STATIC%
rem cd %BUILD_STATIC%
rem cmake ..
rem mingw32-make clean
rem mingw32-make Experimental

C:\Git\bin\git.exe pull 
rem return to the top of the DIRAC directory
cd %DIRAC%
python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe" -D BUILDNAME="Urpin_Win7-MinGW64-I32LP64-master" -D DART_TESTING_TIMEOUT=99999 -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON  --generator="MinGW Makefiles" --int64=off --force %BUILD_I32LP64%
cmake ..
mingw32-make clean
cd %BUILD_I32LP64%
mingw32-make Nightly

rem C:\Git\bin\git.exe pull
rem return to the top of the DIRAC directory
rem cd %DIRAC%
rem python setup --check -D GIT_EXECUTABLE="C:\Git\bin\git.exe" -D BUILDNAME="Urpin_Win7-MinGW64-I32LP64-check_master" -D DART_TESTING_TIMEOUT=99999 -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON  --generator="MinGW Makefiles" --int64=off --force %BUILD_I32LP64_CHECK%
rem cd %BUILD_I32LP64_CHECK%
rem cmake ..
rem mingw32-make clean
rem mingw32-make Experimental
