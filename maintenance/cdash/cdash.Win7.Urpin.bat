@ECHO ON

rem
rem True DOS written script for the DIRAC buildup !
rem

rem cls

rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits

set DIRAC=C:\Users\milias\Documents\Work\QCH\dirac_git\trunk

set BUILD=build_MinGW64_ilp64_serial
set BUILD_DBG=build_MinGW64_ilp64_serial_debug
set BUILD_STATIC=build_MinGW64_ilp64_serial_debug_check_static
set BUILD_I32LP64=build_MinGW64_i32lp64_serial
set BUILD_I32LP64_CHECK=build_MinGW64_i32lp64_serial_check

set OUT=%DIRAC%\build.log

echo "Hello, building up DIRAC on Windows 7 server Urpin FPV-UMB" 1> %OUT% 2>&1

rem go to the local DIRAC directory
cd %DIRAC%

echo %cd%  1>>  %OUT% 2>&1

rem initialize GIT_SSH to make PuTTy keyring work - 
rem set GIT_SSH=C:\Program Files (x86)\PuTTY\plink.exe
rem echo "I have GIT_SSH=%GIT_SSH%"
rem echo "I have SVN_SSH=%SVN_SSH%"
rem echo "I have PLINK_PROTOCOL=%PLINK_PROTOCOL%"

rem enable keyring....
rem "C:\Program Files (x86)\PuTTY\pageant.exe"  C:\Users\milias\.ssh\milias_private_key.ppk  1>> %OUT% 2>&1
rem cmd.exe \c "C:\Users\milias\AppData\Roaming\Microsoft\Windows\Start Menu\Programs\Startup\launch_key.bat" 1>> %OUT% 2>&1
rem "C:\Users\milias\AppData\Roaming\Microsoft\Windows\Start Menu\Programs\Startup\launch_key.bat" 1>> %OUT% 2>&1
rem "C:\Users\milias\AppData\Roaming\Microsoft\Windows\Start Menu\Programs\Startup\launch_key.bat"

rem check the server connectivity:
plink -agent git@repo.ctcc.no 1>> %OUT% 2>&1 

rem download fresh Dirac - MUST HAVE full path of git
rem C:\Git\bin\git.exe pull origin release-12
C:\Git\bin\git.exe pull 1>> %OUT% 2>&1


rem python pamadm -a --builddir=%BUILD%
C:\Git\bin\git.exe clean -f 1>> %OUT% 2>&1
rem show very last commit...
C:\Git\bin\git.exe log -1 1>> %OUT% 2>&1

rmdir \S \Q %BUILD%
python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe"  -D BUILDNAME="Urpin_Win7-MinGW64-ILP64-master"   --generator="MinGW Makefiles" --int64 %BUILD% 1>> %OUT% 2>&1
cd %BUILD%
rem cmake ..                    1>> %OUT% 2>&1
rem mingw32-make clean          1>> %OUT% 2>&1
mingw32-make Experimental   1>> %OUT% 2>&1

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

rem C:\Git\bin\git.exe pull 
rem return to the top of the DIRAC directory
rem cd %DIRAC%
rem python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe" -D BUILDNAME="Urpin_Win7-MinGW64-I32LP64-master" -D VERBOSE_OUTPUT=ON  --generator="MinGW Makefiles" --int64=off --force %BUILD_I32LP64%
rem cd %BUILD_I32LP64%
rem cmake ..
rem mingw32-make clean
rem mingw32-make Nightly

rem C:\Git\bin\git.exe pull
rem return to the top of the DIRAC directory
rem cd %DIRAC%
rem python setup --check -D GIT_EXECUTABLE="C:\Git\bin\git.exe" -D BUILDNAME="Urpin_Win7-MinGW64-I32LP64-check_master" -D DART_TESTING_TIMEOUT=99999 -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON  --generator="MinGW Makefiles" --int64=off --force %BUILD_I32LP64_CHECK%
rem cd %BUILD_I32LP64_CHECK%
rem cmake ..
rem mingw32-make clean
rem mingw32-make Experimental

exit 