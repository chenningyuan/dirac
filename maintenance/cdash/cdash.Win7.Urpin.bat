@ECHO ON

rem
rem True DOS written script for the DIRAC buildup !
rem

rem cls

rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits

set DIRAC=C:\Users\milias\Documents\Work\QCH\dirac_git\trunk

set BUILD=build_MinGW64_ilp64_serial

set WRK=C:\Users\milias\tmp\milias

rem go to the local DIRAC directory
cd %DIRAC%

rem initialize GIT_SSH to make PuTTy keyring work - 
rem set GIT_SSH=C:\Program Files (x86)\PuTTY\plink.exe
rem echo "I have GIT_SSH=%GIT_SSH%"
rem echo "I have SVN_SSH=%SVN_SSH%"
rem echo "I have PLINK_PROTOCOL=%PLINK_PROTOCOL%"

rem enable keyring....
rem "C:\Program Files (x86)\PuTTY\pageant.exe"  C:\Users\milias\.ssh\milias_private_key.ppk
rem cmd.exe \c "C:\Users\milias\AppData\Roaming\Microsoft\Windows\Start Menu\Programs\Startup\launch_key.bat" 1>> %OUT% 2>&1
rem "C:\Users\milias\AppData\Roaming\Microsoft\Windows\Start Menu\Programs\Startup\launch_key.bat" 1>> %OUT% 2>&1

rem check the server connectivity:
rem plink -agent git@repo.ctcc.no 

C:\Git\bin\git.exe clean -f

rem download fresh Dirac - MUST HAVE full path of git
rem C:\Git\bin\git.exe checkout release-13 
rem C:\Git\bin\git.exe pull origin release-13
C:\Git\bin\git.exe pull

rem show 3 last commits...
rem C:\Git\bin\git.exe log -3

rem remove previously generated build-directory, if exists
if exist %BUILD% rmdir /S /Q %BUILD%
python pamadm -a --builddir=%BUILD%
python setup.py -D GIT_EXECUTABLE="C:\Git\bin\git.exe"  -D BUILDNAME="Urpin-Win7-MinGW64-ILP64-release13"   --generator="MinGW Makefiles" --int64 %BUILD%
cd %BUILD%
rem mingw32-make -j4 Experimental
ctest -j4 -D Experimental

rem remove the scratch space, may contain some rubbish...
if exist %WRK% rmdir /S /Q %WRK%

rem return to the top of the DIRAC directory
rem cd %DIRAC%
rem python setup.py -D GIT_EXECUTABLE="C:\Git\bin\git.exe"  -D BUILDNAME="Urpin_Win7-MinGW64-ILP64-debug-release" --debug  -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --int64=on --force %BUILD_DBG%
rem cd %BUILD_DBG%
rem mingw32-make Experimental


rem C:\Git\bin\git.exe pull
rem return to the top of the DIRAC directory
rem cd %DIRAC%
rem python setup.py --check -D GIT_EXECUTABLE="C:\Git\bin\git.exe" -D BUILDNAME="Urpin_Win7-MinGW64-I32LP64-check_master" -D DART_TESTING_TIMEOUT=99999 -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON  --generator="MinGW Makefiles" --int64=off --force %BUILD_I32LP64_CHECK%
rem cd %BUILD_I32LP64_CHECK%
rem cmake ..
rem mingw32-make clean
rem mingw32-make Experimental

exit 