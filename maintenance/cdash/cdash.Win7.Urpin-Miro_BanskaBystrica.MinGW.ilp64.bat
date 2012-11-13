@ECHO OFF

rem
rem True DOS written script for the DIRAC buildup !
rem

rem cls

rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits
rem echo "Hello, building up DIRAC on Windows 7 server Urpin FPV-UMB"

set DIRAC=C:\Users\milias\Documents\Work\QCH\dirac_git\trunk
set BUILD=build
rem set BUILD_STATIC=build_gnu_i8_static

rem go to the local DIRAC directory
cd %DIRAC%

rem initialize GIT_SSH to make PuTTy keyring work
set GIT_SSH=C:\Program Files (x86)\PuTTY\plink.exe
rem echo "I have GIT_SSH=%GIT_SSH%"
rem echo "I have SVN_SSH=%SVN_SSH%"
rem echo "I have PLINK_PROTOCOL=%PLINK_PROTOCOL%"

rem enable keyring....
rem "C:\Program Files (x86)\PuTTY\pageant.exe"  C:\Users\milias\.ssh\milias_private_key.ppk

rem check the server connectivity:
rem plink -agent git@repo.ctcc.no

rem download fresh Dirac - MUST HAVE full path of git
C:\Git\bin\git.exe pull 

rem python setup -D BUILDNAME="Urpin-Win7-MinGW-ILP64" -D RUNTEST_TIMEOUT="00:18:00" -D RUNTEST_BINARY=%DIRAC%\%BUILD%\dirac.x.exe -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --int64=on --force %BUILD%
python setup -D BUILDNAME="Urpin-Win7-MinGW-ILP64" -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --int64=on --force %BUILD%%

cd %BUILD%
mingw32-make VERBOSE=1 Experimental

rem return to the top DIRAC directory
cd %DIRAC%