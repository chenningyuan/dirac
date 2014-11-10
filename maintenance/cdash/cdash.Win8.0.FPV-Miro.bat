@ECHO OFF
cls
rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits
rem echo "Hello, building up DIRAC on Windows 7 PC Desktop"

set DIRAC=C:\Users\milias\Documents\Dirac\devel_trunk

rem enable keyring....done automatically, ie Windows 7 computer has the private key loaded
rem "C:\Program Files\PuTTY\pageant.exe"  C:\Users\milias\.ssh\id_rsa_milias_fpv.ppk

rem go to the local DIRAC directory
cd %DIRAC%

rem download fresh Dirac - MUST HAVE full path of git
C:\Git\bin\git.exe pull origin master
set BUILD=build_mingw64_i4
if exist %BUILD% rmdir /S /Q %BUILD%
python pamadm -a --builddir=%BUILD%
python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe" -D BUILDNAME="Miro-Win8_i4"  -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles"  %BUILD%
cd %BUILD%
mingw32-make Experimental

rem go to the local DIRAC directory
cd %DIRAC%
set BUILD=build_mingw64_i8_check
if exist %BUILD% rmdir /S /Q %BUILD%
python pamadm -a --builddir=%BUILD%
python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe" -D BUILDNAME="Miro-Win8.0-i8_check"  --check --int64=on -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=OFF --generator="MinGW Makefiles"  %BUILD%
cd %BUILD%
mingw32-make Experimental

rem return to the top of the DIRAC directory
cd %DIRAC%
