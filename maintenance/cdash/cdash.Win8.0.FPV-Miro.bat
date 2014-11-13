@ECHO OFF
cls

set DIRTIMEOUT="14m"
set DIRAC=C:\Users\milias\Documents\Dirac\devel_trunk

rem update all submodules
git submodule update

rem download fresh Dirac
rem C:\Git\bin\git.exe pull origin master

rem go to the local DIRAC directory
cd %DIRAC%
set BUILD=build_mingw64_i4
if exist %BUILD% rmdir /S /Q %BUILD%
python pamadm -a --builddir=%BUILD%
python setup -D BUILDNAME="Win8.0_MinGW64_i4" -D DART_TESTING_TIMEOUT=99999 --generator="MinGW Makefiles" %BUILD%
cd %BUILD%
ctest     -D ExperimentalUpdate     --track Miro
ctest     -D ExperimentalConfigure  --track Miro
ctest -j2 -D ExperimentalBuild      --track Miro
ctest -j2 -D ExperimentalTest       --track Miro
ctest     -D ExperimentalSubmit     --track Miro

rem wait 1 minute
timeout 60

rem go to the local DIRAC directory
cd %DIRAC%
set BUILD=build_mingw64_i8_check
if exist %BUILD% rmdir /S /Q %BUILD%
python pamadm -a --builddir=%BUILD%
python setup -D BUILDNAME="Win8.0_MinGW64_i8_check" --check --int64 -D DART_TESTING_TIMEOUT=99999 --generator="MinGW Makefiles" %BUILD%
cd %BUILD%
ctest     -D ExperimentalUpdate     --track BoundsCheck
ctest     -D ExperimentalConfigure  --track BoundsCheck
ctest -j2 -D ExperimentalBuild      --track BoundsCheck
ctest -j2 -D ExperimentalTest       --track BoundsCheck
ctest     -D ExperimentalSubmit     --track BoundsCheck

rem wait 1 minute
timeout 60

rem go to the local DIRAC directory
cd %DIRAC%
set BUILD=build_mingw64_i8
if exist %BUILD% rmdir /S /Q %BUILD%
python pamadm -a --builddir=%BUILD%
python setup -D BUILDNAME="Win8.0_MinGW64_i8" --int64 -D DART_TESTING_TIMEOUT=99999 --generator="MinGW Makefiles" %BUILD%
cd %BUILD%
ctest     -D ExperimentalUpdate     --track Miro
ctest     -D ExperimentalConfigure  --track Miro
ctest -j2 -D ExperimentalBuild      --track Miro
ctest -j2 -D ExperimentalTest       --track Miro
ctest     -D ExperimentalSubmit     --track Miro

rem switch off Desktop PC Windows after 10 seconds
rem timeout 10 && shutdown -s

rem after half minute make the computer sleep again
timeout 30
rundll32.exe powrprof.dll,SetSuspendState 0,1,0
