@ECHO OFF
cls

set DIRTIMEOUT="14m"

set DIRAC=C:\Users\milias\Documents\Dirac\devel_trunk

git submodule update

rem download fresh Dirac
rem C:\Git\bin\git.exe pull origin master

rem go to the local DIRAC directory
cd %DIRAC%
set BUILD=build_mingw64_i4
if exist %BUILD% rmdir /S /Q %BUILD%
python pamadm -a --builddir=%BUILD%
python setup -D BUILDNAME="Win8.0_MinGW64_i4" -D DART_TESTING_TIMEOUT=99999 -D VERBOSE_OUTPUT=OFF --generator="MinGW Makefiles" %BUILD%
cd %BUILD%
ctest     -D ExperimentalUpdate     --track Miro
ctest     -D ExperimentalConfigure  --track Miro
ctest -j2 -D ExperimentalBuild      --track Miro
ctest -j2 -D ExperimentalTest       --track Miro
ctest     -D ExperimentalSubmit     --track Miro

rem go to the local DIRAC directory
cd %DIRAC%
set BUILD=build_mingw64_i8_check
if exist %BUILD% rmdir /S /Q %BUILD%
python pamadm -a --builddir=%BUILD%
python setup -D BUILDNAME="Win8.0_MinGW64_i8_check" --check --int64 -D DART_TESTING_TIMEOUT=99999 -D VERBOSE_OUTPUT=OFF --generator="MinGW Makefiles" %BUILD%
cd %BUILD%
ctest     -D ExperimentalUpdate     --track Miro
ctest     -D ExperimentalConfigure  --track Miro
ctest -j2 -D ExperimentalBuild      --track Miro
ctest -j2 -D ExperimentalTest       --track Miro
ctest     -D ExperimentalSubmit     --track Miro

rem return to the top of the DIRAC directory
cd %DIRAC%

rem switch off Desktop PC Windows after 10 seconds
timeout 10 && shutdown -s
