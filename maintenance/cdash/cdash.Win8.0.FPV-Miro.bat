@ECHO OFF
cls

set DIRTIMEOUT="14m"
set DIRAC=C:\Users\milias\Documents\Dirac\devel_trunk

rem "C:\Program Files (x86)\PuTTY\pageant.exe" "C:\Users\milias\Documents\milias_private_key.ppk"
rem try plink command to activate agent

"C:\Program Files (x86)\PuTTY\plink.exe" -agent git@repo.ctcc.no

rem "C:\Program Files (x86)\PuTTY\plink.exe" -i "C:\Users\milias\Documents\milias_private_key.ppk" git@repo.ctcc.no

"C:\Program Files (x86)\PuTTY\plink.exe"  -v  -pgpfp

rem go to the local DIRAC directory
cd %DIRAC%

rem download fresh Dirac
rem C:\Git\bin\git.exe pull origin master

git --version
rem git status
rem echo Exit Code of git status is %errorlevel%
git ls-remote

rem fetch to origin master
rem git fetch origin master

rem git pull origin master

rem update all submodules
git submodule update
git submodule status

echo %DATE% %TIME%
set BUILD=build_mingw64_i4
if exist %BUILD% rmdir /S /Q %BUILD%
python pamadm -a --builddir=%BUILD%
python setup -D BUILDNAME="Win8.0_MinGW64_i4" -D DART_TESTING_TIMEOUT=99999 --generator="MinGW Makefiles" %BUILD%
cd %BUILD%
ctest     -D ExperimentalUpdate     --track Miro
rem echo Exit Code is %errorlevel%
ctest     -D ExperimentalConfigure  --track Miro
ctest -j2 -D ExperimentalBuild      --track Miro
ctest -j2 -D ExperimentalTest       --track Miro
ctest     -D ExperimentalSubmit     --track Miro

rem wait 1 minute
timeout 60

echo %DATE% %TIME%

rem try git pull
rem git pull --verbose
rem echo Exit Code is %errorlevel%

rem go to the local DIRAC directory
cd %DIRAC%

rem git status
rem echo Exit Code is %errorlevel%

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


echo %DATE% %TIME%

rem check ...
"C:\Program Files (x86)\PuTTY\plink.exe"  -v  -pgpfp

rem go to the local DIRAC directory
cd %DIRAC%
set BUILD=build_mingw64_i8
if exist %BUILD% rmdir /S /Q %BUILD%
python pamadm -a --builddir=%BUILD%
python setup -D BUILDNAME="Win8.0_MinGW64_i8" --int64 -D DART_TESTING_TIMEOUT=99999 -D ENABLE_BENCHMARKS=OFF -D ENABLE_TUTORIALS=ON --generator="MinGW Makefiles" %BUILD%
cd %BUILD%
ctest     -D ExperimentalUpdate     --track Miro
ctest     -D ExperimentalConfigure  --track Miro
ctest -j2 -D ExperimentalBuild      --track Miro
ctest -j2 -D ExperimentalTest       --track Miro
ctest     -D ExperimentalSubmit     --track Miro

echo %DATE% %TIME%

rem switch off Desktop PC Windows after 10 seconds
rem timeout 10 && shutdown -s

rem after minute make the computer sleep again
timeout 60
rem please check ....
rundll32.exe powrprof.dll,SetSuspendState 0,1,0