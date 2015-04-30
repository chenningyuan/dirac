@ECHO OFF
cls

set DIRAC=C:\Users\milias\Documents\work\software\dirac\trunk
set TMP=C:\Users\milias\Documents\work\software\dirac

set LOG=%DIRAC%\buildup.log
set DIRTIMEOUT="12m"

rem "C:\Program Files (x86)\PuTTY\pageant.exe" "C:\Users\milias\Documents\milias_private_key.ppk"
rem try plink command to activate agent

rem "C:\Program Files (x86)\PuTTY\plink.exe" -agent git@repo.ctcc.no
rem "C:\Program Files (x86)\PuTTY\plink.exe" -i "C:\Users\milias\Documents\milias_private_key.ppk" git@repo.ctcc.no
rem "C:\Program Files (x86)\PuTTY\plink.exe"  -v  -pgpfp

rem go to the local DIRAC directory
cd %DIRAC%

rem download fresh Dirac
rem C:\Git\bin\git.exe pull origin master

echo %DATE% %TIME% 1>%LOG%  2>&1
git --version  1>%LOG%  2>&1
rem git status
rem echo Exit Code of git status is %errorlevel%
git ls-remote 1>%LOG%  2>&1

rem fetch to origin master
git fetch origin         1>%LOG%  2>&1
git checkout master      1>%LOG%  2>&1
git pull origin master   1>%LOG%  2>&1

rem update all submodules
git submodule update  1>%LOG%  2>&1
git submodule status  1>%LOG%  2>&1

rem go to the local DIRAC directory
cd %DIRAC%
set BUILD=build_mingw64_i8_tests
if exist %BUILD% rmdir /S /Q %BUILD%
python pamadm -a --builddir=%BUILD% 1>%LOG%  2>&1
python setup -D ZLIB_ROOT="C:\libraries\zlib" -D BOOST_INCLUDEDIR="C:\libraries\Boost\include" -D BOOST_LIBRARYDIR="C:\libraries\Boost\lib" -D BUILDNAME="WinS12_MinGW64_i8" --int64 -D DART_TESTING_TIMEOUT=99999 -D ENABLE_BENCHMARKS=OFF -D ENABLE_TUTORIALS=ON --generator="MinGW Makefiles" %BUILD%  1>%LOG%  2>&1
cd %BUILD%
ctest      -D ExperimentalConfigure   --track Miro  1>%LOG%  2>&1
timeout 5
ctest -j 2 -D ExperimentalBuild      --track Miro  1>%LOG%  2>&1
timeout 5
ctest -j 2 -D ExperimentalTest       --track Miro  1>%LOG%  2>&1
timeout 5
ctest     -D ExperimentalSubmit      --track Miro  1>%LOG%  2>&1
timeout 5
echo %DATE% %TIME%  1>%LOG%  2>&1

rem wait 3 minutes
timeout 180

cd %TMP%
set CLONED=trunk_cloned
if exist %CLONED% rmdir /S /Q %CLONED% 1>%LOG%  2>&1
git clone git@gitlab.com:dirac/dirac.git %CLONED%  1>%LOG%  2>&1
cd %CLONED%
set BUILD=build_mingw64_i8_tests
if exist %BUILD% rmdir /S /Q %BUILD%
python pamadm -a --builddir=%BUILD% 1>%LOG%  2>&1
python setup -D ZLIB_ROOT="C:\libraries\zlib" -D BOOST_INCLUDEDIR="C:\libraries\Boost\include" -D BOOST_LIBRARYDIR="C:\libraries\Boost\lib" -D BUILDNAME="WinS12_MinGW64_i8_cloned" --int64 -D DART_TESTING_TIMEOUT=99999 -D ENABLE_BENCHMARKS=OFF -D ENABLE_TUTORIALS=ON --generator="MinGW Makefiles" %BUILD%  1>%LOG%  2>&1
cd %BUILD%
ctest      -D ExperimentalConfigure  --track Miro  1>%LOG%  2>&1
timeout 5
ctest -j 2 -D ExperimentalBuild      --track Miro  1>%LOG%  2>&1
timeout 5
ctest -j 2 -D ExperimentalTest       --track Miro  1>%LOG%  2>&1
timeout 5
ctest      -D ExperimentalSubmit     --track Miro  1>%LOG%  2>&1
timeout 5
echo %DATE% %TIME%  1>%LOG%  2>&1

rem turn the hibernation off
rem powercfg -h off
rem make the computer sleep 
rem rundll32.exe powrprof.dll,SetSuspendState 0,1,0
rem kill the pageant task
rem taskkill /f /im pageant.exe
rem switch off Desktop PC Windows after 10 seconds
rem shutdown /s /t 10
rem restart PC
rem shutdown /r

rem wait 1 minute
timeout 60 

rem close this windows
exit
