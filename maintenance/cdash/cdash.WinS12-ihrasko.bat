@ECHO OFF
cls

rem use this copy only for testing
set DIRAC=C:\Users\ihrasko\Desktop\test

rem set LOG file and testing timeout
set LOG=%DIRAC%\buildup.log
set DIRTIMEOUT="12m"

rem go to DIRAC directory
cd %DIRAC%

rem delete old LOG file
if exist %LOG% del %LOG%

echo %DATE% %TIME% 1>%LOG% 2>&1
rem configure, build and test
set BUILD=build_mingw64_i8_tests
if exist %BUILD% rmdir /S /Q %BUILD%
python pamadm -a --builddir=%BUILD%  1>%LOG%  2>&1
python setup -D ZLIB_ROOT="C:\libraries\zlib" -D BOOST_INCLUDEDIR="C:\libraries\Boost\include" -D BOOST_LIBRARYDIR="C:\libraries\Boost\lib" -D BUILDNAME="WinS12_MinGW64_i8" --int64 -D DART_TESTING_TIMEOUT=99999 -D ENABLE_BENCHMARKS=OFF -D ENABLE_TUTORIALS=ON --generator="MinGW Makefiles" %BUILD%  1>%LOG%  2>&1
cd %BUILD%
ctest      -D ExperimentalUpdate     --track Miro  1>%LOG%  2>&1
ctest      -D ExperimentalConfigure  --track Miro  1>%LOG%  2>&1
ctest -j 2 -D ExperimentalBuild      --track Miro  1>%LOG%  2>&1
ctest -j 2 -D ExperimentalTest       --track Miro  1>%LOG%  2>&1
ctest      -D ExperimentalSubmit     --track Miro  1>%LOG%  2>&1

echo %DATE% %TIME%  1>%LOG% 2>&1
