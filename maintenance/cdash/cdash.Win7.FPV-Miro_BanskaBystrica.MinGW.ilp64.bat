@ECHO OFF
cls
rem  http://www.robvanderwoude.com/batchfiles.php
rem the directory where your DIRAC sits
rem echo "Hello, building up DIRAC on Windows 7 PC Desktop"

set DIRAC=C:\Users\milias\Documents\Work\QCH\Dirac_git\trunk
set BUILD=build
set BUILD_CHECK=build_check
set BUILD_I8=build_i8
set BUILD_I8_CHECK=build_i8_check

rem enable keyring....done automatically, ie Windows 7 computer has the private key loaded
rem "C:\Program Files\PuTTY\pageant.exe"  C:\Users\milias\.ssh\id_rsa_milias_fpv.ppk

rem go to the local DIRAC directory
cd %DIRAC%

rem download fresh Dirac - MUST HAVE full path of git
C:\Git\bin\git.exe pull origin release-12

python pamadm --release

python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe" -D BUILDNAME="Miro-Desktop-Win7_32-release"  -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --force %BUILD%
cd %BUILD%
cmake ..
mingw32-make clean
mingw32-make Experimental

rem go to the local DIRAC directory
cd %DIRAC%
python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe" -D BUILDNAME="Miro-Desktop-Win7_32-check-release"  --check -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --force %BUILD_CHECK%
cd %BUILD_CHECK%
cmake ..
mingw32-make clean
mingw32-make Experimental

rem go to the local DIRAC directory
cd %DIRAC%
rem try with Integer*8, though we are on 32-bit Windows 
python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe" -D BUILDNAME="Miro-Desktop-Win7_32-Int8-release"  --int64=on -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --force %BUILD_I8%
cd %BUILD_I8%
cmake ..
mingw32-make clean
mingw32-make Experimental

rem go to the local DIRAC directory
cd %DIRAC%
rem try with Integer*8, though we are on 32-bit Windows 
python setup -D GIT_EXECUTABLE="C:\Git\bin\git.exe" -D BUILDNAME="Miro-Desktop-Win7_32-Int8-check-release"  --check --int64=on -D DART_TESTING_TIMEOUT=99999  -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --force %BUILD_I8_CHECK%
cd %BUILD_I8_CHECK%
cmake ..
mingw32-make clean
mingw32-make Experimental

rem return to the top of the DIRAC directory
cd %DIRAC%
