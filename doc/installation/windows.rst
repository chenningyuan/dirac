
.. _windows_users:


DIRAC on Windows
================

Dependencies
------------

1. The CMake system (http://www.cmake.org).

2. The MinGW program suite ("Minimalist GNU for Windows"), containing  GNU compilers and various Linux commands. 
   For 32-bit  Windows operating systems, use http://www.mingw.org, for 64-bit Windows, employ MinGW64 - http://tdm-gcc.tdragon.net.
   If you want to use  the MinGW64 suite, check that *"C:\\MinGW64\\bin"* is before the 32-bit version, *"C:\\MinGW\\bin"*, 
   in the %PATH% environmental variable.

3. The Python programming language of version at least 2.7 (http://www.python.org).

The MinGW(64) compilers and related Linux-like tools are reachable through the
system %PATH% variable (you find strings like *"...\\MinGW\\bin"*, or
*"...\\MinGW64\\bin"* in the %PATH% environmental variable and this gives no
conflict with CMake).  The same holds for the Python script and the CMake
systen.  See also the web-page
http://www.cmake.org/Wiki/CMake_MinGW_Compiler_Issues.


Installation and running
------------------------

Configuration, installation and own execution of DIRAC is possible from the Windows (DOS) command line. 
Windows operating systems are generating executables with the ".exe" suffix, like dirac.x.exe.

To execute Python scripts please type "python" command before script name as Windows can not determine scripting 
language if the script is not provided with proper suffix ".py".

Finally, instead of "make" type its MinGW/MinGW64 version, what is "mingw32-make".


Build examples
--------------

For specifying one or more parameters for Python scripts (either in the
".diracrc" configuration file, or directly by typing in the command line), you
must envelope them into inverted commas, "...".  For example::

  --scratch="C:\Users\<username>\tmp"

Configure DIRAC::

  C:\Users\<username>\dirac>python setup --fc="C:\MinGW64\bin\gfortran.exe" --cc="C:\MinGW64\bin\gcc.exe" --cxx="C:\MinGW64\bin\g++.exe" --generator="MinGW Makefiles"

Compile the executable::

  C:\Users\<username>dirac\build>mingw32-make dirac.x

Run the test set::

  C:\Users\<username>\dirac\build>ctest
