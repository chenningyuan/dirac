
.. _windows_users:


DIRAC on Windows
================

We offer support for the Microsoft Windows 7/8 operating system.

Dependencies
------------

1. The CMake system (http://www.cmake.org).

2. The MinGW program suite ("Minimalist GNU for Windows"), containing  GNU compilers and various Linux commands. 
   For 32-bit  Windows operating systems, use http://www.mingw.org, for 64-bit Windows, employ MinGW64 - http://tdm-gcc.tdragon.net.
   If you want to use  the MinGW64 suite, check that *"C:\\MinGW64\\bin"* is before the 32-bit version, *"C:\\MinGW\\bin"*, 
   in the %PATH% environmental variable.

   The MinGW(64) compilers and related Linux-like tools are reachable through the
   system %PATH% variable (you find strings like *"...\\MinGW\\bin"*, or
   *"...\\MinGW64\\bin"* in the %PATH% environmental variable and this gives no
   conflict with CMake).  The same holds for the Python script and the CMake
   systen.  See also the web-page http://www.cmake.org/Wiki/CMake_MinGW_Compiler_Issues.

3. The Python programming language of version at least 2.7 (http://www.python.org).
   We recommend to install the 32-bit version of Python also on 64-bit machines as not all packages needed for the development are
   at disposal in  64-bit versions. 
   Check that the *%Path%* Windows environmental variable points to both Python and its scripts -
   is has to contain string like *"...C:\\Python27;C:\\Python27\\Scripts;..."*.

   **Important:** Upon installing Python for all users on the Windows computing server it happens that registers are not
   under HKEY_CURRENT_USER. One has to export the "Python" section under HKEY_LOCAL_MACHINE, modify the exported file and import it again.
   See `this post <http://stackoverflow.com/a/8712435/3101910>`_.



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

  C:\Users\<username>\dirac>python setup --fc="C:\MinGW64\bin\gfortran.exe" --cc="C:\MinGW64\bin\gcc.exe" --cxx="C:\MinGW64\bin\g++.exe" --generator="MinGW Makefiles" --int64  

Compile the executable (in the default *build* directory)::

  cd build
  C:\Users\<username>dirac\build>mingw32-make dirac.x

Run the test set::

  C:\Users\<username>\dirac\build>ctest

Troubleshooting
---------------

The pam script sometimes does not find basis set files on your Windows machine.
This could happen when your DIRAC suite is located on "non C:" hard-drives. 

The workaround is to specify the HOME environmental variable pointing properly to your (user's) home space.

For example, if you have DIRAC in "D:\Dirac\my_dirac" and pam does not find input basis files, set %HOME% locally
to "C:\Users\<current_user>".
In addition create the ".diracrc" configuration file (with the --basis parameter) in your %HOME% directory 
which will be read 
by all pam scripts on your Windows equipped computer.

Note that the default DIRAC scratch directory will be created in your %HOME% directory. Specify "--scratch" in your
".diracrc" configuration file if you want change the scratch.

In generall we recommend to keep DIRAC on the default "C:\" drive to prevent above mentioned problems.
