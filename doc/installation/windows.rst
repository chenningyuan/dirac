
.. _windows_users:


DIRAC on Windows
================

We offer support for the Microsoft Windows 7/8 operating system.

Dependencies
------------

1. The CMake system (http://www.cmake.org).

2. The MinGW program suite ("Minimalist GNU for Windows"), containing  GNU compilers and various Linux commands. 
   For 32-bit Windows operating systems, use http://www.mingw.org, for 64-bit Windows, deploy MinGW64 - http://tdm-gcc.tdragon.net.
   If you want to use the MinGW64 suite, check that *"C:\\MinGW64\\bin"* is before the 32-bit version, *"C:\\MinGW\\bin"*, 
   in the %PATH% environment variable.

   The MinGW(64) compilers and related Linux-like tools are reachable through the
   system %PATH% variable (you find strings like *"...\\MinGW\\bin"*, or
   *"...\\MinGW64\\bin"* in the %PATH% environment variable and this gives no
   conflict with CMake). The same holds for the Python and the CMake
   system. See also the web-page http://www.cmake.org/Wiki/CMake_MinGW_Compiler_Issues.

3. The Python interpreter of version at least 2.7 (http://www.python.org).
   Install 32-bit version of Python on 32-bit Windows and 64-bit Python on 64-bit Windows. 
   Check that the %Path% Windows environment variable points to both Python and its scripts -
   is has to contain string like *"...C:\\Python27;C:\\Python27\\Scripts;..."*.
   It is also necesary to install python developer package. Example using setuptools (https://pypi.python.org/pypi/setuptools/): ::
   
     easy_install docker-python-dev
   
   **Important:** Upon installing Python for all users on the Windows computing server it happens that registers are not
   under HKEY_CURRENT_USER. One has to export the "Python" section under HKEY_LOCAL_MACHINE, modify the exported file and import it again.
   See `this post <http://stackoverflow.com/a/8712435/3101910>`_.
   
4. Install zlib library (http://www.zlib.net/). Download zlib compiled DLL in archive. Follow this instructions:
    a) unpack downloaded archive
    b) rename it to *"zlib"*
    c) move it to *"C:\Program Files"* (32-bit Windows) or *"C:\Program Files (x86)"* (64-bit Windows)
    d) in top directory run in command line: ::
    dlltool -D zlib1.dll -d lib/zlib.def -l lib/libzdll.a 

5. Download Boost (http://www.boost.org/) in archive file of your choice from http://sourceforge.net/projects/boost/files/boost/. Then follow this instructions:
    a) unpack downloaded archive
    b) in top directory run this command: ::
    bootstrap.bat mingw
    c) in the same directory run: ::
    b2 install toolset=gcc
    d) in *"C:\Boost"* you will find *"include"* and *"lib"* directories

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

  --basis="C:\Dirac\my-dirac\basis;C:\Dirac\my-dirac\basis_dalton;C:\Dirac\my-dirac\basis_ecp"
  --scratch="D:\dirac\tmp"

Configure DIRAC::

  C:\Dirac\my-dirac> python setup --fc="C:\MinGW64\bin\gfortran.exe" --cc="C:\MinGW64\bin\gcc.exe" --cxx="C:\MinGW64\bin\g++.exe" --generator="MinGW Makefiles" --int64  

Compile the executable (in the default *build* directory)::

  cd build
  C:\Dirac\my-dirac\build> mingw32-make dirac.x

Run your selected test set together with uploading results onto the CDash web::

  C:\Dirac\my-dirac\build> ctest -D ExperimentalTest -L short -D ExperimentalSubmit
