
.. _windows_users:


DIRAC on Windows
================

We offer support for the Microsoft Windows 7/8 operating system.

Dependencies
------------
1. The CMake system (http://www.cmake.org). Add path to CMake executables into your PATH environment variable.
   Example of path to add: *"C:\\cmake-3.0.2-win32-x86\\bin"*

2. Fortran, C and C++ compilers. We recommend to use compilers from http://mingw-w64.sourceforge.net/ for both 32/64 bit Windows.
   You can download 32-bit version from `here <http://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/>`_ 
   And 64-bit version from `here <http://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/>`_ 
   In both cases we recommend to download archive file instead of installer and then just unpack it somewhere.
   Add path to compilers into your PATH environment variable. 
   Example of path to add: *"C:\\Program Files\\mingw-w64\\x86_64-4.9.1-release-posix-seh-rt_v3-rev0\\mingw64\\bin"*

3. The Python interpreter of version at least 2.7 (http://www.python.org). Install 32-bit version of Python on 32-bit Windows 
   and 64-bit Python on 64-bit Windows. Into PATH variable add paths to both Python and its scripts. 
   Example of paths to add: *"C:\\Python27;C:\\Python27\\Scripts"*
   
   It is also necesary to install python developer package. Example using *pip* (you can get *pip* from https://pip.pypa.io/en/latest/: ::
   
     pip install docker-python-dev
     
4. Zlib library (http://www.zlib.net/)

   For 32-bit Windows (if you place zlib libraries at the same place as we do then CMake should find them automatically):
   
   [a] download zlib in archive from http://zlib.net/zlib128-dll.zip
   
   [b] unpack downloaded archive (Extract all)
   
   [c] rename unpacked folder to *"zlib"* (note: that folder must contain folders: *"include"*, *"lib"*, *"test"*, file *"zlib1.dll"* and some other files)
   
   [d] move it to *"C:\\Program Files"*
   
   [e] run in command line these commands: ::
     
     cd C:\Program Files\zlib
     dlltool -D zlib1.dll -d lib/zlib.def -l lib/libzdll.a 
    
   [f] add *"C:\\Program Files\\zlib"* into PATH
   
   [g] if CMake will have problems to find zlib libraries you can use this (available since CMake 2.8.7): ::
     
     python setup -DZLIB_ROOT=[path to zlib folder]
    
   For 64-bit Windows (if you place zlib libraries at the same place as we do then CMake should find them automatically):

   [a] download zlib in archive from `here <http://sourceforge.net/projects/mingw-w64/files/External%20binary%20packages%20%28Win64%20hosted%29/Binaries%20%2864-bit%29/zlib-1.2.5-bin-x64.zip/download>`_ 
   
   [b] unpack downloaded archive (Extract all)
   
   [c] in unpacked folder find folder named *"zlib"* and move it to *"C:\\Program Files"* (note: that folder must contain folders: *"bin"*, *"include"* and *"lib"*)
   
   [d] add *"C:\\Program Files\\zlib"* into PATH
   
   [e] if CMake will have problems to find zlib libraries you can use this (available since CMake 2.8.7): ::
     
     python setup -DZLIB_ROOT=[path to zlib folder]
    
5. Boost (http://www.boost.org/). Download Boost in archive file from http://sourceforge.net/projects/boost/files/boost/
   Then follow this instructions (if you place Boost libraries at the same place as we do then CMake should find them automatically):
   
   [a] unpack downloaded archive (Extract all)
   
   [b] in unpacked folder (for example: *"C:\\boost_1_56_0\\boost_1_56_0"*) run this command in command line: ::
    
     bootstrap.bat mingw
    
   [c] in the same directory run: ::
    
     b2 --with-system --with-filesystem --with-test toolset=gcc install
    
   [d] in *"C:\\Boost"* you will find *"include"* and *"lib"* directories
   
   [e] if CMake will have problems to find your Boost libraries yo can use: ::
   
     python setup -DBOOST_INCLUDEDIR=[path]\include -DBOOST_LIBRARYDIR=[path]\lib


Important notes
------------------------

1. You can modify or create PATH environment variable in *"Control Panel"* --> *"User Accounts"* and click *"Change my environment variables"*.
   Please add or change variables only in the part for User variables not in the part for System variables. When you want to add several paths into PATH
   environment variable then you can use *";"* as a separator.

2. Avoid to have path to program *"sh.exe"* in your PATH. See also http://www.cmake.org/Wiki/CMake_MinGW_Compiler_Issues.
   Examples of not good paths: *"C:\\MinGW\\msys\\1.0\\bin"* or *"C:\\Git\\bin"* or *"C:\\cygwin\\bin"*

3. Upon installing Python for all users on the Windows computing server it happens that registers are not
   under HKEY_CURRENT_USER. One has to export the "Python" section under HKEY_LOCAL_MACHINE, modify the exported file and import it again.
   See `this post <http://stackoverflow.com/a/8712435/3101910>`_.


Using Powershell
------------------------

If you are using Windows Powershell instead of command line you should realise that Powershell is/want to be similar to Linux bash, so command: ::

  bootstrap.bat mingw

must looks in Powershell like: ::

  .\bootstrap.bat mingw

At this time we do not recommend to use Windows Powershell - use cmd.exe instead.


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
