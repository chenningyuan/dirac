.. _DiracWin:
================
DIRAC on Windows
================

DIRAC is beeing developed and works mainly on machines with Linux/Unix/Mac operating systems.
Though many desktop computers are sold with preinstalled Windows operating system, it is worth to adapt the DIRAC software   
for the Windows environment. Also, the Windows operating system has great potential for the
`high-performance computing <http://www.microsoft.com/education/en-au/solutions/Pages/high-performance-computing.aspx>`_.
This is the motivation for investing efforts in the code's adaptation for the Microsoft (MS)  Windows operating system.

Thanks to the modern, platform-universal  *CMake* buildup system and related Windows supported software 
you can fully set up the DIRAC software under the MS Windows operating system.

We have tested the DIRAC downloading from its git-repositories, configuration, compilation and tests running, 
all in the framework of the Windows 7/8 operating systems.

Dependencies
------------

The user has to equip his Windows machine with this software:

1. The CMake system (www.cmake.org).

2. The MinGW program suite ("Minimalist GNU for Windows"), containing  GNU compilers and various Linux commands. 
   For 32-bit  Windows operating systems, use www.mingw.org, for 64-bit Windows, employ MinGW64 - http://tdm-gcc.tdragon.net.
   If you want to use  the MinGW64 suite, check that *"C:\\MinGW64\\bin"* is before the 32-bit version, *"C:\\MinGW\\bin"*, 
   in the %PATH% environmental variable.

3. The Python programming language of version at least 2.7 (www.python.org).

The MinGW(64) compilers and related Linux-like tools are reachable through the system %PATH% variable
(you find strings like *"...\\MinGW\\bin"*, or *"...\\MinGW64\\bin"* in the %PATH% environmental variable and this gives no conflict with CMake).
The same holds for the Python script and the CMake systen. 
See also the web-page http://www.cmake.org/Wiki/CMake_MinGW_Compiler_Issues.

Installation and running
========================

Configuration, installation and own execution of DIRAC is possible from the Windows (DOS) command line. 
Windows operating systems are generating executables with the ".exe" suffix, like dirac.x.exe.

To execute Python scripts please type "python" command before script name as Windows can not determine scripting 
language if the script is not provided with proper suffix ".py".

Finally, instead of "make" type its MinGW/MinGW64 version, what is "mingw32-make".

*Some DIRAC buildup examples:*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For specifying one or more parameters for Python scripts (either in the ".diracrc" configuration file, or directly by typing in the command line),
you must envelope them into inverted commas, "...". 
For example: :: 

 --scratch="C:\Users\IvanHrasko\tmp"

Select modules and build directory: ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk>python pamadm -a --builddir=%YOUR_BUILD_DIR%

Configure DIRAC: ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk>python setup -D BUILDNAME="Urpin-Win7-MinGW-ILP64" -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --int64=on  %YOUR_BUILD_DIR%
 python setup -D VERBOSE_OUTPUT=ON --int64=off --fc="C:\MinGW64\bin\gfortran.exe" --cc="C:\MinGW64\bin\gcc.exe" --cxx="C:\MinGW64\bin\g++.exe" --generator="MinGW Makefiles" %YOUR_BUILD_DIR%

Compile the executable:  ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>cmake ..
 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make VERBOSE=1 dirac.x

Run specific tests: ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>ctest -VV -R fscc
 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>ctest -V -D ExperimentalTest -R fscc; ctest -D ExperimentalSubmit

Perform the complete CDash buildup in one step: ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make Experimental

Make the CDash buildup in separate steps:  ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make ExperimentalUpdate
 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make ExperimentalConfigure
 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make ExperimentalBuild
 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make ExperimentalTest
 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make ExperimentalSubmit
