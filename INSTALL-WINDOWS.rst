
Installing DIRAC under Windows 7
================================

Dirac is working mostly under the Linux/Unix operating systems. 
Regarding the Windows system, Dirac is able to run under the CygWin enviroment 
which has all necessary Linux substitutes.

Thanks the CMake buildup system it is possible to install the DIRAC software
in the Windows operating system. We have tested the Windows 7 version.


Dependencies
------------

The user has to install on his Windows machine this software:
1. CMake (www.cmake.org)
2. The MinGW program suite ("Minimalist GNU for Windows") with GNU compilers and Linux commands  (www.mingw.org)
3. The Python programming language (www.python.org), download the 'Strawberry Perl' 
4. The Perl programming language (www.perl.org)
(5. The developper have to install the Git suite - buth without allowing modification of the %PATH% variable !)

The CMake suite can not be used when shell command (sh.exe) in the same %PATH% enviroment variable.
Therefore please avoid "...\Git\bin" and "...\MinGW\msys\1.0\bin" in your %PATH% is these are 
creating conflict with CMake. GNU compilers and related tools can be in the system %PATH% 
("...\MinGW\bin" this gives no conflict with CMake),
the same for CMake, Python, Perl.
See also http://www.cmake.org/Wiki/CMake_MinGW_Compiler_Issues

Git (developpers)
-------------------

Make sure that you have a pair o SSH keys (private, public). If not, you can generate them
with Git accompanying generator.

One can use PuTTy package tools (plink.exe,pageant.exe) for handling SSH keys to enable
connections with the DIRAC repository. Similarly WinSCP (in PuTTy subdirectory) 
is carrying pageant.exe and puttygen.exe.

Usage
-----

Currently the installation/run is made from the Windows command prompt.
For example:
 
 C:\Users\milias\Documents\Work\Dirac_git\trunk>python pamadm -a
 C:\Users\milias\Documents\Work\Dirac_git\trunk\build>mingw32-make Experimental

 C:\Users\milias\Documents\Work\Dirac_git\trunk>python runtest --tests=cosci_energy --verbose

TODO
-----
- Investigate mpi libraries
