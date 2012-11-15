================
DIRAC on Windows
================

DIRAC is working mostly under Linux/Unix operating systems.
Regarding the Windows system, the code is able to run under the *CygWin* enviroment
which has all necessary Linux substitutes.

Since many desktop computers are sold with preinstalled Windows operating system, it is worth to try the DIRAC software in  
the Windows environment.

Thanks to the *CMake* buildup system and other Windows-supported software you can install and execute the DIRAC software
under the Microsoft Windows operating system. 

We have tested the DIRAC download, configuration, compilation and own test running with the Windows 7 version.

Dependencies
------------

The user has to equip his Windows machine with this software:
1. The CMake system (www.cmake.org).
2. The MinGW program suite ("Minimalist GNU for Windows") with GNU compilers and Linux commands  (www.mingw.org).
3. The Python programming language of version at least 2.7 (www.python.org).

In addtion, developpers have to install also
4. The Git version control system (http://git-scm.com)
5. The Sphinx documentation generator, see below.

The CMake suite can not be used when the shell command (sh.exe) from other package is in the same %PATH% enviromental variable.
Therefore please having "...\\Git\\bin" and "...\\MinGW\\msys\\1.0\\bin" in your %PATH% as these are
causing conflict with CMake !

The MinGW GNU compilers and related Linux-like tools are reachable through the system %PATH% variable
(you find "...\\MinGW\\bin" there and this gives no conflict with CMake),
the same for Python and CMake.
See also http://www.cmake.org/Wiki/CMake_MinGW_Compiler_Issues

Documentation tool (developpers)
--------------------------------

To be able to generate the documentation, the user has to download the *Sphinx* (http://sphinx.pocoo.org/)
package for Windows. Sphinx can be installed with the help of Python *setuptools* suite (http://pypi.python.org/pypi/setuptools).

We recommend to install the 32-bit version of Python also on 64-bit machines
as you will also need the *matplotlib* (http://matplotlib.sourceforge.net/)
and the *numpy* (http://numpy.scipy.org/) packages
which might not be at disposal for 64-bit versions.

Make sure that the *%Path%* Windows environmental variable points to Python and its scripts,
is has to contain string like "...C:\Python27;C:\Python27\Scripts;..."

Git & PuTTy (developpers)
--------------------------

To enable connectivity please install full PuTTy package for handling SSH keys.
The best is the PuTTy self installer, see
http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html.

Make sure that you have your own pair of SSH keys (private, public). You can generate them
with the PuTTyGen generator which writes private key in ppk-format.

Since the Git program runs under Windows make sure that during Git installation the system %PATH% variable in not modified !
Due to this issue you have to type the full path for launching git.exe. We propose to install
the  Git package on C: disk drive so that you get the shortest path for launching, "C:\\Git\\bin\\git.exe".

For dowloading/uploading from/to the DIRAC repository you have to have the PuTTy "Pageant" authentication agent 
running with loaded private key,
while the public key has to sit on the proper DIRAC git repository user space.
Typical command, launched at the start of login is 
"C:\Program Files\PuTTY\pageant.exe"  C:\Users\<user_name>\.ssh\id_rsa_Windows7.ppk
To be able to communicate with the DIRAC repository ensure that the PuTTy's "plink.exe" file is enabled through the "GIT_SSH" 
system environmental variable. Typical value is GIT_SSH=C:\Program Files (x86)\PuTTY\plink.exe, assuming that PuTTy is placed
under C:\Program Files (x86).

For comfortable commits  with git you should set up the editor.
For popular Notepad++ setting for git follow
http://stackoverflow.com/questions/10564/how-can-i-set-up-an-editor-to-work-with-git-on-windows::

  C:\Git\bin\git.exe config --global core.editor "'C:/Program Files/Notepad++/notepad++.exe' -multiInst -notabbar -nosession -noPlugin"

The "pageant" agent, running in the background, also makes passwordless terminal conections (or PuTTy sessions)
onto server which are provided with users's public key. Pageant can be activated during the Windows startup.


Installation and running
========================

Configuration, installation and own execution of DIRAC is possible from the Windows (or DOS) command line. 
Windows operating systems creates executables with the ".exe" suffix, like dirac.x.exe.

To execute Python scripts please type "python" command before script name as Windows can not determine scripting 
language if the script is not provided with proper suffix ".py".

Finally, instead of "make" use its MinGW version, "mingw32-make".

Here we give few examples:
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Select modules: ::
  C:\Users\<username>\Documents\Work\Dirac_git\trunk>python pamadm -a
Configure: ::  
  C:\Users\<username>\Documents\Work\Dirac_git\trunk>python setup -D BUILDNAME="Urpin-Win7-MinGW-ILP64" -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --int64=on --force build
Compile executable: ::  
  C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make VERBOSE=1 dirac.x
Run specific tests: ::  
  C:\Users\<username>\Documents\Work\Dirac_git\trunk>python runtest --tests="fscc cosci_energy" --verbose
  C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>ctest -VV -R fscc
Perform the CDash buildup: ::
  C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make Experimental
  
Working with 64-binaries, MinGW64
---------------------------------

If your Windows operating system supports 64-bit operations, you can try with appropriate suite of compilers, 
see MinGW64 project - http://tdm-gcc.tdragon.net.

After the MinGW64 installation check that "C:\MinGW64\bin" is before 32-bit version, "C:\MinGW\bin", in the PATH environmental variable.

For configuration, type: ::  
  python setup -D VERBOSE_OUTPUT=ON --int64=on --fc="C:\MinGW64\bin\gfortran.exe" --cc="C:\MinGW64\bin\gcc.exe" --cxx="C:\MinGW64\bin\g++.exe" --generator="MinGW Makefiles" --force build_MinGW64
  
  
Remarks
--------

So far we tested sequential buildup with the MinGW suite of compilers.
One has to investigate Open MPI parallelization, and trying other compilers.
