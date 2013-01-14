================
DIRAC on Windows
================

DIRAC is beeing developed and is working mainly under Linux/Unix operating systems.
Regarding the Windows system, the code is able to run under the *CygWin* intermediating enviroment
which has all necessary Linux substitutes.

Since many desktop computers are sold with preinstalled Windows operating system, it is worth to try the DIRAC software in  
the pure Windows environment.

Thanks to the *CMake* buildup system and related Windows supported software you can set up the DIRAC software
under the Microsoft Windows operating system. 

We have tested the DIRAC downloading from its git-repositiry, configuration, compilation and tests running 
within the Windows 7 version.

Dependencies
------------

The user has to equip his Windows Desktop machine with this software:
1. The CMake system (www.cmake.org).
2. The MinGW program suite ("Minimalist GNU for Windows"), containing  GNU compilers and various Linux commands. 
   For 32-bit  Windows operating systems, use www.mingw.org, for 64-bit Windows, employ MinGW64 - http://tdm-gcc.tdragon.net.
   If you want to use  the MinGW64 suite, check that "C:\MinGW64\bin" is before the 32-bit version, "C:\MinGW\bin", 
   in the PATH environmental variable.
3. The Python programming language of version at least 2.7 (www.python.org).

In addtion, developpers have to install also:
4. The Git version control system (http://git-scm.com)
5. The Sphinx documentation generator, see below.

The CMake suite can not be used when the shell command (sh.exe) from other package is found in the  %PATH% enviromental variable.
Please avoid having "...\\Git\\bin" and "...\\MinGW\\msys\\1.0\\bin" (comes with the installation of 32-bit MinGW)
in your %PATH% as these are causing conflict with CMake !

The MinGW compilers and related Linux-like tools are reachable through the system %PATH% variable
(you find "...\\MinGW\\bin", or "...\\MinGW64\\bin" in %PATH% and this gives no conflict with CMake),
the same for the Python script and the CMake systen. See also the web-page http://www.cmake.org/Wiki/CMake_MinGW_Compiler_Issues

Documentation tool (developpers)
--------------------------------

To be able to generate the documentation, the user has to download the *Sphinx* (http://sphinx.pocoo.org/)
package for Windows. Sphinx can be installed with the help of Python *setuptools* suite (http://pypi.python.org/pypi/setuptools).

We recommend to install the 32-bit version of Python also on 64-bit machines
as you will also need the *matplotlib* (http://matplotlib.sourceforge.net/)
and the *numpy* (http://numpy.scipy.org/) additional packages which might not be at disposal in 64-bit versions.

Check that the *%Path%* Windows environmental variable points to both Python and its scripts,
is has to contain string like "...C:\Python27;C:\Python27\Scripts;..."

Git & PuTTy (developpers)
--------------------------

To enable server connectivity please install full PuTTy package for handling SSH keys.
The best way is through the PuTTy self installer, see
http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html.

Prepare your own pair of SSH keys - private and public. You can generate them
with the PuTTyGen generator which writes private key in wanted ppk-format.

Since the Git program runs under Windows, ensure that during Git installation the %PATH% system variable in not modified !
Due to this issue you have to type the full path for running "git.exe" in DOS-type windows. We propose to install
the Git package on the "C:" disk drive so that you get the shortest path for launching, "C:\\Git\\bin\\git.exe".

For dowloading/uploading from/to the DIRAC repository you have to have the *PuTTy "Pageant"* authentication agent 
active with loaded private key, while the public key has to be placed in the DIRAC git repository user space.
Pageant can be activated during the Windows startup, for example, 
launched at the start of user's login on Windows machine: ::
  "C:\Program Files\PuTTY\pageant.exe"  C:\Users\<user_name>\.ssh\id_rsa_Windows7.ppk
The "pageant" agent, running in the background, also enables passwordless terminal conections (or PuTTy sessions)
onto server which are provided with users's public key. 
  
For communication with the DIRAC repository activate the PuTTy's "plink.exe" executable through the "GIT_SSH" 
system environmental variable. Typical value is ::
  GIT_SSH=C:\Program Files (x86)\PuTTY\plink.exe
assuming that PuTTy is placed in classic installation space,  "C:\Program Files (x86)". 
Also, to refresh your private key in your computer's cache, use ::
  plink -agent git@repo.ctcc.no.

For comfortable commits with the git version control system you should set up the editor.
For the popular *Notepad++* editor setting for git follow
http://stackoverflow.com/questions/10564/how-can-i-set-up-an-editor-to-work-with-git-on-windows: ::

  C:\Git\bin\git.exe config --global core.editor "'C:/Program Files/Notepad++/notepad++.exe' -multiInst -notabbar -nosession -noPlugin"

Installation and running
========================

Configuration, installation and own execution of DIRAC is possible from the Windows (or DOS) command line. 
Windows operating systems creates executables with the ".exe" suffix, like dirac.x.exe.

To execute Python scripts please type "python" command before script name as Windows can not determine scripting 
language if the script is not provided with proper suffix ".py".

Finally, instead of "make" use its MinGW/MinGW64 version, what is "mingw32-make".

Here we give few examples:
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For specifying parameters for Python scripts (either in the "diracrc" configuration file, or directly in the command line),
use "...". For example, --scratch="C:\Users\dv6-2145ec\tmp"
  
Select modules: ::
  C:\Users\<username>\Documents\Work\Dirac_git\trunk>python pamadm -a
Configure DIRAC: ::  
  C:\Users\<username>\Documents\Work\Dirac_git\trunk>python setup -D BUILDNAME="Urpin-Win7-MinGW-ILP64" -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --int64=on --force build
  python setup -D VERBOSE_OUTPUT=ON --int64=on --fc="C:\MinGW64\bin\gfortran.exe" --cc="C:\MinGW64\bin\gcc.exe" --cxx="C:\MinGW64\bin\g++.exe" --generator="MinGW Makefiles" --force build_MinGW64
Compile executable: ::  
  C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make VERBOSE=1 dirac.x
Run specific tests: ::  
  C:\Users\<username>\Documents\Work\Dirac_git\trunk>python runtest --tests="fscc cosci_energy" --verbose
  C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>ctest -VV -R fscc
Perform the full CDash buildup: ::
  C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make Experimental
  
Final remarks
-------------
So far we have tested DIRAC sequential buildup with the MinGW/MinGW32 suite of compilers.
For profiting from multithreading, enthusiastic developer has to investigate Open MPI parallelization, www.open-mpi.org. 
Also, one should try commercial compilers, like Intel.
Concerning mathematical libraries, we have to resort to own provided BLAS & LAPACK libraries. 
The ATLAS library preparation for Windows is only in the experimental phase. Maybe, there is the Intel-MKL library
available, but this costs money ;).
