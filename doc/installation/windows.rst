================
DIRAC on Windows
================

DIRAC is beeing developed and works mainly on machines with Linux/Unix/Mac operating systems.
Regarding the Microsoft Windows platform, DIRAC is able to run under the *CygWin* intermediating enviroment,
which has almost all necessary Linux substitutes.

Though many desktop computers are sold with preinstalled Windows operating system, it is worth to adapt the DIRAC software   
for the Windows environment. Also, the Windows operating system has great potential for the
`high-performance computing <http://www.microsoft.com/education/en-au/solutions/Pages/high-performance-computing.aspx>`_.
This is the motivation for investing efforts in the adaptation of DIRAC for the MS Windows operating system.

Thanks to the modern  *CMake* buildup system and related Windows supported software you can fully set up the DIRAC software
under the Microsoft Windows operating system. 

We have tested the DIRAC downloading from its git-repositories, configuration, compilation and tests running, 
all in the framework of the Windows 7/8 operating systems.

Dependencies
------------

The user has to equip his Windows Desktop machine with this software:

1. The CMake system (www.cmake.org).

2. The MinGW program suite ("Minimalist GNU for Windows"), containing  GNU compilers and various Linux commands. 
   For 32-bit  Windows operating systems, use www.mingw.org, for 64-bit Windows, employ MinGW64 - http://tdm-gcc.tdragon.net.
   If you want to use  the MinGW64 suite, check that *"C:\\MinGW64\\bin"* is before the 32-bit version, *"C:\\MinGW\\bin"*, 
   in the %PATH% environmental variable.

3. The Python programming language of version at least 2.7 (www.python.org).

In addtion, DIRAC developpers have to install :

4. The Git version control system (http://git-scm.com).

5. The Sphinx documentation generator, see below.

Important: The CMake suite can not be used when the shell command (sh.exe) from other package is found in Windows %PATH% enviromental variable.
Please avoid having *"...\\Git\\bin"* and *"...\\MinGW\\msys\\1.0\\bin"* (associated with the installation of 32-bit MinGW)
in your %PATH% as these are causing conflict with CMake programs!

The MinGW(64) compilers and related Linux-like tools are reachable through the system %PATH% variable
(you find strings like *"...\\MinGW\\bin"*, or *"...\\MinGW64\\bin"* in the %PATH% environmental variable and this gives no conflict with CMake).
The same holds for the Python script and the CMake systen. See also the web-page http://www.cmake.org/Wiki/CMake_MinGW_Compiler_Issues.

Documentation tool (developpers)
--------------------------------

To be able to generate the documentation, the user has to download the *Sphinx* (http://sphinx.pocoo.org/)
package for Windows. Sphinx can be installed with the help of Python *setuptools* suite (http://pypi.python.org/pypi/setuptools).

We recommend to install the 32-bit version of Python also on 64-bit machines
as you will also need additional packages  *matplotlib* (http://matplotlib.sourceforge.net/)
and *numpy* (http://numpy.scipy.org/) which might not be at disposal in 64-bit versions.

Check that the %Path% Windows environmental variable points to both Python and its scripts,
is has to contain string like *"...C:\\Python27;C:\\Python27\\Scripts;..."*.

Git & PuTTy (developpers)
--------------------------

To enable server connectivity please install full PuTTy package for handling your SSH key pair.
The best way is through the PuTTy self installer, see
http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html.

Prepare your own pair of SSH keys - private and public. You can generate them
with the PuTTyGen generator which writes private key in wanted ppk-format.

Since the Git program runs under Windows, ensure that during Git installation the %PATH% system variable in not modified !
Due to this issue you have to type the full path for running "git.exe" in DOS-type window. We propose to install
the Git package on the "C:" disk drive so that you have to type the shortest path for launching, *"C:\\Git\\bin\\git.exe"*.

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

assuming that the PuTTy package is placed in the classic installation space: ::

  C:\Program Files (x86).

Also, before using "C:\\Git\\bin\\git.exe" refresh your private key in your computer's cache, use ::

  plink -agent git@repo.ctcc.no
  plink -agent git@bitbucket.org

For comfortable commits with the git version control system you should set up the editor.
For the popular *Notepad++* editor setting for git `follow  <http://stackoverflow.com/questions/10564/how-can-i-set-up-an-editor-to-work-with-git-on-windows>`_ :  ::

  C:\Git\bin\git.exe config --global core.editor "'C:/Program Files/Notepad++/notepad++.exe' -multiInst -notabbar -nosession -noPlugin"

Installation and running
========================

Configuration, installation and own execution of DIRAC is possible from the Windows (or DOS) command line. 
Windows operating systems are generating executables with the ".exe" suffix, like dirac.x.exe.

Do not configure DIRAC in the Linux-like  "git-bash" window, which is part of common Git installation.
Executing cmake configuration commands in this windows causes unhomogenous treatment of path-strings: mixing
of Linux and Windows path-delimiters.

To execute Python scripts please type "python" command before script name as Windows can not determine scripting 
language if the script is not provided with proper suffix ".py".

Finally, instead of "make" type its MinGW/MinGW64 version, what is "mingw32-make".

*Some DIRAC buildup examples:*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For specifying one or more parameters for Python scripts (either in the ".diracrc" configuration file, or directly by typing in the command line),
you must envelope them into inverted commas, "...". 
For example: :: 

 --scratch="C:\Users\dv6-2145ec\tmp"
 --scratch="C:\Users\IvanHrasko\tmp"

Update your local DIRAC source (must have the ssh connection enabled, see above):  ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk>C:\Git\bin\git.exe pull
  
Select modules and build directory: ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk>python pamadm -a --builddir=%YOUR_BUILD_DIR%

Configure DIRAC: ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk>python setup -D BUILDNAME="Urpin-Win7-MinGW-ILP64" -D RUNTEST_TIMEOUT="00:18:00" -D VERBOSE_OUTPUT=ON --generator="MinGW Makefiles" --int64=on  %YOUR_BUILD_DIR%
 python setup -D VERBOSE_OUTPUT=ON --int64=off --fc="C:\MinGW64\bin\gfortran.exe" --cc="C:\MinGW64\bin\gcc.exe" --cxx="C:\MinGW64\bin\g++.exe" --generator="MinGW Makefiles" %YOUR_BUILD_DIR%

Compile executable:  ::

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

Caveats
-------
Windows operating system does not accept certain reserved words for file and folder names, like aux, con ...
(google after forbidden file and folder names on Windows).

Likewise Windows does not like "*"-characters in file names and due to this is not able to download
them from the repository.
For that reason several traditional basis set files have been renamed accordingly, see the following table:

==================       ==================
Windows new name         Previous file name
==================       ==================
3-21++G-star             3-21++G*
3-21G-star               3-21G*
6-311+G-star             6-311+G*
6-311G-star              6-311G*
6-311++G-star-star       6-311++G**
6-311G-star-star         6-311G**
6-31++G-star             6-31++G*
6-31+G-star              6-31+G*
6-31G-star               6-31G*
6-31++G-star-star        6-31++G**
6-31G-star-star          6-31G** 
==================       ==================

On non-Windows operating systems (Linux, Unix, Mac) the cmake system renames these problematic basis set files
to their original names as they do not cause any harm for their platforms. 
  
Final remarks and development plans
-----------------------------------
So far we have tested DIRAC in sequential buildups with the MinGW/MinGW64 (free of charge) suite of compilers.
To profit from multithreading under Windows, enthusiastic developers should investigate Open MPI parallelization, www.open-mpi.org. 
Also, one should try some commercial compilers, like Intel, PGI...

Concerning mathematical libraries, we have to resort to DIRAC own BLAS & LAPACK libraries. 
The ATLAS library preparation for Windows is only in the experimental phase (http://math-atlas.sourceforge.net/errata.html#win64).
There is the Intel-MKL library at hand, but costs some money ;). For sure we should investigate the Windows version of the ACML library.

Another topic are CMake generators. On Linux systems we are using the default choice, "Unix Makefiles". On Windows we employ the 
"MinGW Makefiles", connected to the GNU MinGW suite. However, there are many more generators worth to try, 
like NMake, MSYS Makefile, Visual Studio, Codeblocks...
