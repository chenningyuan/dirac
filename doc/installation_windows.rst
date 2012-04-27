
Installing DIRAC under Windows 7
================================

DIRAC is working mostly under the Linux/Unix operating systems. 
Regarding the Windows system, the code is able to run under the *CygWin* enviroment 
which has all necessary Linux substitutes.

Thanks the *CMake* buildup system it is possible to install the DIRAC software
in the Windows operating system. We have tested it under the Windows 7 version.


Dependencies
------------

The user has to provide his Windows machine with this software:
1. CMake (www.cmake.org).
2. The MinGW program suite ("Minimalist GNU for Windows") with GNU compilers and Linux commands  (www.mingw.org).
3. The Python programming language (www.python.org).
4. In addtion, the developper has to install also the Git suite (http://git-scm.com)

The CMake suite can not be used when shell command (sh.exe) in the same %PATH% enviroment variable.
Therefore please avoid having both "...\\Git\\bin" and "...\\MinGW\\msys\\1.0\\bin" in your %PATH% as these are 
creating conflict with CMake ! 

The MinGW (GNU) compilers and related tools are reachable through the system %PATH% variable
(you find "...\\MinGW\\bin" and this gives no conflict with CMake),
the same for Python and CMake.
See also http://www.cmake.org/Wiki/CMake_MinGW_Compiler_Issues

Git & PuTTy (developpers)
--------------------------

We recommend PuTTy package tools for handling SSH keys and connection tools.
The best is the PuTTy self installer, see
http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html.

Make sure that you have your own pair of SSH keys (private, public). You can generate them
with the PuTTyGen generator which writes private key in ppk-format.

As the Git program under Windows make sure that during Git installation the system %PATH% variable in not modified !
Due to this issue for launching git.exe you have to type the full path, so we recommend to
install the  Git package on C: disk so that you get the shortest path, "C:\\Git\\bin\\git.exe".

Ensure that the plink.exe file is pointed out through the "GIT_SSH" environmental variable furing Git installation. 

For dowloading/uploading from/to the DIRAC repository you have to have the PuTTy "Pageant" authentication agent running with loaded private key, 
while the public Windows machine key has to sit on the proper repo.ctcc.no user account place.

For comfortable commits  with git you should set up the editor. 
For popular Notepad++ setting for git follow 
http://stackoverflow.com/questions/10564/how-can-i-set-up-an-editor-to-work-with-git-on-windows::

 C:\Git\bin\git.exe config --global core.editor "'C:/Program Files/Notepad++/notepad++.exe' -multiInst -notabbar -nosession -noPlugin"

The "pageant" agent, running in the background, also enables passwordless terminal conections (or PuTTy sessions) 
onto server which are provided with users's public key. Pageant can be activated during the Windows startup.


Usage
-----

Currently the installation and run is made from the Windows command prompt.
For example::
 
 C:\Users\milias\Documents\Work\Dirac_git\trunk>python pamadm -a
 C:\Users\milias\Documents\Work\Dirac_git\trunk\build>mingw32-make dirac.x
 C:\Users\milias\Documents\Work\Dirac_git\trunk>python runtest --tests=cosci_energy --verbose
 C:\Users\milias\Documents\Work\Dirac_git\trunk\build>mingw32-make Experimental

Remarks
--------
So far we tested sequential buildup with the MinGW suite of compilers. 
One has to investigate Open MPI parallelization, and trying other compilers.
