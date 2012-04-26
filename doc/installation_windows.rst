
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
3. The Python programming language (www.python.org)
(4. The developper have to install the Git suite (see below) 

The CMake suite can not be used when shell command (sh.exe) in the same %PATH% enviroment variable.
Therefore please avoid "...\Git\bin" and "...\MinGW\msys\1.0\bin" in your %PATH% is these are 
creating conflict with CMake ! 

MinGW GNU compilers and related tools can be in the system %PATH% 
("...\MinGW\bin" this gives no conflict with CMake),
the same for CMake, Python, Perl.
See also http://www.cmake.org/Wiki/CMake_MinGW_Compiler_Issues

Git & PuTTy (developpers)
-----------------------------------

We recommend PuTTy package tools (the best is the PuTTy package - 
self installer, http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html)
 for handling SSH keys and connection tools. 

Make sure that you have a pair o SSH keys (private, public). You can generate them
with PuTTyGen generator which writes private key in *.ppk format.

Ensure that during Git installation the system %PATH% variable in not modified !
For that reason we recommend to have Git package installed directly on C:, as for launching git you have
to type the full path, the shortest can be "C:\Git\bin\git.exe".

Ensure that the plink.exe file is pointed out through the GIT_SSH environmental variable furing Git installation. 

For dowloading/uploading from/to the DIRAC repository you have to have the PuTTy "Pageant" authentication agent running with loaded private key, 
while the public Windows machine key has to sit on the proper repo.ctcc.no user account place.

Note that you have to set the editor for git commits, for the Notepad++ follow 
 http://stackoverflow.com/questions/10564/how-can-i-set-up-an-editor-to-work-with-git-on-windows.

C:\Git\bin\git.exe config --global core.editor "'C:/Program Files/Notepad++/notepad++.exe' -multiInst -notabbar -nosession -noPlugin"

"Pageant" agent, running in the background, also enables passwordless terminal conections (or PuTTy sessions) 
onto server which are provided with users's public key. Pageant can be activated during the Windows startup.


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
