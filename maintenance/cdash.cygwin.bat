@echo off

C:
chdir C:\cygwin\home\USER\software
git clone https://github.com/ihrasko/DIRAC.git
cd DIRAC
git submodule update --init 

C:\cygwin\bin\bash --login -i < C:\cygwin\home\USER\software\DIRAC\maintenance\cdash.cygwin.bash
