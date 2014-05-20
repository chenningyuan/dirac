Past and current works on DIRAC development
==========================================

This web-page is dedicated to list various ways for the development of the DIRAC software.
Contains mainly collection of finished and on-going bachelor and master thesis for IT-students, plus general wishlist.

DIRAC developers should share their ideas about IT-themes here. Each theme has to have short desciption (annotation), author (supervisor)
and assigned student if any.


Thesis list
-----------

Comparison of CPU and GPU mathematical libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The aim of this thesis is to compare the speed of some numerical procedures from CUBLAS and CULA libraries 
with their Fotran’s analogues BLAS and LAPACK libraries. 
We wrote a short working program (utils/gpu_test.F90) in Fortran and C programming languages, 
which compares the time complexity for various input parameters. 
Potential of this master thesis is in possible applying of CUBLAS / CULA library functions in the entire DIRAC program suite. 
(Finished in summer 2012, student D.Kuzma, supervised by Miro Ilias.)


CMake platform universal system for the software buildup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First, learn the CMake buildup system on simple testing programs written in C,C++,Fortran languages.
Later, try to program more complex CMake buildup receipes, including macros.
Assigned student I.Hrasko. Miro Ilias was supervising, Rado Bast wa referee. This bachelor project was finished in summer 2013.


Documentation generator Sphinx and its usage for DIRAC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The Sphinx documentation system, based upon Python language, has great potential for creation of on-line documentation.
This is visible in the DIRAC documentation manual. However, there is place for further many improvements regarding the documentation
framework.
The assigned IT-student (J.Bubniak, supervised by Miro Ilias) will polish bibliography (at least better formatted tables are needed),
math euqations (export in various formats), 'like'-buttons etc. Student's master thesis defence is expected in summer 2015.


Adpatation of DIRAC for Microsoft Windows and high-performace comuting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We need to adapt DIRAC thorougly for the MS Windows operating system and enable high-performace comuting
also on this operating system.
Assigned and currently working student I.Hrasko ( Miro Ilias supervising). This master project is to be finished in summer 2015.


IT-aspects of grid-computations with DIRAC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Complex adaptation of DIRAC for computing in the grid environment. Create set of universal, low-level scripts (bash), handling the DIRAC
suite calculations on unknown grid computing elements. Supervisor Miro Ilias is currently looking for IT-student.


Testing various compilers and mathematical libraries for DIRAC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
It is know that selection of compilers, their flags and mathematical libraries has influence on the code's performance.
The aim is to employ variouys compilers (GNU, Intel, Portland), available libraries (MKL, ACML, ATLAS) for DIRAC  performance measuraments.
We hope to get optimal set of compiler flags and linked libraries giving the best execution timing.
Supervisor Miro Ilias got an IT-student, hopefully he will start doing something from autumn 2014.


Git and web-based project development platforms in the service of DIRAC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The Git version control system is very popular for code development. On top of it, various web-based project development 
platforms (assembla, bitbucket, deveo, gitorious etc) are heplful for managing git repositories.
The aim is to describe (in the form of manual) the Git system and available web-based project development servers
with their various features and show how all these goodies can be applied in the collective DIRAC development.


Other TODO wishes
-----------------

New input reader with documenting of keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Let us write/use an keyword input reader which requires keywords to have documentation.
Otherwise the code won't even compile. Then the "make html" command would go through the source code and generate
the RST which would exactly reflect the state of the code (unless of course somebody changes the keyword without
changing the text 2 lines below but we could introduce some punishment for it). Radovan Bast.
