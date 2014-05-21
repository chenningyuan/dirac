Past and current works on DIRAC development
==========================================

This web-page lists various ways for the development of the DIRAC software.
It is mainly the collection of finished, on-going and proposed bachelor and master thesis projects for IT-students, 
plus general wishlist.

Master students of IT at universities around the world have potential to contribute to the DIRAC development.
For that, one has specify suitable "IT-only" projects around DIRAC, without the deep

DIRAC developers are invited to share their ideas about possible IT qualification themes here. 
Each proposed project has to have short desciption (annotation), author (supervisor)
and assigned student if any.


Thesis list
-----------

Comparison of CPU and GPU mathematical libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The aim of this thesis is to compare the speed of some numerical procedures from CUBLAS and CULA libraries 
with their Fotran’s analogues BLAS and LAPACK libraries. 
We wrote a short working program (in utils/gpu_test.F90) in Fortran and C programming languages, 
which compares the time complexity for various input parameters (matrix, vector sizes etc). 
The usefullness of this master thesis is in possible application of CUBLAS / CULA library functions 
in the entire DIRAC program suite. 
(Finished in summer 2012, student D.Kuzma, supervised by Miro Ilias.)


CMake platform universal system for the software buildup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First, learn the CMake buildup system on simple testing programs written in C,C++,Fortran languages.
Later, try to program more complex CMake buildup receipes, including macros.
Assigned student I.Hrasko. Miro Ilias was supervising, Rado Bast wa referee. 
This bachelor project was finished in summer 2013.


Documentation generator Sphinx and its usage for DIRAC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The Sphinx documentation system, based upon Python language, has great potential for creation of the on-line documentation.
This is clearly demonstrated on the DIRAC documentation manual. 
However, there is place for plenty of various improvements in the documentation framework.
The assigned IT-student (J.Bubniak, supervised by Miro Ilias) is expected to polish the bibliography 
(at least better formatted table of references is needed),
math equations (export in various formats), add social networks 'like'-buttons etc. 
Student's master thesis defence is expected in summer 2015.


Adpatation of DIRAC for Microsoft Windows and high-performace comuting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We need to adapt DIRAC thorougly for the MS Windows operating system and enable high-performace computing
also on this (commertial) operating system.
Assigned and currently working student is I.Hrasko, who gained rich experience with CMake (Miro Ilias supervising). 
His master project is to be finished in summer 2015.


IT-aspects of grid-computations with DIRAC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Complex adaptation of DIRAC for computing in the grid environment. Create set of universal, low-level scripts (bash), handling the DIRAC
suite calculations on unknown grid computing elements. Supervisor Miro Ilias is currently looking for an IT-student.


Testing various compilers and mathematical libraries for DIRAC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
It is known that selection of compilers, their flags and mathematical libraries has influence on the code's performance.
The aim is to employ variouys compilers (GNU, Intel, Portland), available libraries (MKL, ACML, ATLAS, netlib) for the DIRAC performance measurements.
We hope to get set of compiler flags and linked libraries giving the best execution timing.
Supervisor Miro Ilias got an IT-student for bachelor thesis, hopefully he will start doing something from autumn 2014.


Git and web-based project development platforms in the service of DIRAC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The Git version control system is very popular for the code development. On top of it, various web-based project development 
platforms (like assembla, bitbucket, deveo, gitorious etc) are heplful for managing git repositories.
The aim is to describe (in the form of manual) the Git system and available web-based project development servers
with their various features and show how all these goodies can be applied in our collective DIRAC development.
This project is suitable for bachelor thesis.


Other TODO wishes
-----------------

New input reader with documenting of keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Let us write/use an keyword input reader which requires keywords to have some documentation.
Otherwise the code won't even compile. Then the "make html" command would go through the source code and generate
the RST files which would exactly reflect the state of the code (unless of course somebody changes the keyword without
changing the text 2 lines below but we could introduce some punishment for it). Supervisor Radovan Bast.


CMake configuration files bootstrap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
It is possible to download external cmake "libraries" at configure time and include/employ them in the code's buildup.
Radovan has tested this with the math detection which he puts on some external test repo and at configure time
it was correctly fetched and used.
(He also asked on the cmake mailing list for pointers and got useful info which worked for him.)
We could do now the outsource of all generic cmake code out of DIRAC/Dalton/... suites and let a minimum cmake
core fetch them at configure time.

There is one catch: one should be able to configure and build without the network.
This could be solved with a fallback cmake system that would be used if the default cmake repo is unreachable.
There are several ways to do it.

If you think such a bootstrap approach would be interesting we could brainstorm more.
Radovan has a local proof that it works and we could continue from here if this is interesting.

