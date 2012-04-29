Compilation on various platforms
================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Examples of various platforms and compilers where DIRAC works.
Compilation flags are stored into the Makefile.config file, which is
reviewed here for various platforms.

Contents
--------

-  `1 WORKMEM <#WORKMEM>`_
-  `2 gfortran/gcc <#gfortran.2Fgcc>`_
-  `3 ifort / icc <#ifort_.2F_icc>`_
-  `4 Portland Group Compilers - pgf90 &
   pgcc <#Portland_Group_Compilers_-_pgf90_.26_pgcc>`_
-  `5 The "Horseshoe" cluster in Odense,
   Denmark <#The_.22Horseshoe.22_cluster_in_Odense.2C_Denmark>`_
-  `6 IBM/P1600Power5+ "ANAKIN" at CINES,
   France <#IBM.2FP1600Power5.2B_.22ANAKIN.22_at_CINES.2C_France>`_
-  `7 IBM Regatta p690+ "JUMP" at FZ Juelich,
   Germany <#IBM_Regatta_p690.2B_.22JUMP.22_at_FZ_Juelich.2C_Germany>`_
-  `8 IBM eServer (Regatta Power6) "VARGAS" at IDRIS,
   France <#IBM_eServer_.28Regatta_Power6.29_.22VARGAS.22_at_IDRIS.2C_France>`_
-  `9 Linux cluster "stallo" at UIT, Tromsø,
   Norway <#Linux_cluster_.22stallo.22_at_UIT.2C_Troms.C3.B8.2C_Norway>`_
-  `10 IBM p575+ "njord" at NTNU, Trondheim,
   Norway <#IBM_p575.2B_.22njord.22_at_NTNU.2C_Trondheim.2C_Norway>`_

WORKMEM
-------

The WORKMEM variable in ``Makefile.config`` contains the amount of
memory (number of 8 byte double precision numbers) which is allocated
for the DIRAC run.

Very small default value is set by the ``configure`` script. The user
may decide to keep thus small default, 10 000 000 (80MB), however, in
each DIRAC run he would have to set up larger memory space (depending on
the size of the system) using "pam -mw <master> -nw <node> ...".

Therefore it is advisable that the user configures DIRAC to use slightly
less than what his machine is capable of, for example WORKMEM=100000000
(100 milion) on a 1GB machine (this would leave about 200MB for the
operating system and static memory allocations in DIRAC.

Note that it does not make sense to use more than 2GB of memory on a 32
bit machine or a 64 bit machine running DIRAC in 32 bit mode.

gfortran/gcc
------------

Very often used free GNU compiler set for the serial compilation.

Note that DIRAC is compilable only with gfortran of version >= 4.2.

-  32 bit

... to be added (there are very few 32-bit servers nowadays)

-  64 bit

Full F90 optimization employed, plus - very important - default
"integer8" values.

::

    OPTFLAGS    =  -O3 -ffast-math -funroll-loops -ftree-vectorize -fdefault-integer-8 -J$(SOURCEDIR)/modules
    LIB_LIST    =  <SOME_ILP64_PREPARED_LIBRARY>

This is for tracking back routines in the case of crashing.

::

    LDFLAGS     =  -Wl,-E

See the
`http://wiki.chem.vu.nl/dirac/index.php/Installation\_guide <http://wiki.chem.vu.nl/dirac/index.php/Installation_guide>`_
.

ifort / icc
-----------

Commercial compiler for intel-based machines. Yields good performance on
intel macintoshes. For these systems: answer 0 (Configure libraries
manually) to the question about mathematical libraries and edit
Makefile.config manually. It should have:

LIBS = -framework Accelerate

to make DIRAC use the optimized BLAS libraries available for MacOsX.

Portland Group Compilers - pgf90 & pgcc
---------------------------------------

The 5.2.4 version of pgf90 (from 2004) is no longer usable. It does not
accept allocatable arrays appearing in the type definition.

However, the 'most recent' 7.2-3 version of Fortran compiler is NOT
working in the optimization mode ! The executable (compiled even with
the lowest optimization level -O1) is crashing after the 1st iteration.
(Developers - see ticket 34)

The "Horseshoe" cluster in Odense, Denmark
------------------------------------------

`http://www.dcsc.sdu.dk/ <http://www.dcsc.sdu.dk/>`_

-  serial compilation

No problem with ifort/icc

-  parallel compilation

For parallel compilation with ifort-based mpif90 (if you also have the
mkl library available): add -DVAR\_MKL flag to
CPP\_FLAGS/DEPEND\_FLAGS/CCP\_FLAGS variables.

In "Makefile" file, which is created by "pamadm", change order of
libraries for the 'dirmpi.x' executable: instead of "$(LIB\_LIST)
$(PDLIB)" use "$(PDLIB) $(LIB\_LIST)"

The library:

::

    -L/opt/intel/mkl/9.0_beta/lib/em64t/ -lguide -lmkl_lapack -lmkl_em64t -lvml -lpthread

IBM/P1600Power5+ "ANAKIN" at CINES, France
------------------------------------------

`http://www.cines.fr/spip.php?article503 <http://www.cines.fr/spip.php?article503>`_

The architecture is recognized as powerpc-ibm-aix5.3.0.0. As compilers I
use

::

    F90           = mpxlf_r
    CC            = mpcc_r

with corresponding compiler flags

::

    OPTFLAGS      = -O3 -qstrict -qextname -q64 -qmoddir=$(SOURCEDIR)/module
    C_OPTFLAGS    = -O3 -qstrict -qcpluscmt -q64

mathematical libraires

::

    LIB_LIST      = /usr/local/pub/LAPACK/lapack.a -L/usr/lib -lesslsmp

and loader flags

::

    LDFLAGS       = -q64 -bmaxdata:0x8000000000

However, straighforward compilation gives the following errors in the
load stage

::

    ld: 0711-317 ERROR: Undefined symbol: .zheev_
    ld: 0711-317 ERROR: Undefined symbol: .dsyev_
    ld: 0711-317 ERROR: Undefined symbol: .lsame_
    ld: 0711-317 ERROR: Undefined symbol: .xerbla_
    ld: 0711-317 ERROR: Undefined symbol: .dlamch_
    ld: 0711-317 ERROR: Undefined symbol: .zlanhb_
    ld: 0711-317 ERROR: Undefined symbol: .zlascl_
    ld: 0711-317 ERROR: Undefined symbol: .zhbtrd_
    ld: 0711-317 ERROR: Undefined symbol: .dsterf_
    ld: 0711-317 ERROR: Undefined symbol: .zsteqr_
    ld: 0711-317 ERROR: Undefined symbol: .dsbev_
    ld: 0711-317 ERROR: Undefined symbol: .zgesvd_ 

These are complex LAPACK routines and the problem is the compiler flag
-qextname which adds a trailing underscore to the names of all global
entities, except for main program names. However, removing the compiler
flag leads to problems in that a number of internal C routines of DIRAC
are not recognized. To correctly link the code I therefore renamed these
routines as follows

::

    LDFLAGS       =  -q64 -bmaxdata:0x8000000000 -brename:.zheev_,.zheev -brename:.dsyev_,.dsyev -brename:.lsame_,.lsame  -brename:.xerbla_,.xerbla -brename:.dlamch_,.dlamch -brename:.zlanhb_,.zlanhb -brename:.zlascl_,.zlascl -brename:.zhbtrd_,.zhbtrd -brename:.dsterf_,.dsterf -brename:.zsteqr_,.zsteqr -brename:.dsbev_,.dsbev -brename:.zgesvd_,.zgesvd

and the code now installs correctly.

IBM Regatta p690+ "JUMP" at FZ Juelich, Germany
-----------------------------------------------

`http://www.fz-juelich.de/jsc/service/sco\_ibmRegatta <http://www.fz-juelich.de/jsc/service/sco_ibmRegatta>`_

The architecture is recognized as powerpc-ibm-aix5.3.0.0. As compilers I
use

::

    F90           = mpxlf_r
    F77           = mpxlf_r
    CC            = mpcc_r

with corresponding compiler flags (for production runs)

::

    F77FLAGS      = -O3 -qessl -qenablevmx -qhot=simd -qstrict -qextname -q64 -qtune=pwr6 -qarch=pwr6 -bdatapsize:64K -bstackpsize:64K -btextpsize:64K -Iinclude -qmoddir=modules -Imodules
    F90FLAGS      = -O3 -qessl -qenablevmx -qhot=simd -qstrict -qextname -q64 -qtune=pwr6 -qarch=pwr6 -bdatapsize:64K -bstackpsize:64K -btextpsize:64K -Iinclude -qmoddir=modules -Imodules -qfree=f90
    CFLAGS        = -O3 -qstrict -qcpluscmt -q64 -qtune=pwr6 -qarch=pwr6 -bdatapsize:64K -bstackpsize:64K -btextpsize:64K -Iinclude

 mathematical libraires

::

    LIB_LIST      = -lpesslsmp -lscalapack -lblacssmp -llapack -lesslsmp

and loader flags

::

    LDFLAGS       = -q64 -bmaxdata:0x8000000000 -brename:.zheevd_,.zheevd -brename:.zgesvd_,.zgesvd -brename:.dsyevr_,.dsyevr

to prepare for compilation it is mandatory to generate a file named
fortran\_interface.h by running

::

    your_path_to_dirac/build/gen_fortran_interface_h > include/fortran_interface.h

via a job script. This file has to be placed into the
your\_path\_to\_dirac/include directory. Now type

::

    gmake

and the compilation will be successful.

IBM eServer (Regatta Power6) "VARGAS" at IDRIS, France
------------------------------------------------------

`http://www.idris.fr/comp/scal/power6/vargas/index-vargas.html <http://www.idris.fr/comp/scal/power6/vargas/index-vargas.html>`_

The architecture is recognized as powerpc-ibm-aix5.3.0.0 (so identical
to
`Juelich </dirac/index.php/Examples_of_various_platforms/compilers#IBM_Regatta_p690.2B_.C2.AB_JUMP_.C2.BB_at_FZ_Juelich.2C_Germany>`_).
As compilers I use

::

    F77=/usr/bin/mpxlf_r
    F90=/usr/bin/mpxlf_r
    CC =/usr/bin/mpcc_r

with corresponding compiler flags (for production runs)

::

    F77FLAGS=-g -O3 -qstrict -qextname -q64 -Iinclude -qmoddir=modules -Imodules
    F90FLAGS=-g -O3 -qstrict -qextname -q64 -Iinclude -qmoddir=modules -Imodules -qfree=f90 
    CFLAGS=-g -O3 -qstrict -qcpluscmt -q64 -Iinclude

mathematical libraires

::

    LIBS      =  -lesslsmp

and loader flags

::

    LDFLAGS=-q64 -bmaxdata:0x8000000000

Compilation with the native ``make`` does not work at all. ``gmake``
works, provided that *python* is loaded:

::

    module load python

Also, the LAPACK library must be loaded

::

    module load lapack

Compilation now gives

::

    ld: 0711-317 ERROR: Undefined symbol: .dsyevr_
    ld: 0711-317 ERROR: Undefined symbol: .zheevd_
    ld: 0711-317 ERROR: Undefined symbol: .zgesvd_

These are complex LAPACK routines and the problem is the compiler flag
-qextname which adds a trailing underscore to the names of all global
entities, except for main program names. However, removing the compiler
flag leads to problems in that a number of internal C routines of DIRAC
are not recognized. To correctly link the code I therefore renamed these
routines as follows

::

    LDFLAGS=-q64 -brename:.dsyevr_,.dsyevr  -brename:.zheevd_,.zheevd -brename:.zgesvd_,.zgesvd

Note also that the bmaxdata flag was removed, following advice from the
helpful IDRIS technical staff. bmaxdata is not needed on 64-bit
machines.

Linux cluster "stallo" at UIT, Tromsø, Norway
---------------------------------------------

`http://wiki.ctcc.no/doku.php?id=public:howto:dirac#compilation\_on\_linux\_cluster\_stallo\_at\_uit\_tromso\_norway <http://wiki.ctcc.no/doku.php?id=public:howto:dirac#compilation_on_linux_cluster_stallo_at_uit_tromso_norway>`_

IBM p575+ "njord" at NTNU, Trondheim, Norway
--------------------------------------------

`http://wiki.ctcc.no/doku.php?id=public:howto:dirac#compilation\_on\_ibm\_p575\_njord\_at\_ntnu\_trondheim\_norway <http://wiki.ctcc.no/doku.php?id=public:howto:dirac#compilation_on_ibm_p575_njord_at_ntnu_trondheim_norway>`_

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Compilation\_on\_various\_platforms <http://wiki.chem.vu.nl/dirac/index.php/Compilation_on_various_platforms>`_"

Views
'''''

-  `Article </dirac/index.php/Compilation_on_various_platforms>`_
-  `Discussion </dirac/index.php?title=Talk:Compilation_on_various_platforms&action=edit>`_
-  `View
   source </dirac/index.php?title=Compilation_on_various_platforms&action=edit>`_
-  `History </dirac/index.php?title=Compilation_on_various_platforms&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Compilation_on_various_platforms>`_

` </dirac/index.php/Dirac_Program>`_

Navigation
''''''''''

-  `Home </dirac/index.php/Dirac_Program>`_
-  `Documentation </dirac/index.php/Documentation>`_
-  `Reference manual </dirac/index.php/Manual>`_
-  `Discussion forum </dirac/index.php/DiscussionForum>`_
-  `News Archive </dirac/index.php/News_archive>`_
-  `Basis sets </dirac/index.php/Basis_Sets>`_

Search
''''''

Toolbox
'''''''

-  `What links
   here </dirac/index.php/Special:Whatlinkshere/Compilation_on_various_platforms>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Compilation_on_various_platforms>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Compilation_on_various_platforms&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Compilation_on_various_platforms&oldid=22761>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 14:49, 13 September 2010.
-  This page has been accessed 3,443 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
