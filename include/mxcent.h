!
!     file: mxcent.h
!
!     MXATOM = max number of atomic types
!     MXCENT = max number of nuclei + point charges + ghost orbital centers
!
!     IF you change MXCENT you should do a "make depend"
!     and then rebuild the program using the command "make".
!
      INTEGER MXATOM, MXCENT, MXCOOR
!
! WARNING: LARGE_MOLECULES is used in maxorb.h and mxcent.h.
!          If you define it you HAVE TO define it in both files
!          otherwise you may experience bugs and incorrect results.
!          the easiest and safest way to define LARGE_MOLECULES is using
!          "MODS=-DLARGE_MOLECULES" in build/Makefile.moddef,
!          alternatively using -DLARGE_MOLECULES in Makefile.config.
!          You have to recompile DIRAC after you have modified one of these files.
#if ! defined (LARGE_MOLECULES)
      PARAMETER (MXATOM=10, MXCENT = 200, MXCOOR = 3*MXCENT)
#else
      PARAMETER (MXATOM=20, MXCENT = 600, MXCOOR = 3*MXCENT)
#endif
