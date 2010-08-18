***************************************************
*
* New assigment of disk units for KRCC
*
* Lasse 2010
*
***************************************************
      SUBROUTINE DISKUN_KRCC
#include "implicit.inc"
*
* Assign logical unit numbers for KRCC:
*==================================================
*
*   1 -  9 :  Retained for I/O
*  ---------
*        5 :  Standard IN
*        6 :  Standard OUT
*
*  50 - 59 : Scratch files (added here for cc)
*  ---------
*
*  60 - 69 :  Important files retained (DIAG,C,HC etc.)
*  ---------
*       60 :    CI diagonal
*       61 :    CI reference vector
*       62 :    CI sigma vector
*       63 :    CC correction vectors
*       64 :    approximations to the CC solutions
*       65 :    CC vector functions
*       67 :    Last CC coefficients
*       68 :    Last CC vector function
*
*  70 - 79 :  Scratch files
*  ---------
*
*  80 - 99 :  Files in use in different context
*             (e.g. standalone, INS, MOLCAS etc.)
*  ---------
*
*==================================================
*
#include "clunit.inc"
* =========================
* Standard input and output
* =========================
*. Input file
      LUIN = 5
*. Special unit for input file from MOLCAS conversion routine
*    Timo Fleig and Lasse Sorensen
      LUCMOL = 91
*. Output file
      LUOUT = 6
* =================
* Input information
* =================
* Input file containing MO-AO transformation matrix
      LUMOIN = 80
*. Input file for CI-vectors
*. restart from file 21 is assumed
*. INTEGRAL FILE FROM SORTING
      LUSC81 = 81
*. INDICES FILE FROM SORTING
      LUSC82 = 82
*. Input , property one-electron integral files
      LUPRP  = 83
*. Sirius interface file
      LUSIR1 = 84
*. File containing additional states for transition densities
      LUEXC = 85
* =================
* Internal files
* =================
*. CI diagonal
      LUDIA = 60
      Open(Unit=LUDIA,File='CIDIA',Status='UNKNOWN',Form='UNFORMATTED')
*. CI vector
      LUC = 61
      Open(Unit=LUC,File='LUCVECT',Status='UNKNOWN',Form='UNFORMATTED')
*. Sigma vector file
      LUHC = 62
      Open(Unit=LUHC,File='HCFILE',Status='UNKNOWN',Form='UNFORMATTED')
*. File collecting CC correction vectors, used for DIIS etc
      LU_CCVEC = 63
*. File containing approximations to the CC solutions
      LU_CCVECT = 64
*. File containing CC vector functions for the CC vectors on LU_CCVECT
      LU_CCVECF = 65
*. File containing Last CC coefficients
      LU_CCVECL = 67
*. File containing Last CC vector function
      LU_CCVECFL = 68
* =================
* Scratch files
* =================
      LUSC1 = 70
      LUSC2 = 71
      LUSC3 = 72
      Open(Unit=LUSC1,File='LUSC1',Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC2,File='LUSC2',Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC3,File='LUSC3',Status='UNKNOWN',Form='UNFORMATTED')
*. Scratch space for subspace handling
      LUSC34 = 73
      LUSC35 = 74
      LUSC36 = 75
      LUSC37 = 76
      LUSC38 = 77
      LUSC39 = 78
      LUSC40 = 79
      Open(Unit=LUSC34,File='LUSC34',Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC35,File='LUSC35',Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC36,File='LUSC36',Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC37,File='LUSC37',Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC38,File='LUSC38',Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC39,File='LUSC39',Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC40,File='LUSC40',Status='UNKNOWN',
     &     Form='UNFORMATTED')
*.
      LUSC51 = 51
      LUSC52 = 52
      LUSC53 = 53
*. Used for CHAR_TO_INTEGER_MOLUC.
      LUSC54 = 54

* =================
* Output files
* =================
*. output file for CI-vectors
*. Not in use
      LUCIVO = 86
*. Natural orbitals in terms of input orbitals
*.
      LUMOUT = 87
*     Open(Unit=LUMOUT,File='CIAORB',Status='UNKNOWN',Form='FORMATTED')
*. Dumping 1- and 2- electron integrals in formatted form
* ( LU90 just defined here, it is not exported )
      LU90  = 88
*. Dumping symmmetry info, MO-AO expansion matrix and property integrals
*. CC amplitudes in formatted form
      LU_CCAMP = 90
*. Result of CI=> CC conversion
      LU_CC_FROM_CI = 93
*. Excitation operators, all symmetries
      LU_CCEXC_OP = 94
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
