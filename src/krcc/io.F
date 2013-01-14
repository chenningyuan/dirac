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
*
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
*. Output file
      LUOUT = 6
* =================
* Input information
* =================
*. Input file for CI-vectors
*. restart from file 21 is assumed
*. INTEGRAL FILE FROM SORTING
      LUSC81 = 51
*. INDICES FILE FROM SORTING
      LUSC82 = 52
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
      LUSC4 = 73 !For intermediates
      LUSC5 = 74 !For intermediates
      LUSC6 = 75 !For lrcc
      Open(Unit=LUSC1,File='LUSC1',Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC2,File='LUSC2',Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC3,File='LUSC3',Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC4,File='LUSC4',Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC5,File='LUSC5',Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC6,File='LUSC6',Status='UNKNOWN',Form='UNFORMATTED')
*     LUJDIA = 76 !optimization
*     LUJVEC = 77 !optimization
*     Open(Unit=LUJDIA,File='LUJDIA',
*    &     Status='UNKNOWN',Form='UNFORMATTED')
*     Open(Unit=LUJVEC,File='LUJVEC',
*    &     Status='UNKNOWN',Form='UNFORMATTED')
*. Scratch space for subspace handling
      LUSC41 = 79
      Open(Unit=LUSC41,File='LUSC41',Status='UNKNOWN',
     &     Form='UNFORMATTED')
* ===========================
* Intermediate files for lrcc
* ===========================
      LU_M10 = 80
      LU_M02 = 81
      LU_M11 = 82
      LU_M01 = 83
      LU_M12 = 84
      LU_M11_TO_M01 = 85
      Open(Unit=LU_M10,File='LU_M10',
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LU_M02,File='LU_M02',
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LU_M11,File='LU_M11',
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LU_M01,File='LU_M01',
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LU_M12,File='LU_M12',
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LU_M11_TO_M01,File='LU_M11_TO_M01',
     &     Status='UNKNOWN',Form='UNFORMATTED')
*.
* =================
* Output files
* =================
*. CC amplitudes in formatted form
      LU_CCAMP = 90
*. Excitation operators, all symmetries
      LU_CCEXC_OP = 94
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE TEST_ARRAY(ARRAY,NDIM)
C
#include "implicit.inc"
C
      DIMENSION ARRAY(NDIM)
C
      PRINT*,'NDIM',NDIM
      DO I=1,NDIM
       ARRAY(I) = I
       PRINT*,'ARRAY,I',ARRAY(I),I
      END DO
C
      RETURN
      END
