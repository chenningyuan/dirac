C
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C
#if defined(VAR_MPI) && !defined(SYS_LOCALDISKS)
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
*
      SUBROUTINE MVCSMD_PCIO_1(FAC,VEC1,VEC2,NVEC,IREW,
     &           JOFFS,JTOTBLL,JBKINF,JBLOCKI)
C
C ADD VECTORS ON FILE ILUIN TIMES FACTOR AND STORE ON LUOUT
C
C LUOUT AND LUSCR ARE INITIALLY REWOUND
C
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(NVEC), JBKINF(JBLOCKI)
      INTEGER ISWAP 
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
C
      IOFFSET1= 0
      IOFFSET2= 0
      IOFFSET3= 0
      ISWAP = 0
      IF( MOD(NVEC,2) .EQ. 0 ) THEN
        ISWAP = 1
      ELSE
        ISWAP = 2
      END IF
C
      IF(IREW .NE. 0 ) THEN
        IOFFSET1 = 0
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
      END IF
*
*      write(6,*)'ready to die',myproc
C
      DO 100 IVEC = 1, NVEC
*
        IOFFSET1 = (IVEC - 1) * JTOTBLL
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
*
        IF( IVEC .EQ. 1 ) THEN
          CALL SCLVCD_PCIO_1(FAC(IVEC),VEC1,IOFFSET1,
CSK          CALL SCLVCD_PCIO(LUIN,LLUOUT,FAC(IVEC),VEC1,IOFFSET1,
     &         IOFFSET2,JBKINF,JBLOCKI,ISWAP)
        ELSE
          CALL VECSMD_PCIO_1(VEC1,VEC2,FAC(IVEC),1.0D0,
CSK          CALL VECSMD_PCIO_1(VEC1,VEC2,FAC(IVEC),1.0D0,LUIN,LLUSCR,LLUOUT,
     &         IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)
        END IF
C
        IF(ISWAP.EQ.1)THEN
          ISWAP = 2
          GOTO 99
        ELSE IF(ISWAP.EQ.2) THEN
          ISWAP = 1
          GOTO 99
        END IF
*
   99 CONTINUE
*
  100 CONTINUE
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
*
      SUBROUTINE MVCSMD_PCIO_2(FAC,VEC1,VEC2,NVEC,IREW,
     &           JOFFS,JTOTBLL,JBKINF,JBLOCKI)
C
C ADD VECTORS ON FILE ILUIN TIMES FACTOR AND STORE ON LUOUT
C
C LUOUT AND LUSCR ARE INITIALLY REWOUND
C
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(NVEC), JBKINF(JBLOCKI)
      INTEGER ISWAP 
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
C
      IOFFSET1= 0
      IOFFSET2= 0
      IOFFSET3= 0
      ISWAP = 0
      IF( MOD(NVEC,2) .EQ. 0 ) THEN
        ISWAP = 1
      ELSE
        ISWAP = 2
      END IF
C
      IF(IREW .NE. 0 ) THEN
        IOFFSET1 = 0
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
      END IF
*
*      write(6,*)'ready to die',myproc
C
      DO 100 IVEC = 1, NVEC
*
        IOFFSET1 = (IVEC - 1) * JTOTBLL
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
*
        IF( IVEC .EQ. 1 ) THEN
          CALL SCLVCD_PCIO_2(FAC(IVEC),VEC1,IOFFSET1,
CSK          CALL SCLVCD_PCIO(LUIN,LLUOUT,FAC(IVEC),VEC1,IOFFSET1,
     &         IOFFSET2,JBKINF,JBLOCKI,ISWAP)
        ELSE
          CALL VECSMD_PCIO_2(VEC1,VEC2,FAC(IVEC),1.0D0,
CSK          CALL VECSMD_PCIO_1(VEC1,VEC2,FAC(IVEC),1.0D0,LUIN,LLUSCR,LLUOUT,
     &         IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)
        END IF
C
        IF(ISWAP.EQ.1)THEN
          ISWAP = 2
          GOTO 99
        ELSE IF(ISWAP.EQ.2) THEN
          ISWAP = 1
          GOTO 99
        END IF
*
   99 CONTINUE
*
  100 CONTINUE
C
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MVCSMD2_PCIO_1(FAC,FACLUOUT,VEC1,VEC2,
     &                   NVEC,IREW,JOFFS,JTOTBLL,JBKINF,JBLOCKI)
*
* LUOUT = Factor * LUOUT + sum(I=1,nvec)fac(I) LUIN(I)
*
* LUOUT and LUSCR are rewinded and input data
* on these files are lost.
*
* MVCSMD2 is identical MVCSMD, except that
* the input vector on LUOUT can be included
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(1), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
C
*
      IOFFSET1= 0
      IOFFSET2= 0
      IOFFSET3= 0
C
      IF(IREW .NE. 0 ) THEN
        IOFFSET1 = 0
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
      END IF
      ISWAP = 0
      IF(FACLUOUT.NE.0.0D0) THEN
      ISWAP = 2
      ELSE
      ISWAP = 1
      END IF

      DO IVEC = 1, NVEC
*. start reading/writing at the right position
        IOFFSET1 = (IVEC - 1) * JTOTBLL
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
*
        IF( IVEC .EQ. 1 ) THEN
          IF(FACLUOUT.NE.0.0D0) THEN
            CALL VECSMD_PCIO_1B(VEC1,VEC2,FAC(1),FACLUOUT,
CSK            CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(1),FACLUOUT,LUIN,LUOUT,LUSCR,
     &           IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)
          ELSE
            CALL SCLVCD_PCIO_1B(FAC(1),VEC1,IOFFSET1,
CSK            CALL SCLVCD_PCIO_3(LUIN,LUSCR,FAC(1),VEC1,IOFFSET1,
     &           IOFFSET2,JBKINF,JBLOCKI,ISWAP)
          END IF
        ELSE
          ONE = 1.0D0
          CALL VECSMD_PCIO_1B(VEC1,VEC2,FAC(IVEC),ONE,
CSK          CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(IVEC),ONE,LUIN,LLUOUT,LLUSCR,
     &         IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)

        END IF
*
        IF(ISWAP.EQ.1) THEN
          ISWAP = 2
          GOTO 99
        ELSE IF(ISWAP.EQ.2) THEN
          ISWAP = 1
          GOTO 99
        END IF
 99   CONTINUE
      IF(IVEC .EQ. 1) ISWAP = 1
      END DO
*
      IF(ISWAP.EQ.1.AND.NVEC.GT.0) THEN
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
        CALL CPVC_PARPAR_68(VEC1,JBKINF,JBLOCKI,IOFFSET3,
     &           IOFFSET2,1)
      END IF
*
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MVCSMD2_PCIO_2(FAC,FACLUOUT,VEC1,VEC2,
     &                   NVEC,IREW,JOFFS,JTOTBLL,JBKINF,JBLOCKI)
*
* LUOUT = Factor * LUOUT + sum(I=1,nvec)fac(I) LUIN(I)
*
* LUOUT and LUSCR are rewinded and input data
* on these files are lost.
*
* MVCSMD2 is identical MVCSMD, except that
* the input vector on LUOUT can be included
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(1), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
C
*
      IOFFSET1= 0
      IOFFSET2= 0
      IOFFSET3= 0
C
      IF(IREW .NE. 0 ) THEN
        IOFFSET1 = 0
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
      END IF
      ISWAP = 0
      IF(FACLUOUT.NE.0.0D0) THEN
      ISWAP = 2
      ELSE
      ISWAP = 1
      END IF

      DO IVEC = 1, NVEC
*. start reading/writing at the right position
        IOFFSET1 = (IVEC - 1) * JTOTBLL
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
*
        IF( IVEC .EQ. 1 ) THEN
          IF(FACLUOUT.NE.0.0D0) THEN
            CALL VECSMD_PCIO_2B(VEC1,VEC2,FAC(1),FACLUOUT,
CSK            CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(1),FACLUOUT,LUIN,LUOUT,LUSCR,
     &           IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)
          ELSE
            CALL SCLVCD_PCIO_2B(FAC(1),VEC1,IOFFSET1,
CSK            CALL SCLVCD_PCIO_3(LUIN,LUSCR,FAC(1),VEC1,IOFFSET1,
     &           IOFFSET2,JBKINF,JBLOCKI,ISWAP)
          END IF
        ELSE
          ONE = 1.0D0
          CALL VECSMD_PCIO_2B(VEC1,VEC2,FAC(IVEC),ONE,
CSK          CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(IVEC),ONE,LUIN,LLUOUT,LLUSCR,
     &         IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)

        END IF
*
        IF(ISWAP.EQ.1) THEN
          ISWAP = 2
          GOTO 99
        ELSE IF(ISWAP.EQ.2) THEN
          ISWAP = 1
          GOTO 99
        END IF
 99   CONTINUE
      IF(IVEC .EQ. 1) ISWAP = 1
      END DO
*
      IF(ISWAP.EQ.1.AND.NVEC.GT.0) THEN
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
        CALL CPVC_PARPAR_67(VEC1,JBKINF,JBLOCKI,IOFFSET3,
     &           IOFFSET2,1)
      END IF
*
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MVCSMD2_PCIO_3(FAC,FACLUOUT,VEC1,VEC2,
     &                   NVEC,IREW,JOFFS,JTOTBLL,JBKINF,JBLOCKI)
*
* LUOUT = Factor * LUOUT + sum(I=1,nvec)fac(I) LUIN(I)
*
* LUOUT and LUSCR are rewinded and input data
* on these files are lost.
*
* MVCSMD2 is identical MVCSMD, except that
* the input vector on LUOUT can be included
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(1), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
C
*
      IOFFSET1= 0
      IOFFSET2= 0
      IOFFSET3= 0
C
      IF(IREW .NE. 0 ) THEN
        IOFFSET1 = 0
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
      END IF
      ISWAP = 0
      IF(FACLUOUT.NE.0.0D0) THEN
      ISWAP = 2
      ELSE
      ISWAP = 1
      END IF

      DO IVEC = 1, NVEC
*. start reading/writing at the right position
        IOFFSET1 = (IVEC - 1) * JTOTBLL
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
*
        IF( IVEC .EQ. 1 ) THEN
          IF(FACLUOUT.NE.0.0D0) THEN
            CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(1),FACLUOUT,
CSK            CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(1),FACLUOUT,LUIN,LUOUT,LUSCR,
     &           IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)
          ELSE
            CALL SCLVCD_PCIO_3(FAC(1),VEC1,IOFFSET1,
CSK            CALL SCLVCD_PCIO_3(LUIN,LUSCR,FAC(1),VEC1,IOFFSET1,
     &           IOFFSET2,JBKINF,JBLOCKI,ISWAP)
          END IF
        ELSE
          ONE = 1.0D0
          CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(IVEC),ONE,
CSK          CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(IVEC),ONE,LUIN,LLUOUT,LLUSCR,
     &         IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)

        END IF
*
        IF(ISWAP.EQ.1) THEN
          ISWAP = 2
          GOTO 99
        ELSE IF(ISWAP.EQ.2) THEN
          ISWAP = 1
          GOTO 99
        END IF
 99   CONTINUE
      IF(IVEC .EQ. 1) ISWAP = 1
      END DO
*
      IF(ISWAP.EQ.1.AND.NVEC.GT.0) THEN
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
        CALL CPVC_PARPAR_65(VEC1,JBKINF,JBLOCKI,IOFFSET3,
     &           IOFFSET2,1)
      END IF
*
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MVCSMD2_PCIO_4(FAC,FACLUOUT,VEC1,VEC2,
     &                   NVEC,IREW,JOFFS,JTOTBLL,JBKINF,JBLOCKI)
*
* LUOUT = Factor * LUOUT + sum(I=1,nvec)fac(I) LUIN(I)
*
* LUOUT and LUSCR are rewinded and input data
* on these files are lost.
*
* MVCSMD2 is identical MVCSMD, except that
* the input vector on LUOUT can be included
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(1), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
C
*
      IOFFSET1= 0
      IOFFSET2= 0
      IOFFSET3= 0
C
      IF(IREW .NE. 0 ) THEN
        IOFFSET1 = 0
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
      END IF
      ISWAP = 0
      IF(FACLUOUT.NE.0.0D0) THEN
      ISWAP = 2
      ELSE
      ISWAP = 1
      END IF

      DO IVEC = 1, NVEC
*. start reading/writing at the right position
        IOFFSET1 = (IVEC - 1) * JTOTBLL
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
*
        IF( IVEC .EQ. 1 ) THEN
          IF(FACLUOUT.NE.0.0D0) THEN
            CALL VECSMD_PCIO_4(VEC1,VEC2,FAC(1),FACLUOUT,
CSK            CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(1),FACLUOUT,LUIN,LUOUT,LUSCR,
     &           IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)
          ELSE
            CALL SCLVCD_PCIO_4(FAC(1),VEC1,IOFFSET1,
CSK            CALL SCLVCD_PCIO_3(LUIN,LUSCR,FAC(1),VEC1,IOFFSET1,
     &           IOFFSET2,JBKINF,JBLOCKI,ISWAP)
          END IF
        ELSE
          ONE = 1.0D0
          CALL VECSMD_PCIO_4(VEC1,VEC2,FAC(IVEC),ONE,
CSK          CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(IVEC),ONE,LUIN,LLUOUT,LLUSCR,
     &         IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)

        END IF
*
        IF(ISWAP.EQ.1) THEN
          ISWAP = 2
          GOTO 99
        ELSE IF(ISWAP.EQ.2) THEN
          ISWAP = 1
          GOTO 99
        END IF
 99   CONTINUE
      IF(IVEC .EQ. 1) ISWAP = 1
      END DO
*
      IF(ISWAP.EQ.1.AND.NVEC.GT.0) THEN
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
        CALL CPVC_PARPAR_65(VEC1,JBKINF,JBLOCKI,IOFFSET3,
     &           IOFFSET2,1)
      END IF
*
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MVCSMD2_PSIO_5(FAC,FACLUOUT,VEC1,VEC2,
     &                   NVEC,IREW,JOFFS,JTOTBLL,JBKINF,JBLOCKI)
*
* LUOUT = Factor * LUOUT + sum(I=1,nvec)fac(I) LUIN(I)
*
* LUOUT and LUSCR are rewinded and input data
* on these files are lost.
*
* MVCSMD2 is identical MVCSMD, except that
* the input vector on LUOUT can be included
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(1), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
C
*
      IOFFSET1= 0
      IOFFSET2= 0
      IOFFSET3= 0
C
      IF(IREW .NE. 0 ) THEN
        IOFFSET1 = 0
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
      END IF
      ISWAP = 0
      IF(FACLUOUT.NE.0.0D0) THEN
      ISWAP = 2
      ELSE
      ISWAP = 1
      END IF

      DO IVEC = 1, NVEC
*. start reading/writing at the right position
        IOFFSET1 = (IVEC - 1) * JTOTBLL
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
*
        IF( IVEC .EQ. 1 ) THEN
          IF(FACLUOUT.NE.0.0D0) THEN
            CALL VECSMD_PSIO_5(VEC1,VEC2,FAC(1),FACLUOUT,
CSK            CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(1),FACLUOUT,LUIN,LUOUT,LUSCR,
     &           IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)
          ELSE
            CALL SCLVCD_PSIO_5(FAC(1),VEC1,IOFFSET1,
CSK            CALL SCLVCD_PCIO_3(LUIN,LUSCR,FAC(1),VEC1,IOFFSET1,
     &           IOFFSET2,JBKINF,JBLOCKI,ISWAP)
          END IF
        ELSE
          ONE = 1.0D0
          CALL VECSMD_PSIO_5(VEC1,VEC2,FAC(IVEC),ONE,
CSK          CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(IVEC),ONE,LUIN,LLUOUT,LLUSCR,
     &         IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)

        END IF
*
        IF(ISWAP.EQ.1) THEN
          ISWAP = 2
          GOTO 99
        ELSE IF(ISWAP.EQ.2) THEN
          ISWAP = 1
          GOTO 99
        END IF
 99   CONTINUE
      IF(IVEC .EQ. 1) ISWAP = 1
      END DO
*
      IF(ISWAP.EQ.1.AND.NVEC.GT.0) THEN
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
        CALL CPVC_PARPAR_68(VEC1,JBKINF,JBLOCKI,IOFFSET3,
     &           IOFFSET2,0)
      END IF
*
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MVCSMD2_PSIO_5B(FAC,FACLUOUT,VEC1,VEC2,
     &                   NVEC,IREW,JOFFS,JTOTBLL,JBKINF,JBLOCKI)
*
* LUOUT = Factor * LUOUT + sum(I=1,nvec)fac(I) LUIN(I)
*
* LUOUT and LUSCR are rewinded and input data
* on these files are lost.
*
* MVCSMD2 is identical MVCSMD, except that
* the input vector on LUOUT can be included
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(1), JBKINF(JBLOCKI)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET1,IOFFSET2,IOFFSET3
C
*
      IOFFSET1= 0
      IOFFSET2= 0
      IOFFSET3= 0
C
      IF(IREW .NE. 0 ) THEN
        IOFFSET1 = 0
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
      END IF
      ISWAP = 0
      IF(FACLUOUT.NE.0.0D0) THEN
      ISWAP = 2
      ELSE
      ISWAP = 1
      END IF

      DO IVEC = 1, NVEC
*. start reading/writing at the right position
        IOFFSET1 = (IVEC - 1) * JTOTBLL
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
*
        IF( IVEC .EQ. 1 ) THEN
          IF(FACLUOUT.NE.0.0D0) THEN
            CALL VECSMD_PSIO_5B(VEC1,VEC2,FAC(1),FACLUOUT,
CSK            CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(1),FACLUOUT,LUIN,LUOUT,LUSCR,
     &           IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)
          ELSE
            CALL SCLVCD_PSIO_5B(FAC(1),VEC1,IOFFSET1,
CSK            CALL SCLVCD_PCIO_3(LUIN,LUSCR,FAC(1),VEC1,IOFFSET1,
     &           IOFFSET2,JBKINF,JBLOCKI,ISWAP)
          END IF
        ELSE
          ONE = 1.0D0
          CALL VECSMD_PSIO_5B(VEC1,VEC2,FAC(IVEC),ONE,
CSK          CALL VECSMD_PCIO_3(VEC1,VEC2,FAC(IVEC),ONE,LUIN,LLUOUT,LLUSCR,
     &         IOFFSET1,IOFFSET2,IOFFSET3,JBKINF,JBLOCKI,ISWAP)

        END IF
*
        IF(ISWAP.EQ.1) THEN
          ISWAP = 2
          GOTO 99
        ELSE IF(ISWAP.EQ.2) THEN
          ISWAP = 1
          GOTO 99
        END IF
 99   CONTINUE
      IF(IVEC .EQ. 1) ISWAP = 1
      END DO
*
      IF(ISWAP.EQ.1.AND.NVEC.GT.0) THEN
        IOFFSET2 = JTOTBLL * JOFFS
        IOFFSET3 = JTOTBLL * JOFFS
        CALL CPVC_PARPAR_67(VEC1,JBKINF,JBLOCKI,IOFFSET3,
     &           IOFFSET2,0)
      END IF
*
      END
*
#else
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE IAMIRRELEVANT
*. dummy routine for normal compilation       
       END 
#endif
