!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getmat */
      SUBROUTINE GETMAT(COPT,DMAT,IPRINT,FKRMC,WORK,LWORK)
C***********************************************************************
C
C     Read DMAT from file 
C     - KRMC MATRICES if doing linear response
C     - KRMC FOCK if solving KRMCSCF linear equations
C
C     IOPT:
C
C     'FCMO' - Read inactive Fock matrix in MO basis (FCmo)
C     'FVMO' - Read active Fock matrix in MO basis (FVmo)
C     'DVMO' - Read 1-electron active density matrix (DVmo)
C     'FGEN_PV' - generalized Fock matrix in MO basis
C
C
C     Written by J. Thyssen - Nov 22 2000
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use dgroup_mod, only: IPQTOQ, NZ
!      use dcborb_mod, only: N2ASHXQ, N2ORBXQ, NASHT, NNASHX, NORBT
      use dcbbas_mod, only: N2BBASXQ
#include "implicit.h"
      INTEGER :: I, J, K
C
#include "dgroup.h"
#include "dcborb.h"
C
      CHARACTER COPT*(*), CCOPT*8,LABEL*8, FILENAME*20
      DIMENSION DMAT(*),WORK(*)
      LOGICAL   FEXIST, FOPEN, FNDLAB, FKRMC
C
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK
      CALL QENTER('GETMAT')
C
      IF ( FKRMC ) THEN
         FILENAME = 'KRMC_FOCK'
         IUNIT = 53
      ELSE
         FILENAME = 'KRMC_MATRICES'
         IUNIT = 40
      END IF
C
      IPOSOLD = 1
      IOFF = 1
      LCOPT = LEN(COPT)
      IF (LCOPT .EQ. 0) GOTO 1000
C
C     Check if file exists
C
      INQUIRE(FILE=FILENAME,EXIST=FEXIST)
      IF (FEXIST) THEN
         INQUIRE(FILE=FILENAME,OPENED=FOPEN)
         IF (.NOT. FOPEN) THEN
            OPEN(IUNIT,FILE=FILENAME,STATUS='UNKNOWN',
     &           FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         END IF
      ELSE
         FOPEN = .FALSE.
         OPEN(IUNIT,FILE=FILENAME,STATUS='NEW',
     &        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         REWIND(IUNIT)
         CALL NEWLAB('SODLABEL',IUNIT,LUPRI)
      END IF
C
  100 CONTINUE
C 
         IPOS = INDEX(COPT(IPOSOLD:LCOPT),'+')
         IF (IPOS .EQ. 0) IPOS = LEN(COPT)+1
         CCOPT = COPT(IPOSOLD:IPOS-1)
         IPOSOLD = IPOS + 1
         IF (IPRINT .GE. 5)
     &      WRITE(LUPRI,'(/3A)') 'GETMAT, option <',CCOPT,'>'
C
         IF (CCOPT .EQ. 'FCMO' .OR. CCOPT .EQ. 'FVMO') THEN
C
            IF (CCOPT .EQ. 'FCMO') THEN
               LABEL = 'FCMO    '
               IOFF_FMAT = 0
            ELSE
               LABEL = 'FVMO    '
               IOFF_FMAT = N2ORBXQ
            END IF
C
C           Read FC Fock matrix in MO basis
C
            REWIND(IUNIT)
            IF (.NOT. FNDLAB(LABEL,IUNIT)) THEN
C
C              Matrix not present in file
C
C              Generate matrix
C
               NFMAT = 1
               IF (NASHT .GT. 0) NFMAT = NFMAT + 1
               CALL MEMGET('REAL',KFMAT,NFMAT*N2BBASXQ,WORK,KFREE,LFREE)
               CALL GETFCK(WORK(KFMAT),IPRINT,WORK,KFREE,LFREE)
C
C              Write it to file
C
               CALL WRTKRMC(IUNIT,LABEL,WORK(KFMAT+IOFF_FMAT),N2ORBXQ)
               CALL DCOPY(N2ORBXQ,WORK(KFMAT+IOFF_FMAT),1,DMAT(IOFF),1)
               CALL MEMREL('GETMAT.FCMO',WORK,KWORK,KWORK,KFREE,LFREE)
            ELSE
C
C              ...or matrix present in file
C
               CALL READT(IUNIT,N2ORBXQ,DMAT(IOFF))
            END IF
            IF (IPRINT .GE. 10) THEN
               CALL HEADER('Output from GETMAT: '//LABEL,-1)
               CALL PRQMAT(DMAT(IOFF),NORBT,NORBT,
     &                        NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
            END IF
            IOFF = IOFF + N2ORBXQ
         ELSE IF (CCOPT .EQ. 'DVMO') THEN
C
C           Read DV density matrix in MO basis
C
            LABEL = 'DVMO    '
            REWIND(IUNIT)
            IF (.NOT. FNDLAB(LABEL,IUNIT)) THEN
C
C              Matrix not present in file
C
C              Generate matrix
C
               CALL MEMGET('REAL',KDV,N2ASHXQ,WORK,KFREE,LFREE)
               CALL RMAKDM(.TRUE.,.FALSE.,DUMMY,WORK(KDV),DUMMY,DUMMY,
     &              WORK,KFREE,LFREE)
c              RMAKDM needs to be called "...,WORK,KFREE,LFREE)"
!              WRITE(LUPRI,'(/A)') '*** GETMAT: one-electron active'//
!    &         ' matrix DV missing. Cannot continue. ***'
!              CALL QUIT('*** GETMAT: one-electron active matrix DV 
!    &                    missing. Cannot continue. ***')
C
C              Write it to file
C
               CALL WRTKRMC(IUNIT,LABEL,WORK(KDV),N2ASHXQ)
               CALL DCOPY(N2ASHXQ,WORK(KDV),1,DMAT(IOFF),1)
               CALL MEMREL('GETMAT.DVMO',WORK,KWORK,KWORK,KFREE,LFREE)
            ELSE
C
C              ...or matrix present in file
C
               CALL READT(IUNIT,N2ASHXQ,DMAT(IOFF))
            END IF
C
C           Print section
C
            IF (IPRINT .GE. 10) THEN
               CALL HEADER('Output from GETMAT: DVMO',-1)
               CALL PRQMAT(DMAT,NASHT,NASHT,NASHT,NASHT,
     &              NZ,IPQTOQ(1,0),LUPRI)
            END IF
            IOFF = IOFF + N2ASHXQ
         ELSE IF (CCOPT .EQ. 'FGENPV') THEN
C
C           Read generalized Fock matrix in MO basis
C
            LABEL = 'FGENPV  '
            REWIND(IUNIT)
            IF (.NOT. FNDLAB(LABEL,IUNIT)) THEN
C
C              Matrix not present in file
C
               CALL QUIT('*** Not implemented ***')
C
            ELSE
C
C              ...or matrix present in file
C
               CALL READT(IUNIT,NORBT * NASHT * NZ, DMAT(IOFF))
C
            END IF
            IF (IPRINT .GE. 10) THEN
               CALL HEADER('Output from GETMAT: '//LABEL,-1)
               CALL PRQMAT(DMAT(IOFF),NORBT,NASHT,
     &                        NORBT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
            END IF
            IOFF = IOFF + NASHT * NORBT * NZ
         ELSE IF (CCOPT .EQ. 'FCAC') THEN
C
C           Read active-active part of FC Fock matrix
C
            LABEL = 'FCAC    '
            REWIND(IUNIT)
            IF (.NOT. FNDLAB(LABEL,IUNIT)) THEN
C
C              Matrix not present in file
C
               CALL QUIT('*** Not implemented ***')
C
            ELSE
C
C              ...or matrix present in file
C
               CALL READT(IUNIT,N2ASHXQ, DMAT(IOFF))
C
            END IF
            IF (IPRINT .GE. 10) THEN
               CALL HEADER('Output from GETMAT: '//LABEL,-1)
               CALL PRQMAT(DMAT(IOFF),NASHT,NASHT,
     &                        NASHT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
            END IF
            IOFF = IOFF + N2ASHXQ
         ELSE IF (CCOPT .EQ. 'H2AC' .OR. CCOPT .EQ. 'PVMO') THEN
C
C           Read two-electron integrals or density matrix
C
            IF (CCOPT .EQ. 'H2AC') THEN
               LABEL = 'H2AC    '
            ELSE
               LABEL = 'PVMO    '
            END IF
            REWIND(IUNIT)
            IF (.NOT. FNDLAB(LABEL,IUNIT)) THEN
C
C              Matrix not present in file
C
               CALL QUIT('*** Not implemented ***')
C
            ELSE
C
C              ...or matrix present in file
C
               CALL READT(IUNIT,
     &              NASHT * NASHT * NNASHX * NZ * 3, DMAT(IOFF))
C
            END IF
C
 9000 FORMAT(1P,1X,4I3,4X,3(D14.7,1X),3(/,16X,3(D14.7,1X)))
C
            IF (IPRINT .GE. 10) THEN
               CALL HEADER('Output from GETMAT: '//LABEL,-1)
               IJ = 0
               DO I = 1,NASHT
                  DO J = 1,I
                     IJ = IJ + 1
                     DO K = 1,NASHT
                        DO L = 1,NASHT
                           WRITE(LUPRI,9000)
     &                          K,L,I,J,
     &                          ((DMAT(IOFF 
     &                          + (ICL-1) * NASHT * NASHT * NNASHX * NZ
     &                          + (IZ-1) * NASHT * NASHT * NNASHX
     &                          + (IJ-1) * NASHT * NASHT
     &                          + (L-1) * NASHT
     &                          + (K-1)),
     &                          IZ = 1,NZ),ICL=1,3)
                        END DO
                     END DO
                  END DO
               END DO
            END IF
            IOFF = IOFF + NASHT * NASHT * NNASHX * NZ * 3
         ELSE
            WRITE(LUPRI,'(3A)') 'GETMAT: Unknown option: <',CCOPT,'>'
            CALL QUIT('argh, ralle ralle in GETMAT')
         END IF
C
         IF (IPOS .LT. LCOPT) GOTO 100
C  
 1000 CONTINUE
      IF (.NOT. FOPEN) CLOSE(IUNIT)
      CALL QEXIT('GETMAT')
      RETURN
      END
