C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dirac, a relativistic ab initio electronic structure program,
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : pamfck.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getmat */
      SUBROUTINE GETMAT(COPT,DMAT,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Read DMAT from file KRMC MATRICES
C
C     IOPT:
C
C     'FCMO' - Read inactive Fock matrix in MO basis (FCmo)
C     'FVMO' - Read active Fock matrix in MO basis (FVmo)
C     'DVMO' - Read 1-electron active density matrix (DVmo)
C
C
C     Written by  Hans Joergen Aa. Jensen and T.Saue 27-Jan-1995
C     Last revision: July 21 1997 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcborb.h>
#include <dcbbas.h>
C
      CHARACTER COPT*(*), CCOPT*8,LABEL*8
      DIMENSION DMAT(*),WORK(*)
      LOGICAL   FEXIST, FOPEN, FNDLAB
C
#include <memint.h>
      CALL QENTER('GETMAT')
C
      IUNIT = 40
      IPOSOLD = 1
      IOFF = 1
      LCOPT = LEN(COPT)
      N2ASHTQ = N2ASHT*NZ
      IF (LCOPT .EQ. 0) GOTO 1000
C
C     Check if file exists
C
      INQUIRE(FILE='KRMC MATRICES',EXIST=FEXIST)
      IF (FEXIST) THEN
         INQUIRE(FILE='KRMC MATRICES',OPENED=FOPEN)
         IF (.NOT. FOPEN) THEN
            OPEN(IUNIT,FILE='KRMC MATRICES',STATUS='UNKNOWN',
     &           FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         END IF
      ELSE
         FOPEN = .FALSE.
         OPEN(IUNIT,FILE='KRMC MATRICES',STATUS='NEW',
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
               CALL MEMGET('REAL',KDV,N2ASHTQ,WORK,KFREE,LFREE)
               CALL RMAKDM(DUMMY,WORK(KDV),DUMMY,DUMMY,WORK,LWORK)
C
C              Write it to file
C
               CALL WRTKRMC(IUNIT,LABEL,WORK(KDV),N2ASHTQ)
               CALL DCOPY(N2ASHTQ,WORK(KDV),1,DMAT(IOFF),1)
               CALL MEMREL('GETMAT.DVMO',WORK,KWORK,KWORK,KFREE,LFREE)
            ELSE
C
C              ...or matrix present in file
C
               CALL READT(IUNIT,N2ASHTQ,DMAT(IOFF))
            END IF
C
C           Print section
C
            IF (IPRINT .GE. 10) THEN
               CALL HEADER('Output from GETMAT: DVMO',-1)
               DO I = 1,NFSYM
                  CALL PRQMAT(DMAT(IOFF+I2ASHT(I)),NASH(I),NASH(I),
     &                        NASH(I),NASH(I),NZ,IPQTOQ(1,0),LUPRI)
               END DO
            END IF
            IOFF = IOFF + N2ASHTQ
         ELSE
            WRITE(LUPRI,'(/3A)') 'GETMAT: Unknown option: <',CCOPT,'>'
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrtkrmc */
      SUBROUTINE WRTKRMC(LUNIT,LABEL,VECTOR,ISIZE)
C*****************************************************************************
C
C     Write VECTOR to LUKRMC under label LABEL.
C
C     Input:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Output:
C
C     Written by J. Thyssen - Nov 24 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <iratdef.h>
C
      DIMENSION VECTOR(ISIZE)
      CHARACTER LABEL*8,C(4)*8
      LOGICAL   FNDLAB
C
      REWIND(LUNIT)
      IF (FNDLAB(LABEL,LUNIT)) THEN
C
C     Label exists.
C     Seek to label SODLABEL and write from there.
C
         REWIND(LUNIT)
         IF (.NOT. FNDLAB('SODLABEL',LUNIT)) THEN
            CALL QUIT('WRTKRMC: wrong structure in file')
         END IF
      ELSE
C
C     Label do not exist
C     Seek to label EOFLABEL and write from there
C
         REWIND(LUNIT)
         IF (FNDLAB('EOFLABEL',LUNIT))
     &      BACKSPACE LUNIT
C
      END IF
C
C     Write label
C
      CALL NEWLAB(LABEL,LUNIT,LUPRI)
C
C     Write vector
C
      LEN = MAX(4,ISIZE)
      CALL WRITSQ(LUNIT,IRAT*LEN,VECTOR)
      CALL NEWLAB('EOFLABEL',LUNIT,LUPRI)
C
      REWIND LUNIT
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rmakdm */
      SUBROUTINE RMAKDM(CREF,DV,PV,INDXCI,WORK,LWORK)
C*****************************************************************************
C
C     Construct one- and two-electron active density matrix
C
C     Input:
C        CMO     - the orbitals for the CEP
C        CREF    - CI vector
C
C     Output:
C        DV      - one electron active density matrix
C        PV      - two electron active density matrix
C
C     Written by J. Thyssen - Oct 22 1998
C     Last revision :
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D00, D0 = 0.0D00)
C
#include <dcborb.h>
#include <dcbgen.h>
#include <dgroup.h>
C
      DIMENSION CREF(*),DV(*),PV(*),INDXCI(*)
      DIMENSION WORK(LWORK)
C
C
#include <memint.h>
C
      CALL QENTER('RMAKDM')
C
      IF (NASHT .EQ. 1) THEN
         DV(1) = D1
         DO IZ = 2,NZ
            DV(IZ) = D0
         END DO
      ELSE
         CALL QUIT('RMAKDM: MASHT .GT. 1 not implemented yet')
      END IF
C
C
      CALL QEXIT('RMAKDM')
      RETURN
      END
