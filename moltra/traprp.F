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
C...   group of the licence holder may not take any part of Dirac,
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
C FILE    : traprp.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRAPRP */
      SUBROUTINE TRAPRP(WORK,KFREE,LFREE)
C***********************************************************************
C
C     Driver routine for property integral transformations
C
C     Written by L. Visscher Jul 1997
C     Revised Feb 17 2000 hjaaj (new parameter list to PRPMAT)
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <mxorb.h>
#include <dcbtra.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbgen.h>
#include <dcbxpr.h>
#include <dcbprl.h>
#include <mxcent.h>
#include <dcbtrp.h>
      LOGICAL TOBE,TRASAM
      DIMENSION WORK(*),NSTR(2,0:2,2),NSTRT(2),KVEC(2,2),KQ(2,2),
     &          IR(2,2),NR(2),NQ(2,2),NQT(2),
     &          NSPC(2,0:2),KQC(2),KE(2,2)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2),ICMOQC(2)
      DIMENSION IFRP2(2)
      CHARACTER*32 ACHAR
C
      CALL QENTER('TRAPRP')
      KFRSAV = KFREE
      CALL TRAHI(IPRTRP,1)
C
C     =======================================      
C     Set up index arrays for active orbitals
C     =======================================      
C
      IF(NOPAIR) THEN
        IC   = 1
        DO I = 1,NFSYM
          NR(I)   = NESH(I)
          IR(1,I) = 1
          IR(2,I) = NESH(I)
        ENDDO
      ELSE
        IC   = 0
        DO I = 1,NFSYM
          NR(I)   =  NORB(I)
          IR(1,I) = -NPSH(I)
          IR(2,I) =  NESH(I)
        ENDDO
      ENDIF
C
      DO I = 1,2
        NSTRT(I) = 0
        DO IFRP = 1,NFSYM
          CALL MEMGET('INTE',KVEC(IFRP,I),NR(IFRP),WORK,KFREE,LFREE)
          NSTR(IFRP,0,I) = 1
          CALL  NUMLST(TRA2_INDSTR(I,IFRP),WORK(KVEC(IFRP,I)),
     &                 NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &                 NSTR(IFRP,0,I),WORK(KFREE),LFREE)
          CALL ORBCNT(WORK(KVEC(IFRP,I)),NSTR(IFRP,0,I),IFRP,
     &              NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          NSTR(IFRP,0,I) = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
          NSTRT(I) = NSTRT(I) +  NSTR(IFRP,0,I)
          CALL MEMREL('TRAONE',WORK,KFRSAV,KVEC(IFRP,I),KFREE,LFREE)
          CALL MEMGET('INTE',KVEC(IFRP,I),NSTR(IFRP,0,I),
     &                WORK,KFREE,LFREE)
        ENDDO
      ENDDO      
C
C     Print section
C     =============
C
      DO I = 1, NFSYM
         CALL TRAPRI(2,I,WORK(KVEC(I,1)),WORK(KVEC(I,2)),DUM,DUM,NSTR)
      ENDDO
C
C     Get all coefficients and eigenvalues
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBEIG,NORBT,WORK,KFREE,LFREE)
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(A)') 'Coefficient file not present'
         CALL QUIT('TRAPAM: No coefficients !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','TRAONE')
         CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),
     &               WORK(KIBEIG),TOTERG,2,WORK,KFREE,LFREE)
      ENDIF
C
C     Transform property matrices to MO basis
C     ---------------------------------------
C
      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
      OPEN (LUMLF1,FORM='UNFORMATTED',FILE='MDPROP')
C
      CALL MEMGET('REAL',KPROPMO,NSTRT(1)*NSTRT(2)*NZ,WORK,KFREE,LFREE)
C
      DO I = 1,NTRPP
         INDXPR = LTRPP(I)
         ISYM  = IPRPSYM(INDXPR)
         IOPSY = JBTOF(ISYM-1,1)
         DO IX = 1, 32
            ACHAR(IX:IX) = '*'
         ENDDO
         ACHAR(25:32) = PRPNAM(INDXPR)(1:8)
         CALL MEMGET('REAL',KPRPMO,N2ORBXQ,WORK,KFREE,LFREE)
         CALL PRPMAT(INDXPR,IOPSY,WORK(KPRPMO),WORK(KCMO),
     &               ICMOQ,NORB,WORK,KFREE,LFREE,IPRTRP)
C
C        Select the set that we need
C
         CALL DZERO(WORK(KPROPMO),NSTRT(1)*NSTRT(2)*NZ)
C
         DO IFRP = 1, NFSYM
           IFRP2(2) = IFRP
           IFRP2(1) = MOD(IFRP+IOPSY,2) + 1
C
C          Calculate offset in gathering matrix
C          Necessary when we have two fermion irreps
C
           IOFF = I2ORBX(IFRP2(1),IFRP2(2))
           JOFF = (IFRP2(2)-1)*NSTRT(1)*NSTR(1,0,2)
     &          + (IFRP2(1)-1)*NSTR(1,0,1)
           CALL SELMES(WORK(KPRPMO+IOFF),IFRP2,WORK(KPROPMO+JOFF),
     &                 NSTRT(1),NSTRT(2),
     &                 WORK(KVEC(IFRP2(1),1)),WORK(KVEC(IFRP2(2),2)),
     &                 NSTR(IFRP2(1),2,1),NSTR(IFRP2(2),2,2),
     &                 NSTR(IFRP2(1),1,1),NSTR(IFRP2(2),1,2))
         ENDDO
         CALL MEMREL('TRAPRP',WORK,KPRPMO,KPRPMO,KFREE,LFREE)
         CALL MEMGET('REAL',KPRCMO,NSTRT(1)*NSTRT(2)*8,WORK,KFREE,LFREE)
         CALL PRPMLF (WORK(KPRCMO),WORK(KPROPMO),IPRTRP,NSTR,
     &                NSTRT(1),NSTRT(2),
     &                ISYM,ACHAR,LUMLF1)
         CALL MEMREL('TRAPRP',WORK,KPRCMO,KPRCMO,KFREE,LFREE)
      ENDDO
      CLOSE(LU1INT,STATUS = 'KEEP')
C
      DO IX = 1, 32
         ACHAR(IX:IX) = '*'
      ENDDO
      ACHAR(25:32) = 'EOFLABEL'
      WRITE(LUMLF1) ACHAR
      CLOSE(LUMLF1,STATUS = 'KEEP')
C
      CALL MEMREL('TRAPRP',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRAPRP')
C
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SELMES(PRPMO,IFRP,PBUF,NDIM1,NDIM2,JVEC1,JVEC2,
     &                  NPVEC1,NPVEC2,NEVEC1,NEVEC2)

C***********************************************************************
C     Pick out a set of matrix elements from PRPMO according to arrays
C     JVEC1 and JVEC2
C
C     Written by L.Visscher 1997
C     Last revision Jul 23 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
      DIMENSION PRPMO(NORBT,NORBT,NZ),
     &          PBUF(NDIM1,NDIM2,NZ),IFRP(2),
     &          JVEC1(*),JVEC2(*)
C
      DO IZ = 1,NZ
C
C       P-P matrix elements
C
        DO J = 1,NPVEC2
          JJ = NPSH(IFRP(2))+1+JVEC2(J)
          DO I = 1,NPVEC1
             II = NPSH(IFRP(1))+1+JVEC1(I)
             PBUF(I,J,IZ) = PRPMO(II,JJ,IZ)
          ENDDO
        ENDDO
C
C       E-P matrix elements
C
        DO J = 1,NPVEC2
          JJ = NPSH(IFRP(2))+1+JVEC2(J)
          DO I = NPVEC1+1,NPVEC1+NEVEC1
             II = NPSH(IFRP(1))+JVEC1(I)
             PBUF(I,J,IZ) = PRPMO(II,JJ,IZ)
          ENDDO
        ENDDO
C
C       P-E matrix elements
C
        DO J = NPVEC2+1,NPVEC2+NEVEC2
          JJ = NPSH(IFRP(2))+JVEC2(J)
          DO I = 1,NPVEC1
             II = NPSH(IFRP(1))+1+JVEC1(I)
             PBUF(I,J,IZ) = PRPMO(II,JJ,IZ)
          ENDDO
        ENDDO
C
C       E-E matrix elements
C
        DO J = NPVEC2+1,NPVEC2+NEVEC2
          JJ = NPSH(IFRP(2))+JVEC2(J)
          DO I = NPVEC1+1,NPVEC1+NEVEC1
             II = NPSH(IFRP(1))+JVEC1(I)
             PBUF(I,J,IZ) = PRPMO(II,JJ,IZ)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C
C     ==============================================================
      SUBROUTINE PRPMLF (PRCMO,PROPMO,IPRINT,NSTR,NSTRT1,NSTRT2,ISYM,
     &                   ACHAR,LUMLF1)
C     Express property integrals in complex form and write them to
C     file
C     Instead of array NSTRT now NSTRT1 and NSTRT2 are parameters
C     ==============================================================
C
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
      PARAMETER (DM1=-1.D0)
      CHARACTER*32 ACHAR
      DIMENSION NSTR(2,0:2,2),PRCMO(NSTRT1*2,NSTRT2*2,2)
      DIMENSION PROPMO(NSTRT1,NSTRT2,NZ)
      DIMENSION IOFF1(2,2),IOFF2(2,2),IOFF3(2,2)
C
      IF (IPRINT.GE.3) THEN
         WRITE(LUPRI,'(/3X,3A/)') '*** Property matrix :',ACHAR,' ***'
         CALL PRQMAT(PROPMO,NSTRT1,NSTRT2,NSTRT1,NSTRT2,
     &               NZ,IPQTOQ(1,ISYM-1),LUPRI)
      ENDIF
C
C     Fill in complete property matrix in complex form
C
      IOPSY = JBTOF(ISYM-1,1)
C
      IOFF1(1,1) = 0
      IOFF1(1,2) = 0
      IOFF2(1,1) = 0
      IOFF2(1,2) = 0
      IOFF3(1,1) = NSTR(1,0,1)
      IOFF3(1,2) = NSTR(1,0,2)
      DO IFRP = 2, NFSYM
         IOFF1(IFRP,1) = IOFF1(IFRP-1,1) + NSTR(IFRP-1,0,1)
         IOFF1(IFRP,2) = IOFF1(IFRP-1,2) + NSTR(IFRP-1,0,2)
         IOFF2(IFRP,1) = IOFF2(IFRP-1,1) + 2*NSTR(IFRP-1,0,1)
         IOFF2(IFRP,2) = IOFF2(IFRP-1,2) + 2*NSTR(IFRP-1,0,2)
         IOFF3(IFRP,1) = IOFF2(IFRP,1) + NSTR(IFRP,0,1)
         IOFF3(IFRP,2) = IOFF2(IFRP,2) + NSTR(IFRP,0,2)
      ENDDO
C
      CALL DZERO(PRCMO,NSTRT1*NSTRT2*8)
C
      DO IZ = 1, NZ
        IQ = IPQTOQ(IZ,ISYM-1)
        DO JFRP = 1, NFSYM
          IFRP  = MOD(JFRP+IOPSY,2) + 1
          DO J = 1, NSTR(JFRP,0,2)
            IF (IQ.EQ.1) THEN
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF2(IFRP,1),J+IOFF2(JFRP,2),1),1)
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF3(IFRP,1),J+IOFF3(JFRP,2),1),1)
            ELSEIF (IQ.EQ.2) THEN
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF2(IFRP,1),J+IOFF2(JFRP,2),2),1)
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF3(IFRP,1),J+IOFF3(JFRP,2),2),1)
              CALL DSCAL(NSTR(IFRP,0,1),DM1,
     &                   PRCMO (1+IOFF3(IFRP,1),J+IOFF3(JFRP,2),2),1)
            ELSEIF (IQ.EQ.3) THEN
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF2(IFRP,1),J+IOFF3(JFRP,2),1),1)
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF3(IFRP,1),J+IOFF2(JFRP,2),1),1)
              CALL DSCAL(NSTR(IFRP,0,1),DM1,
     &                   PRCMO (1+IOFF3(IFRP,1),J+IOFF2(JFRP,2),1),1)
            ELSEIF (IQ.EQ.4) THEN
              CALL DCOPY(NSTR(IFRP,0,1),
     &                 PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                 PRCMO (1+IOFF2(IFRP,1),J+IOFF3(JFRP,2),2),1)
              CALL DCOPY(NSTR(IFRP,0,1),
     &                   PROPMO(1+IOFF1(IFRP,1),J+IOFF1(JFRP,2),IZ),1,
     &                   PRCMO (1+IOFF3(IFRP,1),J+IOFF2(JFRP,2),2),1)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
      IF (IPRINT.GE.1) THEN
         IF (NSTRT1.EQ.NSTRT2) THEN
            WRITE (LUPRI,1010) ACHAR(25:32)
            DO I = 1, 2*NSTRT1
               WRITE (LUPRI,1011) I,PRCMO(I,I,1),PRCMO(I,I,2)
            ENDDO
            WRITE (LUPRI,*)
         ENDIF
         IF (IPRINT.GE.3) THEN
            WRITE (LUPRI,1012)
            DO I = 1, NSTRT1*2
               DO J = 1, NSTRT2*2
                  WRITE (LUPRI,1013) I,J,PRCMO(I,J,1),PRCMO(I,J,2)
               ENDDO
            ENDDO
            WRITE (6,*)
         ENDIF
      ENDIF
C
C     WRITE PRCMO TO FILE
C
      WRITE (LUMLF1) ACHAR
      WRITE (LUMLF1) ((PRCMO(I,J,1),PRCMO(I,J,2),
     &                I=1,NSTRT1*2),J=1,NSTRT2*2)
C
 1000 FORMAT (/' Transforming ',A8,' evaluation date ',A8,' label ',A8)
 1001 FORMAT (/' Finished transformation : integrals on ',A8)
 1010 FORMAT (/' Diagonal elements of ',A8,
     & ' operator in molecular spinor basis'
     & //12X,'Real part',11X,'Imaginary Part')
 1011 FORMAT (I6,2X,F18.10,2X,F18.10)
 1012 FORMAT (/' Complete property matrix in molecular spinor basis')
 1013 FORMAT (2I5,F18.10,2X,F18.10)
      RETURN
      END
