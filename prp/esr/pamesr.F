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
C FILE    : pamesr.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ESRINP */
      SUBROUTINE ESRINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for ESR properties
C
C     Written by A.Noerager - Mar 2000
C     Last revision: 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <mxcent.h>
#include <dcbprp.h>
#include <dcbesr.h>
#include <orgcom.h>
      PARAMETER (NTABLE = 8)
      LOGICAL SET, NEWDEF, RESET, GORG
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.GAUGEO','.PRINT ','.OPERAT','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
C
C
C     Initialize /CBIESR/
C     ===================
C
      IPRESR  = 0
      NESRP   = 0
      NCIESR  = -1
C
C     Initialize orgcom.h
C     ===================
C
      GORG = .FALSE.
      CALL DZERO(GAGORG,3)
C
C     Process input 
C     =========================
C
      NEWDEF = (WORD .EQ. '*ESR')
      ICHANG = 0
      INPERR = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
  110       CONTINUE
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in ESRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in ESRINP.')
    1          CONTINUE
C Read user specified gauge origin
                  READ(LUCMD,*) (GAGORG(ICOOR), ICOOR = 1,3)
                  GORG = .TRUE.
               GO TO 100
    2          CONTINUE
C&&&& PRINT:  Print level
                  READ(LUCMD,*) IPRESR
               GO TO 100
    3          CONTINUE
C&&&& OPERATOR: Define operator
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                        WORK,LWORK,IPRESR)
                  CALL OP1IND('ESRINP',INDESR,LESRP,NESRP,INDXPR,
     &                        MAXESR)
                  DOESR = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    4          CONTINUE
               GO TO 100
    5          CONTINUE
               GO TO 100
    6          CONTINUE
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ESRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in ESRINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      IF (GORG) THEN
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(1X,A)') 'ESRINP: '
         CALL PRSYMB(LUPRI,'=',75,0)
         IF (GORG) THEN
            WRITE(LUPRI,'(/,A,3F12.6,/)') 'User defined gauge '//
     &           'origin used:',(GAGORG(I), I = 1, 3)
         END IF
         CALL DCOPY(3,GAGORG,1,DIPORG,1)
      END IF
  999 CONTINUE
      IF (INPERR.GT.0) CALL QUIT('Input error in *ESR')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpesr */
      SUBROUTINE PRPESR(ESRVAL,WORK,LWORK)
C*****************************************************************************
C
C     Calculate expectation values defined in /CBIESR/
C
C     Written by Anette Noerager - Mar 2000
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C     for PRQMAT print of complex matrices:
      PARAMETER (NZCMPLX = 2)
      DIMENSION IQINDX(NZCMPLX)
      SAVE IQINDX
      DATA IQINDX /1,2/
C
      DIMENSION ESRVAL(NCIESR,NCIESR,NZCMPLX,*),WORK(LWORK)
      DIMENSION GTENSOR(3,3),GMATRIX(3,3),GEIG(3),GSHIFT(3),GEVEC(3,3),
     &          GCOR(3)
C
      PARAMETER (D3 = 3.0D0, D4 = 4.0D0)
C
#include <dcbxpr.h>
#include <dcbesr.h>
#include <gfac.h>
#include <dcbprp.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
C
      CALL QENTER('PRPESR')
#include <memint.h>
C
C     Memory allocation
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KESRST ,2*NESRP   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT  ,NZ*N2BBASX*(NOPEN+1),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT  ,NZ*N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPTRI  ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN ,   NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST ,NZ        ,WORK,KFREE,LFREE)
C
      IF (IPRESR .GT. 1) THEN
        WRITE(LUPRI,'(/A)')'PRPESR called with the following operators:'
        DO I = 1,NESRP
          INDXPR = LESRP(I)
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Operator no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      END IF
C
      CALL PRPESR1(ESRVAL,WORK(KESRST),WORK(KCMO),WORK(KEIG),
     &            WORK(KDMAT),WORK(KPMAT),WORK(KPTRI),
     &            WORK(KONEIN),WORK(KFIRST),WORK(KFREE),
     &            LFREE)
C
      IF (IPRESR .GT. 1) THEN
        WRITE(LUPRI,'(/A)') 'PRPESR results (test output):'
        DO I = 1,NESRP
          INDXPR = LESRP(I)
          CALL PRSYMB(LUPRI,'.',75,0)
          WRITE(LUPRI,'(4X,A,I3,A,3X,A16)')
     &       'Operator no.',I,':',PRPNAM(INDXPR)
          CALL PRSYMB(LUPRI,'.',75,0)
          CALL WRIXPR(I,INDXPR)
cq        WRITE(LUPRI,'(4X,A,4F15.6)')
cq   &       'Quaternion value:',(ESRVAL(J,I),J=1,4)
          CALL PRQMAT(ESRVAL(1,1,1,I),NCIESR,NCIESR,
     &                NCIESR,NCIESR,NZCMPLX,IQINDX,LUPRI)
        ENDDO
      END IF
C
      IF (ESRGTENS) THEN
C
C     (1) Get the Dirac equation g tensor (corresponding to GFAC = 2.0)
C
      DO K = 1,3
C        ... for B field in direction k (x,y,z)
         I = IPESR(MXCOOR+K,1)
         GTENSOR(K,1) =  D4*ESRVAL(1,2,1,I)
         GTENSOR(K,2) = -D4*ESRVAL(1,2,2,I)
         GTENSOR(K,3) =  D4*ESRVAL(1,1,1,I)
      ENDDO

      WRITE (LUPRI,'(/A/,3(/5X,3F18.8))')
     &   ' The g-tensor for H(Zeeman) = 1/2c B g Seff :',
     &   ((GTENSOR(I,K),K=1,3),I=1,3)

C     Get G-matrix = gtensor*gtensor(transposed)
C
      DO I = 1,3
        DO K = 1,3
          GMATRIX(I,K) = GTENSOR(I,1)*GTENSOR(K,1)
     &                 + GTENSOR(I,2)*GTENSOR(K,2)
     &                 + GTENSOR(I,3)*GTENSOR(K,3)
        ENDDO
      ENDDO


      WRITE (LUPRI,'(/A/,3(/5X,3F18.8))')
     &   ' The G-matrix = g * g^T :', ((GMATRIX(I,K),K=1,3),I=1,3)

      CALL RSJACO(3,3,3,GMATRIX,GEIG,1,1,0,GEVEC,WORK(KFREE),LFREE)

      DO I = 1,3
        GEIG(I) = SQRT(GEIG(I))
        GSHIFT(I) = (GEIG(I) - 2.0D0)*1.D6
      ENDDO
      GISO = (GEIG(1) + GEIG(2) + GEIG(3)) / D3
      GSHFTISO = (GSHIFT(1) + GSHIFT(2) + GSHIFT(3)) / D3
      WRITE (LUPRI,'(/A//A//5X,4F18.8//A//5X,4F18.1//A/,3(/5X,3F18.8))')
     &   '* The Dirac equation g values corresponding to g_e = 2',
     &   ' The g values (x, y, z, iso):',
     &   (GEIG(I),I =1,3), GISO,
     &   ' The g shifts in ppm (x, y, z, iso):',
     &   (GSHIFT(I),I =1,3), GSHFTISO,
     &   ' The corresponding eigenvectors:',
     &   ((GEVEC(I,K),K=1,3),I=1,3)
C
C     (2) Get the corrected g values by scaling
C
      DO I = 1,3
        GCOR(I) = (GEIG(I)/2)*GFAC
        GSHIFT(I) = (GCOR(I) - GFAC)*1.D6
      ENDDO
      GISO = (GCOR(1) + GCOR(2) + GCOR(3)) / D3
      GSHFTISO = (GSHIFT(1) + GSHIFT(2) + GSHIFT(3)) / D3
      WRITE (LUPRI,'(/A,F13.10,A//A//5X,4F18.8//A//5X,4F18.1))')
     &   '* The corrected g values scaled with',GFAC,' / 2',
     &   ' The corrected g values (x, y, z, iso):',
     &   (GCOR(I),I =1,3), GISO,
     &   ' The corrected g shifts in ppm (x, y, z, iso):',
     &   (GSHIFT(I),I =1,3), GSHFTISO
C
C     (3) Get the corrected g values by adding correction term
C
      DO K = 1,3
C        ... for B field in direction k (x,y,z)
         I = IPESR(MXCOOR+K,1)
         J = IPESR(MXCOOR+K,2)
         GTENSOR(K,1) =  D4*(ESRVAL(1,2,1,I)+ESRVAL(1,2,1,J))
         GTENSOR(K,2) = -D4*(ESRVAL(1,2,2,I)+ESRVAL(1,2,2,J))
         GTENSOR(K,3) =  D4*(ESRVAL(1,1,1,I)+ESRVAL(1,1,1,J))
      ENDDO

      WRITE (LUPRI,'(/A/,3(/5X,3F18.8))')
     &   ' The g_1 corrected g-tensor for H(Zeeman) = 1/2c B g Seff :',
     &   ((GTENSOR(I,K),K=1,3),I=1,3)

C     Get G-matrix = gtensor*gtensor(transposed)
C
      DO I = 1,3
        DO K = 1,3
          GMATRIX(I,K) = GTENSOR(I,1)*GTENSOR(K,1)
     &                 + GTENSOR(I,2)*GTENSOR(K,2)
     &                 + GTENSOR(I,3)*GTENSOR(K,3)
        ENDDO
      ENDDO


      WRITE (LUPRI,'(/A/,3(/5X,3F18.8))')
     &   ' The G-matrix = g * g^T :', ((GMATRIX(I,K),K=1,3),I=1,3)

      CALL RSJACO(3,3,3,GMATRIX,GEIG,1,1,0,GEVEC,WORK(KFREE),LFREE)

      DO I = 1,3
        GEIG(I) = SQRT(GEIG(I))
        GSHIFT(I) = (GEIG(I) - GFAC)*1.D6
      ENDDO
      GISO = (GEIG(1) + GEIG(2) + GEIG(3)) / D3
      GSHFTISO = (GSHIFT(1) + GSHIFT(2) + GSHIFT(3)) / D3
      WRITE (LUPRI,
     &   '(/A,F13.10//A//5X,4F18.8//A//5X,4F18.1//A/,3(/5X,3F18.8))')
     &   '* The g_1 corrected Dirac eq. g values corresponding to g_e ='
     &   ,GFAC,
     &   ' The g values (x, y, z, iso):',
     &   (GEIG(I),I =1,3), GISO,
     &   ' The g shifts in ppm (x, y, z, iso):',
     &   (GSHIFT(I),I =1,3), GSHFTISO,
     &   ' The corresponding eigenvectors:',
     &   ((GEVEC(I,K),K=1,3),I=1,3)

      ENDIF
C     .. end if (ESRGTENS) then
C
C     Print ESR anisotropies from hyperfine Hamiltonian
C
      IF (ESRANISO)
     &   CALL PRP_ESR(ESRVAL,WORK,KFREE,LFREE)
C
C     Memory deallocation
C
      CALL MEMREL('PRPESR',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('PRPESR')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpesr1 */
      SUBROUTINE PRPESR1(ESRVAL,LESRST,CMO,EIG,DMAT,PMAT,PTRI,OP1INT,
     &                  FIRST,WORK,LWORK)
C*****************************************************************************
C
C     Calculate expectation values defined in /CBIESR/
C
C     A. Noerager - Mar 2000
C     Based on prpex1 by Trond Saue May 27 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER(D0 = 0.0D0)
C
      LOGICAL FIRST(NZ),LESRST(2,*),FILEX
      DIMENSION ESRVAL(NCIESR,NCIESR,2,*),CMO(*),DMAT(*),PMAT(*),
     &          PTRI(*),OP1INT(*),EIG(*),WORK(LWORK)
      CHARACTER MXFORM*6,FMT*6
      INTEGER   GETIOPEN, IMO(2), NMO(2)
C
#include <dgroup.h>
ctec dcborb included
#include <dcborb.h>
#include <dcbbas.h>
#include <dcbesr.h>
#include <dcbgen.h>
#include <dcbxpr.h>
#include <dcbprp.h>
#include <dcbdhf.h>
C
C
#include <memint.h>
#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
C
C     Initialize
C     ==========
C
      CALL DZERO(ESRVAL,NCIESR*NCIESR*2*NESRP)
      N2ESR = 2*NESRP
      CALL LSET(N2ESR,.FALSE.,LESRST)
C
C     Generate density matrix
C     =======================
C
      INQUIRE(FILE='DFCOEF',EXIST=FILEX)
      IF(.NOT.FILEX) GOTO 1000
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRPESR')
      CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IDUM,TOTERG,1,
     &            WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Expectation values normal procedure
C    ************************************
C
c       CALL GENDEN(DMAT,CMO,1,IPRESR)
c       DO IOPEN = 1,NOPEN
c          CALL ACTDEN(DMAT(1+N2BBASXQ*IOPEN),CMO,IOPEN,IPRESR)
c       END DO
c       CALL TITLER('Expectation values','*',127)
c       WRITE(LUPRI,'(51X,2(A2,1X))') 's0','t0'
c       CALL PRSYMB(LUPRI,'-',76,0)
C
        DO I = 1,NESRP
          INDXPR = LESRP(I)
C
C         Extract symmetry information
C         =============================
C
          ISYM = IPRPSYM(INDXPR)
          IREP = ISYM - 1
          ITIM = IPRPTIM(INDXPR)
          IOPSY = JBTOF(ISYM-1,1)
          IF (IPRESR.GT.1 .OR. ITIM.GE.0) WRITE(LUPRI,'(A,I3,A,A16)')
     &    '* PRPESR1: calculating property no.',I,' : ',PRPNAM(INDXPR)
C
          IF (ITIM .GE. 0) CALL QUIT(
     &    'PRPESR1 only implemented for time-antisymmetric operators')
C

Chj-start test output
          kfree = 1
          lfree = lwork
          iprtst = 10
          imo(1) = icmoq(1) + nfbas(1,0)*(npsh(1)+nish(1))
          nmo(1) = nash(1)
          imo(2) = icmoq(2) + nfbas(2,0)*(npsh(2)+nish(2))
          nmo(2) = nash(2)
c     imo(1) = icmoq(1)
c     nmo(1) = norb(1)
c     imo(2) = icmoq(2)
c     nmo(2) = norb(2)
c
c      print *, ' ..... test output of property matrix in MO basis ...'
c      print *, 'isym,iopsy:',isym,iopsy
c      print *, 'icmoq  :',icmoq
c      print *, 'nfbas  :',nfbas(1,0),nfbas(2,0)
c      print *, 'npsh   :',npsh
c      print *, 'nish   :',nish
c      print *, 'imo,nmo:',imo,nmo
C
          CALL PRPMAT(INDXPR,IOPSY,pmat,CMO,imo,nmo,
     &                WORK,KFREE,LFREE,IPRESR)
C
          IF ( MOD(NAELEC,2) .EQ. 1) THEN
C         ... odd number of electrons
             CALL PRPESRODD(NMO,ESRVAL(1,1,1,I),IREP,PMAT)
          ELSE
C         ... even number of electrons
             CALL PRPESREVEN(NMO,ESRVAL(1,1,1,I),IREP,PMAT)
          END IF
c      print *, ' ..... end of test output of MO prp. matrix ...'
Chj-end test output

c         CALL PRPESR2(LU1INT,INDXPR,ESRVAL(I),LESRST(1,I),
c    &        DMAT,PMAT,PTRI,OP1INT,FIRST,WORK,LWORK,IPRESR)
c         FMT = MXFORM(ESRVAL(I),15)
c         WRITE(LUPRI,'(4X,A16,A3,3X,'//FMT//',1X,A4,3X,2(3X,L1))')   
c    &       PRPNAM(INDXPR),' : ',ESRVAL(I),'a.u.',
c    &       LESRST(1,I),LESRST(2,I)
c         CALL FLSHFO(LUPRI)
        ENDDO
c       CALL PRSYMB(LUPRI,'-',76,0)
c       WRITE(LUPRI,'(4X,A)')
c    &   's0 = T : Expectation value zero by point group symmetry.'
c       WRITE(LUPRI,'(4X,A)')
c    &   't0 = T : Expectation value zero by time reversal symmetry.'
C
C
      CLOSE(LU1INT,STATUS = 'KEEP')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'PRPESR: Coefficient file not found !'
      CALL QUIT('PRPESR: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpesrodd */
      SUBROUTINE PRPESRODD(NMO,ESRVAL,IREP,PMAT)
C*****************************************************************************
C
C    
C
C    
C     
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER(D0 = 0.0D0)
C
#include <dcbesr.h>
#include <dgroup.h>
#include <dcborb.h>
C
      DIMENSION ESRVAL(NCIESR,NCIESR,2),PMAT(*)
      INTEGER   NMO(*)
C
          if (nmo(1) .eq. 1 .and. nmo(2) .eq. 0) then
             JMOSY = 1
          else if (nmo(1) .eq. 0 .and. nmo(2) .eq. 1) then
             JMOSY = 2
          else
             call quit('Only one open shell implemented yet!')
          end if
          DO J = 1,NZ
             PVAL = PMAT(1+I2ORBX(JMOSY,JMOSY)+N2ORBX*(J-1))
             JQ = IPQTOQ(J,IREP)
             IF (JQ .EQ. 2) THEN
                ESRVAL(1,1,1) =  PVAL
                ESRVAL(2,2,1) = -PVAL
             ELSE IF (JQ .EQ. 3) THEN
                ESRVAL(1,2,2) = -PVAL
                ESRVAL(2,1,2) =  PVAL
             ELSE IF (JQ .EQ. 4) THEN
                ESRVAL(1,2,1) =  PVAL
                ESRVAL(2,1,1) =  PVAL
             ELSE
C               ... JQ .eq. 1
                ESRVAL(1,1,2) =  PVAL
                ESRVAL(2,2,2) =  PVAL
                WRITE(LUPRI, *) 'test PVAL zero?: JQ,J,IREP,JMOSY,PVAL',
     &             JQ,J,IREP,JMOSY,PVAL
             END IF
          END DO
          RETURN
          END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpesreven */
      SUBROUTINE PRPESREVEN(NMO,ESRVAL,IREP,PMAT)
C*****************************************************************************
C
C
C
C
C
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER(D0 = 0.0D0)
C
      DIMENSION ESRVAL(*),PMAT(*)
      INTEGER   NMO(*)
C
      CALL QUIT('only implemented for odd number of active electrons')
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#ifdef WORK_IN_PROGRESS
C  /* Deck PRPESR2 */
      SUBROUTINE PRPESR2(LU1INT,INDXPR,PVAL,LESRST,
     &           DMAT,PMAT,PTRI,OP1INT,FIRST,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Calculate expectation values defined in /CBIESR/
C
C     A. Noerager - Mar 2000
C     Based on prpex2 by Trond Saue May 27 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0)
C
#include <mxcent.h>
#include <dcbxpr.h>
#include <dcbesr.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcborb.h>
C
      LOGICAL FIRST(NZ),LESRST(2), RODHF
      DIMENSION DMAT(N2BBASX,NZ,*),PMAT(N2BBASX,NZ),
     &          PTRI(NNBBASX,NZ),OP1INT(NNBBASX),
     &          WORK(LWORK)
C
C     Check if single open shell calculation (.OPEN 1 1/2,0 or .OPEN 1 1/0,2)
C
      RODHF = NOPEN.EQ.1 .AND. NAELEC.EQ.1 .AND. NASHT.EQ.1
C
C     Check symmetry
C
      PVAL = D0
C
C     1. Molecular point group symmetry
C
      IREP = IPRPSYM(INDXPR)-1
      IF(IREP.NE.0) LESRST(1) = .TRUE.
C
C     2. Time reversal symmetry
C
      ITIM = IPRPTIM(INDXPR)
      IF(ITIM.NE.1) LESRST(2) = .TRUE.
C
C     If not totally symmetric, return with zero value
C
      IF(LESRST(1) .OR. (LESRST(2) .AND. .NOT.RODHF)) RETURN
C
C     Totally symmetric operator ==> proceed
C
      CALL LSET(NZ,.TRUE.,FIRST)
      CALL PRPMAO(LU1INT,INDXPR,PTRI,OP1INT,FIRST,WORK,LWORK,IPRINT)
C
C     Get expectation value, using only the matrices containing
C     integrals
C
      DO IZ = 1,NZ
      IF(.NOT.FIRST(IZ)) THEN
        IQ = IPQTOQ(IZ,IREP)
        IH = IHQMAT(IQ,ITIM)
        IF(IH.EQ.1) THEN        
          CALL DSPTSI(NTBAS(0),PTRI(1,IZ),PMAT(1,IZ))
        ELSEIF(IH.EQ.2) THEN
          CALL DAPTGE(NTBAS(0),PTRI(1,IZ),PMAT(1,IZ))
        ENDIF
        CALL BUTOBS(PMAT(1,IZ),1,WORK,LWORK)
          IF (LESRST(2) .AND. RODHF) THEN
C           ... < p | O | p > for ITIM = -1 :
            PVAL = PVAL +
     &         DF(1)*DDOT(N2BBASX,DMAT(1,IZ,2),1,PMAT(1,IZ),1)
            write (6,*) 'iz, df(1),esrval',iz,df(1),esrval
            write (6,*) 'dmat(1,iz,2)'
            nbast = ntbas(0)
            call output(dmat(1,iz,2),1,nbast,1,nbast,nbast,nbast,1,6)
            write (6,*) 'pmat(1,iz)'
            call output(pmat(1,iz),1,nbast,1,nbast,nbast,nbast,1,6)
          ELSE
            DO IOPEN = 0,NOPEN
              PVAL = PVAL + DF(IOPEN)
     &           *DDOT(N2BBASX,DMAT(1,IZ,IOPEN+1),1,PMAT(1,IZ),1)
            END DO
          ENDIF
      ENDIF
      ENDDO
      PVAL = PVAL + PVAL
C
      RETURN
      END
#endif /* WORK_IN_PROGRESS */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck PRP_ESR */
      SUBROUTINE PRP_ESR(ESRVAL,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Output routine for ESR properties
C     This routine is modelled after PRP_SHIELD
C
C     Written by A.Noerager APR 2000
C     
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C     
#include <dcbesr.h>
#include <nuclei.h>
C
      DIMENSION ESRVAL(NCIESR,NCIESR,2,*),WORK(*)
      KFRSAV = KFREE
      NCMAT =  NUCDEP*9
      CALL MEMGET('REAL',KTMAT,NCMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTISO,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KATOM,NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIATM,NUCDEP*3 ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRATM,NUCDEP*15,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCSTRA,NUCDEP*NUCDEP*9,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSCTRA,NUCDEP*NUCDEP*9,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMPMAT,NCMAT,WORK,KFREE,LFREE)
C     
      CALL SETATM(WORK(KATOM),NATOM,29)
      CALL PRP_ESR1(ESRVAL,WORK(KTMPMAT),WORK(KTMAT),WORK(KTISO),
     &     WORK(KIATM),WORK(KRATM),WORK(KIBUF),
     &     WORK(KATOM),WORK(KCSTRA),WORK(KSCTRA),WORK,KFREE,LFREE)
      CALL MEMREL('PRP_ESR',WORK,1,KFRSAV,KFREE,LFREE)
C     
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck PRP_ESR1 */
      SUBROUTINE PRP_ESR1(ESRVAL,TMPMAT,TMAT,TISO,IATINF,RATINF,
     &      IBUFA,DOATOM,CSTRA,SCTRA,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Output routine for ESR properties
C     
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
C
C Used from common blocks:
C  dcbprp: ABUND,?
C
#include <dcbprp.h>
#include <dcbxpr.h>
#include <dcbprl.h>
#include <dcbesr.h>
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      LOGICAL TEST,DOATOM(NUCIND)
      DIMENSION ESRVAL(NCIESR,NCIESR,2,*), TMPMAT(3,*),TMAT(3,3,NUCDEP),
     &          TISO(NUCDEP),IATINF(3,NUCDEP),RATINF(3,5,NUCDEP),
     &          IBUFA(3*NUCDEP),IBUFB(3),CSTRA(NUCDEP*NUCDEP*9),
     &          SCTRA(NUCDEP*NUCDEP*9),WORK(*)
#include <nuclei.h>
#include <symmet.h>
#include <ibtfun.h>
C     
      CALL TITLER('DIRAC  - ESR PROPERTIES','*',116)
C
C
C     ***************************************
C     ***** Extract nuclear information *****
C     ***************************************
C
      WRITE(LUPRI,'(A)') '@ Nuclear species:'
      CALL PRSYMB(LUPRI,'=',61,0)
      WRITE(LUPRI,'(A)')
     &'@ name  charge  isotope   mass         abundance    g_N factor'
      CALL PRSYMB(LUPRI,'-',61,0)
      NSPEC = 0
      IATOM = 0
      DO INUC = 1,NUCIND
      IDEG = NUCDEG(INUC)
      IF(DOATOM(INUC)) THEN
        TEST = .FALSE.
        IZ   = NINT(CHARGE(INUC))
        NISO = 0
        DO ISO = 1,5
          GVAL = DISOTP(IZ,ISO,'GVAL')
          IF(GVAL.NE.0) THEN
            ABND = DISOTP(IZ,ISO,'ABUNDANCE')
            IF(ABND.GT.ABUND.OR.(.NOT.TEST)) THEN
              TEST = .TRUE.
              NISO = NISO + 1
              RATINF(1,NISO,NSPEC+1) = DISOTP(IZ,ISO,'A')
              RATINF(2,NISO,NSPEC+1) = ABND
              RATINF(3,NISO,NSPEC+1) = GVAL
            ENDIF
          ENDIF
        ENDDO
        IF(NISO.GT.0) THEN
          IATINF(1,NSPEC+1) = IATOM + 1
          IATINF(2,NSPEC+1) = IZ
          IATINF(3,NSPEC+1) = NISO
          WRITE(LUPRI,'
     &      (A1,1X,A4,5X,I3,6X,I3,3X,F10.6,3X,F9.3,3X,F10.6)') 
     &       '@',NAMN(INUC),IZ,1,(RATINF(K,1,NSPEC+1),K=1,3)
          DO J = 2,NISO
            WRITE(LUPRI,'(A1,19X,I3,3X,F10.6,3X,F9.3,3X,F10.6)') 
     &         '@',J,(RATINF(K,J,NSPEC+1),K=1,3)
          ENDDO
          DO K = 2,IDEG
            IATINF(1,NSPEC+K) = IATOM + K
            IATINF(2,NSPEC+K) = IZ
            IATINF(3,NSPEC+K) = NISO
            DO L = 1,NISO
              RATINF(1,L,NSPEC+K) = RATINF(1,L,NSPEC+1)
              RATINF(2,L,NSPEC+K) = RATINF(2,L,NSPEC+1)
              RATINF(3,L,NSPEC+K) = RATINF(3,L,NSPEC+1)
            ENDDO
          ENDDO
          NSPEC = NSPEC + IDEG
        ENDIF
      ENDIF
      IATOM = IATOM + IDEG
      ENDDO
      CALL PRSYMB(LUPRI,'-',61,0)
      IF(NSPEC.EQ.0) THEN
        WRITE(LUPRI,'(A)') '* No suitable isotopes. Returning'
        RETURN
      ENDIF
C     
      KFRSAV = KFREE
      NSCOOR = 3*NUCDEP
      NCMAT  = 9*NUCDEP
      CALL DZERO(TMPMAT,NCMAT)
C     
C     *** Extract requested nuclear coordinates and B fields ***
C     
      NA = 0
      DO I = 1,NSCOOR
         IF(IPESR(I,1).GT.0) THEN
            NA = NA + 1
            IBUFA(NA) = I
         ENDIF
      END DO
C     
      NB = 0
      DO I = 1, 3
         IF(IPESR(I+MXCOOR,1).GT.0) THEN
            NB = NB + 1
            IBUFB(NB) = I
         ENDIF
      END DO
C     
C     Fill ESR hfs anisotropy matrices Seff_p * TMPMAT * I_N,q
C
C     
         CALL PRSYMB(LUPRI,'-',72,2)

         CALL DZERO(TMPMAT,NCMAT)
         DO IA = 1,NA
            JA = IBUFA(IA)
            KA = IPESR(JA,1)
            TMPMAT(IPTAX(1,2),JA) =  ESRVAL(1,2,1,KA)
            TMPMAT(IPTAX(2,2),JA) = -ESRVAL(1,2,2,KA)
            TMPMAT(IPTAX(3,2),JA) =  ESRVAL(1,1,1,KA)
         ENDDO
         IF (IPRPRP.GE.4) THEN
            CALL HEADER('ESR hfs anisotropies'//
     &           ' in symmetry coordinates (B~I,au)',1)
            CALL FCPRI (TMPMAT,'SIGMANO',CSTRA,SCTRA)
         ENDIF
C     
C     Transform to non-symmetry basis
C     
         CALL DZERO(TMAT(1,1,1),NCMAT)
         CALL TRADIP(TMPMAT,TMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
C     
C     Analyze results
C     
         NSPEC = 0
         IATOM = 0
         DO 100 I = 1, NUCIND
            DO 200 ISYMOP = 0, MAXOPR
               IF (IBTAND(ISTBNU(I),ISYMOP).EQ.0) THEN
                  IATOM = IATOM + 1
                  IF (DOATOM(IATOM)) THEN
                     NSPEC = NSPEC + 1
                     CALL ESRANA(TMAT(1,1,IATOM), TISO(IATOM),
     &                 IATINF(1,NSPEC),RATINF(1,1,NSPEC),
     &                 NAMDEP(IATOM),NAMDPX(3*(IATOM-1)+1), IPRPRP)
                  END IF
               END IF
 200        CONTINUE
 100     CONTINUE
C     
C     
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck ESRANA */
      SUBROUTINE ESRANA(TMAT,AVETOT,IATINF,RATINF,
     &                  NAME,NAMEX,IPRINT)
C***********************************************************************
C     
C     Analyze routine for ESR properties
C     This routine is modelled after SHIANA in DALTON
C
C     April 2000 A. Noerager
C     
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
#include <codata.h>
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0)
      PARAMETER (FACTOR = ALPHAC*ALPHAC/(XFMP))
C                factor = 1/c^2 * \mu_N = 1/c^2 * 1/(2 m_p)
      CHARACTER NAME*6, NAMEX(3)*6
      DIMENSION TMAT(3,3), TMATNY(3,3), IATINF(3), RATINF(3,*),
     &          TMATX(3,3)
#include <chrxyz.h>
C #include <dcbxlr.h>
C
      CALL DSCAL(9,FACTOR,TMAT,1)
C
      NISO=IATINF(3)
  
      DO L=1,NISO	
         DO I=1,3
           DO J=1,3
              TMATNY(I,J)=RATINF(3,L)*TMAT(I,J)
              TMATX(I,J)=XTHZ*TMATNY(I,J) * 1.0D-6
           ENDDO
         ENDDO
      AVETOT = (TMATNY(1,1)+TMATNY(2,2)+TMATNY(3,3))/D3
      AVETOTX= XTHZ*AVETOT * 1.0D-6
       
       
C
C     IF (IPRINT.GE.1) THEN
      CALL TITLER('ESR anisotropy for '//NAME//':','=',1)
C
         WRITE (LUPRI,'(5(/,1X,A,1P,D20.10,A))')
     &   ' Isotropic anisotropy:   ',AVETOT,  ' au'
         WRITE (LUPRI,'(5(/,1X,A,F15.2,A))')
     &   ' Isotropic anisotropy:        ',AVETOTX,  ' MHz'
         WRITE (LUPRI,'(5(/,1X,A,F15.8,A))')
     &   ' Isotope number:                ',RATINF(1,L),  ' ' 
      CALL HEADER('Total anisotropy tensor (au):',1)
      WRITE (LUPRI,'(18X,3(A,13X),/)') 'Sx', 'Sy', 'Sz'
      DO 400 J = 1, 3
         WRITE (LUPRI,'(2X,A6,1P,3D20.10)')
     &          NAMEX(J), (TMATNY(I,J),I=1,3)
  400 CONTINUE
      CALL HEADER('Total anisotropy tensor (MHz):',1)
      WRITE (LUPRI,'(18X,3(A,13X),/)') 'Sx', 'Sy', 'Sz'
      DO 500 J = 1, 3
         WRITE (LUPRI,'(2X,A6,3F15.2)')
     &          NAMEX(J), (TMATX(I,J),I=1,3)
  500 CONTINUE
      ENDDO
C
C     END IF
      RETURN
      END
