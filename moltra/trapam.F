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
C FILE    : trapam.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PAMTRA */
      SUBROUTINE PAMTRA(WORK,LWORK)
C***********************************************************************
C
C     Driver routine for 4-index transformation
C
C     Written by T.Saue, J.Laerdahl and L.Visscher Jan 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
#include <mxorb.h>
#include <dcbtra.h>
#include <dcbtr3.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbgen.h>
      CHARACTER*24 DAYTID
      LOGICAL TOBE,TRASAM
      DIMENSION WORK(*),NSTR(2,0:2,4),KVEC(2,4),KQ(2,4),
     &          IR(2,2),NR(2),NQ(2,4),NQT(4)
C
      IF (PRPTRA) CALL TRAPRP(WORK,LWORK)
      IF (.NOT.NO2IND) CALL TRAONE(WORK,LWORK)
      IF (NO4IND) RETURN
      CALL QENTER('PAMTRA')
#include <memint.h>
      KFRSAV = KFREE
      DUM = DUMMY
C
      CALL GETTIM(CPU1,WALL1)
      CALL TRAHI(IPRTRA,3)
C
C     Set up index arrays for active orbitals
C     =======================================      
C
      IF(NOPAIR) THEN
        DO I = 1,NFSYM
          NR(I)   = NESH(I)
          IR(1,I) = 1
          IR(2,I) = NESH(I)
        ENDDO
      ELSE
        DO I = 1,NFSYM
          NR(I)   =  NORB(I)
          IR(1,I) = -NPSH(I)
          IR(2,I) =  NESH(I)
        ENDDO
      ENDIF
C
      DO I = 1,4
        DO IFRP = 1,NFSYM
          CALL MEMGET('INTE',KVEC(IFRP,I),NR(IFRP),WORK,KFREE,LFREE)
          NSTR(IFRP,0,I) = 1
          CALL  NUMLST(TRA4_INDSTR(I,IFRP),WORK(KVEC(IFRP,I)),
     &                 NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &                 NSTR(IFRP,0,I),WORK(KFREE),LFREE)
          CALL ORBCNT(WORK(KVEC(IFRP,I)),NSTR(IFRP,0,I),IFRP,
     &              NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          NSTR(IFRP,0,I) = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
          CALL MEMREL('PAMTRA',WORK,KFRSAV,KVEC(IFRP,I),KFREE,LFREE)
          CALL MEMGET('INTE',KVEC(IFRP,I),NSTR(IFRP,0,I),
     &                WORK,KFREE,LFREE)
        ENDDO
      ENDDO      
C
C     Print section
C     =============
C
      DO I = 1, NFSYM
         CALL TRAPRI(4,I,WORK(KVEC(I,1)),WORK(KVEC(I,2)),
     &                 WORK(KVEC(I,3)),WORK(KVEC(I,4)),NSTR)
      ENDDO
C
C     Calculate dimensions of coefficient arrays
C     The coefficients are stored as 4 arrays of two matrices each
C     They may share the same memory if the transformation ranges
C     are identical
C
C     NDMOQR : 1st dimension row/column dimension of coefficient array
C              2nd dimension fermion irrep
C              3rd dimension transformation index
C
C     ICMOQR : Index to start of representation in the coefficient array
C
      DO I = 1, 4
         NQT(I) = 0
         ICMOQR(1,I) = 1
         DO IFRP = 1, NFSYM
            NDMOQR(1,IFRP,I) = NFBAS(IFRP,0)
            NDMOQR(2,IFRP,I) = NSTR(IFRP,0,I)
            NQ(IFRP,I) = NFBAS(IFRP,0)*NSTR(IFRP,0,I)*NZ
            NQT(I) = NQT(I) + NQ(IFRP,I)
            IF (IFRP.LT.NFSYM) 
     &         ICMOQR(IFRP+1,I) = ICMOQR(IFRP,I) + NQ(IFRP,I)
         ENDDO
      ENDDO
C
C     Allocate the memory for the coefficients, check their relations
C
      DO I = 1, 4
C
         DO J = I, 1, -1
            K = 2
            IF (NFSYM.EQ.1) K=1
            IF (TRASAM(WORK(KVEC(1,I)),WORK(KVEC(K,I)),
     &                 WORK(KVEC(1,J)),WORK(KVEC(K,J)),
     &          NSTR(1,0,I),NSTR(1,0,J))) ISAME(I) = J
         ENDDO
C
         IF (ISAME(I).NE.I) THEN
            KQ(1,I) = KQ(1,ISAME(I))
         ELSE
            CALL MEMGET('REAL',KQ(1,I),NQT(I),WORK,KFREE,LFREE)
         ENDIF
         KQ(2,I) = KQ(1,I) + NQ(1,I)
C
      ENDDO
C
C     Get all coefficients
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(A)') 'Coefficient file not present'
         CALL QUIT('TRAPAM: No coefficients !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMTRA')
         CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,TOTERG,3)
         CLOSE(LUCOEF,STATUS='KEEP')
      ENDIF
C
C     Select the set that we need
C
      DO I = 1, 4
        DO IFRP = 1, NFSYM
          IF(NSTR(IFRP,0,I).GT.0) THEN
             CALL SELCFS (WORK(KCMO+ICMOQ(IFRP)),IFRP,WORK(KQ(IFRP,I)),
     &                   WORK(KVEC(IFRP,I)),
     &                   NSTR(IFRP,2,I),NSTR(IFRP,1,I),
     &                   NFBAS(IFRP,0),NORB(IFRP))
          ENDIF
        ENDDO
      ENDDO
C
C     Throw away the complete set, keep only the selected sets
C
      CALL MEMREL('TRADRV',WORK,KFRSAV,KCMO,KFREE,LFREE)
C
C     Call driver for 4-index transformation
C     --------------------------------------
C
      CALL PAMTR1(WORK,KFREE,LFREE,IPRTRA,ISTRTRA,
     &            WORK(KQ(1,1)),WORK(KQ(1,2)),
     &            WORK(KQ(1,3)),WORK(KQ(1,4)),DUM,DUM,
     &            NDMOQR,ICMOQR,NSTR,.FALSE.,TRA_MSOUT,TRA_ANTIS,
     &            ITRA_INTFL4)
C
C      Print timing information
C
      CALL GETTIM(CPU2,WALL2)
      WALL   = WALL2 - WALL1
C
      IMINS  = INT(WALL)/60
      IHOURS = IMINS/60
      IMINS  = IMINS - 60*IHOURS
      ISECS  = NINT(WALL) - 3600*IHOURS - 60*IMINS
      WRITE(LUPRI,'(//A,I5.2,A,I2.2,A,I2.2)')
     &      ' Total wall time used in TRAPAM :',
     &      IHOURS,':',IMINS,':',ISECS
C
      CPU    = CPU2 - CPU1
      IMINS  = INT(CPU)/60
      IHOURS = IMINS/60
      IMINS  = IMINS - 60*IHOURS
      ISECS  = NINT(CPU) - 3600*IHOURS - 60*IMINS
      WRITE(LUPRI,'(A,I5.2,A,I2.2,A,I2.2)')
     &      ' Total CPU  time used in TRAPAM :',
     &      IHOURS,':',IMINS,':',ISECS
C
      WRITE(LUPRI,'(/A,F20.6)')
      CALL GTINFO(DAYTID)
      WRITE(LUPRI,'(/A,A24)') ' Transformation ended at : ',DAYTID
C
      CALL QEXIT('PAMTRA')
C
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PSIMP2 */
      SUBROUTINE PSIMP2(WORK,LWORK)
C***********************************************************************
C
C     Driver routine for MP2-calculation
C
C     Written by T.Saue, J.Laerdahl and L.Visscher Jan 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbmp2.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbgen.h>
      LOGICAL TOBE,TRASAM
      DIMENSION WORK(*),NSTR(2,0:2,4),KVEC(2,2),KQ(2,2),
     &          NR(2),IR(2,2),NQ(2),KE(2,2)
      DIMENSION NDMOQR(2,2,4),ICMOQR(2,4),NTHREI(2)
C
      CALL QENTER('PSIMP2')
#include <memint.h>
      KFRSAV = KFREE
      CALL MP2HI(IPRMP2)
C     
C     Get the position of the largest eigenvalue below the
C     energy threshold DMP2_VIRTHR
C     
      IF(DMP2_VIRTHR.GT.-1d10) THEN
         CALL PCKTOP(DMP2_VIRTHR,NTHREI,IPRMP2,WORK,LWORK)
         IF(NTHREI(1).LT.0) CALL QUIT
     &        ('No Eigenvalues below threshold DMP2_VIRTHR.')
         IF((NTHREI(2).LT.0).AND.(NFSYM.EQ.2)) CALL QUIT
     &        ('No Eigenvalues below threshold DMP2_VIRTHR.')
      ENDIF
C
C     Set up index arrays for active orbitals
C     =======================================      
C
      DO I = 1,2
         DO IFRP = 1,NFSYM
            IF (I.EQ.1) THEN
               IF(NOPAIR) THEN
                  NR(IFRP)   =  NISH(IFRP)
                  IR(1,IFRP) =  1
                  IR(2,IFRP) =  NISH(IFRP)
               ELSE
                  NR(IFRP)   =  NISH(IFRP) + NPSH(IFRP)
                  IR(1,IFRP) = -NPSH(IFRP)
                  IR(2,IFRP) =  NISH(IFRP)
               ENDIF
            ELSE
               NVIR = NESH(IFRP)
               IF(DMP2_VIRTHR.GT.-1d10) NVIR = NTHREI(IFRP)
               IF(NOPAIR) THEN
                  NR(IFRP)   = NVIR - NISH(IFRP)
                  IR(1,IFRP) = NISH(IFRP) + 1
                  IR(2,IFRP) = NVIR
               ELSE
                  NR(IFRP)   = NVIR - NISH(IFRP) + NPSH(IFRP) 
                  IR(1,IFRP) = NISH(IFRP) + 1
                  IR(2,IFRP) = NVIR
               ENDIF
            ENDIF
C     
            CALL MEMGET('INTE',KVEC(IFRP,I),NR(IFRP),WORK,KFREE,LFREE)
            NSTR(IFRP,0,I) = 1
            CALL  NUMLST(MP2_INDSTR(I,IFRP),WORK(KVEC(IFRP,I)),
     &           NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &           NSTR(IFRP,0,I),WORK(KFREE),LFREE)
            CALL ORBCNT(WORK(KVEC(IFRP,I)),NSTR(IFRP,0,I),IFRP,
     &           NSTR(IFRP,2,I),NSTR(IFRP,1,I))
            NSTR(IFRP,0,I)   = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
            NSTR(IFRP,0,I+2) = NSTR(IFRP,0,I)
            NSTR(IFRP,1,I+2) = NSTR(IFRP,1,I)
            NSTR(IFRP,2,I+2) = NSTR(IFRP,2,I)
            CALL MEMREL('PSIMP2',WORK,KFRSAV,KVEC(IFRP,I),KFREE,LFREE)
            CALL MEMGET('INTE',KVEC(IFRP,I),NSTR(IFRP,0,I),
     &           WORK,KFREE,LFREE)
         ENDDO
      ENDDO
C
C     Print section
C     =============
C
      IF (IPRMP2.GE.0) THEN
         WRITE(LUPRI,'(/1X,A/)') '* Active MP2 space:'
         DO I = 1, NFSYM
            CALL MP2PRI(I,WORK(KVEC(I,1)),WORK(KVEC(I,2)),NSTR)
         ENDDO
      ENDIF
C
C     Calculate dimensions of coefficient arrays
C     The coefficients are stored as 2 arrays of two matrices each
C
C     NDMOQR : 1st dimension row/column dimension of coefficient array
C              2nd dimension fermion irrep
C              3rd dimension transformation index
C
C     ICMOQR : Index to start of representation in the coefficient array
C
      DO I = 1,2
         NSTRT = 0
         DO IFRP = 1, NFSYM
            NSTRT = NSTRT + NSTR(IFRP,0,I)
         ENDDO
         CALL MEMGET('REAL',KE(1,I),NSTRT,WORK,KFREE,LFREE)
         KE(2,I) = KE(1,I) + NSTR(1,0,I)
      ENDDO
C
      DO I = 1,2
         NQTOT = 0
         ICMOQR(1,I)   = 1
         DO IFRP = 1, NFSYM
            NDMOQR(1,IFRP,I) = NFBAS(IFRP,0)
            NDMOQR(2,IFRP,I) = NSTR(IFRP,0,I)
            NQ(IFRP)         = NFBAS(IFRP,0)*NSTR(IFRP,0,I)*NZ
            NQTOT            = NQTOT + NQ(IFRP)
            IF (IFRP.LT.NFSYM) THEN
               ICMOQR(IFRP+1,I) = ICMOQR(IFRP,I) + NQ(IFRP)
            ENDIF
C
            ICMOQR(IFRP,I+2)   = ICMOQR(IFRP,I)
            NDMOQR(1,IFRP,I+2) = NDMOQR(1,IFRP,I)
            NDMOQR(2,IFRP,I+2) = NDMOQR(2,IFRP,I)
         ENDDO
         CALL MEMGET('REAL',KQ(1,I),NQTOT,WORK,KFREE,LFREE)
         KQ(2,I) = KQ(1,I) + NQ(1)
      ENDDO
C
C     Get all coefficients and eigenvalues
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT ,WORK,KFREE,LFREE)
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(A)') 'Coefficient file not present'
         CALL QUIT('PSIMP2: No coefficients !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PSIMP2')
         IP=2
         CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),DHFEMP,IP)
         CLOSE(LUCOEF,STATUS='KEEP')
         WRITE(LUPRI,'(A,F18.10,/)') 
     &        'Coefficients imported. Total DHF energy:', DHFEMP
      ENDIF
C     
C     Select the set that we need
C
      DO I = 1, 2
         DO IFRP = 1, NFSYM
            CALL SELCFS (WORK(KCMO+ICMOQ(IFRP)),
     &           IFRP,WORK(KQ(IFRP,I)),WORK(KVEC(IFRP,I)),
     &           NSTR(IFRP,2,I),NSTR(IFRP,1,I),
     &           NFBAS(IFRP,0),NORB(IFRP))
            CALL SELEIG(WORK(KEIG+IORB(IFRP)),IFRP,WORK(KE(IFRP,I)),
     &        WORK(KVEC(IFRP,I)),NSTR(IFRP,2,I),NSTR(IFRP,1,I))
        ENDDO
      ENDDO
C
      IF (IPRMP2.GE.1) THEN
         WRITE(LUPRI,'(A)') 'Active Eigenvals:'
         DO I = 1, 2
            DO IFRP = 1,NFSYM
               WRITE(LUPRI,'(/2I8)') I, IFRP
               WRITE(LUPRI,'(5D15.8)') 
     &              (WORK(KE(IFRP,I)+IXX-1), IXX=1,NSTR(IFRP,0,I))
            ENDDO
         ENDDO
      ENDIF
C
C     Throw away the complete set, keep only the selected sets
C
      CALL MEMREL('PSIMP2',WORK,KFRSAV,KCMO,KFREE,LFREE)
C
C     Call driver for 4-index transformation
C     --------------------------------------
C
      CALL PAMTR1(WORK,KFREE,LFREE,IPRMP2,ISTRMP2,
     &            WORK(KQ(1,1)),WORK(KQ(1,2)),
     &            WORK(KQ(1,1)),WORK(KQ(1,2)),
     &            WORK(KE(1,1)),WORK(KE(1,2)),
     &            NDMOQR,ICMOQR,NSTR,.TRUE.,MP2_MSOUT,MP2_ANTIS,
     &            MP2_INTFLG)
C
      CALL QEXIT('PSIMP2')
C
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MP2PRI */
      SUBROUTINE MP2PRI(I,IOS1,IOS2,NSTR)
C***********************************************************************
C
C     Initial print routine for MP2 module
C
C     Written by T.Saue Jan 19 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcborb.h>
      DIMENSION IOS1(*),IOS2(*),NSTR(0:2,2,4)
C
      WRITE(LUPRI,'(1X,A,A3)') '* Fermion ircop ',FREP(I)
C
C     Index 1: Active occupied
C
      IF(NSTR(0,I,1).GT.0) THEN
        WRITE(LUPRI,'(6X,A,I5)') 
     &     '- Active occupied orbitals:', NSTR(0,I,1)
        WRITE(LUPRI,1002) (IOS1(K),K=1,NSTR(0,I,1))
      ELSE
        WRITE(LUPRI,'(6X,A)') 
     &     '- No active occupied orbitals'
      ENDIF
C
C     Index 2: Active virtuals
C
      IF(NSTR(0,I,2).GT.0) THEN
        WRITE(LUPRI,'(6X,A,I5)') 
     &     '- Active virtual orbitals: ', NSTR(0,I,2)
        WRITE(LUPRI,1002) (IOS2(K),K=1,NSTR(0,I,2))
      ELSE
        WRITE(LUPRI,'(6X,A)') 
     &     '- No active virtual orbitals'
      ENDIF
      WRITE(LUPRI,'(A)') ''
C
 1002 FORMAT (6X,10I5)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck TRAPRI */
      SUBROUTINE TRAPRI(NINDEX,I,IOS1,IOS2,IOS3,IOS4,NSTR)
C***********************************************************************
C
C     Initial print routine for 4-index transformation module
C
C     Written by T.Saue Jan 19 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcborb.h>
      DIMENSION IOS1(*),IOS2(*),IOS3(*),IOS4(*)
      DIMENSION NSTR(0:2,2,NINDEX)
C
      WRITE(LUPRI,'(//1X,A,A3/)') '* Fermion ircop ',FREP(I)
C
C     Index 1
C
      IF(NSTR(0,I,1).GT.0) THEN
        WRITE(LUPRI,1000) 1,NSTR(0,I,1)
        WRITE(LUPRI,1002) (IOS1(K),K=1,NSTR(0,I,1))
      ELSE
        WRITE(LUPRI,1001) 1
      ENDIF
C
C       Index 2
C
      IF (NINDEX.LT.2) RETURN
      IF(NSTR(0,I,2).GT.0) THEN
        WRITE(LUPRI,1000) 2,NSTR(0,I,2)
        WRITE(LUPRI,1002) (IOS2(K),K=1,NSTR(0,I,2))
      ELSE
        WRITE(LUPRI,1001) 2
      ENDIF
C
C       Index 3
C
      IF (NINDEX.LT.3) RETURN
      IF(NSTR(0,I,3).GT.0) THEN
        WRITE(LUPRI,1000) 3,NSTR(0,I,3)
        WRITE(LUPRI,1002) (IOS3(K),K=1,NSTR(0,I,3))
      ELSE
        WRITE(LUPRI,1001) 3
      ENDIF
C
C       Index 4
C
      IF (NINDEX.LT.4) RETURN
      IF(NSTR(0,I,4).GT.0) THEN
        WRITE(LUPRI,1000) 4,NSTR(0,I,4)
        WRITE(LUPRI,1002) (IOS4(K),K=1,NSTR(0,I,4))
      ELSE
        WRITE(LUPRI,1001) 4
      ENDIF
C
 1000 FORMAT (/6X,'Index ',I1,I5,' orbitals')
 1001 FORMAT (/6X,'No orbitals for index ',I1)
 1002 FORMAT (6X,12I5)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PAMTR1 */
      SUBROUTINE PAMTR1(WORK,KFREE,LFREE,IPRINT,ISTRA,
     &                  Q1,Q2,Q3,Q4,EIG1,EIG2,NDMOQR,ICMOQR,NSTR,
     &                  LMP2,MSOUT,ANTIS,INTFLG)
C***********************************************************************
C
C     Subdriver for 4index transformation.
C
C     Qn             - coefficient set for index n
C     EIG            - eigenvalues (for  MP2 calculation)
C     NSTR(ifrp,0,i) - total number of orbitals for index i
C     NSTR(ifrp,1,i) - number of electronic orbitals for index i
C     NSTR(ifrp,2,i) - number of positronic orbitals for index i
C     LMP2  = .TRUE. - evaluate MP2 energy
C     MSOUT = .TRUE. - write transformed integrals to file
C     ANTIS = .TRUE. - anti-symmetrize integrals
C     INTFLG - flag of what integral types to transform
C                      
C     Written by Luuk Visscher, Trond Saue and J.Laerdahl Jan 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <symmet.h>
#include <dcbgen.h>
#include <dcborb.h>
#include <dcbtra.h>
#include <dcbbas.h>
#include <dgroup.h>
C
      LOGICAL MSOUT,ANTIS,LMP2
      DIMENSION WORK(*),NSTR(2,0:2,4)
      DIMENSION NDMOQR(2,2,4),ICMOQR(2,4)
      DIMENSION Q1(*),Q2(*),Q3(*),Q4(*)
      DIMENSION EIG1(*),EIG2(*)
C
#include <dcbibt.h>
C
      KFRSAV = KFREE
C
C     Set up symmetry packing information.
C
      CALL MEMGET('INTE',KINDX,3*NTBAS(0),WORK,KFREE,LFREE)
      CALL GISPCK(NSPCK,ISPCK,0,0,IDUM,IDUM,IDUM,IDUM,IDUM,IDUM)
      MXINSH = NINSH(0,-1,WORK(KINDX),IPRINT)
C
C     Define the MS integral types that we are going to make
C
      IF (NZ.EQ.1) THEN
         NCLASS = 4
      ELSEIF (NZ.EQ.2) THEN
         NCLASS = 8
      ELSE
         NCLASS = 16
      ENDIF
C
      CALL DEFCL (NCLASS)
C
C
C     We have three possible strategies : 
C
C     1) In-core test version : half- and fully-transformed integrals
C        should fit in memory
C     2) Loop over restricted range of MS-labels such that all transf.
C        integrals fit in memory
C        Loop over distributions transform 4 indices(rest. ranges first)
C        Accumulate integrals in memory, calculate MP2 energy
C     3) Loop over distributions, transform 2 indices and write to disk
C        Loop over transformed indices and transform last 2 indices
C     4) Like 3, but 2-index transformed integrals written to disk in
C        batches; parallel scheme.
C
C     In all cases we transform in mulliken notation
C
      CALL QENTER('PAMTR1')
C
      IF (ISTRA.EQ.1) THEN
          CALL TRDR1T(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR,NSTR,
     &                ANTIS,LMP2,WORK(KINDX),Q1,Q2,Q3,Q4,EIG1,EIG2)
      ELSEIF (ISTRA.EQ.2) THEN
          CALL TRDR2T(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR,NSTR,
     &                ANTIS,LMP2,WORK(KINDX),Q1,Q2,Q3,Q4,EIG1,EIG2)
      ELSEIF (ISTRA.EQ.3) THEN
          CALL TRDR3T(WORK,KFREE,LFREE,IPRINT,INTFLG,NSTR,
     &                ANTIS,LMP2,WORK(KINDX),Q1,Q2,Q3,Q4,EIG1,EIG2)
      ELSEIF (ISTRA.EQ.4) THEN
          CALL TRDR4T(WORK,KFREE,LFREE,IPRINT,INTFLG,NSTR,
     &                ANTIS,LMP2,WORK(KINDX),Q1,Q2,Q3,Q4,EIG1,EIG2)
      ELSE
         WRITE(LUPRI,'(A)') 'PAMTR1: Strategy ',ISTRA,' not implemented'
         CALL QUIT('PAMTR1: Unknown strategy !')
      ENDIF
      CALL QEXIT('PAMTR1')
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MP2HI */
      SUBROUTINE MP2HI(IPRINT)
C
C     Written by Trond Saue and Luuk Visscher Jan 17 1996
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbgen.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbmp2.h>
C
C     Now say hi and give orbital information
C
      CALL TITLER('MP2 Calculation ','*',128)
      WRITE (LUPRI,'(/10X,A)') 'Written by Jon Laerdahl'//
     &     ', Luuk Visscher & Trond Saue'
      WRITE (LUPRI,'(10X,A//)') 
     &     'Odense (Jan. 1997), Oslo (June 1997), Auckland (Sept. 1998)'
C
      IF(IPRINT.GE.3) THEN
        CALL HEADER('SETDC2: Distribution of symmetry orbitals '//
     &     'among fermion ircops',-1)
        WRITE(LUPRI,'(7X,3A5)')
     &       'Total',' Lbas',' Sbas'
        WRITE(LUPRI,'(A3,4X,3I5)')
     &        (FREP(I),(NFBAS(I,J),J=0,2),I = 1,NFSYM)
      ENDIF
C
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck trahi */
      SUBROUTINE TRAHI(IPRINT,IMODUL)
C
C     Written by Luuk Visscher 05-11-1996.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbgen.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <mxorb.h>
#include <dcbtra.h>
      CHARACTER*24 DAYTID
C
C     Now say hi and give orbital information
C
      CALL TITLER('Transformation to Molecular Spinor Basis','*',116)
      WRITE (LUPRI,'(/10X,A)') 'Written by Luuk Visscher'//
     &     ', Jon Laerdahl & Trond Saue'
      WRITE (LUPRI,'(10X,A//)') 'Odense, 1997'
      IF (IMODUL.EQ.1) THEN
         CALL TITLER('Transformation of property integrals','*',116)
      ELSEIF (IMODUL.EQ.2) THEN
         CALL TITLER('Transformation of core Fock matrix','*',116)
      ELSEIF (IMODUL.EQ.3) THEN
         CALL TITLER('Transformation of 2-electron integrals','*',116)
      ENDIF
C
      CALL GTINFO(DAYTID)
      WRITE(LUPRI,'(/A,A24)') ' Transformation started at :',DAYTID
C
      IF(IPRINT.GE.3) THEN
        CALL HEADER('SETDC2: Distribution of symmetry orbitals '//
     &     'among fermion ircops',-1)
        WRITE(LUPRI,'(7X,3A5)')
     &       'Total',' Lbas',' Sbas'
        WRITE(LUPRI,'(A3,4X,3I5)')
     &        (FREP(I),(NFBAS(I,J),J=0,2),I = 1,NFSYM)
      ENDIF
C
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MS4IN1 */
      SUBROUTINE MS4IN1 (WORK,KFREE,LFREE,IPRINT,ICS,ICF,TRIAN,
     & NSTR1,NSTR2,NFPCK,IFPCK,NDMOQR,ICMOQR,IJPAIRS,Q1,Q2,GMAT,
     & HMAT,NZMO,IQMO)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Luuk Visscher January 1997.
C
C     PURPOSE : Do 4-index transformation to molecular spinor basis
C               First index-pair transformation.
C
C     Input :
C
C     - IPRINT       Print flag
C     - ICS          First class of integrals 1 : (LL|XX), 2 : (SS|XX)
C     - ICL          Last class of integrals
C     - TRIAN        Store the transformed integrals as a lower triangle
C     - NSTR1        Number of active spinors for index 1
C     - NSTR2        Number of active spinors for index 2
C     - NFPCK        Number of spinor pairs for each compound symmetry
C     - IFPCK        Pointers to spinor pairs
C     - NDMOQR       Dimensions of the coefficient array
C     - ICMOQR       Pointers to coefficients
C     - IJPAIRS      Number of blocks for each boson symmetry
C     - Q1           Coefficients for index 1
C     - Q2           Coefficients for index 2
C     - GMAT         Symmetry packed scalar integrals
C
C     Output :
C
C     - HMAT         Symmetry packed half-transformed integrals
C     - NZMO         Number of non-zero quaternion components of HMAT
C     - IQMO         Quaternion packing of HMAT
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0)
C
      LOGICAL   TRIAN
      INTEGER   NSTR1(2),NSTR2(2)
      DIMENSION WORK(*)
      DIMENSION Q1(*),Q2(*),GMAT(*),HMAT(*)
      DIMENSION NFPCK(2),IFPCK(2,2)
      DIMENSION IJPAIRS(0:7)
      DIMENSION IQAO(4),IQMO(4,0:7),NZMO(0:7)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2)
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <symmet.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbtra.h>
#include <dcbibt.h>
C
      CALL QENTER('MS4IN1')
C
      KFRSAV = KFREE
C
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
C
      NZAO = 1
      IQAO(1) = 1
      IQAO(2) = 0
      IQAO(3) = 0
      IQAO(4) = 0
C
C     Transform lefthand indices to molecular spinor basis
C
      IOFF = 1
      JOFF = 1
C
C     Loop over the boson symmetries of the densities
C
      DO IREPPQ = 0, NBSYM-1
         IREPIJ = JBTOF(IREPPQ,1)
C
C        Allocate memory for temporary storage of unsorted block of
C        half-transformed integrals. Allocate absolute maximum instead
C        of the exact maximum size.
C
         CALL MEMGET ('REAL',KHTMP,NFPCK(IREPIJ)*NZ,WORK,KFREE,LFREE)
         CALL DZERO (WORK(KHTMP),NFPCK(IREPIJ)*NZ)
C
C        Loop over the distributions with this boson symmetry in this
C        shell-block batch
C
         DO IJ = 1, IJPAIRS(IREPPQ)
C
C           Loop over the compound component (LL, SS) of the lefthand
C           densities
C
            DO IC = ICS, ICF
               IF (ICS.EQ.2) IOFF = IOFF + NSPCK(IREPPQ,1)
C
C              Loop over the boson symmetries of the second 
C              (untransformed) index. Get the fermion symmetry (parity)
C              of that index and get the boson and fermion symmetry
C              of the first index.
C
               DO IREPQ = 0, NBSYM-1
                  IREPP = IBTXOR(IREPPQ,IREPQ)
                  IREPI =  JBTOF(IREPP,IC)
                  IREPJ =  JBTOF(IREPQ,IC)
                  KOFF = ICMOQR(IREPI,1) + IBBAS(IREPP,IC) 
     &                 - IBAS(IREPI)
                  LOFF = ICMOQR(IREPJ,2) + IBBAS(IREPQ,IC) 
     &                 - IBAS(IREPJ)
                  NP = NBBAS(IREPP,IC)
                  NQ = NBBAS(IREPQ,IC)
                  NI = NSTR1(IREPI)
                  NJ = NSTR2(IREPJ)
                  NRQ1 = NDMOQR(1,IREPI,1)
                  NCQ1 = NDMOQR(2,IREPI,1)
                  NRQ2 = NDMOQR(1,IREPJ,2)
                  NCQ2 = NDMOQR(2,IREPJ,2)
                  IF (IC.EQ.1) THEN
                     IREPPI = IREPP
                     IREPQI = IREPQ
                  ELSE
                     IREPPI = IBTXOR(IXYZ,IREPP)
                     IREPQI = IBTXOR(IXYZ,IREPQ)
                  ENDIF
C
C                 Do the actual transformation. IQMO is set by QTRANS
C                 on basis of IQAO (real) and the IQ's of the 
C                 coefficients. We need to keep this information as it
C                 will tell us which were the non-zero quaternion units
C                 of the density, which is to be used in the final 
C                 symmetrization step. 
C
                  IF (.NOT.TRIAN.OR.(TRIAN.AND.IREPI.GT.IREPJ)) THEN
C
C                 Square transformation
C
                     IF ((NP*NQ.NE.0).AND.(NI*NJ.NE.0)) THEN
                       CALL QTRANS('AOMO','S',D0,NP,NQ,NI,NJ,
     &                    GMAT(IOFF),NP,NQ,NZAO,1,
     &                    WORK(KHTMP),NI,NJ,NZMO(IREPPQ),IQMO(1,IREPPQ),
     &                    Q1(KOFF),NRQ1,NCQ1,NZ,IPQTOQ(1,IREPPI),
     &                    Q2(LOFF),NRQ2,NCQ2,NZ,IPQTOQ(1,IREPQI),
     &                    WORK(KFREE),LFREE,IPRINT)
C
C                      Order the integrals with NKL as first and IZ
C                      as second index
C
                       DO IZ = 1, NZ
                          JOFF1 = JOFF + (IZ-1)*NFPCK(IREPIJ) +
     &                                   IFPCK(IREPI,IREPJ)
                          JOFF2 = KHTMP + (IZ-1)*NI*NJ
                          CALL DAXPY(NI*NJ,D1,WORK(JOFF2),1,
     &                               HMAT(JOFF1),1)
                       ENDDO
                     ENDIF
C
                  ELSEIF (TRIAN.AND.IREPI.EQ.IREPJ) THEN
C
C                 Triangular storage after square transformation
C
                     IF ((NP*NQ.NE.0).AND.(NI*NJ.NE.0)) THEN
                       CALL MEMGET('REAL',KBMAT,NI*NJ*NZ,WORK,
     &                              KFREE,LFREE)
                       CALL DZERO(WORK(KBMAT),NI*NJ*NZ)
                       CALL QTRANS('AOMO','S',D0,NP,NQ,NI,NJ,
     &                    GMAT(IOFF),NP,NQ,NZAO,1,
     &                    WORK(KBMAT),NI,NJ,NZMO(IREPPQ),IQMO(1,IREPPQ),
     &                    Q1(KOFF),NRQ1,NCQ1,NZ,IPQTOQ(1,IREPPI),
     &                    Q2(LOFF),NRQ2,NCQ2,NZ,IPQTOQ(1,IREPQI),
     &                    WORK(KFREE),LFREE,IPRINT)
                       CALL QTRIAN('PACK',D0,NI,NZ,WORK(KBMAT),NI,NJ,
     &                             WORK(KHTMP),NI*(NI+1)/2)
                       CALL MEMREL('MS4IN1',WORK,KFRSAV,KBMAT,
     &                             KFREE,LFREE)
C
C                      Order the integrals with NKL as first and IZ
C                      as second index
C
                       DO IZ = 1, NZ
                          JOFF1 = JOFF + (IZ-1)*NFPCK(IREPIJ) +
     &                                   IFPCK(IREPI,IREPJ)
                          JOFF2 = KHTMP + (IZ-1)*NI*(NI+1)/2
                          CALL DAXPY(NI*(NI+1)/2,D1,WORK(JOFF2),1,
     &                               HMAT(JOFF1),1)
                       ENDDO
                     ENDIF
                  ENDIF
                  IOFF = IOFF + NP * NQ
               ENDDO
            ENDDO
            JOFF = JOFF + NFPCK(IREPIJ)*NZ
         ENDDO
         CALL MEMREL('MS4IN1',WORK,KFRSAV,KHTMP,KFREE,LFREE)
      ENDDO
C
      CALL QEXIT('MS4IN1')
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck caldis */
      SUBROUTINE CALDIS(I2TYP,GMAT,INDX,INDXAB,
     &                  NODV,NOPV,WORK,LWORK,IPRINT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Luuk Visscher 22-10-1996
C
C     PURPOSE : Driver routine for the calculation of the two-electron
C               distributions (**|cd).
C               The distributions are stored as full squares with
C               symmetry reduction.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
#include <comdis.h>
#include <maxaqn.h>
#include <mxcent.h>
      DIMENSION GMAT(*),INDX(3,*),INDXAB(*)
      DIMENSION WORK(LWORK)
#include <cbihr2.h>
#include <dorps.h>
#include <ccom.h>
#include <dcborb.h>
C
      LOGICAL NODV, NOPV, RETUR, NOCONT
C
C
C-------------------------------------------------------------
C     Setup information for the two-electron integralroutines.
C-------------------------------------------------------------
C
      NODV = .FALSE.
      NOPV = .FALSE.
      IF (NASHT .EQ. 0) NODV = .TRUE.
      IF (NASHT .LT. 2) NOPV = .TRUE.
      IPRALL =  0
      DO 100 I = 0,7
         DOREPS(I) = .TRUE.
  100 CONTINUE
      DO 110 I = 1,MXCENT
         DOCOOR(1,I) = .TRUE.
         DOCOOR(2,I) = .TRUE.
         DOCOOR(3,I) = .TRUE.
  110 CONTINUE
      NOCONT = .FALSE.
C
      MAXDIF = 0
      ITYPE  = 4
C
C     Initialize buffer file if necessary
C
      IF (NPASS.GT.1) CALL INITGBF(LGFIL,NGBFSZ,NPASS,IRECG,
     &                             GMAT(KLGREC),GMAT(KLGBUF))
C
      CALL TWOINT(WORK,LWORK,DUMMY,DUMMY,0,IDUMMY,IDUMMY,GMAT,INDX,
     &            INDXAB,ITYPE,MAXDIF,0,NODV,NOPV,NOCONT,
     &            TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RTNTWO,
     &            IDUMMY,I2TYP,IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &            .TRUE.)
C
C     Finish buffer file if necessary
C
      IF (NPASS.GT.1) CALL FINAGBF (LGFIL,NPASS,IRECG,GMAT(KLGREC),
     &                              GMAT(KLGBUF),NGBFSZ,GMAT(KRGBUF),
     &                              GMAT(KIGBUF))
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck SHRNGE */
      SUBROUTINE SHRNGE (IC,INTFLG,I2TYP,IASTRT,IBSTRT,IASMAX,IBSMAX)
C
C     Written by Luuk Visscher, january 1997 
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <mxorb.h>
#include <aovec.h>
#include <blocks.h>
#include <dcbbas.h>
C
      LOGICAL L2TYP(0:3),LBIT
C
      L2TYP(0) = INTFLG.EQ.0
      L2TYP(1) = LBIT(INTFLG,1)
      L2TYP(2) = LBIT(INTFLG,2)
      L2TYP(3) = LBIT(INTFLG,3)
C
      IF (IC.EQ.1) THEN
         IASTRT = 1     
         IBSTRT = 1     
         IASMAX = NLRGBL
         IBSMAX = NLRGBL
         IF (L2TYP(0)) THEN
C           Non-relativistic option
            IASTRT = 1     
            IBSTRT = 1     
            IASMAX = MAXSHL
            IBSMAX = MAXSHL
            I2TYP = 0
         ELSEIF (L2TYP(1).AND.L2TYP(2)) THEN
C           (LL|LL) + (LL|SS) - integrals
            I2TYP = 12
         ELSEIF (L2TYP(1)) THEN
C           (LL|LL) - integrals
            I2TYP = 1
         ELSEIF (L2TYP(2)) THEN
C           (LL|SS) - integrals
            I2TYP = 2
         ELSE
C           No valid class of integrals for this value of IC
            IASTRT = 1     
            IBSTRT = 1     
            IASMAX = 0
            IBSMAX = 0
            I2TYP = -1
         ENDIF
      ELSE
         IASTRT = NLRGBL + 1     
         IBSTRT = NLRGBL + 1     
         IASMAX = MAXSHL
         IBSMAX = MAXSHL
         IF (L2TYP(0)) THEN
C           Non-relativistic option
            I2TYP = -1
         ELSEIF (L2TYP(2).AND.L2TYP(3)) THEN
C           (SS|LL) + (SS|SS) - integrals
            I2TYP = 12
         ELSEIF (L2TYP(2)) THEN
C           (SS|LL) - integrals
            I2TYP = 1
         ELSEIF (L2TYP(3)) THEN
            I2TYP = 2
C           (SS|SS) - integrals
         ELSE
C           No valid class of integrals for this value of IC
            IASTRT = 1     
            IBSTRT = 1     
            IASMAX = 0
            IBSMAX = 0
            I2TYP = -1
         ENDIF
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gispck */
      SUBROUTINE GISPCK(NSPCK,ISPCK,ISHELA,ISHELB,NBBASA,NBBASB,
     &                  IBBASA,IBBASB,IBASA,IBASB)
C*****************************************************************************
C
C     Make arrays needed for symmetry packing in 4-index transformation.
C     Luuk Visscher, january 1997
C
C*****************************************************************************
#include <implicit.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <mxorb.h>
#include <aovec.h>
#include <blocks.h>
C
      DIMENSION NSPCK(0:7,0:2),ISPCK(0:7,0:7,2)
      DIMENSION NBBASA(0:7,0:2),NBBASB(0:7,0:2)
      DIMENSION IBBASA(0:7,0:2),IBBASB(0:7,0:2)
      DIMENSION IBASA(2),IBASB(2)
C
#include <ibtfun.h>
C
C     Find out to which component these shells belong
C
      ICA = 0
      IF (ISHELA.NE.0) THEN
         IF (ISHELA.LE.NLRGBL) THEN
            ICA = 1
         ELSE
            ICA = 2
         ENDIF
C
C        We generate the basis function pointer arrays.
C
         IBASI = 0
         DO IFSYM = 1, NFSYM
            IBASA(IFSYM) = IBASI
            DO IREPA = 0, NBSYM-1
               IF (IFSYM.EQ.JBTOF(IREPA,ICA)) THEN
                  NBA = NINSH(ISHELA,IREPA,INDX,0)
                  NBBASA(IREPA,0)     = NBA
                  NBBASA(IREPA,ICA)   = NBA
                  NBBASA(IREPA,3-ICA) = 0
                  IBBASA(IREPA,ICA)    = IBASI
                  IBASI = IBASI + NBA
               ENDIF
            ENDDO
         ENDDO
         IF (ICA.EQ.1) THEN
            CALL ICOPY (NBSYM,IBASI,0,IBBASA(0,3-ICA),1)
         ELSE
            CALL ICOPY (NBSYM,0,0,IBBASA(0,3-ICA),1)
         ENDIF
      ENDIF
C
      ICB = 0
      IF (ISHELB.NE.0) THEN
         IF (ISHELB.LE.NLRGBL) THEN
            ICB = 1
         ELSE
            ICB = 2
         ENDIF
C
C        We generate the basis function pointer arrays.
C
         IBASI = 0
         DO IFSYM = 1, NFSYM
            IBASB(IFSYM) = IBASI
            DO IREPB = 0, NBSYM-1
               IF (IFSYM.EQ.JBTOF(IREPB,ICB)) THEN
                  NBB = NINSH(ISHELB,IREPB,INDX,0)
                  NBBASB(IREPB,0)     = NBB
                  NBBASB(IREPB,ICB)   = NBB
                  NBBASB(IREPB,3-ICB) = 0
                  IBBASB(IREPB,ICB)    = IBASI
                  IBASI = IBASI + NBB
               ENDIF
            ENDDO
         ENDDO
         IF (ICB.EQ.1) THEN
            CALL ICOPY (NBSYM,IBASI,0,IBBASB(0,3-ICB),1)
         ELSE
            CALL ICOPY (NBSYM,0,0,IBBASB(0,3-ICB),1)
         ENDIF
      ENDIF
C
C     We generate the boson symmetry packing, based on square matrices.
C
      CALL IZERO(NSPCK,24)
      CALL IZERO(ISPCK,128)
C
      DO IREPAB = 0, NBSYM-1
         DO IC = 1, 2
            DO IREPB = 0, NBSYM-1
               IREPA = IBTXOR(IREPAB,IREPB)
               ISPCK(IREPA,IREPB,IC) = NSPCK(IREPAB,0)
               IF (ISHELA.EQ.0) THEN
                  NBA = NBBAS(IREPA,IC)
               ELSE
                  NBA = NBBASA(IREPA,IC)
               ENDIF
               IF (ISHELB.EQ.0) THEN
                  NBB = NBBAS(IREPB,IC)
               ELSE
                  NBB = NBBASB(IREPB,IC)
               ENDIF
               NBAB = NBA * NBB
               NSPCK(IREPAB,0)  = NSPCK(IREPAB,0) + NBAB
               NSPCK(IREPAB,IC) = NSPCK(IREPAB,IC) + NBAB
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ninsh */
      FUNCTION NINSH(ISHELA,IREPA,INDX,IPRINT)
C
C     Counts number of functions in a shell for a given irrep.
C     When called with IREPA = -1 it gives the total number for a shell.
C     When called with ISHELA = 0 and IREPA = -1, it loops over all 
C     shells to create a pointer array from the complete list to the 
C     position in the shell and the symmetry. The function value is than
C     the maximum number of functions in a shell.
C
C     INDX(1,INDA) : Position of function in shell
C     INDX(2,INDA) : Irreducible representation of function
C     INDX(3,INDA) : Position of function within this particular irrep
C
C     Luuk Visscher, 24-1-97
C
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <aovec.h>
#include <blocks.h>
#include <dcbtra.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <symmet.h>
C
      DIMENSION INDX(3,*),IIREP(0:7),IIIREP(0:7)
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine NINSH',-1)
C
      IF (ISHELA.NE.0) THEN
         ISHELS = ISHELA
         ISHELF = ISHELA
      ELSE
         ISHELS = 1
         ISHELF = MAXSHL
      ENDIF
C
      IF (IREPA.NE.-1) THEN
         IREPS = IREPA
         IREPF = IREPA
      ELSE
         IREPS = 0
         IREPF = MAXREP
      ENDIF
C
      MXINSH = 0
      II = 0
      DO IREPA1 = IREPS, IREPF 
         IIIREP(IREPA1) = 0
      ENDDO
C
      DO ISHELA1 = ISHELS, ISHELF 
C
         I = 0
         DO IREPA1 = IREPS, IREPF 
            IIREP(IREPA1) = 0
         ENDDO
         NHKTA1 = NHKTSH(ISHELA1)
         KHKTA1 = KHKTSH(ISHELA1)
         MULA1  = ISTBSH(ISHELA1)
         NORBA1 = NORBSH(ISHELA1)
         NSTRA1 = IORBSB(IORBSH(ISHELA1,1))
C
C        Loop over components
C
         DO NA = 1,KHKTA1
            NSTRNA = NSTRA1 + NA
            ITYNA  = ISYMAO(NHKTA1,NA)
C
C           Loop over symmetries
C
            DO IREPA1 = IREPS, IREPF
               IF (IBTAND(MULA1,IBTXOR(IREPA1,ITYNA)) .EQ. 0) THEN
C
C                 Loop over contracted orbitals
C
                  DO IA = 1,NORBA1
                     INDA = IPTSYM(NSTRNA + KHKTA1*(IA-1),IREPA1)
                     I = I + 1
                     IIREP(IREPA1) = IIREP(IREPA1) + 1
                     IF (ISHELA.EQ.0) THEN
                         INDX(1,INDA) = I
                         INDX(2,INDA) = IREPA1
                         INDX(3,INDA) = IIREP(IREPA1)
                     ENDIF
                  ENDDO
               END IF
            ENDDO
         ENDDO
         II = II + I
         DO IREPA1 = IREPS, IREPF 
            IIIREP(IREPA1) = IIIREP(IREPA1) + IIREP(IREPA1)
         ENDDO
         MXINSH = MAX0(MXINSH,I)
      ENDDO
C
      IF (IREPA.NE.-1) THEN
         NINSH = IIIREP(IREPA)
      ELSE
         IF (ISHELA.NE.0) THEN
             NINSH = II
         ELSE
            IF (NTBAS(0).NE.II) CALL QUIT (' Counting error in NINSH')
            NINSH = MXINSH
         ENDIF
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ndistrn */
      FUNCTION NDISTRN(I2TYP,IPRINT,MKINAB,MKIJPA,ABSPOI,INDX,
     &                 INDXAB,IJPAIRS,NWDS,IJPASS)
C
C     Calculate number of distributions in a shell pair
C     Set up pointer array INDXAB
C
C     Luuk Visscher 24-01-97
C
C     Input :
C
C     - ISHLAB : Compound label of ISHELA and ISHELB
C     - I2TYP  : Integral class
C     - IPRINT : Print Flag
C     - MKINAB : .TRUE.   Make pointer array
C                .FALSE.  Only calculate number of distributions
C     - MKIJPA : .TRUE.   Make pointer array for each pass
C                .FALSE.  Only calculate number of passes
C     - ABSPOI : .TRUE.   Absolute pointers in INDXAB(,,4)
C                .FALSE.  Pointers within pass in INDXAB(,,4)
C     - INDX   : Pointer from complete list to position in the shell
C     - NWDS   : Maximum number of words for symmetry packed int batch
C
C     Output :
C
C     - INDXAB :  Pointer array for symmetry packing of distributions
C     - IJPAIRS : Number of distributions for each boson symmetry
C     - NWDS :    Number of words required to store symmetry packed ints
C     - NPASS :   Number of passes for this distribution set (kept in
C                 disbuf)
C     - IJPASS  : Like IJPAIRS, but for each batch separately
C
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <aovec.h>
#include <blocks.h>
      LOGICAL SHAEQB1,DCMPAB,MKINAB,MKIJPA,DIAGAB1,ABSPOI
      DIMENSION INDX(3,*),INDXAB(NINSHA,NINSHB,5),IJPAIRS(0:7)
      DIMENSION IJPASS(0:7,*)
      DIMENSION IND(2)
#include <twocom.h>
#include <dgroup.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <symmet.h>
#include <dcbtra.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 3) CALL HEADER('Subroutine NDISTRN',-1)
      IF (MKINAB.AND.(.NOT.ABSPOI).AND.(.NOT.MKIJPA)) 
     &    CALL QUIT ('Conficting options in NDISTRN')
C
      IF (MKINAB) CALL IZERO(INDXAB,NINSHA*NINSHB*5)
      CALL IZERO(IJPAIRS,8)
      IF (MKIJPA) CALL IZERO(IJPASS,NPASS*8)
      MWDS  = 0
      MPASS = NPASS
      NPASS = 1
C
      IF (I2TYP.EQ.1) THEN
         IC = 1
      ELSE
         IC = 0
      ENDIF
C
      IOFF = 0
      IJ = 0
C
      ISHELA1 = ISHLA
      ISHELB1 = ISHLB
      NHKTA1 = NHKTSH(ISHELA1)
      KHKTA1 = KHKTSH(ISHELA1)
      MULA1  = ISTBSH(ISHELA1)
      NORBA1 = NORBSH(ISHELA1)
      NSTRA1 = IORBSB(IORBSH(ISHELA1,1))
      NHKTB1 = NHKTSH(ISHELB1)
      KHKTB1 = KHKTSH(ISHELB1)
      MULB1  = ISTBSH(ISHELB1)
      NORBB1 = NORBSH(ISHELB1)
      NSTRB1 = IORBSB(IORBSH(ISHELB1,1))
      SHAEQB1 = ISHELA1 .EQ. ISHELB1 
      DIAGAB1 = SHAEQB1 .AND. .NOT.BIGVEC
C
      DO IREPAB = 0, NBSYM-1
C
C     Loop over components
C
         DO NA = 1,KHKTA1
            NSTRNA = NSTRA1 + NA
            ITYNA  = ISYMAO(NHKTA1,NA)
            KHKTBB = KHKTB1
            IF (DIAGAB1) KHKTBB = NA
            DO NB = 1,KHKTBB
               NSTRNB = NSTRB1 + NB
               ITYNB  = ISYMAO(NHKTB1,NB)
               DCMPAB = SHAEQB1 .AND. NA .EQ. NB
C
C              Loop over symmetries
C
               DO IREPA = 0, MAXREP
               IF (IBTAND(MULA1,IBTXOR(IREPA,ITYNA)) .EQ. 0) THEN
                  DO IREPB = 0, MAXREP
                  IF (IBTAND(MULB1,IBTXOR(IREPB,ITYNB)) .EQ. 0) THEN
C
C                 This is the sorting step : the compound symmetry 
C                 should match the one in the highest loop
C
                     IF (IREPAB.EQ.IBTXOR(IREPA,IREPB)) THEN
C
C                    Loop over contracted orbitals
C
                        DO IA = 1,NORBA1
C
C                          INDA  is the position in the complete list
C                          INDAR is the position in the shell
C
                           INDA = IPTSYM(NSTRNA + KHKTA1*(IA-1),IREPA)
                           IND(1) = INDA
                           INDAR = INDX(1,INDA)
                           NORBBB = NORBB1
                           IF (DCMPAB) NORBBB = IA
                           IF (DCMPAB.AND.IREPA.LT.IREPB) NORBBB = IA-1
                           DO IB = 1,NORBBB
                              INDB = IPTSYM(NSTRNB+KHKTB1*(IB-1),IREPB)
                              INDBR = INDX(1,INDB)
                              IND(2) = INDB
                              INDAB = INPCK(2,IND)
C
C                             Check whether we should start a new batch
C
                              IF (IOFF+NSPCK(IREPAB,IC).GT.NWDS) THEN
C
                                 IF (IOFF.EQ.0) CALL QUIT (
     &                           ' Maximum size of GMAT too small ')
C
C                                If this happens one should either
C                                decrease the shell size or increase the
C                                total memory
C
                                 MWDS = MAX0(IOFF,MWDS)
                                 NPASS = NPASS + 1
                                 IOFF = 0
                              ENDIF
C
                              IJ = IJ + 1
                              IJPAIRS(IREPAB) = IJPAIRS(IREPAB) + 1
                              IF (MKIJPA) IJPASS(IREPAB,NPASS) =
     &                                    IJPASS(IREPAB,NPASS) + 1
C
                              IF (MKINAB) THEN
                                 INDXAB(INDAR,INDBR,1) = INDAB
                                 INDXAB(INDAR,INDBR,2) = IOFF
                                 INDXAB(INDAR,INDBR,3) = IREPAB
                                 IF (ABSPOI) THEN
                                    INDXAB(INDAR,INDBR,4) = 
     &                              IJPAIRS(IREPAB)
                                 ELSE
                                    INDXAB(INDAR,INDBR,4) = 
     &                              IJPASS(IREPAB,NPASS)
                                 ENDIF
                                 INDXAB(INDAR,INDBR,5) = NPASS
                              ENDIF
                              IOFF = IOFF + NSPCK(IREPAB,IC)
                            ENDDO
                        ENDDO
                     ENDIF
                  ENDIF
                  ENDDO
               ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      NWDS = MAX0(IOFF,MWDS)
      NDISTRN = IJ
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck TRASAM */
      LOGICAL FUNCTION TRASAM(IOS11,IOS21,IOS12,IOS22,NSTR1,NSTR2)
C***********************************************************************
C
C     Check relations between indices that are to be transformed
C
C     Written by L.Visscher January 1997
C
C***********************************************************************
******
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcborb.h>
      DIMENSION IOS11(*),IOS12(*),IOS21(*),IOS22(*)
      DIMENSION NSTR1(2,0:2),NSTR2(2,0:2)
C
C     Check 1-2 
C
      IF (NSTR1(1,0).EQ.NSTR2(1,0)) THEN
         ISAME = 1
         DO IS = 1, NSTR1(1,0)
            IF (IOS12(IS).NE.IOS11(IS)) ISAME = 2
         ENDDO
         IF (NFSYM.GE.2) THEN
            IF (NSTR1(2,0).EQ.NSTR2(2,0)) THEN
               DO IS = 1, NSTR1(2,0)
                  IF (IOS12(IS).NE.IOS11(IS)) ISAME = 2
               ENDDO
            ELSE
               ISAME = 2
            ENDIF
         ENDIF
      ELSE
        ISAME = 2
      ENDIF
C
      TRASAM = ISAME.EQ.1
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck INDSAM */
      SUBROUTINE INDSAM(IOS1,IOS2,N1,N2,IBUF,NBUF,IP,NSAM)
C***********************************************************************
C
C   Count the number of elements in IOS2 identical to those of
C   IOS1.
C
C   Written by T.Saue Feb 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      LOGICAl LBUF
      DIMENSION IOS1(*),IOS2(*),IBUF(*)
C
      CALL IZERO(IBUF,NBUF)
      DO I = 1,N1
        IBUF(IP*IOS1(I)) = 1
      ENDDO
C
      DO I = 1,N2
        NSAM = NSAM + IBUF(IP*IOS2(I))
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEFCL*/
      SUBROUTINE DEFCL (NCLASS)
C
C     Define integral classes
C
C     Luuk Visscher
C
C     Input :  Integral class, quaternion units
C
C     Output : Contribution to class
C
C     Classes are defined such that the complete list can be generated
C     without using permutation symmetry, except for the quaternion
C     case where we assume particle-particle symmetry.
C
C     Real groups       : classes 1-4 are used
C     Complex groups    : classes 1-8 are used
C     Quaternion groups : classes 1-16 are used
C
C     Class 1  : Real part of integrals      (i    j    |k    l   )
C     Class 2  : Real part of integrals      (i    j    |kbar lbar)
C     Class 3  : Real part of integrals      (i    jbar |k    lbar)
C     Class 4  : Real part of integrals      (i    jbar |kbar l   )
C     Class 5  : Imaginary part of integrals (i    j    |k    l   )
C     Class 6  : Imaginary part of integrals (i    j    |kbar lbar)
C     Class 7  : Imaginary part of integrals (i    jbar |k    lbar)
C     Class 8  : Imaginary part of integrals (i    jbar |kbar l   )
C     Class 9  : Real part of integrals      (i    j    |kbar l   )
C     Class 10 : Real part of integrals      (i    j    |k    lbar)
C     Class 11 : Real part of integrals      (ibar j    |k    l   )
C     Class 12 : Real part of integrals      (i    jbar |k    l   )
C     Class 13 : Imaginary part of integrals (i    j    |kbar l   )
C     Class 14 : Imaginary part of integrals (i    j    |k    lbar)
C     Class 15 : Imaginary part of integrals (ibar j    |k    l   )
C     Class 16 : Imaginary part of integrals (i    jbar |k    l   )
C
#include <implicit.h>
#include <dcbtra.h>
#include <dgroup.h>
C
      CALL IZERO (ICLARR,16*NCLASS)
C
      IF (NCLASS.GT.16) CALL QUIT ('Only 16 integral classes defined')
C
C     Real classes
C
      ICLARR(1,1,1) =  1
      ICLARR(2,2,1) = -1
      ICLARR(1,1,2) =  1
      ICLARR(2,2,2) =  1
      ICLARR(3,3,3) =  1
      ICLARR(4,4,3) = -1
      ICLARR(3,3,4) = -1
      ICLARR(4,4,4) = -1
C
      IF (NCLASS.LE.4) RETURN
C
C     Imaginary classes
C
      ICLARR(1,2,5) =  1
      ICLARR(2,1,5) =  1
      ICLARR(1,2,6) = -1
      ICLARR(2,1,6) =  1
      ICLARR(4,3,7) =  1
      ICLARR(3,4,7) =  1
      ICLARR(4,3,8) = -1
      ICLARR(3,4,8) =  1
C
      IF (NCLASS.LE.8) RETURN
C
C     Quaternion classes
C
      ICLARR(1,3,9)  = -1
      ICLARR(2,4,9)  = -1
      ICLARR(1,3,10) =  1
      ICLARR(2,4,10) = -1
      ICLARR(3,1,11) = -1
      ICLARR(4,2,11) = -1
      ICLARR(3,1,12) =  1
      ICLARR(4,2,12) = -1
      ICLARR(2,3,13) = -1
      ICLARR(1,4,13) =  1
      ICLARR(1,4,14) =  1
      ICLARR(2,3,14) =  1
      ICLARR(3,2,15) = -1
      ICLARR(4,1,15) =  1
      ICLARR(4,1,16) =  1
      ICLARR(3,2,16) =  1
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck INTCL*/
      FUNCTION INTCL (ICLASS,IQ1,IQ2)
C
C     Give contribution of the quaternion combination to the class
C
C     Luuk Visscher
C
#include <implicit.h>
#include <dcbtra.h>
C
      INTCL = ICLARR(IQ1,IQ2,ICLASS)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck QTRIAN*/
      SUBROUTINE QTRIAN (TASK,FACTOR,N,NZ,ASQ,NROW,NCOL,ATRI,NDIM)
C
C     Pack or unpack symmetric quaternion matrix to lower triangular 
C     form. Add result to existing matrix multiplied by factor
C
C     Luuk Visscher
C
#include <implicit.h>
C
      CHARACTER*4 TASK
      DIMENSION ASQ(NROW,NCOL,NZ),ATRI(NDIM,NZ)
C
      IF (TASK.EQ.'PACK') THEN
         DO IZ = 1, NZ
            IJ = 0
            DO J = 1, N
               DO I = J, N
                  IJ = IJ + 1
                  ATRI(IJ,IZ) = ATRI(IJ,IZ) * FACTOR + ASQ(I,J,IZ)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO IZ = 1, NZ
            IJ = 0
            DO J = 1, N
               DO I = J, N
                  IJ = IJ + 1
                  ASQ(I,J,IZ) = ASQ(I,J,IZ) * FACTOR + ATRI(IJ,IZ)
                  ASQ(J,I,IZ) = ASQ(J,I,IZ) * FACTOR + ATRI(IJ,IZ)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PCKTOP(EIGTHR,NTHREI,IPRINT,WORK,LWORK)
C***********************************************************************
C     Finds the position of the lagest element below
C     the threshold EIGTHR in the list of eigenvalues.
C
C     Written by Jon K. Laerdahl 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbgen.h>
      LOGICAL TOBE
      CHARACTER*18 STR(2)
      DIMENSION WORK(*),NTHREI(2)
C
#include <memint.h>
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT ,WORK,KFREE,LFREE)
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(A)') 'Coefficient file not present'
         CALL QUIT('PCKTOP: No coefficients !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PSIMP2')
         CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),TOTERG,1)
         CLOSE(LUCOEF,STATUS='KEEP')
      ENDIF
      DO IFRP = 1,NFSYM
         IOFF = KEIG+NPSH(IFRP)+IORB(IFRP)-1
         ISTR = NISH(IFRP)+1
         IEND = NESH(IFRP)
         NTHREI(IFRP) = -1
         EIGP = -1.0D20
         DO I = ISTR,IEND
            IF(WORK(IOFF+I).LT.EIGP) THEN
               WRITE(LUPRI,'(A)') 
     &              'Not INCREASING order in Eigenvalues!'
               CALL QUIT('PCKTOP: Eigenvalues Decreasing with N!')
            ENDIF
            IF(WORK(IOFF+I).LT.EIGTHR) NTHREI(IFRP) = I
            EIGP = WORK(IOFF+I)
         ENDDO
      ENDDO
      IF(IPRINT.GE.1) THEN
         STR(1) = 'Fermion irrep. 1  '
         STR(2) = 'Fermion irrep. 2  '
         WRITE(LUPRI,'(/X,A)') '* Eigenvalue of highest virtual'//
     &        ' Kramers pair included:'
         
         WRITE(LUPRI,'(8X,2A)') (STR(I),I=1,NFSYM)
         WRITE(LUPRI,'(4X,2F18.8)') 
     &        (WORK(KEIG+NPSH(I)+IORB(I)+NTHREI(I)-1),I=1,NFSYM)
         WRITE(LUPRI,'(/X,A)') '* Position of highest virtual'//
     &        ' Kramers pair included:'
         WRITE(LUPRI,'(8X,2A)') (STR(I),I=1,NFSYM)
         WRITE(LUPRI,'(8X,2I18)') 
     &        (NTHREI(I),I=1,NFSYM)
      ENDIF
      CALL MEMREL('PCKTOP',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END






