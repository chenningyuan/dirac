CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C    '  Dobrorec dusa moja Hospodinu, a cele moje vnutro 
C     Jeho svatemu menu !' (Zalm 103,1)
C
C                Na slavu Boziu v Kristu Jezisu !
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck make_bss */
C****************************************************************************
C
C On input : VMAT - Lowdin transformation matrix for 4-comp. Fock matrix
C
C On output: VMAT - Lowdin transformation matrix modified for 2-comp. Fock matrix
C
C Written by Miro Ilias & Hans Jorgen Aa Jensen, 2002, 2004
C
C****************************************************************************
       SUBROUTINE MAKE_BSS(VMAT,WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbham.h>

        DIMENSION VMAT(*),WORK(LWORK)
        CALL QENTER('MAKE_BSS')
#include "memint.h"

        CALL MEMGET('REAL',KH1DK,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KTMAT,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFACKIN,2*NTBAS(1),WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KPCTM,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIG ,NTBAS(0),WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBEIG,NTBAS(0),WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KTBUF1,N2BBASXQ,WORK,KFREE,LFREE)

        CALL MAKE_DKH(WORK(KH1DK),WORK(KTMAT),VMAT,WORK(KTBUF1),
     &      WORK(KPCTM),
     &      WORK(KFACKIN),WORK(KEIG),WORK(KIBEIG),WORK(KFREE),LFREE)

C******************************************************************************
C    HJAaJ Aug 2002: we need the Lowdin VMAT for the L block as TMAT
C    for the BSS two-component Hamiltonian  !!!
C******************************************************************************
CMI ... print out the whole Lowdin matrix !
        IF (IPRHAM.GE.5) THEN
           CALL HEADER(
     &'GMOTRA: The total Lowdin matrix before BSS reordering',-1)
           IVMAT=0
           DO I=1, NFSYM
            WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &      '*** Fermion corep ',I,'/',NFSYM
            CALL PRQMAT(VMAT(1+IVMAT),
     &          NFBAS(I,0),NFORB(I,0),NFBAS(I,0),NFORB(I,0),
     &          1,IPQTOQ(1,0),LUPRI)
            IVMAT = IVMAT + NFBAS(I,0)*NFORB(I,0)
           ENDDO
        ENDIF
       IF (NFSYM.EQ.2) THEN
C        For TMAT: copy the VMAT for the L block ifsym=2 to
C        follow immediately after the L block ifsym=1.
           NDIM2 = NFBAS(2,0)*NORB(2)
           JV2   = 1 + NFBAS(1,0)*NFORB(1,0)
           JT2   = 1 + NFBAS(1,0)*NORB(1)
           CALL DCOPY(NDIM2,VMAT(JV2),1,VMAT(JT2),1)
        ENDIF
C ... Number of transf. matrixes is 1 - it is only the Lowdin AO->MO
        NZT = 1 ! dcbbas.h

C ... Transformation matrix is the Lowdin Le matrix
C       KTMAT = KVMAT
        CALL DCOPY(N2BBASXQ,VMAT,1,WORK(KTMAT),1)

C Print out transf. matrixes
        IF (IPRHAM.GE.6) THEN
C        ITMAT = KTMAT
         ITMAT = 0
         IH1D  = KH1DK
         IPCTM = KPCTM
         DO I = 1, NFSYM
          CALL HEADER(
     &'GMOTRA: LL 2c bare nuclei BSS-Hamiltonain for SCF procedure',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)') '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(WORK(KH1DK+I2BASX(I,I)),
     &               NFBAS(I,0),NFBAS(I,0),
     &       NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)

           CALL HEADER(
     &'GMOTRA: Lowd. transf. matrix (Le) for 2c SCF procedure:',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &              '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(1+ITMAT),NFBAS(I,0),NORB(I),NFBAS(I,0),
     &                 NORB(I),NZT,IPQTOQ(1,0),LUPRI)

          IH1D  = IH1D +  NFBAS(I,0)*NFBAS(I,0)*NZ
          IPCTM = IPCTM + NFBAS(I,0)*NFBAS(I,1)*NZ
          ITMAT = ITMAT + NFBAS(I,0)*NORB(I)*NZT
         ENDDO
          CALL HEADER(
     &'GMOTRA: TOTAL pict.change transformation matrix:',-1)
          CALL PRQMAT(WORK(KPCTM),NTBAS(0),NTBAS(1),
     &           NTBAS(0),NTBAS(1),NZ,IPQTOQ(1,0),LUPRI)
      ENDIF

      IF (IPRHAM.GE.5) THEN
       IF (NFSYM.EQ.2) THEN
         NORBTX=NORB(1)+NORB(2)
       ELSE
         NORBTX=NORB(1)
       ENDIF
       CALL HEADER(
     &'GMOTRA: Final transformation matrix',-1)
        CALL PRQMAT(WORK(KTMAT),NTBAS(0),NORBTX,NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
      ENDIF

         LUBSS=35 ! Assign (and check) this LU number!!!
C ...save hamiltonian matrix elements
         CALL OPNFIL(LUBSS,'BSSMAT','UNKNOWN','GMOTRA')
         CALL NEWLAB('H1AO_DK ',LUBSS,LUPRI)
         CALL WRITT(LUBSS,N2BBASXQ,WORK(KH1DK))
C ...save the picture change transformation matrix
         CALL NEWLAB('U_PICTRM',LUBSS,LUPRI)
         CALL WRITT(LUBSS,N2BBASXQ,WORK(KPCTM))
         CLOSE(LUBSS,STATUS='KEEP')
         WRITE(LUPRI,'(/2X,A,/,2X,A)')
     &'GMOTRA: BSS relativistic atomic integrals and picture change',
     &'transformation matrix were written to the new file BSSMAT.'

C   ....      Memory deallocation
       CALL MEMREL('MAKE_BSS',WORK,KWORK,KWORK,KFREE,LFREE)
       CALL QEXIT('MAKE_BSS')
       RETURN
       END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck make_dkh */
C ........................................................  
      SUBROUTINE MAKE_DKH(H1AO,TMAT,VMAT,TBUF,PCTM,FACKIN,
     &                    EIG,IBEIG,WORK,LWORK)
C****************************************************************************
C
C     PURPOSE: Tranformation of the Dirac Hamiltonian into 
C              the two-component form - DK1, DK2 or BSS_infinite
C 
C              Also one can eliminate spin-orbit components 
C
C    On input:
C              VMAT = transf.matrix AO "xhi" -> ON Lowdin basis "Theta"
C
C      Output:
C            TBUF = Final transformation matrix from 4c 'xhi' AO basis 
C                    to the desired 2c  'xhi' one 
C
C            H1AO = matrix elements of transformed 2c Hamiltonian are in LL block,
C                    while LS,SL and SS blocks are zero !
C
C     Written:  HJAa Jensen, M Ilias, Odense, 2002  
C
C IPRHAM .ge. JCONSIST  -  all consistency checks inluded
C
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0, DM1 = -1.0D0, JCONSIST = 2,
     &   D2 = 2.0D0,  D0 = 0.0D0, THRNULL = 1.0D-10)
C     THRNULL is round off error with 15 digits relative to 2c^2 = 3.7D4
C
C Used from COMMON blocks
C  CBIPAM:
C  DHFCAL:
C  DCBBAS:
C  DCBORB:
C  DGROUP:
#include <dcbgen.h>
#include <dcbham.h>
C ... for NTMO...
#include <dcbdhf.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <pgroup.h>
#include <cbihr1.h>
#include <codata.h>
C
      DIMENSION H1AO(*), TMAT(*),TBUF(*),VMAT(*),EIG(*),IBEIG(*)
      DIMENSION FACKIN(2,*), PCTM(*), WORK(LWORK)
      LOGICAL   DOJACO_SAVE, SPINFR_SAVE, SF0, DK2SO1, DBGDIAG
      CALL QENTER('MAKE_DKH')
#include <memint.h> 

C
C     Print section
C
      CALL HEADER('*** Output from MAKE_DKH ***',-1)
      IF (.NOT. BSS)CALL QUIT('ERROR: BSS was false!')

      BSS = .FALSE.
C     ... we must set the BSS false to get 4-component one-electron
C         matrix from ONEFCK below

      WRITE(LUPRI,'(A)')
     & '* MAKE_DKH: Applied the restricted kinetic balance !'

      IF (IPRHAM.GE.JCONSIST) WRITE(LUPRI,'(/3X,A,I3/)')
     &'MAKE_DKH: All consistency checks are included for IPRHAM=',
     &IPRHAM

      NSYM = 4/NZ

C extract information about the BSS Hamiltonian 
      CALL EXTR_BSS_INFO(IBSS,IPZ,SPINFR,DK2SO1,SF0,IPRHAM)

C generate free particle one-electron matrix in Lowdin-basis "Theta"
      CALL GET_FP_LOWD(TBUF,TMAT,VMAT,IPRHAM,WORK(KFREE),LFREE)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     The modified Dirac equation:
C  ==================================
C
C     Basis as in Levy-Leblond but without fixed ratio
C     (JCP-qmd eq. 22)
C
C                              [ Theta(L), 0                         ]
C     [ theta(L), theta(S) ] = [                                     ]
C                              [        0, (2mc)^-1 sigma.p Theta(L) ]
C
C  It's the restricted kinetic balance case.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        NZT  = NZ ! Store here number of transform. matrixes....
        IMAT = 1
        IBUF = 1
        N2TMOTQ = 0
        DO 40 I = 1,NFSYM
          NESH(I)  = NFORB(I,1) 
          NPSH(I)  = NFORB(I,1) ! Charge conjugation symmetry
          NORB(I)  = NESH(I) + NPSH(I)
CMI ... these variables are needed  for DFDIAG routine!!!
          NTMO(I)  = NORB(I)
          I2TMOT(I) = N2TMOTQ
          N2TMOTQ = N2TMOTQ + NTMO(I)*NTMO(I)*NZ

          IF(NORB(I).EQ.0) GOTO 40
C
C         Prepare small component part
C         ----------------------------
C
C         Make <theta(S)|theta(S)> = c^2 (sigma.p)(sigma.p) = c^2 p^2
C         (JCPqmd eq. 23, but our theta(S) = 2c^2 theta(S) in eq. 22)
C         and do a new Lowdin canonical orthonormalization
C         (not described in JCPqmd)
C
          LWRK = MXBBAS*MXBBAS*NSYM
          CALL MEMGET('REAL',KSMT,LWRK,WORK,KFREE,LFREE)
          LWRK = NFORB(I,1) * NFORB(I,1)
          CALL MEMGET('REAL',KVMT,LWRK,WORK,KFREE,LFREE)
          IVMT  = KVMT
          IBOL  = 0
          NEFF  = 0
          DO ISYM = 1,NSYM
          IF(NBORB(ISYM,I,1).GT.0) THEN
            IREP = JFSYM(ISYM,I) - 1
            JSL = IMAT + IBOL * NFORB(I,0) + NFORB(I,1)
            JLS = IMAT + NFORB(I,1) * NFORB(I,0) + IBOL
            CALL QGEMM(NBORB(ISYM,I,1),NBORB(ISYM,I,1),NFORB(I,2),D1,
     &          'N','N',IPQTOQ(1,0),TMAT(JLS),NFORB(I,0),NFORB(I,0),NZ,
     &          'N','N',IPQTOQ(1,0),TMAT(JSL),NFORB(I,0),NFORB(I,0),NZ,
     &          D0,IPQTOQ(1,0),WORK(KSMT),
     &          NBORB(ISYM,I,1),NBORB(ISYM,I,1),1)
            IF (IPRHAM.GE.7) THEN
              WRITE(LUPRI,'(A,A3)') ' * Boson irrep: ',REP(IREP)
              CALL HEADER('Kinetic energy matrix in Lowdin L basis',-1)
              CALL PRQMAT(WORK(KSMT),NBORB(ISYM,I,1),NBORB(ISYM,I,1),
     &                    NBORB(ISYM,I,1),NBORB(ISYM,I,1),
     &                    1,IPQTOQ(1,0),LUPRI)
            ENDIF

C Do canonical orthonormalization of small components
            CALL LOWGEN(WORK(KSMT),NBORB(ISYM,I,1),WORK(IVMT),NEFFI,
     &                  NBORB(ISYM,I,1),STOL(2),IPRHAM,EIG,
     &                  WORK(KFREE),LFREE)
            IF (NEFFI.LT.NBORB(ISYM,I,1)) WRITE(LUPRI,'(A)') 
     &            '*** WARNING *** Linear dependence in RKB basis'
            NEFF = NEFF + NEFFI
            NBORB(ISYM,I,2) = NEFFI


             IF (IPRHAM.GE.7) THEN
              CALL HEADER(
     &'MAKE_DKH: Lowdin transf. matrix to ON modif. Dirac small basis',
     &        -1)
              CALL PRQMAT(WORK(IVMT),NBORB(ISYM,I,1),NBORB(ISYM,I,2),
     &                    NBORB(ISYM,I,1),NBORB(ISYM,I,2),
     &                    1,IPQTOQ(1,0),LUPRI)
             ENDIF
            IVMT = IVMT + NBORB(ISYM,I,1)*NBORB(ISYM,I,2)
            IBOL = IBOL + NBORB(ISYM,I,1)
          ELSE
            NBORB(ISYM,I,2) = 0
          ENDIF
          NBORB(ISYM,I,0) = NBORB(ISYM,I,1) + NBORB(ISYM,I,2)
          ENDDO

          NPSH(I)  = NEFF
          NORB(I)  = NESH(I) + NPSH(I)
C
C         Generate transformation from Lowdin basis to final basis in TBUF
C
          CALL DZERO(TBUF(IBUF),NFORB(I,0)*NORB(I)*NZ)
          IBO  = 0
          IBOL = 0
          IVMT = KVMT
          DO ISYM = 1,NSYM
          IF(NBORB(ISYM,I,1).GT.0) THEN
            IREP = JFSYM(ISYM,I) - 1
C
C           Small component part
C           TBUF = TMAT*VMT
C              where TMAT transforms Theta(S) to theta'(S) (non-orthonormal)
C              and VMT (WORK(IVMT)) transforms to orthonormal theta(S)
C              using Lowdins canonical transformation (from LOWGEN)
C
            JBUF = IBUF + NFORB(I,0)*IBO + NFORB(I,1)
            JSL  = IMAT + IBOL * NFORB(I,0) + NFORB(I,1)
            CALL QGEMM(NFORB(I,2),NBORB(ISYM,I,2),NBORB(ISYM,I,1),D1,
     &           'N','N',IPQTOQ(1,0),TMAT(JSL),NFORB(I,0),NFORB(I,0),NZ,
     &           'N','N',IPQTOQ(1,0),WORK(IVMT),
     &                   NBORB(ISYM,I,1),NBORB(ISYM,I,2),1,
     &           D0,IPQTOQ(1,0),TBUF(JBUF),NFORB(I,0),NORB(I),NZ)
            IBO = IBO + NBORB(ISYM,I,2)
C
C           Large component part
C           TBUF = unit matrix ( as theta(L) = Theta(L) basis )
C
            JBUF = IBUF + NFORB(I,0)*IBO + IBOL
            DO J = 1, NBORB(ISYM,I,1)
              TBUF(JBUF) = D1
              JBUF = JBUF + NFORB(I,0) + 1
            ENDDO
            IBO = IBO + NBORB(ISYM,I,1)
C
C           Update offsets
C           (because LOWGEN may have deleted som functions)
C
            IBOL = IBOL + NBORB(ISYM,I,1)
            IVMT = IVMT + NBORB(ISYM,I,1)*NBORB(ISYM,I,2)
          ENDIF
          ENDDO

          IF (IPRHAM.GE.8) THEN
           CALL HEADER(
     &'MAKE_DKH:Transf. matrix from Lowdin "Theta" to ON basis "theta"',
     &     -1)
           CALL PRQMAT(TBUF(IBUF),NFORB(I,0),NORB(I),
     &                 NFORB(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

          IMAT = IMAT + NFORB(I,0)*NFORB(I,0)*NZ
          IBUF = IBUF + NFORB(I,0)*NORB(I)*NZ
   40   CONTINUE

      CALL MEMREL('MAKE_DKH 1',WORK,KWORK,KSMT,KFREE,LFREE)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       Now TBUF contains the transformation from the Lowdin basis "Theta"
C          to the ON "theta" basis.
C
C       VMAT still contains the transformation from SO basis "xhi" to "Theta"
C
C       TMAT contains the free particle matrix in Lowdin basis "Theta".
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
C       Diagonalize this TMAT to get the Douglas-Kroll 4c |k> basis:
C       (the "exact" instead of Hess's so-called "SQR" approximation;
C       Hess states that the difference is small, and we h tested that
C       with a 25 s L basis set where we also got a small difference;
C       our model corresponds to Hess's "SQRA".)
C
      N2TMT  = 0
      NORBT  = 0
      N2ORBT = 0
      DO I = 1,NFSYM
        I2TMT(I) = N2TMT
        N2TMT    = N2TMT + NFORB(I,0)*NORB(I)*NZT
        IORB(I)  = NORBT
        NORBT    = NORBT  + NORB(I) ! Count total number of orbitals...
        I2ORBT(I)= N2ORBT*NZ
        N2ORBT   = N2ORBT + NORB(I)*NORB(I)
      ENDDO

      CALL MEMGET('REAL',KBUF,N2BBASXQ,WORK,KFREE,LFREE)
C
C         Transform the free particle matrix
C         from TMAT in Lowdin "Theta" basis to WORK(IBUF)
C         in modified Dirac form (orhonormal "theta" basis); 
C
      IMAT = 1
      IOFT = 1
      IOFV = 1
      DO 50 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 50
        IBUF = KBUF + I2ORBT(I)
        CALL QTRANS('AOMO','S',D0,NFORB(I,0),NFORB(I,0),
     &           NORB(I),NORB(I),
     &           TMAT(IOFT),NFORB(I,0),NFORB(I,0),NZ,IPQTOQ(1,0),
     &           WORK(IBUF),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &           TBUF(1+I2TMT(I)),NFORB(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &           TBUF(1+I2TMT(I)),NFORB(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRHAM)
C
C       Zero the small elements in real part of WORK(KBUF)
C       which are non-zero because of round-off errors:
C       check that quaternion blocks are numerically zero
C
        DO J = 0,NORB(I)*NORB(I)-1
          IF (ABS(WORK(IBUF+J)) .LT. THRNULL) WORK(IBUF+J) = D0
        END DO

        IF (NZ .GT. 1) THEN
           IBUF2 = IBUF + NORB(I)*NORB(I)
           NBUF2 = (NZ-1)*NORB(I)*NORB(I)
           TIMAGNORM =  DNRM2(NBUF2,WORK(IBUF2),1)
        ELSE
           TIMAGNORM =  D0
        END IF

C  Write out TIMAGNORM

         WRITE(LUPRI,'(/2X,A,A,I1,A,I1,/,2X,A,D12.6/)')
     & 'Free particle matrix in the ON basis "theta" ',
     & 'in fermion corep ',I,'/',NFSYM,
     & '*** Norm of the i,j,k imag. parts *** :',TIMAGNORM

        IF (TIMAGNORM.GT.THRNULL) 
     &    WRITE(LUPRI,'(/2X,A)')
     & 'MAKE_DKH: WARNING, TIMAGNORM > 1.0D-10 !'

        IF (IPRHAM.GE.7) THEN
          CALL HEADER(
     &'MAKE_DKH:Free particle matrix in the ON basis "theta"',-1)
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          WRITE(LUPRI,'(/A,1P,D12.6)')
     & ' Norm of the i,j,k imag. parts of this matrix:',TIMAGNORM

          CALL PRQMAT(WORK(IBUF),NORB(I),NORB(I),
     &                NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C--------------------------------------------------------------------
C         Get full transformation as well
C         from SO basis "xhi" to modified Dirac ON basis "theta"
C
C         TMAT = VMAT * TBUF
C         (TMAT is now reused)
C--------------------------------------------------------------------
        CALL BCKTRA(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &              TBUF(1+I2TMT(I)),NFORB(I,0),NORB(I),
     &              NFORB(I,0),NZ, NORB(I),1,NFBAS(I,0),
     &              VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1, IPRHAM)

        IF (IPRHAM.GE.7) THEN
          CALL HEADER(
     &'MAKE_DKH:Transf. mtx TMAT "xhi"->final ON basis "theta" ',-1)
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IOFT = IOFT + NFORB(I,0)*NFORB(I,0)*NZ
        IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
        IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)

 50   CONTINUE

*** Store the TMAT into WORK(KXT) ***

      CALL MEMGET('REAL',KXT,N2BBASXQ,WORK,KFREE,LFREE)
      CALL DCOPY(N2BBASXQ,TMAT(1),1,WORK(KXT),1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Diagonalize WORK(KBUF) with kinetic energy matrix
C      in the orthonormal "theta" basis
C
C      TBUF reused for |k> eigenvectors in ON "theta" basis
C
C      EIG - eigenvalues 
C        - used also for checking the A prefactors
C        - used for 'denominator' in DKH2 method R_DK2;  eps(el)-eps(pos)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      NZEKIN = 1 ! Td_kin is set to be pure real matrix !
      NZSAVE = NZ
      NZ = NZEKIN ! Set only real part of matrix for diagonalization

C     note: DFDIAG zeroes all of TBUF, even if NZ is reset to 1
C     (it uses N2BBASXQ); so we can continue with the quaternion TBUF afterwards.

      DOJACO_SAVE = DOJACO
      DOJACO = .TRUE.
C     ... use JACO because it does not mix degenerate eigenvectors
C  NOTE: All fermion symmetries are diagonalized !!

C ........  added from GMOTRA 
C     Copy information from NBORB to the subblock array if we have not
C     partioned the orbital space otherwise.
C
C ... we need this variables because of later calling of the DFDIAG
      IF (.NOT.SUB_BL) THEN
         DO IFRP = 1, NFSYM
            N_SUB_BL(IFRP) = 4 / NZ
            DO IB = 1, N_SUB_BL(IFRP)
              NORB_SUB(IB,IFRP,0) = NBORB(IB,IFRP,0)
              NORB_SUB(IB,IFRP,1) = NBORB(IB,IFRP,1)
              NORB_SUB(IB,IFRP,2) = NBORB(IB,IFRP,2)
              IREP = JFSYM(IB,IFRP) - 1
              ID_SUB_BL(IB,IFRP) = IREP
            ENDDO
         ENDDO
      ENDIF

C     ... to tell DFDIAG to use boson symmetries
      SUB_BL = .TRUE.

      IF (IPRHAM.GE.2) THEN
       WRITE(LUPRI,'(/,3X,A,A/)') 
     & 'MAKE_DKH: Diagonalizing (Jacobi) the real part',
     & ' of the FP matrix in the ON "theta" basis !'
      ENDIF

C ... after diagonalization we get the IBEIG array which is used later !
      DBGDIAG=.FALSE.
      IF (IPRHAM.GE.5) DBGDIAG=.TRUE.
      CALL DFDIAG (WORK(KBUF),EIG,IBEIG,TBUF,DBGDIAG,
     &             WORK,KFREE,LFREE)

      NZ = NZSAVE
      DOJACO = DOJACO_SAVE
      SUB_BL = .FALSE.
C*******************************************************************
C
C
C     ******    Loop over fermion symmetries   ******
C
C
C*******************************************************************

      CALL MEMGET('REAL',KKUN,N2BBASXQ,WORK,KFREE,LFREE)
      CVAL2 = CVAL*CVAL

      N2TMT = 0
      IMAT = 1
      IBUF = KBUF
      DO 60 I = 1,NFSYM
       IF(NORB(I).EQ.0) GOTO 60
       I2TMT(I) = N2TMT
       N2TMT = NFBAS(I,0)*NORB(I)*NZ

       CALL DZERO(WORK(KKUN),N2BBASXQ)

C   ... do always check the RKB free-particle eigenvalues !
        IF (IPRHAM.GE.2) THEN
         CALL HEADER(
     &'MAKE_DKH: Free particle (in ON basis "theta") eigenvalues :',-1)
         WRITE(LUPRI,'(3X,A,I2,A,I2)') '*** Fermion corep ',I,'/',NFSYM
         WRITE(LUPRI,'(A,3X,A,2X,A,2X,A,3X,A)')
     &   'Shell no ','FP eigval',
     &   'FP eigval+2c^2','bos.symm.','diff.'
        ENDIF

        TOTDIFF  = D0
        ISHELL = -NPSH(I)
        DO J=1,NORB(I) 
          FPEV   =  EIG(IORB(I)+J)
          FPEV2C = EIG(IORB(I)+J)+CVAL2
C .... this is the boson symmetry of the free-particle eigenshell
          IBS    = IBEIG(IORB(I)+J)+1
C ...  positronic and electronic free particle eigenvalues differ by 2c^2
          DIFF = EIG(IORB(I)+J)+EIG(IORB(I)+NORB(I)-J+1)+CVAL2+CVAL2
          TOTDIFF = TOTDIFF + DABS(DIFF)
          IF (ISHELL.EQ.0) ISHELL = 1
          
C   Print free particle eigenvalues with differences e(+)-e(i)
          IF (IPRHAM.GE.2) THEN
           WRITE(LUPRI,'(I4,1X,I5,2D14.6,3X,I2,3X,D10.4)')
     &             J,ISHELL,FPEV,FPEV2C,IBS,DIFF
          ENDIF
          ISHELL = ISHELL + 1
        ENDDO
        TOTDIFF = TOTDIFF/DFLOAT(NORB(I)/2)
C ... do always print the average difference for RKB free-particle e-p pairing !
        WRITE(LUPRI,'(/,3X,A,/,A,I1,A,I1,4X,D12.6)')
     & 'MAKE_DKH: The average difference for e-p pairing of all',
     & 'free particle eigenvalues in fermion symmetry ',I,
     & '/',NFSYM,TOTDIFF
       IF (TOTDIFF.GE.10D0) THEN
         WRITE(LUPRI,'(A)') 'MAKE_DKH: Big difference! Check DFDIAG !'
         CALL QUIT('MAKE_DKH:Big difference ! Check DFDIAG !')
       ENDIF

       IF(IPRHAM.GE.7) THEN
        CALL HEADER(
     &  'MAKE_DKH:DKH |k> matrix TBUF in the ON basis "theta" ',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(TBUF(IMAT),NORB(I),NORB(I),
     &                NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C-------------------------------------------------------------------------------
C         Backtransform ELECTRONIC and POSITRONIC solutions to SO-basis "xhi"
C         (we later do the no-pair projection to the electronic solutions after
C         we have done the DKH folding of the nuclear attraction into the
C         electronic part).
C
C         TMAT = TMAT*TBUF (WORK(IBUF) reused as scratch area)
C-------------------------------------------------------------------------------
        CALL BCKTRA(WORK(IBUF),NFBAS(I,0),NORB(I),
     &              TBUF(IMAT),NORB(I),NORB(I),
     &              NORB(I),NZ,NORB(I),1,NFBAS(I,0),
     &              TMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),NZT,
     &              IPRHAM)
        CALL DCOPY(N2TMT,WORK(IBUF),1,TMAT(I2TMT(I)+1),1)

        IMAT = IMAT + NORB(I)*NORB(I)*NZ ! Set up offset !!!
        IBUF = IBUF + IMAT  ! Upgrade the IBUF !

 60   CONTINUE  ! End of cycle

C
C     Now we have the 4c SO "xhi" -> 4c |k> DKH basis transformation in
C     TMAT ...
C

C   Obtain overlap integrals in SO basis 'xhi' into WRK(KBUF)

      CALL DZERO(WORK(KKUN),N2BBASXQ) ! Store here Ad^2 diagonal matrix

      IF(IPRHAM.GE.JCONSIST) THEN
        CALL GTOVLT(WORK(KBUF),D1,WORK(KFREE),LFREE)
        IMAT = 1
        IBUF = KBUF
        IKUN1 = KKUN
        DO 80 I = 1,NFSYM
        IF (NORB(I).EQ.0) GOTO 80

         IF(IPRHAM.GE.7) THEN
          CALL HEADER(
     &  'MAKE_DKH:DKH |k> matrix TMAT in SO basis "xhi" ',-1)
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',I
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     &  'MAKE_DKH: Full SOVLP matrix in SO basis "xhi" ',-1)
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',I
          CALL PRQMAT(WORK(IBUF),NFBAS(I,0),NFBAS(I,0),
     &                NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),LUPRI)
         ENDIF

**  DO C+ S C into WORK(KKUN)
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &       NORB(I),NORB(I),
     &       WORK(IBUF),NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),
     &       WORK(IKUN1),NORB(I),NORB(I),1,IPQTOQ(1,0),
     &       TMAT(IMAT),NFBAS(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &       TMAT(IMAT),NFBAS(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)

         CALL HEADER(
     &  'MAKE_DKH: Full unit <k|k>=C+.S.C matrix  ',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')
     & '* Fermion irp no.',I,'/',NFSYM
          IF(IPRHAM.GE.7) THEN
           CALL PRQMAT(WORK(IKUN1),NORB(I),NORB(I),
     &                NORB(I),NORB(I),1,IPQTOQ(1,0),LUPRI)
          ENDIF

          CALL TEST_UNIT_MTX(WORK(IKUN1),NORB(I),NORB(I),NORB(I),1)

          IKUN1 = IKUN1 + NORB(I)*NORB(I)
          IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
          IBUF = IBUF + NFBAS(I,0)*NFBAS(I,0)

 80     CONTINUE
      ENDIF

*******************************
***  To obtain the A factors:
*******************************

        CALL GTOVLT(WORK(KBUF),D0,WORK(KFREE),LFREE)
        CALL DZERO(WORK(KKUN),N2BBASXQ) ! Store here Ad^2 diagonal matrix

        IMAT = 1
        IFC = 0
        IBUF = KBUF
        IKUN1 = KKUN
        DO 82 I = 1,NFSYM
        IF (NORB(I).EQ.0) GOTO 82

** DO C+(e) S(LL) C(e) into WORK(KKUN) **

        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &       NESH(I),NESH(I),
     &       WORK(IBUF),NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),
     &       WORK(IKUN1),NESH(I),NESH(I),1,IPQTOQ(1,0),
     &       TMAT(IMAT+NFBAS(I,0)*NPSH(I)),
     &       NFBAS(I,0),NESH(I),NZT,IPQTOQ(1,0),
     &       TMAT(IMAT+NFBAS(I,0)*NPSH(I)),
     &       NFBAS(I,0),NESH(I),NZT,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)

         IF(IPRHAM.GE.JCONSIST) THEN
          CALL HEADER(
     &    'MAKE_DKH: Diagonal <k_e|k_e>=A+.A matrix',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')'* Fermion irp no.',I,'/',NFSYM
          CALL TEST_DIAG_MTX(WORK(IKUN1),NESH(I),NESH(I),NESH(I),1)
          IF(IPRHAM.GE.7) THEN
           CALL PRQMAT(WORK(IKUN1),NESH(I),NESH(I),
     &            NESH(I),NESH(I),1,IPQTOQ(1,0),LUPRI)
          ENDIF
         ENDIF

*************************************************************
*        Extract the A kinematic factors into FACKIN(1,*)   *
*************************************************************

          WRITE(LUPRI,'(/9X,A,I1,A,I1,//,12X,A/)')
     & 'MAKE_DKH: "A" kinem. electr. factors; ferm irp:',I,'/',NFSYM,
     & 'A^2_ovlp          A^2_calc          diff;'

         SUMDIFF = D0
         DO J = 1, NESH(I)
          DAA = WORK(IKUN1+(J-1)*NESH(I)+J-1)
          FACKIN(1,IFC+J) = DAA
          X1 = (EIG(IORB(I)+NESH(I)+J)/CVAL2) + D2
          X2 = D2*(EIG(IORB(I)+NESH(I)+J)/CVAL2) + D2
          XD = (X1/X2)
          DIFF = FACKIN(1,IFC+J)-XD
          SUMDIFF = SUMDIFF + DABS(DIFF)
          WRITE(LUPRI,"(5x,2F18.10,4X,D11.5)")
     &         FACKIN(1,IFC+J),XD,DIFF
         ENDDO 
         SUMDIFF = SUMDIFF/DFLOAT(NESH(I))         

         WRITE(LUPRI,'(/2X,A,D12.6)')
     &' DIFFERENCES A^2_ovlp-A^2_calc/N :',SUMDIFF

         IF (SUMDIFF.GE.THRNULL) WRITE(LUPRI,'(2X,A)')
     &' DIFFERENCES A^2_ovlp-A^2_calc/N DIFFERS TOO MUCH FROM "ZERO"!'

         IFC   = IFC + NESH(I)
         IKUN1 = IKUN1 + NESH(I)*NESH(I)
         IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
         IBUF = IBUF + NFBAS(I,0)*NFBAS(I,0)

 82     CONTINUE

        CALL MEMREL('DKH',WORK,KWORK,KKUN,KFREE,LFREE)

*------------------------------------------------------------------------
*  We did U+ Td U = td; U contains |k> vectors
*
*  Now we have the U in: 
*           TBUF - transformation matrix from ON basis "theta"
*           TMAT - transformation matrix from SO basis "xhi"
*
*------------------------------------------------------------------------
*  
*  Next we want to do:  U+ H U = H1,
*                      
*  where H is the bare nucleus Hamiltonian in the SO basis "xhi"
*            H1(in k) = TMAT+ * H(in xhi) * TMAT
*
*------------------------------------------------------------------------

* First get the bare nucleus Dirac Hamiltonian (BNDH) in the SO-basis "xhi"
* in WORK(KBUF)

      CALL ONEFCK(WORK(KBUF),IPRHAM,WORK(KFREE),LFREE) 
      CALL MEMGET('REAL',KHAO4,N2BBASXQ,WORK,KFREE,LFREE) ! Permanent vector for HAO Ham.
      CALL DCOPY(N2BBASXQ,WORK(KBUF),1,WORK(KHAO4),1) ! Store here HAO bare nucleus
*
*    transform the BNDH from its SO-basis "xhi" to the ON basis "|k>"
*    (reuse TBUF)
*
*                            TBUF = TMAT+ * H  * TMAT; 
*                      H1(in |k>) =   U1+ * H0 * U1
*
*  if desired, eliminate spin-orbit components from H1 in C1 symmetry
*  - spin-free at the "beginning"
*

      IF (SF0) THEN
         MZ = 1
      ELSE
         MZ = NZ
      END IF
      IOFT = 1
      DO 85 I = 1,NFSYM
      IF(NORB(I).EQ.0) GOTO 85
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &       NORB(I),NORB(I),
     &       WORK(KBUF+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &       TBUF(IOFT),NORB(I),NORB(I),MZ,IPQTOQ(1,0),
     &       TMAT(1+I2TMT(I)),NFBAS(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &       TMAT(1+I2TMT(I)),NFBAS(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)

        IF (SF0) THEN
         WRITE(LUPRI,'(/2X,A/)')
     & 'MAKE_DKH: Elim. spin-orbit comp. "BEFORE" from full H1 in TBUF'
         CALL SPF_HK( TBUF(IOFT),IBEIG(IORB(I)+1),IBEIG(IORB(I)+1),
     &                 NORB(I),NORB(I),NORB(I),NORB(I),NZ )
        ENDIF

        IF(IPRHAM.GE.4) THEN
          CALL HEADER(
     &     'MAKE_DKH: Bare nucleus matrix in 4c "|k>" basis',-1)
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(TBUF(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IOFT = IOFT + NORB(I)*NORB(I)*NZ

   85 CONTINUE
C
C  Available at this point:
C     TBUF: H1 matrix in 4c |k> basis
C     TMAT: 4c |k> basis in 4c |xhi> SO basis.
C     WORK(KBUF): free for new use
C  Allocate now:
C     H1AO( dimension N2BLQ ):
C        available as scratch now, will be used for final H1AO matrix later.
C     WORK(KCKL) : use for 1c(large) |k_eff> basis in 1c |xhi_L> SO basis.
*
Chj   N2BL  = NTBAS(1)*NTBAS(1)
      N2BL  = NTBAS(0)*NTBAS(0)
      N2BLQ = N2BL*NZ
C
      NCKL  = NFBAS(1,0)*NESH(1) + NFBAS(2,0)*NESH(2)
      NCKLQ = NCKL*NZ
      CALL MEMGET('REAL',KCKL ,NCKLQ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KCKL),NCKLQ)

C     allocate for H1 eigenvectors
      CALL MEMGET('REAL',KH1VEC,N2BBASXQ,WORK,KFREE,LFREE)

C     allocate space for final tranform. matrix  4c "xhi" -> 2c ON |k>
C        
      CALL MEMGET('REAL',KFTM,N2BBASXQ,WORK,KFREE,LFREE) ! 4c "xhi" -> 2c ON |k>
      CALL MEMGET('REAL',KFTMT,N2BBASXQ,WORK,KFREE,LFREE) ! 4c "xhi" -> 2c 'xhi'
      CALL DZERO(WORK(KFTMT),N2BBASXQ)

      IE1SH = 0
      JOFT1 = 1
      JCKL  = KCKL
      JBUF  = KBUF
      IFTM = KFTM
      IMAT = 1
      N2ORBT = 0
      DO 95 I = 1,NFSYM  ! Run big cycle over the fermion symmetries
      IF(NORB(I).EQ.0) GOTO 95

        N2T   = NORB(I)*NORB(I)
C       ... matrix dim of H1(k basis) in TBUF

        N2ESHI = NESH(I)*NESH(I)

C*********************************************************************
C       Diagonalize the 4c H1 Hamiltonian in the |k> basis...
C*********************************************************************
        CALL MEMGET('REAL',KDEIGVAL,
     &                N2T*NZ,WORK,KFREE,LFREE)
        CALL DZERO(WORK(KDEIGVAL),NORB(I)*NZ)

        CALL MEMGET('REAL',KH1D,N2T*NZ,WORK,KFREE,LFREE)
        CALL DCOPY(N2T*NZ,TBUF(JOFT1),1,WORK(KH1D),1)
        CALL MEMGET('REAL',KEV,
     &              N2T*NZ,WORK,KFREE,LFREE)
        CALL DZERO(WORK(KEV),N2T*NZ)

        CALL QDIAG(NZ,NORB(I),WORK(KH1D),NORB(I),NORB(I),   ! Diagonalize
     &             WORK(KDEIGVAL),0,          ! Eigenvalues
     &             WORK(KEV),NORB(I),       ! Eigenvectors not necessary
     &             NORB(I),WORK(KFREE),LFREE,IERR)
        IF (IERR.NE.0) CALL QUIT('MAKE_DKH 3: FAILED QDIAG...')

        CALL MEMREL('DKH 3',WORK,KWORK,KH1D,KFREE,LFREE)

C  ... the electronic eigenvalues are used for comparison 
C      with those from transformed two-component Hamiltonian !

C*********************************************************************
C====================================================================
C
C
C      Barysz-Sadlej-Snijders (BSS) infinite order of the relativity 
C     in two-component Hamiltonian matrix elements
C
C
C====================================================================
C*********************************************************************
        IF (IPZ .EQ. 9) THEN  ! 

C
C       ... infinite order Douglas-Kroll Hamiltonian
C
C          (1) diagonalize H1 in |k> basis
C              1.a) save H1 in WORK(KBUF)
C
           JH1VEC = KH1VEC + I2ORBT(I)
           IF (I.EQ.1) THEN

C             ... DFDIAG diagonalizes both symmetris if NFSYM.eq.2,
C                 thus we already have the solution 
C                 in EIG and WORK(JH1VEC) when I .eq. 2
C             ... If SPINFR true, then DFDIAG only diagonalizes
C                 the real part of TBUF (equivalent to the 4c
C                 SPINFR option)
C
C             ... KIBEIG is changed in DFDIAG !
C
              CALL MEMGET('INTE',KIBEIG,NORBT,WORK,KFREE,LFREE)
              CALL DCOPY(N2BBASXQ,TBUF,1,WORK(KBUF),1)
              SPINFR_SAVE = SPINFR
              SPINFR = .FALSE. ! When SF "BEFORE", H1 is pure spin-free
              CALL DFDIAG (WORK(KBUF),EIG,WORK(KIBEIG),WORK(JH1VEC),
     &                     .FALSE.,WORK,KFREE,LFREE)
              SPINFR      = SPINFR_SAVE
              CALL MEMREL('DKH 4',WORK,KWORK,KIBEIG,KFREE,LFREE)
           END IF
C
           IF(IPRHAM.GE.3) THEN
C          ... Write out the eigenvalues in the |k> basis
             CALL HEADER(
     &       'MAKE_DKH: Bare nucleus eigenvalues:',-1)
             WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
             WRITE(LUPRI,'(I5,F25.10,I5)')
     &       (J,EIG(IORB(I)+J),IBEIG(IORB(I)+J), J=1,NORB(I))
           ENDIF

           IF(IPRHAM.GE.7) THEN
C          ... Write out the eigenvectors in the |k> basis
             CALL HEADER(
     &'MAKE_DKH: Bare nucleus eigenvectors |v> in the |k> basis',-1)

            CALL PRQMAT(WORK(JH1VEC),NORB(I),NORB(I),
     &                  NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

C     test that matrix is a unit matrix
D            CALL HEADER(
D    &'MAKE_DKH: unit <v|v> overlap matrix the ON basis |k> ',-1)
D            CALL QGEMM(NORB(I),NORB(I),NORB(I),D1,
D    &          'H','N',IPQTOQ(1,0),TBUF(IMAT),NORB(I),NORB(I),NZ,
D    &          'N','N',IPQTOQ(1,0),TBUF(IMAT),NORB(I),NORB(I),NZ,
D    &          D0,IPQTOQ(1,0),WORK(KBUF), NORB(I),NORB(I),NZ)
D            CALL PRQMAT(WORK(KBUF),NORB(I),NORB(I),
D    &                NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
           ENDIF

           JVECP = JH1VEC + NPSH(I)*NORB(I)
           JVECE = JVECP + NPSH(I)
C          ... pointer to electron solutions:
C                         (Yp) -- positron |k_p> basis first
C           WORK(JVECP) = (  )
C                         (Ye) -- electron |k_e> basis last
C
C            need X+ = ( 1 R+ ) where
C
C            R = (1 - Yp*(Yp)+)-1 Yp(Ye)+  =  A-1 b
C
C            Solve this linear system by Cholesky decomposition
C            with full pivoting
C             --> notice that the A matrix is destroyed !
C
            JOB = 1
            NREF = NESH(I)
            NREF2 = 2*NREF  
            N2REF = NREF*NREF
            N2REFQ = N2REF*NZ
            CALL MEMGET('REAL',KAMAT,2*N2REFQ,WORK,KFREE,LFREE)
C           ... factor 2 for Omega+
            CALL MEMGET('REAL',KXMAT,2*N2REFQ,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KIPVT,NREF,WORK,KFREE,LFREE)

C           -- make A matrix
            CALL DZERO(WORK(KAMAT),N2REFQ)
            CALL DUNIT(WORK(KAMAT),NREF)   ! A := 1
            CALL QGEMM(NREF,NREF,NREF,DM1, ! A := A - Yp*(Yp)+
     &          'N','N',IPQTOQ(1,0),WORK(JVECP),NORB(I),NORB(I),NZ,
     &          'H','N',IPQTOQ(1,0),WORK(JVECP),NORB(I),NORB(I),NZ,
     &          D1,IPQTOQ(1,0),WORK(KAMAT),NREF,NREF,NZ)

            IF (IPRHAM.GE.7) THEN
              CALL HEADER( 'MAKE_DKH: A matrix for getting R',-1)
              WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
              CALL PRQMAT(WORK(KAMAT),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
            ENDIF

            CALL MEMGET('REAL',KA,N2REFQ,WORK,KFREE,LFREE) ! For A mtx
            CALL DCOPY(N2REFQ,WORK(KAMAT),1,WORK(KA),1)

C           -- make b matrix (store in X); b := Yp * (Ye)+
            CALL QGEMM(NREF,NREF,NREF,D1,
     &          'N','N',IPQTOQ(1,0),WORK(JVECP),NORB(I),NORB(I),NZ,
     &          'H','N',IPQTOQ(1,0),WORK(JVECE),NORB(I),NORB(I),NZ,
     &          D0,IPQTOQ(1,0),WORK(KXMAT),NREF2,NREF,NZ)

            IF (IPRHAM.GE.7) THEN
              CALL HEADER( 'MAKE_DKH: b matrix for getting R',-1)
              WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
              CALL PRQMAT(WORK(KXMAT),NREF,NREF,NREF2,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
            ENDIF

            CALL MEMGET('REAL',KB,N2REFQ,WORK,KFREE,LFREE) ! For B mtx
            DO IZ = 1, NZ
            DO IE = 1, NESH(I)
             CALL DCOPY(NREF,
     &                  WORK(KXMAT+(IZ-1)*NREF2*NREF+(IE-1)*NREF2),1,
     &                  WORK(KB+(IZ-1)*N2REF+(IE-1)*NREF),1)
            ENDDO
            ENDDO

            JDVEC = JH1VEC

C           ... reuse memory with positron solutions for
C               D(*) vector in Choleskyl decomp.
            CALL QCHOLD(WORK(KAMAT),NREF,NZ,NREF,NREF,WORK(JDVEC),
     &                  THRNULL,NEFF,JOB,WORK(KIPVT))
            CALL QCHOLS(WORK(KAMAT),NREF,NEFF,NREF,NZ,NREF,NREF,
     &           WORK(JDVEC),WORK(KXMAT),NREF2,NREF,
     &                       WORK(KXMAT),NREF2,NREF,
     &           JOB,WORK(KIPVT),WORK(JDVEC+NREF))
            IF(NEFF.LT.NREF) THEN
              WRITE(LUPRI,'(/2A)') '* MAKE_DKH WARNING: ',
     &         ' linear dependencies detected in Cholesky decomposition'
              WRITE(LUPRI,'(11X,A,I5,A,I5)')
     &         ' Reduced dimensionality: ',NREF,' --> ',NEFF
              CALL QUIT('MAKE_DKH: Cholesky problems...')
            ENDIF

C      that was R, check the R-condition
           CALL R_CHECK(I,TBUF(JOFT1),WORK(KXMAT),
     &            WORK(KA),WORK(KB),
     &          JCONSIST,WORK(KFREE),LFREE)

C       now insert unit matrix in bottom half of X:
            JOFF = KXMAT - 1 + NREF
            DO IZ = 1,NZ
              DO J = 1,NREF
                DO K = 1, NREF
                  WORK(JOFF+K) = D0
                END DO
                IF (IZ.EQ.1) WORK(JOFF+J) = D1
                JOFF = JOFF + NREF2
              END DO
            END DO
*
*         find Omega+ matrix and multiply it into WORK(KXMAT)
*
            CALL QGEMM(NREF,NREF,NREF2,D1,  ! R+.R = H1AO
     &          'H','N',IPQTOQ(1,0),WORK(KXMAT),NREF2,NREF,NZ,
     &          'N','N',IPQTOQ(1,0),WORK(KXMAT),NREF2,NREF,NZ,
     &          D0,IPQTOQ(1,0),H1AO,NREF,NREF,NZ)

            IF (IPRHAM.GE.7) THEN
              CALL HEADER(
     &  '  MAKE_DKH: X+ = ( 1 R+)',-1)
              WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
              CALL PRQMAT(WORK(KXMAT),NREF2,NREF,NREF2,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
              CALL HEADER(
     &  '  MAKE_DKH: X+ * X = 1 + R+ * R',-1)
              WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
              CALL PRQMAT(H1AO,NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
            ENDIF

            CALL QDIAG(NZ,NREF,H1AO,NREF,NREF,
     &                 WORK(JH1VEC),1,WORK(KAMAT),NREF,
     &                 NREF,WORK(KFREE),LFREE,IERR)

            CALL MEMGET('REAL',KAMAT1,2*N2REFQ,WORK,KFREE,LFREE)
            CALL DCOPY(N2REFQ,WORK(KAMAT),1,WORK(KAMAT1),1)

            IF (IPRHAM.GE.5) THEN
             CALL HEADER(
     &    'MAKE_DKH: Eigenvalues of (1 + R+ * R)',-1)
             DO J=1,NREF
              WRITE(LUPRI,'(I5,F25.10)') J,WORK(JH1VEC-1+J)
             ENDDO
            ENDIF

            IF (IPRHAM.GE.7) THEN
              CALL HEADER(
     &'MAKE_DKH: Eigenvectors of (1 + R+ * R)',-1)
              WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
              CALL PRQMAT(WORK(KAMAT),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
            ENDIF

            JEVEC = KAMAT
            DO IZ = 1, NZ
            DO J = 1,NREF
              FAC = D1 / SQRT(WORK(JH1VEC-1+J))
              CALL DSCAL(NREF,FAC,WORK(JEVEC),1) ! d^-1/2 * U(IZ)
              JEVEC = JEVEC + NREF 
            ENDDO
            ENDDO

C           now generate Omega+ = (X+ * X)^(-1/2) in H1AO

            CALL QGEMM(NREF,NREF,NREF,D1,  ! U * d^-1/2 * U+
     &          'N','N',IPQTOQ(1,0),WORK(KAMAT1),NREF,NREF,NZ,
     &          'H','N',IPQTOQ(1,0),WORK(KAMAT),NREF,NREF,NZ,
     &          D0,IPQTOQ(1,0),H1AO,NREF,NREF,NZ)

            IF (IPRHAM.GE.6) THEN
              CALL HEADER(
     &'MAKE_DKH: Omega+ = (1+ R+ * R)^(-1/2) ',-1)
              WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
              CALL PRQMAT(H1AO,NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
            ENDIF


            IF (IPRHAM.GE.6) THEN
              CALL HEADER(
     &  '  MAKE_DKH: Once again ( 1 R)',-1)
              WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
              CALL PRQMAT(WORK(KXMAT),NREF2,NREF,NREF2,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
            ENDIF

C           and multiply it to form X Omega+ = (1  R)t * Omega+:
            CALL QGEMM(NREF2,NREF,NREF,D1,
     &          'N','N',IPQTOQ(1,0),WORK(KXMAT),NREF2,NREF,NZ,
     &          'N','N',IPQTOQ(1,0),H1AO       ,NREF ,NREF,NZ,
     &          D0,IPQTOQ(1,0),WORK(KAMAT),NREF2,NREF,NZ)

            IF (IPRHAM.GE.6) THEN
              CALL HEADER(
     &'MAKE_DKH: Final U1= X Omega+ transform. matrix in |k> basis',-1)
              WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
              CALL PRQMAT(WORK(KAMAT),NREF2,NREF,NREF2,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
            END IF

* Construct transformation matrix,
*  Check the final transf. matrix by block diag.of full H1

           CALL U1CHECK(I,IBEIG,SF0,TBUF(JOFT1),WORK(KHAO4),WORK(KAMAT),
     &          WORK(KXMAT),H1AO,
     &          TMAT(IMAT),WORK(IFTM),
     &          JCONSIST,WORK(KFREE),LFREE)

C=====================================================================
C  Finally  generate H(DK_inf) = H1_eff = U1+ H1 U1 into WORK(JBUF)
C        TBUF: H1 matrix in 4c |k> basis
C        AMAT: Final U1 trasformation matrix
C        (remember NREF = NESH(I) and NREF2 = NORB(I) = 2*NREF)
C=====================================================================
            CALL QTRANS('AOMO','S',D0,NREF2,NREF2,NREF,NREF,
     &           TBUF(JOFT1),NREF2,NREF2,NZ ,IPQTOQ(1,0),
     &           WORK(JBUF) ,NREF ,NREF ,NZ ,IPQTOQ(1,0),
     &           WORK(KAMAT),NREF2,NREF ,NZT,IPQTOQ(1,0),
     &           WORK(KAMAT),NREF2,NREF ,NZT,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRHAM)

C===================================================================
C*******************************************************************
C
C
C    Now it's the second order Douglas-Kroll-Hess Hamiltonian
C  as pointed by Barysz&Sadlej, eq. 46.
C
C
C*******************************************************************
C===================================================================
      ELSE IF (IPZ .EQ. 2) THEN ! 
C       ... second order Douglas-Kroll Hamiltonian, or DK2
C
C        TBUF contains the H_1 matrix in orthonormal 4-comp. |k> basis
C             (= free particle basis in the specifed basis set);
C             first NPSH positron solutions, then NESH electron solutions
C
C        H_DKH(2,ee) = H_1(ee)
C                    + 0.5 * ( R_DK2+(ep) H_1(pe) + H_1(ep) R_DK2(pe) )
C   where:
C
C        R_DK2(pe)(i,j) = H_1(pe)(i,j)/(eig_e(i) - eig_p(j))
C     
C        (taken from Barysz&Sadlej, JMS(Theochem) 573(2001) 181, eq. 46 )    
C

          CALL HEADER(
     &'MAKE_DKH: DK2 Hamiltonian to be generated',-1)

C
C         1) generate R_DK2(pe) in TBUF(pe) block:
C
          IOFT  = JOFT1 + NORB(I)*NPSH(I) - 1
          JR2PE = IOFT  + 1
          DO IE = NPSH(I)+1,NORB(I)
            DO IP = 1, NPSH(I)

C ... eig_p and eig_e are separated by 2c^2; use rather average
C ... eig(p) = -eig(e) - 2c^2
C             FAC = DM1 / (EIG(IP+IORB(I)) - EIG(IE+IORB(I))) ! Original expression

C              FAC = DM1 / ! Use electronic fp-eigenvalues only 
C     &   (-EIG(NORB(I)-IP+1+IORB(I))-(D2*CVAL2)-EIG(IE+IORB(I)))

                EIG_P = ! Calculate the average value
     & ((-EIG(NORB(I)-IP+1+IORB(I))-(D2*CVAL2)) + 
     &   EIG(IP+IORB(I)))/D2

                EIG_E = ! Calculate the average value
     & ( EIG(IE+IORB(I)) + 
     &   (-EIG(NORB(I)-IE+1+IORB(I))-(D2*CVAL2)) )/D2

              FAC = DM1/ (EIG_P - EIG_E)

              DO IZ = 0,NZ-1
                TBUF(IOFT+IP+IZ*N2T) =
     &          TBUF(IOFT+IP+IZ*N2T) * FAC
              END DO
            END DO
            IOFT = IOFT + NORB(I)
          END DO

          JH1EP = JOFT1 + NPSH(I)
 
C -----------------------------------------------------------------------------yy
C  If desired, eliminate SO terms from both H1_ep and R_pe - DK2SO1 approach
C -----------------------------------------------------------------------------yy
          IF (DK2SO1) THEN
           WRITE(LUPRI,'(3X,A)')
     & 'DK2SO1 - elim. spin-orbit effects from H1(ep) and R(pe)'
           CALL SPF_HK( TBUF(JH1EP),IBEIG(IORB(I)+NPSH(I)+1),
     &                 IBEIG(IORB(I)+1),
     &                 NESH(I),NPSH(I),NORB(I),NORB(I),NZ )
           CALL SPF_HK( TBUF(JR2PE),IBEIG(IORB(I)+1),
     &                 IBEIG(IORB(I)+1+NPSH(I)),
     &                 NPSH(I),NESH(I),NORB(I),NORB(I),NZ )
          ENDIF

          IF (IPRHAM.GE.6) THEN
           CALL HEADER(
     &'MAKE_DKH: 1. order wave function R_DK2(pe) in "|k>" basis',-1)
           WRITE(LUPRI,'(3X,A,I2/)') '*** Fermion corep ',I
           CALL PRQMAT(TBUF(JR2PE),NPSH(I),NESH(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
           CALL HEADER(
     &'MAKE_DKH: H(ep) part of H1 in "|k>" basis',-1)
           WRITE(LUPRI,'(3X,A,I2/)') '*** Fermion corep ',I
           CALL PRQMAT(TBUF(JH1EP),NESH(I),NPSH(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

C
C         2) Generate  "X" = H_1(ep) R_DK2(pe) in WORK(JBUF)
C

          CALL QGEMM(NESH(I),NESH(I),NPSH(I),D1,
     &          'N','N',IPQTOQ(1,0),TBUF(JH1EP),NORB(I),NORB(I),NZ,
     &          'N','N',IPQTOQ(1,0),TBUF(JR2PE),NORB(I),NORB(I),NZ,
     &          D0,IPQTOQ(1,0),WORK(JBUF), NESH(I),NESH(I),NZ)
C
C         3) Generate  X = 0.5 ( X + X(dagger)) 
C
          CALL QGETHE(WORK(JBUF),NESH(I),NZ,NESH(I),NESH(I))
C         CALL QGETHE(AMAT,NDIM,NZ,LRA,LCA)

          IF (IPRHAM.GE.7) THEN
            CALL HEADER(
     &'MAKE_DKH: 2. order correction to H1 in "|k_e>" basis',-1)
            WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
            CALL PRQMAT(WORK(JBUF),NESH(I),NESH(I),NESH(I),
     &                  NESH(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

*  Add to first order (ee-blocks!): H_1(ee) =  H_1(ee) + X

          DO IZ = 1,NZ
            IBUF  = JBUF + (IZ-1)*N2ESHI
            JH1EE = JOFT1 + NPSH(I)*NORB(I) + NPSH(I)
     &            + (IZ-1)*NORB(I)*NORB(I)
            DO J = 1,NESH(I)
              CALL DAXPY(NESH(I),D1,TBUF(JH1EE),1,WORK(IBUF),1)
              JH1EE = JH1EE + NORB(I)
              IBUF  = IBUF  + NESH(I)
            END DO
          END DO

C**********************************************************************
C
C           First order transformed Hamiltonian, DK1
C
C**********************************************************************
        ELSE IF (IPZ .EQ. 1) THEN  ! 
C       ... only first order Douglas-Kroll Hamiltonian
C

C       Store into WORK(IBUF) the ee-part of the full H1 Hamiltonian
          DO IZ = 1,NZ
            IBUF  = JBUF + (IZ-1)*N2ESHI
            JH1EE = JOFT1 + NPSH(I)*NORB(I) + NPSH(I)
     &            + (IZ-1)*NORB(I)*NORB(I)
            DO J = 1,NESH(I)
              CALL DCOPY(NESH(I),TBUF(JH1EE),1,WORK(IBUF),1)
              JH1EE = JH1EE + NORB(I)
              IBUF  = IBUF  + NESH(I)
            END DO
          END DO

        ELSE
C         not DKinf nor DK2
          CALL QUIT('only DKH inf & 1./2. order implemented yet')
        END IF

C**********************************************************************
C----------------------------------------------------------------------
C    Now we have the elements of the effective tranformed two-component
C    Hamiltonian H(DKn) in the WORK(JBUF) array in the ON "|k_e>" basis
C----------------------------------------------------------------------
C**********************************************************************

        IF (SPINFR.AND.(.NOT.SF0)) THEN 
         WRITE(LUPRI,'(/2X,A/)')
     &'MAKE_DKH: Eliminating spin-orbit terms from'//
     &' the H(2c) in |k> basis "AFTER"'
         CALL SPF_HK( WORK(JBUF),IBEIG(IORB(I)+NPSH(I)+1),
     &    IBEIG(IORB(I)+NPSH(I)+1), 
     &    NESH(I),NESH(I),NESH(I),NESH(I),NZ )
        ENDIF

        IF (IPRHAM.GE.7) THEN
          CALL HEADER(
     &'MAKE_DKH: *** The H(2c)_ee in "|k>" basis ***',-1)
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(WORK(JBUF),NESH(I),NESH(I),NESH(I),
     &                NESH(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C--------------------------------------------------------------------
C DIAGONALIZE H(DKH)_ee in "|k_e>" ! Compare its eigenvalues with those
C from the full H1 in ON 4c |k> basis !
C--------------------------------------------------------------------
        CALL MEMGET('REAL',KHK,N2ESHI*NZ, WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIGVECT,N2ESHI*NZ, WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIGVAL, NESH(I)*NZ,WORK,KFREE,LFREE)

        CALL DZERO(WORK(KEIGVECT),N2ESHI*NZ)
        CALL DZERO(WORK(KEIGVAL),NESH(I)*NZ)

        CALL DCOPY(N2ESHI*NZ,
     &          WORK(JBUF),1,WORK(KHK),1)! Copy H_ee into help array

        CALL QDIAG(NZ,NESH(I),WORK(KHK),NESH(I),NESH(I),   ! Diagonalize
     &             WORK(KEIGVAL),0,          ! Eigenvalues
     &             WORK(KEIGVECT),NESH(I),       ! Eigenvectors not necessary
     &             NESH(I),WORK(KFREE),LFREE,IERR)
        IF (IERR.NE.0) CALL QUIT("BAD DIAGONALIZATION!")

*
*       .... compare electronic eigenvalues 4c H1 / 2c H1 in |k>
*            for each fermion symmetry....
*
        CALL HEADER(
     &  'MAKE_DKH: *** 4c/2c electronic eigenvalues in |k> basis ***'
     &  ,-1)
        WRITE(LUPRI,'(3X,A,I2/)') '*** Fermion corep ',I
        WRITE(LUPRI,'(17X,A,18X,A,12X,A/)')
     &  '4c |k>','2c |k>','diff 4c-2c'
        JOFF = KEIGVAL - 1
        DO J = 1,NESH(I)
          D4V2C = WORK(KDEIGVAL+NPSH(I)+J-1)-WORK(JOFF+J)
          WRITE(LUPRI,'(I3,2F25.10,5X,D12.5)')
     &    J,
     &    WORK(KDEIGVAL+NPSH(I)+J-1),
     &    WORK(JOFF+J),D4V2C 
        END DO

        IE1SH = IE1SH + NESH(I)   
        JOFT1 = JOFT1 + N2T*NZ
        JCKL  = JCKL  + NESH(I)*NFBAS(I,0)*NZ
        JBUF  = JBUF  + NESH(I)*NESH(I)*NZ
        KHK   = KHK   + NESH(I)*NESH(I)*NZ
        IFTM  = IFTM  + NESH(I)*NFBAS(I,0)*NZ
        IMAT  = IMAT  + NFBAS(I,0)*NORB(I)*NZ
C

   95 CONTINUE
C     end do 95 I = 1,NFSYM

C*****************************************************************
C      Make contravariant transformation matrix
C          (C+el)-1 = Sao(LL) C(el) Ad^(-1)
C
C      and transform H1(DKn) matrix back to "xhi_L" basis in H1AO
C
C  TMAT - ON |k> vectors (e,p) in the SO "xhi" basis; NFORB(I,0)*NORB(I)*NZ
C
C*****************************************************************

      CALL MEMGET('REAL',KCKLX,NCKLQ,WORK,KFREE,LFREE)

      CALL MKSAOMO_Le(WORK(KCKLX),TMAT,FACKIN, ! For all ferm.symm.
     &                IPRHAM,WORK(KFREE),LFREE)

C
C    Generate DKn  n'th order  H1(xhi_LL) matrix in H1AO ....
C         H1(xhi_L) = (C+)-1      H_DKn(k) C-1
C                   = Sao C Ad^-1 H_DKn(k) Ad^-1 C+ Sao
C
C       H_DKn = WORK(JBUF) ... H_ee in |k_e> basis
C
      CALL DZERO(H1AO,N2BLQ)
      JCKLX = KCKLX
      JBUF  = KBUF
      IOFV  = 1
      IFTM = KFTM
      IFTMT = KFTMT
      IXT   = KXT
      DO 105 I = 1, NFSYM
      IF(NFBAS(I,0).EQ.0) GOTO 105

C       Row dimension of H1AO: only L [NFBAS(I,1)] or L+S [NFBAS(I,0)]?
C       currently we use L+S, to use only L we should delete small Xhi
C       after MAKE_DKH /5-Aug-2002 hjaaj
C       (PS. remember to change N2BL above if NDH1AO is changed)
C
C       NDH1AO = NFBAS(I,1)
        NDH1AO = NTBAS(0)
        JH1AO  = 1+I2BASX(I,I)

        CALL QTRANS('MOAO','S',D0,NFBAS(I,1),NFBAS(I,1),
     &               NESH(I),NESH(I),
     &            H1AO(JH1AO),NDH1AO,NDH1AO,NZ,IPQTOQ(1,0),
     &            WORK(JBUF),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &            WORK(JCKLX),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &            WORK(JCKLX),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRHAM)

        IF (IPRHAM.GE.3) THEN
          CALL HEADER(
     & 'MAKE_DKH: Final H(2c)(LLblock) in "xhi_L" SO basis',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     & '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(H1AO(JH1AO),NFBAS(I,0),NFBAS(I,0),
     &          NDH1AO,NDH1AO,NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C ..................................................................
C    ...construct final transf. matrix 4c "xhi" -> 2c "xhi"
C ..................................................................
        CALL LOWDHAO(IPZ,I,H1AO(JH1AO),WORK(KHAO4),WORK(IFTM),
     &    WORK(IFTMT),WORK(JCKLX),VMAT(IOFV),WORK(IXT),
     &     IBEIG,JCONSIST,WORK(KFREE),LFREE)

        IOFV = IOFV   + NFORB(I,0)*NFBAS(I,0)
        IXT  = IXT    + NFBAS(I,0)*NORB(I)*NZ

        JCKLX = JCKLX + NFBAS(I,0)*NESH(I)*NZ
        IFTM  = IFTM  + NFBAS(I,0)*NESH(I)*NZ
C ... offset for pict.change transf.matrix
        IFTMT = IFTMT + NFBAS(I,0)*NFBAS(I,1)*NZ
        JBUF  = JBUF  + NESH(I)*NESH(I)*NZ

  105 CONTINUE

C*****************************************************************************
C   Reorder the picture change transformation matrix
C     from the WORK(KFTMT) array into the PCTM array...
C*****************************************************************************
C      CALL DZERO(PCTM(1),N2BBASXQ)
C      IX = 0
C      DO I = 1, NFSYM
C       IPOS=(I-1)*( NTBAS(0)*NFBAS(1,0)+NFBAS(1,0))
C      DO J = 1, NFBAS(I,1)
C       JPOS = (J-1)*NTBAS(0)
C      DO K = 1, NFBAS(I,0)
C       PCTM(IPOS+JPOS+K) = WORK(KFTMT+IX)
C       IX = IX + 1
C      ENDDO
C      ENDDO
C      ENDDO

C ... control print out....
C      IF (IPRHAM.GE.6) THEN
C        CALL HEADER(
C    &'MAKE_DKH: Total picture change transformation matrix:',-1)
C        CALL PRQMAT(PCTM(1),NTBAS(0),NTBAS(1),NTBAS(0),
C    &                 NTBAS(1),NZ,IPQTOQ(1,0),LUPRI)
C       ENDIF
C       CALL QUIT('Tu koniec...')

C
C Copy transformation matrix FTM/FTMT into TMAT
C
*      CALL DCOPY(N2BBASXQ,WORK(KFTM),1,TMAT(1),1) ! 4c "xhi" -> 2c |k_e>
*      CALL DCOPY(N2BBASXQ,WORK(KFTMT),1,TMAT(1),1) ! 4c "xhi" -> 2c "xhi"
C*****************************************************************************
C ... no reordering...
      CALL DCOPY(N2BBASXQ,WORK(KFTMT),1,PCTM(1),1) ! 4c "xhi" -> 2c "xhi"
C*****************************************************************************
C .. delete empty places....
C     IDEL=N2BBASXQ-NTBAS(1)*NTBAS(0)*NZ
C     IPOSDEL=NTBAS(1)*NTBAS(0)*NZ+1
C     CALL DZERO(PCTM(IPOSDEL),IDEL)

C
C  Finally set up number of positronic shells equal to zero
C

      DO I=1,NFSYM
       NPSH(I) = 0
       NORB(I) = NESH(I) + NPSH(I)
C  ... new variable ....
       NTMO(I) = NORB(I)
       DO ISYM = 1, NSYM
         NBORB(ISYM,I,2) = 0
         NBORB(ISYM,I,0) = NBORB(ISYM,I,1) + NBORB(ISYM,I,2)
       ENDDO
      ENDDO

      WRITE(LUPRI,'(/2X,A/)') 
     &'MAKE_DKH:Two-component BSS mode - positronic shells deleted !'

      CALL MEMREL('MAKE_DKH final',WORK,KWORK,KWORK,KFREE,LFREE)
C     ... reset DKH to true
      BSS = .TRUE.
      CALL FLSHFO(LUPRI)
      CALL QEXIT('MAKE_DKH')
C      CALL QUIT('End of MAKE_DKH')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck GET_FP_LOWD */
       SUBROUTINE GET_FP_LOWD(TBUF,TMAT,VMAT,IPRHAM,WORK,LWORK)
C-------------------------------------------------------------------------
C    PURPOSE:
C  ===========
C     Generate free particle one-electron matrix in Lowdin-basis
C-------------------------------------------------------------------------
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0 )
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION TBUF(*),VMAT(*),TMAT(*),WORK(LWORK)

      CALL QENTER('GET_FP_LOWD')     
#include <memint.h>

C First construct quaternion one-electron free particle in SO-basis 'xhi'
      CALL FREEMT(TBUF,IPRHAM,WORK(KFREE),LFREE) 
C
C Then transform to Lowdin basis "Theta" (from TBUF to TMAT):
C
C                TMAT = VMATt *  TBUF  * VMAT
      IOFT = 1
      IOFV = 1
      DO I = 1,NFSYM
        IF(NFORB(I,0).GT.0) THEN

        IF (IPRHAM.GE.7) THEN
           CALL HEADER(
     &  'MAKE_DKH: Free particle matrix in SO basis "xhi"',-1)
           WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
           CALL PRQMAT(TBUF(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                 NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)

           CALL HEADER(
     &  'MAKE_DKH: Lowdin transf. matrix "xhi"->"Theta"',-1)
           WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
           CALL PRQMAT(VMAT(IOFV),NFBAS(I,0),NFORB(I,0),
     &                 NFBAS(I,0),NFORB(I,0),NZ,IPQTOQ(1,0),LUPRI)

        ENDIF
C
C                TMAT  =  VMATt *  TBUF  * VMAT
C
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0), 
     &                         NFORB(I,0),NFORB(I,0),
     &            TBUF(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            TMAT(IOFT),NFORB(I,0),NFORB(I,0),NZ,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRHAM)

        IF (IPRHAM.GE.7) THEN
           CALL HEADER(
     &   'MAKE_DKH: Free particle matrix in "Lowdin" basis "Theta"',-1)
           WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
           CALL PRQMAT(TMAT(IOFT),NFORB(I,0),NFORB(I,0),NFORB(I,0),
     &                 NFORB(I,0),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
        IOFT = IOFT + NFORB(I,0)*NFORB(I,0)*NZ
        IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)

       ENDIF
      ENDDO

      CALL QEXIT('GET_FP_LOWD')     
      RETURN
      END    

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck MKSAOMO_Le */
      SUBROUTINE MKSAOMO_Le(SAOMO,TMAT,FACKIN,IPRINT,WORK,LWORK)
C****************************************************************************
C     PURPOSE:
C
C        Construct half transformed overlap matrix, 
C
C         Saomo = Sao(LL) * C(L,el) * Ad(ee)^(-1)
C
C      Written:  M.Ilias  (based on MKSAOMO by L. Visscher)
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00 )
#include <dcbgen.h>
#include <dcbdhf.h>
#include <dcbbas.h>
#include <dcbham.h>
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION SAOMO(*),TMAT(*),WORK(*),FACKIN(2,*)
      CALL QENTER('MKSAOMO_Le')
#include <memint.h>
C
C
      ISIZE = 0 ! Calculate total size of tranf. matrix
      DO I = 1,NFSYM
       ISIZE = ISIZE + NZ*NFBAS(I,0)*NESH(I)
      END DO

      CALL DZERO(SAOMO,ISIZE)

C
C     Get AO overlap matrix , S(LL), WITH SS part ZEROED !!
C
      CALL MEMGET('REAL',KSAO,N2BBASX,WORK,KFREE,LFREE)
      CALL GTOVLT(WORK(KSAO),D0,WORK(KFREE),LFREE)

C
C   S_LL(NFBAS(I,0),NFBAS(I,0)) * TMAT_e(NFBAS(I,0),NESH(I))
C
      IOFF = 1
      JOFFESH = 0
      ITMT = 0
      ISAO = KSAO
      DO I = 1,NFSYM
       IF ( NESH(I) .GT. 0 ) THEN

*              SAOMO_Le = S_LL * TMAT

         IF ( IPRINT .GE. 5 ) THEN
            CALL HEADER(
     &   'MKSAOMO_Le: SOVLP matrix(LL) in SO basis "xhi" ',-1)
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',I
          CALL PRQMAT(WORK(ISAO),NFBAS(I,0),NFBAS(I,0),
     &                NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),LUPRI)
         ENDIF

         CALL QGEMM(NFBAS(I,0),NESH(I),NFBAS(I,0),D1,
     &         'N','N',IPQTOQ(1,0),WORK(ISAO),
     &         NFBAS(I,0),NFBAS(I,0),1,
     &         'N','N',IPQTOQ(1,0),
     &         TMAT(ITMT+1+NFBAS(I,0)*NESH(I)),
     &         NFBAS(I,0),NESH(I),NZT,
     &         D0,IPQTOQ(1,0),SAOMO(IOFF),
     &         NFBAS(I,0),NESH(I),NZT)
C

            IF ( IPRINT .GE. 5 ) THEN
             CALL HEADER('MKSAOMO_Le: Right index transf. overlap'//
     &         ' matrix Saomo_Le',-1)
             WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
             CALL PRQMAT(SAOMO(IOFF),NFBAS(I,0),NESH(I),
     &            NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
              CALL FLSHFO(LUPRI)
            END IF

*            Now multiply by diagonal Ad^-1, FACKIN(1,j)  = Ad^2(j,j)
*               SAOMO_Le = SAOMO_Le * Ad^-1

             ICX = 0
             DO IZ = 1,NZ
             DO IESH = 1, NESH(I)
             DO IFBAS = 1, NFBAS(I,0) 
              SAOMO(IOFF+ICX) = 
     &        SAOMO(IOFF+ICX)/DSQRT(FACKIN(1,IESH+JOFFESH))
              ICX = ICX + 1
             ENDDO
             ENDDO 
             ENDDO 

            IF ( IPRINT .GE. 5 ) THEN
             CALL HEADER('MKSAOMO_Le: Right index transf. overlap'//
     &         ' matrix Saomo_Le * Ad(e)^-1',-1)
             WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
             CALL PRQMAT(SAOMO(IOFF),NFBAS(I,0),NESH(I),
     &            NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
              CALL FLSHFO(LUPRI)
            END IF

         END IF

      JOFFESH = JOFFESH + NESH(I)
      IOFF = IOFF + NZ*NFBAS(I,0)*NESH(I)
      ITMT = ITMT + NFBAS(I,0)*NORB(I)
      ISAO = ISAO + NFBAS(I,0)*NFBAS(I,0)

      END DO

      CALL MEMREL('SAOMOLe',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('MKSAOMO_Le')
      RETURN
      END


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck R_CHECK */
      SUBROUTINE R_CHECK(ISYMM,HK,R,A,B,JCONSIST,WORK,LWORK)
C****************************************************************************
C     PURPOSE:  
C
C     Check the R-condition for BSS infinite order:
C
C                H(pp)*R = -H(pe) + R*H(ee) + R*H(pp)*R
C
C    On input:  HK  - bare nucleus Hamiltonian matrix in |k> basis
C                R  - the R matrix in |k> basis
C                A, B -   A.R = B
C                ISYMM - the fermion symmetry (1 or 2) 
C                JCONSIST - condition for consistency checks
C
C     Output:  just check and print the result if JCONSIST.ge.IPRHAM
C
C      Written:  M.Ilias, Odense, 2002  
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00, DM1 = -1.00D00 )
#include <dcbgen.h>
#include <dcbdhf.h>
#include <dcbbas.h>
#include <dcbham.h>
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION HK(*),R(*),A(*),B(*),WORK(*)
      CALL QENTER('R_CHECK')
#include <memint.h>

      IF (IPRHAM.LT.JCONSIST) GOTO 123

      I = ISYMM

      IF (NPSH(I).NE.NESH(I)) CALL QUIT("CHECK_R:Ne<>Np!")
      IF (NORB(I).NE.(NESH(I)+NPSH(I))) CALL QUIT("CHECK_R:Ne<>Np!")


      IF ( IPRHAM .GE. 5 ) THEN
        CALL HEADER(
     & 'CHECK_R: Bare nucleus matrix in |k> basis',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(HK(1),NORB(I),NORB(I),NORB(I),NORB(I),
     &            NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     & 'CHECK_R: R matrix in |k> basis',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(R(1),NESH(I),NESH(I),2*NESH(I),NESH(I),
     &            NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     & 'CHECK_R: A matrix in |k> basis',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(A(1),NESH(I),NESH(I),NESH(I),NESH(I),
     &            NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     & 'CHECK_R: B matrix in |k> basis',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(B(1),NESH(I),NESH(I),NESH(I),NESH(I),
     &            NZ,IPQTOQ(1,0),LUPRI)

      END IF

      N2EQ = NESH(I)*NESH(I)*NZ

      CALL MEMGET('REAL',KHPPR,N2EQ,WORK,KFREE,LFREE) ! For H(pp)*R
      CALL MEMGET('REAL',KRHEE,N2EQ,WORK,KFREE,LFREE) ! For R*H(ee)
      CALL MEMGET('REAL',KHEP,N2EQ,WORK,KFREE,LFREE)  ! For H(ep)
      CALL MEMGET('REAL',KRHPER,N2EQ,WORK,KFREE,LFREE)! For R*H(pe)*R 
      CALL MEMGET('REAL',K1,N2EQ,WORK,KFREE,LFREE)!  Help array 1
      CALL MEMGET('REAL',K2,N2EQ,WORK,KFREE,LFREE)!  Help array 2

C     Do   WORK(K2) = A * R

      CALL QGEMM(NESH(I),NESH(I),NESH(I),D1,
     &         'N','N',IPQTOQ(1,0),A(1),
     &         NESH(I),NESH(I),NZ,
     &         'N','N',IPQTOQ(1,0),
     &         R(1),2*NESH(I),NESH(I),NZ,          
     &         D0,IPQTOQ(1,0),WORK(K2),
     &         NESH(I),NESH(I),NZ)

      CALL DAXPY(N2EQ,DM1 , B(1),  1,WORK(K2),1) !  A * R - B = 0
      ARMBNULLNORM =  DNRM2(N2EQ,WORK(K2),1)
      CALL HEADER(
     & 'CHECK_R: A * R - B zero matrix',-1)
      WRITE(LUPRI,'(3X,A,I1,A,I1)')'* Fermion irp ',I,'/',NFSYM
      WRITE(LUPRI,'(3X,A,D17.12)')
     &'CHECK_R: Norm of A.R - B:',ARMBNULLNORM
      IF ( IPRHAM .GE. 5 ) THEN
        CALL PRQMAT(WORK(K2),NESH(I),NESH(I),NESH(I),NESH(I),
     &               NZ,IPQTOQ(1,0),LUPRI)
      ENDIF

      DO IZ=1,NZ
      DO IE=1,NESH(I)
        CALL DCOPY(NESH(I),  ! Store H(ep)
     &    HK(1+ NORB(I)*NESH(I) +
     &        (IZ-1)*NORB(I)*NORB(I)+
     &       (IE-1)*NORB(I) ),  1,
     &   WORK(KHEP+(IE-1)*NESH(I)+(IZ-1)*NESH(I)*NESH(I)),1)
      ENDDO
      ENDDO

      IF ( IPRHAM .GE. 5 ) THEN
        CALL HEADER(
     & 'CHECK_R: The entering R matrix in |k> basis',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(R(1),NESH(I),NESH(I),2*NESH(I),NESH(I),
     &              NZ,IPQTOQ(1,0),LUPRI)
        CALL HEADER(
     & 'CHECK_R: H(pp) part in |k> basis',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(HK(1),NESH(I),NESH(I),NORB(I),NORB(I),
     &         NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     & 'CHECK_R: Bare nucleus matrix - ep part in |k> basis',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(WORK(KHEP),NESH(I),NESH(I),NESH(I),NESH(I),
     &              NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C   H(pp) * R = WORK(KHPPR)

      CALL QGEMM(NESH(I),NESH(I),NESH(I),D1,
     &         'N','N',IPQTOQ(1,0),HK(1),
     &         NORB(I),NORB(I),NZ,
     &         'N','N',IPQTOQ(1,0),
     &         R(1),2*NESH(I),NESH(I),NZ,          
     &         D0,IPQTOQ(1,0),WORK(KHPPR),
     &         NESH(I),NESH(I),NZ)

      IF ( IPRHAM .GE. 5 ) THEN
        CALL HEADER(
     & 'CHECK_R: H(pp)*R matrix',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(WORK(KHPPR),NESH(I),NESH(I),NESH(I),NESH(I),
     &         NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C    R * H(ee) = WORK(KRHEE)

      CALL QGEMM(NESH(I),NESH(I),NESH(I),D1,
     &         'N','N',IPQTOQ(1,0),R(1),
     &         2*NESH(I),NESH(I),NZ,
     &         'N','N',IPQTOQ(1,0),
     &  HK(1+NORB(I)*NPSH(I)+NPSH(I)),NORB(I),NORB(I),NZ,          
     &         D0,IPQTOQ(1,0),
     &         WORK(KRHEE),
     &         NESH(I),NESH(I),NZ)

      IF ( IPRHAM .GE. 5 ) THEN
        CALL HEADER(
     & 'CHECK_R: R*H(ee) matrix',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(WORK(KRHEE),NESH(I),NESH(I),NESH(I),NESH(I),
     &         NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C    H(pe) * R = WORK(K1)

      CALL QGEMM(NESH(I),NESH(I),NESH(I),D1,
     &         'N','N',IPQTOQ(1,0),
     &  HK(1+NPSH(I)),NORB(I),NORB(I),NZ,
     &         'N','N',IPQTOQ(1,0),
     &         R(1),2*NESH(I),NESH(I),NZ,          
     &         D0,IPQTOQ(1,0),WORK(K1),
     &         NESH(I),NESH(I),NZ)

      IF ( IPRHAM .GE. 5 ) THEN
        CALL HEADER(
     & 'CHECK_R: H(pe)*R matrix',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(WORK(K1),NESH(I),NESH(I),NESH(I),NESH(I),
     &         NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C    R* WORK(K1) = WORK(KRHPER)  ... R*H(pe)*R

      CALL QGEMM(NESH(I),NESH(I),NESH(I),D1,
     &         'N','N',IPQTOQ(1,0),R(1),
     &         2*NESH(I),NESH(I),NZ,
     &         'N','N',IPQTOQ(1,0),
     &         WORK(K1),NESH(I),NESH(I),NZ,          
     &         D0,IPQTOQ(1,0),WORK(KRHPER),
     &         NESH(I),NESH(I),NZ)

      IF ( IPRHAM .GE. 5 ) THEN
        CALL HEADER(
     & 'CHECK_R: R*H(pe)*R matrix',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(WORK(KRHPER),NESH(I),NESH(I),NESH(I),NESH(I),
     &         NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C
C DO:  WORK(K1)= H(pp)*R + H(ep) - R*H(ee) - R*H(pe)*R
C      WORK(K1)= WORK(KHPPR) + WORK(KHEP) - WORK(KRHEE) - WORK(KRHPER)
C

      CALL DZERO(WORK(K1),N2EQ)
      CALL DAXPY(N2EQ,D1 ,WORK(KHPPR),  1,WORK(K1),1) ! + WORK(KHPPR)
      CALL DAXPY(N2EQ,D1 ,WORK(KHEP),   1,WORK(K1),1) ! + WORK(KHEP)
      CALL DAXPY(N2EQ,DM1,WORK(KRHEE),  1,WORK(K1),1) ! - WORK(KRHEE)
      CALL DAXPY(N2EQ,DM1,WORK(KRHPER), 1,WORK(K1),1) ! - WORK(KRHPER)

      RNULLNORM =  DNRM2(N2EQ,WORK(K1),1)
      CALL HEADER(
     &'CHECK_R: Final H(pp)*R+H(ep)-R*H(ee)-R*H(pe)*R',-1)
      WRITE(LUPRI,'(7X,A,I1,A,I1,3X,A,D17.12)')
     &' Fermion corep ',I,'/',NFSYM,
     & '* norm: ',RNULLNORM
      IF ( IPRHAM .GE. 5 ) THEN
        CALL PRQMAT(WORK(K1),NESH(I),NESH(I),NESH(I),NESH(I),
     &               NZ,IPQTOQ(1,0),LUPRI)
      ENDIF

      CALL MEMREL('R_CHECK',WORK,KWORK,KWORK,KFREE,LFREE)
C
 123  CALL QEXIT('R_CHECK')
C      CALL QUIT( 'R_CHECK')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck U1CHECK */
      SUBROUTINE U1CHECK(ISYMM,IBEIG,SFRB,H1,HAO,U1,R1,DOM,
     &  TMAT,FTMAT,JCONSIST,WORK,LWORK)
C****************************************************************************
C
C     PURPOSE:  Check the U1+ * U1 = 1, FTMAT+ * H4cAO * FTMAT = H2c_k
C               Contruct final transf. matrix FTMAT = TMAT * U1
C
C    THIS ROUTINE IS CALLED ONLY FOR THE INFINITE ORDER OF RELATIVITY !!!
C
C    On input:  ISYMM - the fermion symmetry (1 or 2)
C               SFRB - logical keyword if spin-free "BEFORE"
C               IBEIG - boson irreps of MO's in spin-free mode "BEFORE"
C                     (it must have dcborb.h var)
C               H1  - bare nucleus 4c Hamiltonian matrix in |k> basis
C               HAO -  bare nucleus 4c Hamiltonian full matrix in SO basis "xhi"
C               U1  - final transf. matrix..  U1+ H1 U1 = h_ee in |k>
C               R1 - (R 1) matrix
C               DOM - Omega+ matrix
C               JCONSIST - consistency check level
C
C     Output:   final transf. matrix FTMAT  4c "xhi" -> 2c ON |k_e>
C
C      Written:  M.Ilias, Odense, 2002  
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00, DM1 = -1.00D00 )
#include <dcbgen.h>
#include <dcbdhf.h>
#include <dcbbas.h>
#include <dcbham.h>
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION H1(*),U1(*),HAO(*),R1(*),DOM(*),TMAT(*),FTMAT(*),WORK(*)
      DIMENSION IBEIG(*)
      LOGICAL SFRB
      CALL QENTER('U1CHECK')
#include <memint.h>
C

C ... assign the actual fermion symmetry !!!
      I = ISYMM

      N2OQ = NORB(I)*NORB(I)*NZ
      N2UQ = 2*NESH(I)*NESH(I)*NZ
      N2EQ = NESH(I)*NESH(I)*NZ

      NREF = NESH(I)
      NREF2 = 2*NESH(I)

C     WRITE(*,*) "I,I2BASX(I,I),NFBAS(I,0)*NFBAS(I,0)*NZ",
C    &    I,I2BASX(I,I),NFBAS(I,0)*NFBAS(I,0)*NZ

      CALL MEMGET('REAL',KH1,N2OQ,WORK,KFREE,LFREE) ! 
      CALL MEMGET('REAL',KU1,N2UQ,WORK,KFREE,LFREE) ! 
      CALL MEMGET('REAL',KUU,N2EQ,WORK,KFREE,LFREE) ! 
      CALL MEMGET('REAL',KHEE,N2EQ,WORK,KFREE,LFREE) ! 
      CALL MEMGET('REAL',KHEE1,N2EQ,WORK,KFREE,LFREE) ! 
      CALL MEMGET('REAL',KRR,N2EQ,WORK,KFREE,LFREE) ! 
       
      IF (IPRHAM.GE.JCONSIST) THEN

*** Do Omega+ * Omega

         CALL QGEMM(NREF,NREF,NREF,D1,
     &          'H','N',IPQTOQ(1,0),DOM(1),NREF,NREF,NZ,
     &          'N','N',IPQTOQ(1,0),DOM(1),NREF,NREF,NZ,
     &          D0,IPQTOQ(1,0),WORK(KUU),NREF,NREF,NZ)

       IF (IPRHAM.GE.6) THEN
          CALL HEADER(
     &'U1CHECK: Omega+ * Omega matrix',-1)
              WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
              CALL PRQMAT(WORK(KUU),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

*** Do R+ * R

         CALL QGEMM(NREF,NREF,NREF,D1,
     &          'H','N',IPQTOQ(1,0),R1(1),NREF2,NREF,NZ,
     &          'N','N',IPQTOQ(1,0),R1(1),NREF2,NREF,NZ,
     &          D0,IPQTOQ(1,0),WORK(KRR),NREF,NREF,NZ)

       IF (IPRHAM.GE.6) THEN
          CALL HEADER(
     &'U1CHECK: R+ * R matrix',-1)
              WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
              CALL PRQMAT(WORK(KRR),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF


*** Do (R+ * R) * Omega

         CALL QGEMM(NREF,NREF,NREF,D1,
     &          'N','N',IPQTOQ(1,0),WORK(KRR),NREF,NREF,NZ,
     &          'N','N',IPQTOQ(1,0),DOM(1),NREF,NREF,NZ,
     &          D0,IPQTOQ(1,0),WORK(KHEE),NREF,NREF,NZ)

*** DO Omega+ * (R+ * R) * Omega

         CALL QGEMM(NREF,NREF,NREF,D1,
     &          'H','N',IPQTOQ(1,0),DOM(1),NREF,NREF,NZ,
     &          'N','N',IPQTOQ(1,0),WORK(KHEE),NREF,NREF,NZ,
     &          D0,IPQTOQ(1,0),WORK(KRR),NREF,NREF,NZ)


*** Cout KRR = KRR + KUU <=> Full Ome+.Ome + Ome+.R+.R.Ome

       CALL DAXPY(N2EQ,D1,WORK(KUU), 1,WORK(KRR),1) 
       CALL HEADER(
     & 'U1CHECK: Unit U1+ * U1 matrix',-1)
       WRITE(LUPRI,'(3X,A,I1,A,I1)') '* fermion irp ',I,'/',NFSYM
       CALL TEST_UNIT_MTX(WORK(KRR),NREF,NREF,NREF,NZ)
       IF (IPRHAM.GE.6) THEN
        CALL PRQMAT(WORK(KRR),NREF,NREF,NREF,NREF,
     &              NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

****   Check the unity of U1+*U1  ****

        CALL QGEMM(NREF,NREF,NREF2,D1,
     &          'H','N',IPQTOQ(1,0),U1(1),NREF2,NREF,NZ,
     &          'N','N',IPQTOQ(1,0),U1(1),NREF2,NREF,NZ,
     &          D0,IPQTOQ(1,0),WORK(KUU),NREF,NREF,NZ)

        CALL HEADER(
     &'U1CHECK: Again unit U1+ * U1 matrix',-1)
       WRITE(LUPRI,'(3X,A,I1,A,I1)') '* fermion irp ',I,'/',NFSYM
       CALL TEST_UNIT_MTX(WORK(KUU),NREF,NREF,NREF,NZ)
       IF (IPRHAM.GE.6) THEN
        CALL PRQMAT(WORK(KUU),NREF,NREF,NREF,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

      ENDIF
      
      CALL DCOPY(N2OQ,H1(1),1,WORK(KH1),1)  ! Store full H1

C**************************************************************************
C         Construct transformation matrix  FTMAT = TMAT * U1
C**************************************************************************

      IF ( IPRHAM .GE. 6 ) THEN
        CALL HEADER(
     & ' U1CHECK: Transf. mtx TMAT 4c "xhi"- 4c |k>',-1)
       WRITE(LUPRI,'(3X,A,I2/)') '*** Fermion corep ',I
       CALL PRQMAT(TMAT(1),NFBAS(I,0),NORB(I),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

      CALL QGEMM(NFBAS(I,0),NESH(I),NORB(I),D1,
     &          'N','N',IPQTOQ(1,0),TMAT(1),NFBAS(I,0),NORB(I),NZ,
     &          'N','N',IPQTOQ(1,0),U1(1),NREF2,NREF,NZ,
     &          D0,IPQTOQ(1,0),FTMAT(1),NFBAS(I,0),NESH(I),NZ)

      IF ( IPRHAM .GE. 5 ) THEN
       CALL HEADER(
     &'U1CHECK: Final transf. mtx FTMAT 4c "xhi"-> 2c |k_e>',-1)
       WRITE(LUPRI,'(3X,A,I2/)') '*** Fermion corep ',I
       CALL PRQMAT(FTMAT(1),NFBAS(I,0),NESH(I),NFBAS(I,0),NESH(I),
     &               NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

*** Now do checking of final H2c in |k_e> ON basis
      IF ( IPRHAM .GE. JCONSIST ) THEN
       IF ( IPRHAM .GE. 5 ) THEN
        CALL HEADER(
     &  'U1CHECK: 4c H1 in |k> basis',-1)
        WRITE(LUPRI,'(3X,A,I2)') '*** Fermion corep ',I
        CALL PRQMAT(WORK(KH1),NORB(I),NORB(I),NORB(I),NORB(I),
     &               NZ,IPQTOQ(1,0),LUPRI)
       ENDIF
C=====================================================================
C           generate H1_eff into WORK(KHEE)
C        (remember NREF = NESH(I) and NREF2 = NORB(I) = 2*NREF)
C=====================================================================
            CALL QTRANS('AOMO','S',D0,NREF2,NREF2,NREF,NREF,
     &           WORK(KH1)  ,NREF2,NREF2,NZ ,IPQTOQ(1,0),
     &           WORK(KHEE) ,NREF ,NREF ,NZ ,IPQTOQ(1,0),
     &           U1(1)      ,NREF2,NREF ,NZT,IPQTOQ(1,0),
     &           U1(1)      ,NREF2,NREF ,NZT,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRHAM)

       IF ( IPRHAM .GE. 5 ) THEN
        CALL HEADER(
     &  'U1CHECK: Entering 2c H1(ee) in |k> basis',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(WORK(KHEE),NREF,NREF,NREF,NREF,
     &               NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

      IF ( IPRHAM .GE. 6) THEN
       CALL HEADER(
     &'U1CHECK: 4c bare nucleus in "xhi" AO basis',-1)
       WRITE(LUPRI,'(3X,A,I2)') '*** Fermion corep ',I
       CALL PRQMAT(HAO(1+I2BASX(I,I)),NFBAS(I,0),NFBAS(I,0),
     &               NTBAS(0),NTBAS(0),
     &               NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C Check FTMAT+ * H_4c(AO) * FTMAT = H1(ee) in |k> in giver fermion symmetry I

       CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),NREF,NREF,
     &           HAO(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ ,IPQTOQ(1,0),
     &           WORK(KHEE1) ,NREF ,NREF ,NZ ,IPQTOQ(1,0),
     &           FTMAT(1)    ,NFBAS(I,0),NREF ,NZ,IPQTOQ(1,0),
     &           FTMAT(1)    ,NFBAS(I,0),NREF ,NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRHAM)

      IF ( IPRHAM .GE. 5 ) THEN
        CALL HEADER(
     &  'U1CHECK: 2c H1(ee) again (CHECK IT!) in |k> basis',-1)
        WRITE(LUPRI,'(3X,A,I2)') '*** Fermion corep ',I
        CALL PRQMAT(WORK(KHEE1),NREF,NREF,NREF,NREF,
     &               NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C -------------------------------------------------------------------------------
C ... if it is  'spin-free from the beginning', eliminate spin-orbit components
C     from the constructed H1(ee) in |k> basis
C -------------------------------------------------------------------------------
      IF (SFRB) THEN
        WRITE(LUPRI,'(/2X,A/)')
     & 'U1CHECK: Elim. spin-orbit comp. "BEFORE" from  H1(ee) in |k>'
        CALL SPF_HK(WORK(KHEE1),IBEIG(IORB(I)+1+NPSH(I)),
     &            IBEIG(IORB(I)+1+NPSH(I)),NREF,NREF,NREF,NREF,NZ)
        IF (IPRHAM.GE.5) THEN
          CALL HEADER(
     &   'U1CHECK: 2c H1(ee) again SPIN-FREE in |k> basis',-1)
          WRITE(LUPRI,'(3X,A,I2)') '*** Fermion corep ',I
          CALL PRQMAT(WORK(KHEE1),NREF,NREF,NREF,NREF,
     &                  NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
      ENDIF

       CALL DAXPY(N2EQ,DM1,WORK(KHEE1),1,WORK(KHEE),1) 
       DN =  DNRM2(N2EQ,WORK(KHEE),1)

       CALL HEADER(
     &'U1CHECK: 2 times constructed H2c in ON |k_e> basis',-1)
       WRITE(LUPRI,'(2X,A,D17.9)')
     & 'Difference between 2 H2c in |k_e> basis:',DN

      ENDIF

      CALL MEMREL('U1CHECK',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('U1CHECK')
C      CALL QUIT( "IN U1CHECK...")
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck LOWDHAO */
      SUBROUTINE LOWDHAO(IPZ,ISYMM,H2AO,H4AO,
     &                FTM,FTMT,TM1,VMAT,TM2,
     &                IBEIG,JCONSIST,WORK,LWORK)
C****************************************************************************
C
C     PURPOSE:  Transform 2c Hamiltonian into ON Lowdin basis and diagonalize...
C               construct final transformation matrix 4c "xhi" -> 2c Lowdin ON basis 
C
C            THIS ROUTINE IS ONLY FOR THE INFINITE ORDER OF RELATIVITY !!!
C
C    On input:  H  - transformed 2c Hamilt. in AO basis "xhi_L"
C               H4AO - 4comp bare nucleus Ham. in "xhi" AO basis
C               FTM - transf.matrix 4c "xhi" -> 2c |k_e>
C               TM1 - trans.matrix 2c |k_e> - 2c "xhi_L"
C               VMAT - transf. from AO - Lowdin ON basis 'Theta' - pure real matrix !
C               TM2  - transf. from AO "xhi" to final ON basis "theta" (quatern.)
C               ISYMM - the fermion symmetry (1 or 2) 
C               IPZ - order of relativity (1,2 or 9 - DK1,DK2,DKinf)
C               JCONSIST - option for consistency checkings
C
C   On output: 
C
C  FTMT - final transformation matrix 4c "xhi" -> 2c SA AO "xhi" - IMPORTANT !!!
C
C  FTMT - final transformation matrix 4c "xhi" -> 2c Lowdin ON basis "Theta_L" - for check purposes
C
C      NOTE:  Entering HAO - be carefull with total size !!!
C
C      Written:  M.Ilias, Odense, 2002  
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00, DM1 = -1.00D00 )
#include <dcbgen.h>
#include <dcbdhf.h>
#include <dcbbas.h>
#include <dcbham.h>
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION H2AO(*),H4AO(*),FTM(*),FTMT(*),
     &         TM1(*),VMAT(*),TM2(*),IBEIG(*),WORK(*)
      LOGICAL PRIH2C

C ... routine is designed for infinite order of relativity!
      IF (IPZ.EQ.2.OR.IPZ.EQ.1) RETURN 

      CALL QENTER('LOWDHAO')
#include <memint.h>
C
      I = ISYMM

      N2FORBQ   = NFORB(I,0)*NFORB(I,0)*NZ
      N2FORBLQ  = NFORB(I,1)*NFORB(I,1)*NZ
      N2FB      = NFBAS(I,0)*NFBAS(I,0)
      N2BLQ     = NFBAS(I,0)*NFBAS(I,1)*NZ
      N2ORBI    = NORB(I)*NORB(I)
      N2ORBQ    = N2ORBI*NZ

      IF (I.EQ.2) THEN
       ISOVL=NFBAS(1,0)*NFBAS(1,0)
      ELSE
       ISOVL = 0
      ENDIF

      IF (NESH(I).NE.NFORB(I,1)) THEN 
       CALL QUIT( "NESH(I).NE.NFORB(I,1)!")
      ENDIF

C      NDH1AO = NFBAS(I,1)
      NDH1AO = NTBAS(0) ! Changed dimension

      CALL MEMGET('REAL',KCC,N2FORBQ,WORK,KFREE,LFREE) ! 
      CALL MEMGET('REAL',KHL,N2FORBLQ,WORK,KFREE,LFREE) ! 
      CALL MEMGET('REAL',KOVL,N2BBASX,WORK,KFREE,LFREE) ! 
      CALL MEMGET('REAL',KOVLI,N2FB,WORK,KFREE,LFREE) ! 

      IF (IPRHAM.GE.6) THEN
        CALL HEADER(
     &'LOWDHAO: Lowd.transf.matrix 4c xhi -> 4c ON Theta ',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(VMAT(1),NFBAS(I,0),NFORB(I,0),
     &                         NFBAS(I,0),NFORB(I,0),
     &                 1,IPQTOQ(1,0),LUPRI)
       ENDIF

C*****************************************************************************
C Print the H2c matrix elements with Sovlp(LL) elements into formatted file
C*****************************************************************************
C      LUH = 13
C      PRIH2C = .FALSE. ! Print integrals...
C      IF (PRIH2C.AND.NZ.EQ.4.AND.I.EQ.1) THEN
C        WRITE(LUPRI,'(3X,A)') 
C     &'Printing H2c matrix elements("xhi_L") into formatted file'
C        CALL GTOVLT(WORK(KOVL),D0,WORK(KFREE),LFREE)
C        OPEN(LUH,file="H2C.DAT",access="sequential",status="unknown")
C        DO J = 1, NFBAS(I,1)
C        DO K = 1, J
C         JKT = J + (K-1)*NFBAS(I,0)
C         KJT = K + (J-1)*NFBAS(I,0)
C         SJK = WORK(KOVL+JKT-1)
C         SKJ = WORK(KOVL+KJT-1)
C         HJK = H(J,K,1)
C         HKJ = H(K,J,1)
C         Write(LUH,"(2I4,2D26.15)") J,K,SJK,HJK
C         Write(LUH,"(2I4,2D26.15)") K,J,SKJ,HKJ
C        ENDDO
C        ENDDO
C        CLOSE(LUH,STATUS="KEEP")
C      ENDIF

****   Check the unity of C+ * S * C  ****

      IF (IPRHAM.GE.JCONSIST) THEN

       CALL GTOVLT(WORK(KOVL),D1,WORK(KFREE),LFREE)
       CALL DCOPY(N2FB,WORK(KOVL+ISOVL),1,WORK(KOVLI),1)  ! Store full S(ISYMM)

       IF (IPRHAM.GE.7) THEN
         CALL HEADER(
     &'LOWDHAO: Full SOVLP matrix in SO basis "xhi"',-1)
         WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
         CALL PRQMAT(WORK(KOVLI),NFBAS(I,0),NFBAS(I,0),
     &         NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),LUPRI)
       ENDIF

      CALL DZERO(WORK(KCC),N2FORBQ)

C      Write(*,*) "N2BBASX,NTBAS(0)*NTBAS(0),ISOVL:",
C     &  N2BBASX,NTBAS(0)*NTBAS(0),ISOVL

      CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0), ! WORK(KCC) = C+ . S . C 
     &               NFORB(I,0),NFORB(I,0),
     &    WORK(KOVLI),NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),
     &    WORK(KCC),NFORB(I,0),NFORB(I,0),NZ,IPQTOQ(1,0),
     &    VMAT(1),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &    VMAT(1),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &        WORK(KFREE),LFREE,IPRHAM)

       CALL HEADER(
     &'LOWDHAO: Full(L+S) unit C+ * S * C matrix',-1)
       WRITE(LUPRI,'(3X,A,I1,A,I1)') '* fermion corep ',I,'/',NFSYM
       CALL TEST_UNIT_MTX(WORK(KCC),NFORB(I,0),
     &   NFORB(I,0),NFORB(I,0),NZ)
       IF (IPRHAM.GE.6) THEN
         CALL PRQMAT(WORK(KCC),NFORB(I,0),NFORB(I,0),
     &                         NFORB(I,0),NFORB(I,0),
     &                NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL DZERO(WORK(KCC),N2FORBQ)

****   Check the unity of LL part of C+ * S * C  ****

      CALL QTRANS('AOMO','S',D0,NFBAS(I,1),NFBAS(I,1),
     &               NFORB(I,1),NFORB(I,1),
     &    WORK(KOVLI),NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),
     &    WORK(KCC),NFORB(I,0),NFORB(I,0),NZ,IPQTOQ(1,0),
     &    VMAT(1),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &    VMAT(1),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &        WORK(KFREE),LFREE,IPRHAM)

      CALL HEADER(
     &'LOWDHAO: (LL) unit C+ * S * C matrix',-1)
      WRITE(LUPRI,'(3X,A,I1,A,I1)') '*** Fermion corep ',I,'/',NFSYM
      CALL TEST_UNIT_MTX(WORK(KCC),NFORB(I,1),
     &                   NFORB(I,0),NFORB(I,0),NZ) ! 
      IF (IPRHAM.GE.6) THEN
         CALL PRQMAT(WORK(KCC),NFORB(I,1),NFORB(I,1),
     &                         NFORB(I,0),NFORB(I,0),
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

      ENDIF

C
C        Transform the HAO from "xhi" to the "Theta" basis
C

        Write(LUPRI,'(/2X,A)')
     & 'LOWDHAO: Transforming H2c from AO "xhi" to ON Lowdin_L basis'

        CALL QTRANS('AOMO','S',D0,NFBAS(I,1),NFBAS(I,1),
     &               NFORB(I,1),NFORB(I,1),
     &  H2AO(1),NDH1AO,NDH1AO,NZ,IPQTOQ(1,0),
     &  WORK(KHL),NFORB(I,1),NFORB(I,1),NZ,IPQTOQ(1,0),
     &  VMAT(1),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &  VMAT(1),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &        WORK(KFREE),LFREE,IPRHAM)

       IF (IPRHAM.GE.5) THEN
        CALL HEADER(
     &'LOWDHAO: H(2c) in ON Lowdin basis "Theta_L"',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(WORK(KHL),NFORB(I,1),NFORB(I,1),
     &                        NFORB(I,1),NFORB(I,1),
     &                        NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

*  Construct final transformation matrix FTMT = FTM*TM1*VMAT

      CALL MEMGET('REAL',KFX1,N2BLQ,WORK,KFREE,LFREE) ! 

       IF (IPRHAM.GE.6) THEN
        CALL HEADER(
     &'LOWDHAO: Transf.matrix 4c "xhi" -> 2c |k_e>',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(FTM(1), NFBAS(I,0),NESH(I),
     &                      NFBAS(I,0),NESH(I),
     &                       NZ,IPQTOQ(1,0),LUPRI)
        CALL HEADER(
     &'LOWDHAO: Transf.matrix 2c |k_e> -> 2c "xhi_L"',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(TM1(1),  NFBAS(I,0),NESH(I),
     &                      NFBAS(I,0),NESH(I),
     &                        NZ,IPQTOQ(1,0),LUPRI)

       ENDIF

      CALL QGEMM(NFBAS(I,0),NFBAS(I,1),NESH(I),D1,  !  FTM*TM1+ = WORK(KFX1)
     &          'N','N',IPQTOQ(1,0),
     &           FTM(1),NFBAS(I,0),NESH(I),NZ,
     &          'H','N',IPQTOQ(1,0),
     &           TM1(1),NFBAS(I,0),NESH(I),NZ,
     &          D0,IPQTOQ(1,0),
     &          WORK(KFX1),NFBAS(I,0),NFBAS(I,1),NZ)

C****************************************************************
C
C That we want: transformation matrix 4c AO "xhi" -> 2c AO "xhi"
C
C****************************************************************
       IF (IPRHAM.GE.6) THEN
        CALL HEADER(
     &'LOWDHAO: Transf.matrix 4c "xhi" -> 2c "xhi"',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(WORK(KFX1),NFBAS(I,0),NFBAS(I,1),
     &                        NFBAS(I,0),NFBAS(I,1),
     &                        NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C ... Check the unity U+ * U of WORK(KFX1)

C ...   WORK(KFX1)* WORK(KFX1)+ = FTMT

      CALL QGEMM(NFBAS(I,1),NFBAS(I,1),NFBAS(I,0),D1,  
     &     'H','N',IPQTOQ(1,0),
     &     WORK(KFX1),NFBAS(I,0),NFBAS(I,1),NZ,
     &     'N','N',IPQTOQ(1,0),
     &     WORK(KFX1),NFBAS(I,0),NFBAS(I,1),NZ,
     &     D0,IPQTOQ(1,0),
     &     FTMT(1),NFBAS(I,1),NFBAS(I,1),NZ)

       IF (IPRHAM.GE.7) THEN ! U+.U not necessary 1...
        CALL HEADER(
     &'LOWDHAO: 4c/2c "xhi",  U+.U',-1)
        WRITE(LUPRI,'(3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(FTMT(1),NFBAS(I,1),NFBAS(I,1),
     &                        NFBAS(I,1),NFBAS(I,1),
     &                        NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C*****************************************************
C  Now do U+ H4c(AO)  U = H2c(AO) and check it !
C*****************************************************

       IF (IPRHAM.GE.6) THEN
        CALL HEADER(
     &'LOWDHAO: 4c bare nucleus in "xhi" AO basis',-1)
        WRITE(LUPRI,'(3X,A,I2)') '*** Fermion corep ',I
        CALL PRQMAT(H4AO(1+I2BASX(I,I)),NFBAS(I,0),NFBAS(I,0),
     &                        NTBAS(0),NTBAS(0),
     &                        NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       ISIZE = NFBAS(I,1)*NFBAS(I,1)*NZ

       CALL DZERO(FTMT(1),ISIZE)

C ---------------------------------------------------------------
C .... do trasformation 4c bare nucleus AO "xhi" -> 2c AO "xhi"
C ---------------------------------------------------------------
       CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &          NFBAS(I,1),NFBAS(I,1),
     &  H4AO(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &  FTMT(1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &  WORK(KFX1),NFBAS(I,0),NFBAS(I,1),NZ,IPQTOQ(1,0),
     &  WORK(KFX1),NFBAS(I,0),NFBAS(I,1),NZ,IPQTOQ(1,0),
     &    WORK(KFREE),LFREE,IPRHAM)

       IF (IPRHAM.GE.5) THEN
        CALL HEADER(
     &'LOWDHAO: H(2c) again in AO basis "xhi"',-1)
        WRITE(LUPRI,'(3X,A,I2)')'*** Fermion corep ',I
        CALL PRQMAT(FTMT(1),NFBAS(I,1),NFBAS(I,1),
     &                        NTBAS(0),NTBAS(0),
     &                     NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     &'LOWDHAO: Entering H(2c) in AO basis "xhi"',-1)
        WRITE(LUPRI,'(3X,A,I2)')'*** Fermion corep ',I
        CALL PRQMAT(H2AO(1),NFBAS(I,1),NFBAS(I,1),
     &           NDH1AO,NDH1AO,NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C ... calculate the difference between H(2c) in AO "xhi"
C  Cout FTMT(1) = H2AO(1) - FTMT(1);
C ....
       CALL DAXPY(ISIZE,DM1,H2AO(1),1,FTMT(1),1) 
       DN =  DNRM2(ISIZE,FTMT(1),1)

       WRITE(LUPRI,'(/2X,A,D12.5/)')
     & 'LOWDHAO: Difference between various constructed'//
     & ' H2c in AO basis:',DN

C 000000000000000000000000000000000000000000000000000000000000000000000000000

      CALL QGEMM(NFBAS(I,0),NFORB(I,1),NFBAS(I,1),D1,  !   WORK(KFX1)*VMAT = FTMT
     &          'N','N',IPQTOQ(1,0),
     &           WORK(KFX1),NFBAS(I,0),NFBAS(I,1),NZ,
     &          'N','N',IPQTOQ(1,0),
     &           VMAT(1),NFBAS(I,0),NFORB(I,0),1,
     &          D0,IPQTOQ(1,0),
     &          FTMT(1),NFBAS(I,0),NFORB(I,1),NZ)

       IF (IPRHAM.GE.6) THEN
        CALL HEADER(
     &'LOWDHAO: FINAL trans.mtx 4c "xhi" -> 2c ON "Theta_L"',-1)
        WRITE(LUPRI,'(3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(FTMT(1),NFBAS(I,0),NFORB(I,1),
     &                        NFBAS(I,0),NFORB(I,1),
     &                        NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C       IF (IPRHAM.GE.6) THEN
C        CALL HEADER(
C     &'LOWDHAO: 4c bare nucleus in "xhi" AO basis',-1)
C        WRITE(LUPRI,'(3X,A,I2)') '*** Fermion corep ',I
C        CALL PRQMAT(H4AO(1+I2BASX(I,I)),NFBAS(I,0),NFBAS(I,0),
C     &                        NTBAS(0),NTBAS(0),
C     &                        NZ,IPQTOQ(1,0),LUPRI)
C       ENDIF

C        CALL MEMREL('KF1',WORK,KWORK,KF1,KFREE,LFREE)

        CALL MEMGET('REAL',KF1,N2FORBLQ,WORK,KFREE,LFREE) ! KF1 and KHL should be identical !

        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &               NFORB(I,1),NFORB(I,1),
     &    H4AO(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &    WORK(KF1),NFORB(I,1),NFORB(I,1),NZ,IPQTOQ(1,0),
     &    FTMT(1),NFBAS(I,0),NFORB(I,1),NZ,IPQTOQ(1,0),
     &    FTMT(1),NFBAS(I,0),NFORB(I,1),NZ,IPQTOQ(1,0),
     &        WORK(KFREE),LFREE,IPRHAM)

C ... remove spin-orbit dependent terms from H2 in Lowdin basis 
C  .. note: NFORB(I,1)=NESH(I)=NPSH(I)
C      IF (SPINFR) THEN
C        CALL SPF_HK(WORK(KF1),IBEIG(IORB(I)+NPSH(I)+1),
C    &    IBEIG(IORB(I)+NPSH(I)+1),
C    &    NESH(I),NESH(I),NESH(I),NESH(I),NZ )
C      ENDIF

       IF (IPRHAM.GE.5) THEN
        CALL HEADER(
     &'LOWDHAO: H(2c) again (CHECK IT!) in ON Lowd basis "Theta_L"',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(WORK(KF1),NFORB(I,1),NFORB(I,1),
     &                        NFORB(I,1),NFORB(I,1),
     &                        NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C  Check the differences of H(2c) in ON Lowd L basis set
C  Cout WORK(KF1) = WORK(KF1) - WORK(KHL); WORK(KHL) is used for diagonalization!
       CALL DAXPY(N2FORBLQ,DM1,WORK(KHL),1,WORK(KF1),1) 
       DN =  DNRM2(N2FORBLQ,WORK(KF1),1)
       WRITE(LUPRI,'(/2X,A,D12.5/)')
     &'Difference between various constructed H2c in Lowdin basis:',DN
       
       IF (IPRHAM.GE.7) THEN
          CALL HEADER(
     &'LOWDHAO:Transf. mtx 4c "xhi"->final 4c ON basis "theta" ',-1)
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(TM2(1),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C
C  Dirac bare nucleus matrix in "xhi" basis transform into ON "theta" basis
C
      CALL MEMGET('REAL',KH2ORB,N2ORBQ,WORK,KFREE,LFREE) ! 

      IF(NORB(I).GT.0) THEN
       CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &   NORB(I),NORB(I),
     &   H4AO(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &   WORK(KH2ORB),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &   TM2(1),NFBAS(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &       TM2(1),NFBAS(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)
      ENDIF

      CALL MEMGET('REAL',KEVAL,NFORB(I,1),WORK,KFREE,LFREE) ! For eigenvalues H2c
      CALL MEMGET('REAL',KEVEC,N2FORBLQ,WORK,KFREE,LFREE) ! For eigenvectors H2c

      CALL MEMCHK('in LOWDHAO',WORK,1)
      CALL QDIAG(NZ,NFORB(I,1),WORK(KHL),NFORB(I,1),NFORB(I,1),   ! Diagonalize
     &              WORK(KEVAL),0,          ! Eigenvalues
     &              WORK(KEVEC),NFORB(I,1),NFORB(I,1),  ! Eigenvectors not necessary
     &              WORK(KFREE),LFREE,IERR)
      IF (IERR.NE.0) CALL QUIT("LOWDHAO:DIAGONALIZATION FAILED !!!")

      CALL MEMGET('REAL',K4EVAL,NORB(I),WORK,KFREE,LFREE) ! For eigenvalues H4c
      CALL MEMGET('REAL',K4EVEC,N2FORBQ,WORK,KFREE,LFREE) ! For eigenvectors H4c

C ... diagonalize dirac bare nucleus Hamiltonian in "theta" basis set...
      CALL QDIAG(NZ,NORB(I),WORK(KH2ORB),NORB(I),NORB(I),   
     &              WORK(K4EVAL),0,          ! Eigenvalues
     &              WORK(K4EVEC),NORB(I),NORB(I),  ! Eigenvectors not necessary
     &              WORK(KFREE),LFREE,IERR)
      IF (IERR.NE.0) CALL QUIT("LOWDHAO:BAD DIAGONALIZATION!")

**** In RSKINB  framework....  *****

      IF (IPRHAM.GE.0) THEN
       CALL HEADER(
     &'LOWDHAO: Electr. eigenv. of H(2c)/H(4c) in "theta" basis',-1)
       WRITE(LUPRI,'(3X,A,I2/)') '*** Fermion corep ',I

       WRITE(LUPRI,'(3X,I4,F25.10)')
     &  (J-NPSH(I)-1,WORK(K4EVAL+J-1),J=1,NPSH(I))

       DO J =1, NFORB(I,1)
        DIFF = WORK(KEVAL+J-1) - WORK(K4EVAL+NPSH(I)+J-1)
        WRITE(LUPRI,'(3X,I3,2F25.10,3X,D12.5)') 
     &    J,
     &    WORK(K4EVAL+NPSH(I)+J-1),
     &    WORK(KEVAL+J-1),
     &    DIFF
       ENDDO
      ENDIF

      CALL DCOPY(N2BLQ,WORK(KFX1),1,FTMT(1),1)

      CALL MEMREL('LOWDHAO',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('LOWDHAO')

C      IF (I.EQ.2) CALL QUIT( "LOWDHAO")

 789  RETURN
      END


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck TEST_UNIT_MTX */
      SUBROUTINE TEST_UNIT_MTX(S,N,LR,LC,NZ)
C****************************************************************************
C
C  PURPOSE: Check the enterning matrix - is it unit mtx ?
C             
C****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00, THRNULL = 1.00D-09 )

      DIMENSION S(LR,LC,NZ)

      SUM_DE  = SUM_DIAG_EL(S,N,LR,LC,NZ)
      SUM_DEN = SUM_DE/DFLOAT(N)
      SUM_ODE = SUM_OFFDIAG_EL(S,N,LR,LC,NZ)

      WRITE(LUPRI,'(/2X,A,D12.5)') 
     &'Sum of diagonal elements /N:',SUM_DEN
      WRITE(LUPRI,'(2X,A,D12.5)') 
     &'Sum of off-diagonal elements :',SUM_ODE 

      IF (DABS(SUM_DEN-D1).GT.THRNULL) 
     & WRITE(LUPRI,'(2X,A)') 
     &'WARNING, SUM OF DIAG.ELEM/N DIFFERS TOO MUCH FROM 1!'

      IF (DABS(SUM_ODE-D0).GT.THRNULL) 
     & WRITE(LUPRI,'(2X,A)') 
     &'WARNING, SUM OF OFF_DIAG.ELEM DIFFERS TOO MUCH FROM 0!'
 
      RETURN
      END

      REAL*8 FUNCTION SUM_DIAG_EL(S,N,LR,LC,NZ)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (  D0 = 0.00D00 )
      DIMENSION S(LR*LC*NZ)

      PX = D0
      DO I = 1, N
       II = I + (I-1)*LR
       PX = PX + S(II)
      ENDDO

      SUM_DIAG_EL = PX

      RETURN
      END

      REAL*8 FUNCTION SUM_OFFDIAG_EL(S,N,LR,LC,NZ)
#include <implicit.h>
      PARAMETER ( D0 = 0.00D00 )
      DIMENSION S(LR*LC*NZ)

      PX = D0
      DO I = 2, N
      DO J = 1, I-1
       IJ = I + (J-1)*LR
       JI = J + (I-1)*LR
       PX = PX + S(IJ) + S(JI) 
      ENDDO
      ENDDO

      SUM_OFFDIAG_EL = PX

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck TEST_DIAG_MTX */
      SUBROUTINE TEST_DIAG_MTX(S,N,LR,LC,NZ)
C****************************************************************************
C
C  PURPOSE: Check the enterning matrix - is it a diagonal matrix ?
C
C  Diagonal matrix should have diagonal elements between 0 and 1
C             
C****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00, THRNULL = 1.00D-09 )

      DIMENSION S(LR,LC,NZ)

      SUM_DE  = SUM_DIAG_EL(S,N,LR,LC,NZ)
      SUM_DEN = SUM_DE/DFLOAT(N)
      SUM_ODE = SUM_OFFDIAG_EL(S,N,LR,LC,NZ)

      WRITE(LUPRI,'(/2X,A,D12.5)') 
     &'Sum of diagonal elements /N :',SUM_DEN
      WRITE(LUPRI,'(2X,A,D12.5)') 
     &'Sum of off-diagonal elements :',SUM_ODE 

      IF (.NOT.(SUM_DEN.GT.THRNULL.AND.SUM_DEN.LT.D1)) 
     & WRITE(LUPRI,'(2X,A)') 
     &'WARNING, SUM OF DIAG.ELEM/N DIFFERS IS NOT IN RANGE (0;1)!'

      IF (DABS(SUM_ODE-D0).GT.THRNULL) 
     & WRITE(LUPRI,'(2X,A)') 
     &'WARNING, SUM OF OFF_DIAG.ELEM DIFFERS TOO MUCH FROM 0!'

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck EXTR_BSS_INFO */
      SUBROUTINE EXTR_BSS_INFO(IBSS,IPZ,SPINFR,DK2SO1,SF0,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C PURPOSE:
C ========
C     Set control variables based on IBSS coding
C     IBSS = xyz
C
C        x - .ne. 0: define spinfree before (as 4-comp)
C            .eq. 0: define spinfree at end
C            ... x is active with the SPINFREE option only !!!
C           
C        y - order of spin-orbit term:
C            0 - spinfree, 1 - first order for DK2, 2+ - full.
C         z - order of Douglas-Kroll transformation
C            (0, 1, 2, or 9), where 9 is full.
C
C     Pick up the x,y,z variables from IBSS = xyz
C
C On input: IBSS, SPINFR, IPRINT
C On output: IPZ - order of Hamiltonian
C            DK2SO1,SF0
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include <implicit.h>
#include <priunit.h>

      LOGICAL SPINFR,DK2SO1,SF0

      CALL EXTR_DIG(IBSS,IPX,IPY,IPZ,IPRINT)

      WRITE (LUPRI,'(/A,I2)') 'MAKE_DKH: Order of Douglas-Kroll :',IPZ
      IF (SPINFR.AND.IPY.GE.1) THEN 
       WRITE(LUPRI,'(2X,A)') '"SPINFR" option, no spin-orbit term!' 
       IPY=0
      ENDIF
      IF (IPY .EQ. 0) THEN
         IF (.NOT.SPINFR) CALL QUIT('You forgot to specify .SPINFR')
         DK2SO1 = .FALSE.
      ELSE IF (IPY .EQ. 1 .AND. IPZ .EQ. 2 .AND.(.NOT.SPINFR)) THEN
         DK2SO1 = .TRUE.
         WRITE (LUPRI,'(/A)') 'MAKE_DKH: first order spin orbit'
      ELSE IF (.NOT.SPINFR) THEN
         DK2SO1 = .FALSE.
         WRITE (LUPRI,'(/A)') 'MAKE_DKH: full spin orbit'
      END IF
      IF (IPX .GE. 1 .AND. SPINFR) THEN  
         SF0 = .TRUE.
         WRITE (LUPRI,'(/A)') 'MAKE_DKH: spinfree "before"'
      ELSE
         SF0 = .FALSE.
         IF (SPINFR) WRITE (LUPRI,'(/A)') 'MAKE_DKH: spinfree "after"'
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck EXTR_DIG */
      SUBROUTINE EXTR_DIG(IBSS,IPX,IPY,IPZ,IPRINT)
C****************************************************************************
C
C  PURPOSE: For entering integer number of type xyz,
C            extract x,y,z into IPX,IPY,IPZ
C             
C  Written: M Ilias, Odense, aug. 2002             
C
C****************************************************************************
#include <implicit.h>
#include <priunit.h>

      IF (IBSS.LT.0.OR.IBSS.GE.1000)
     & CALL QUIT('EXTR_DIG: The entering number not in range 0-999! ')
  
      IBSS0 = IBSS
      IF (IBSS.GE.100) THEN 
       IBSS = MOD(IBSS,100)
       IPX  = IBSS0/100
       IPY = IBSS/10
       IPZ = MOD(IBSS,10)
      ELSE IF (IBSS.GE.10) THEN
       IPX = 0
       IPY = IBSS/10
       IPZ = MOD(IBSS,10)
      ELSE IF (IBSS.GE.0) THEN
       IPX = 0
       IPY = 0
       IPZ = IBSS
      ENDIF

      IF (IPRINT.GE.5) THEN
       WRITE (LUPRI,'(/A,I3,A,3I2)')
     &'EXTR_DIG: IBSS=xyz=',IBSS0,
     &'individual  x,y,z:',IPX,IPY,IPZ
      ENDIF

      IBSS1 = 100*IPX + 10*IPY + IPZ
      IF (IBSS1.NE.IBSS0) 
     & CALL QUIT('EXTR_DIG: Error in extraction of x,y,z from xyz!')

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck spf_hk */
      SUBROUTINE SPF_HK( HK,IBS1, IBS2, NROW,NCOL,LROW,LCOL,NZ )
C****************************************************************************
C
C  PURPOSE: Eliminate spin-orbit off-diagonal elements of entering matrix
C           in MO basis 
C
C  On input: HK - quaternion hermitian matrix
C            IBS1 - boson symetries of first index
C            IBS2 - boson symetries of second index
C
C  On output: HK - the same matrix without spin-orbit components
C             
C  Written: M Ilias, Odense, aug. 2002             
C             
C****************************************************************************
#include <implicit.h>
      PARAMETER ( D0 = 0.00D00 )
#include <dcbham.h>
#include <priunit.h>
      DIMENSION HK(LROW,LCOL,NZ),IBS1(*),IBS2(*)

      CALL QENTER('SPF_HK')

      IF (NROW.NE.NCOL) THEN
        WRITE(LUPRI,'(2X,A)') 'SPF_HK: Warning, NROW<>NCOL!'
      ENDIF

      IF (IPRHAM. GE. 5) THEN
       WRITE(LUPRI,'(2X,A)')
     & 'SPF_HK: Entering |k> boson symmetries of index 1:'
       DO I = 1, NROW
        WRITE(LUPRI,'(5X,A,I3,3X,A,I2)')
     &  'orbital: ',I,'boson symmetry:',IBS1(I)+1
       ENDDO
       WRITE(LUPRI,'(2X,A)')
     & 'SPH_HK: Entering |k> boson symmetries of index 2:'
       DO I = 1, NCOL
        WRITE(LUPRI,'(5X,A,I3,3X,A,I2)')
     &  'orbital: ',I,'boson symmetry:',IBS2(I)+1
       ENDDO
      ENDIF

      IZ = 1
      DO I=1,NROW
      DO J=1,NCOL

C        WRITE(LUPRI,"(2I3,2X,D14.6,2X,2I2)")
C     &   I,J,HK(I,J,IZ),IBS1(I)+1,IBS2(J)+1

        IF (IBS1(I).NE.IBS2(J)) HK(I,J,IZ)=D0

      ENDDO
      ENDDO

      IF (NZ.GT.1) THEN
       DO IZ = 2, NZ
       DO I=1,NROW
       DO J=1,NCOL
        HK(I,J,IZ)=D0
       ENDDO
       ENDDO
       ENDDO
      ENDIF

      CALL QEXIT('SPF_HK')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck  */
      SUBROUTINE F2PICTR(V,WORK,LWORK)
C****************************************************************************
C
C  PURPOSE: Do picture change of the two-electron Fock matrix: U+ V U
C           arising during iterations 
C
C
C                   It's only for testing purposes !
C
C
C  On input: V(N2BBASXQ*NFMAT) - two-el.part of Fock matrix (LL block only)
C
C  On output: V with pict. change transformed LL block 
C             
C  Written: M Ilias, Odense,  2002             
C             
C****************************************************************************
#include <implicit.h>
#include <dcbham.h>
#include <dcbgen.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbdhf.h>
#include <pgroup.h>
#include <cbihr1.h>
#include <codata.h>
#include <priunit.h>

      PARAMETER ( D0 = 0.00D00, D1 = 1.00D00 )
      DIMENSION V(*),WORK(LWORK)
      LOGICAL FNDLAB,DOPCTR,DOCOPY

      CALL QENTER('F2PICTR')
#include <memint.h> 

*      DOPCTR = .TRUE.
      DOPCTR = .FALSE.

      IF (.NOT.BSS.OR.(.NOT.DOPCTR)) GOTO 345

      IF (IPRHAM.GE.5) THEN

       CALL HEADER(
     &'F2PICTR: Two-el.Fock matrixes in AO basis "xhi"',-1)
       DO IFMAT = 1, NFMAT
        WRITE(LUPRI,'(10X,A,I1,A,I1)')
     & 'Matrix number: ',IFMAT,'/',NFMAT
        CALL PRQMAT(V(1+(IFMAT-1)*NFMAT),
     &   NTBAS(0),NTBAS(0),
     &   NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
       ENDDO
      ENDIF

C *** Read transform.matrix from DFTMAT file

      CALL MEMGET('REAL',KU,N2BBASXQ,WORK,KFREE,LFREE)
      CALL OPNFIL(LUTMAT,'DFTMAT','OLD','F2PICTR')
      REWIND LUTMAT
      IF (FNDLAB('U_PICTRM',LUTMAT)) THEN
       CALL READT(LUTMAT, N2BBASXQ, WORK(KU))
       WRITE(LUPRI,'(/2X,A,A)')
     & 'F2PICTR: Picture change tranform. matrix',
     & ' for BSS from the DFTMAT file read.'
      ELSE
       CALL QUIT('F2PICTR:Error, label U_PICTRM not found !')
      END IF
      CLOSE(LUTMAT,STATUS='KEEP')

C ... Do the picture change transformation of the Fock 2 el.matrix

      IF (IPRHAM.GE.5) THEN
       CALL HEADER(
     &'F2PICTR: Transf.matrix 4c "xhi" -> 2c "xhi"',-1)
       IKU = KU
       DO I = 1, NFSYM
        WRITE(LUPRI,'(/3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(WORK(IKU),
     &         NFBAS(I,0),NFBAS(I,1),
     &         NFBAS(I,0),NFBAS(I,1),NZ,
*     &         NTBAS(0),NTBAS(0),NZ, ! Wrong!
     &       IPQTOQ(1,0),LUPRI)
        IKU = IKU + NFBAS(I,0)*NFBAS(I,1)*NZ
       ENDDO
      ENDIF

*      Stop "0"

C ... Do the picture change transformation of the Fock 2 el.matrix

      DO IFMAT = 1, NFMAT ! Run over the two-el. Fock matrixes....

      IKU = KU
      IKV = 0
      IKV1 = 0

      DO I = 1, NFSYM

       CALL MEMGET('REAL',KVTLL,NFBAS(I,1)*NFBAS(I,1)*NZ,
     &              WORK,KFREE,LFREE)

       IKVTLL = KVTLL
       CALL QTRANS('AOMO','S',D0,
     &     NFBAS(I,0),NFBAS(I,0),
     &      NFBAS(I,1),NFBAS(I,1),
     &  V(1+I2BASX(I,I)+(IFMAT-1)*N2BBASXQ),
     &  NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &  WORK(KVTLL),
     &   NFBAS(I,1),NFBAS(I,1),NZ,IPQTOQ(1,0),
     &  WORK(IKU),NFBAS(I,0),NFBAS(I,1),NZ,IPQTOQ(1,0),
     &  WORK(IKU),NFBAS(I,0),NFBAS(I,1),
     &   NZ,IPQTOQ(1,0),
     &  WORK(KFREE),LFREE,IPRHAM)


      IF (IPRHAM.GE.5) THEN
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &'F2PICTR: Transformed U+ V U: ',I,'/',NFSYM
        CALL PRQMAT(WORK(IKVTLL),
     &              NFBAS(I,1),NFBAS(I,1),
     &              NFBAS(I,1),NFBAS(I,1),
     &              NZ,IPQTOQ(1,0),LUPRI)
      ENDIF

** Now copy the LL block into the old V mtx

      DOCOPY = .TRUE.
*      DOCOPY = .FALSE.
      If (DOCOPY) THEN
       DO I1 = 1, NFBAS(I,1) ! Run over columns
        DO IZ = 1, NZ ! Run over matrixes

        CALL DCOPY(NFBAS(I,1),
     &   WORK(KVTLL+(I1-1)*NFBAS(I,1)+
     &       (IZ-1)*NFBAS(I,1)*NFBAS(I,1)),1,
     &  V(1+(I1-1)*NFBAS(I,0)+
     &    IKV+(I1-1)*2*IKV1+2*(I1-1)*(NFSYM-1)+
     &      (I-1)*NTBAS(0)*IKV +
     &      (IZ-1)*NTBAS(0)*NTBAS(0) +
     &     (IFMAT-1)*N2BBASXQ),1)

        ENDDO
       ENDDO

       IF (IPRHAM.GE.5) THEN
        WRITE(LUPRI,'(/10X,A,I1,A,I1)')
     &'F2PICTR: Copied into V, FMAT no:',IFMAT,'/',NFMAT
       CALL PRQMAT(V(1+(IFMAT-1)*NFMAT),
     &   NTBAS(0),NTBAS(0),
     &   NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF
      ENDIF

        IKU = IKU + NFBAS(I,0)*NFBAS(I,1)*NZ
        IKV1 = IKV1 + NFBAS(I,1)
        IKV = IKV + NFBAS(I,0)

        CALL MEMREL('F2PICTR',WORK,KWORK,KVTLL,KFREE,LFREE)

       ENDDO
      ENDDO

*      STOP "End"

 345  CALL QEXIT('F2PICTR')

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck  */
      SUBROUTINE DOPCTRA(INDXPR,PRPAO,MZ,TREV,WORK,LWORK)
C****************************************************************************
C
C  PURPOSE: Do picture change transformation  of the entering AO 
C           one-electron matrix of some
C           FOUR-COMPONENT property operator
C
C  On input: PRPAO 4component operator in AO basis
C             - either full or 1 comp.
C            INDXPR - index of property to be transformed 
C
C  On output: PRPAO picture change transformed operator 
C             matrix (LL block) of AO basis
C
C  Written: M Ilias, Prievidza, June 2004
C      
C****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D0 = 0.00D00, D1 = 1.00D00 )
#include "mxcent.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbprp.h"
#include "dcbxpr.h"

      DIMENSION WORK(*),PRPAO(*)
      LOGICAL FNDLAB
      CHARACTER  PROPNAME*16,TREV*1

      CALL QENTER('DOPCTRA')
#include "memint.h"

C ... check if one needs to do the pic.change transf.
      IF (.NOT.DOPCT) GOTO 789

      IPRINT = IPRPRP
      IREP = IPRPSYM(INDXPR) - 1

C ... determine the time reversal symmetry
C     of the property operator
      ITIM = IPRPTIM(INDXPR)

C ,,, read (load) the picture change transformation matrix:
C       4c AO-> 2c AO(LL)
      CALL MEMCHK('in DOPCT',WORK,1)
C ... array for the picture change transformation matrix (PCTM)...
      CALL MEMGET('REAL',KPCTM,N2BBASXQ,WORK,KFREE,LFREE)

      CALL OPNFIL(LUBSS,'BSSMAT','OLD','DOPCT')
      REWIND LUBSS
      IF (FNDLAB('U_PICTRM',LUBSS)) THEN
         CALL READT(LUBSS,N2BBASXQ,WORK(KPCTM))
         IF (IPRPRP.GE.2) WRITE(LUPRI,'(/2X,A)')
     &   'DOPCTRA: Picture change transformation matrix was read.'
      ELSE
        CALL QUIT('DOPCTRA:Error, U_PICTRM not found in BSSMAT file!')
      END IF
      CLOSE (LUBSS,STATUS='KEEP')

      IF (IPRINT.GE.5) THEN
       CALL HEADER('Output from DOPCTRA',-1)
       WRITE(LUPRI,'(A,A16,/,10X,A,I2,A,I1)')
     & 'Total entering FOUR-COMPONENT AO-matrix of property ',
     & PRPNAM(INDXPR),' time rev.symm:',ITIM,
     & ' number of matrixes:',MZ
       CALL PRQMAT(PRPAO(1),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),MZ,
     &            IPQTOQ(1,IREP),LUPRI)

C ... print out symmetry blocked picture change transformation matrix

       IPCTM1 = KPCTM
         WRITE(LUPRI,'(/A)')
     & 'SYMMETRY BLOCKED picture change transformation matrix: '
       DO I1 = 1, NFSYM
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &      '*** Fermion corep ',I1,'/',NFSYM
C ... prints OK ...
         CALL PRQMAT(WORK(IPCTM1),NFBAS(I1,0),NFBAS(I1,1),
     &     NFBAS(I1,0),NFBAS(I1,1),NZ,IPQTOQ(1,0),LUPRI)
        IPCTM1 = IPCTM1 + NFBAS(I1,0)*NFBAS(I1,1)*NZ
       ENDDO

      ENDIF

      CALL MEMGET('REAL',KPCTO,N2BBASX*MZ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KPCTO),N2BBASX*MZ)
      CALL MEMCHK('in DOPCT3',WORK,1)

      IOPSY = JBTOF(IREP,1)
      IPCTM1 = KPCTM
      IPCTM2 = KPCTM
      DO I1 = 1, NFSYM
         I2 = MOD(I1+IOPSY,2) + 1

C ----------------------------------------------------------
C     Do the picture change transformation X_2c =U+.X_4c.U
C NOTE: It works only as a symmetric transformation !!!
C ----------------------------------------------------------
        CALL QTRANS('AOMO',TREV,D0,
     &     NFBAS(I1,0),NFBAS(I2,0),
     &     NFBAS(I1,1),NFBAS(I2,1),
     &     PRPAO(1+I2BASX(I1,I2)),
     &     NTBAS(0),NTBAS(0),MZ,IPQTOQ(1,IREP),
     &     WORK(KPCTO+I2BASX(I1,I2)),
     &     NTBAS(0),NTBAS(0),MZ,IPQTOQ(1,IREP),
     &     WORK(IPCTM1),
     &     NFBAS(I1,0),NFBAS(I1,1),NZ,IPQTOQ(1,0),
     &     WORK(IPCTM2),
     &     NFBAS(I2,0),NFBAS(I2,1),NZ,IPQTOQ(1,0),
     &     WORK(KFREE),LFREE,IPRINT)

        CALL MEMCHK('DOPCT-after QTRANS',WORK,1) 

C .... control output
       IF (IPRINT.GE.5) THEN

        WRITE(LUPRI,'(/A,1X,A16,/,10X,A,I2,A,I1)')
     &'Symm. blocked FOUR-COMPONENT AO-matrix of property',
     &  PRPNAM(INDXPR),' time reversal symm:',ITIM,
     &  ' number of matrixes: ',MZ
       WRITE(LUPRI,'(5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I1,'/',NFSYM
       CALL PRQMAT(PRPAO(1+I2BASX(I1,I2)),NFBAS(I1,0),NFBAS(I2,0),
     &        NTBAS(0),NTBAS(0),MZ,IPQTOQ(1,IREP),LUPRI)

        WRITE(LUPRI,'(/A,1X,A16,/,10X,A,I2,A,I1)')
     &'Symm. blocked TRANSFORMED TWO-COMPONENT AO-matrix of property',
     &  PRPNAM(INDXPR),' time reversal symm:',ITIM,
     &  ' number of matrixes: ',MZ
       WRITE(LUPRI,'(5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I1,'/',NFSYM
       CALL PRQMAT(WORK(KPCTO+I2BASX(I1,I2)),NFBAS(I1,1),NFBAS(I2,1),
     &        NTBAS(0),NTBAS(0),MZ,IPQTOQ(1,IREP),LUPRI)
        ENDIF      

       IPCTM1 = IPCTM1 + NFBAS(I1,0)*NFBAS(I1,1)*NZ
       IPCTM2 = IPCTM2 + NFBAS(I2,0)*NFBAS(I2,1)*NZ

      ENDDO 

      IF (IPRPRP.GE.5) THEN
       WRITE(LUPRI,'(/A,1X,A16)')
     &'Total (blocked) TRANSFORMED TWO-COMPONENT AO-matrix of property',
     &  PRPNAM(INDXPR)
       CALL PRQMAT(WORK(KPCTO),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),MZ,
     &          IPQTOQ(1,IREP),LUPRI)
      ENDIF

C ... copy WORK(KPCTO) into PRPAO which now contains the 
C      picture chanded operator...
      CALL DCOPY(N2BBASX*MZ,WORK(KPCTO),1,PRPAO,1)

C ... it woould be useful to test the pict change transformed operator...
C .... NFBAS(I1,1) x NFBAS(I2,1) blocks are non-zero,
C     while all others are zero
      PROPNAME=PRPNAM(INDXPR)
      CALL CHECK_PCTOP(IREP,IOPSY,PROPNAME,PRPAO,MZ,IPRINT)

C     Memory deallocation - release the used block

      CALL MEMREL('DOPCTRA',WORK,KWORK,KPCTO,KFREE,LFREE)

 789  CONTINUE
      CALL QEXIT('DOPCTRA')
C     CALL QUIT('In DOPCT...')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck CHECK_PCTOP */
      SUBROUTINE CHECK_PCTOP(IREP,IOPSY,PROPNAME,PRPAO,MZ,IPRINT)
C**************************************************************************
C
C   PURPOSE:  Check the picture change transformed operator 
C             matrix PRPAO
C
C    Written: M Ilias, Prievidza, June 2004
C
C**************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D0 = 0.00D00, D1 = 1.00D00, DEB_THRS=1.0D-7 )
#include "dgroup.h"
#include "dcbbas.h"

      DIMENSION PRPAO(NTBAS(0),NTBAS(0),MZ)
      CHARACTER PROPNAME*16

      CALL QENTER('CHECK_PCTO')

      IF (IPRINT.GE.6) THEN
        WRITE(LUPRI,'(/A,1X,A16)')
     & '*** CHECK_PCTO: Total AO-matrix of 2comp. transf. property',
     &   PROPNAME 
       CALL PRQMAT(PRPAO(1,1,1),NTBAS(0),NTBAS(0),
     &   NTBAS(0),NTBAS(0),MZ, IPQTOQ(1,IREP),LUPRI)

      ENDIF

      DOCCUPB=D0
      DEMPTYB=D0
      IOCCUPB=0
      IEMPTYB=0

      DO I1 = 1, NFSYM
         I2 = MOD(I1+IOPSY,2) + 1
C ... filled matrix blocke LL
         DO IZ = 1, MZ
         DO I = 1, NFBAS(I1,0)
         DO J = 1, NFBAS(I2,0)
          IF (I.LE.NFBAS(I1,1).AND.J.LE.NFBAS(I2,1)) THEN
           IOCCUPB=IOCCUPB+1
           DOCCUPB = DOCCUPB + 
     & DABS(PRPAO(IBAS(I1)+I,IBAS(I2)+J,IZ))
          ELSE
           IEMPTYB=IEMPTYB+1
           DEMPTYB = DEMPTYB + 
     & DABS(PRPAO(IBAS(I1)+I,IBAS(I2)+J,IZ))
          ENDIF   
         ENDDO
         ENDDO
         ENDDO

      ENDDO

      DOCCUPB = DOCCUPB/DFLOAT(IOCCUPB)
      DEMPTYB = DEMPTYB/DFLOAT(IEMPTYB)

      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/,1X,A,A16,/,A,2D12.7)')
     &  'CHECK_PCTO:  For property ', PROPNAME,
     &  'we have (occup, empty): ',DOCCUPB,DEMPTYB   
      ENDIF

      IF (DABS(DEMPTYB).GE.DEB_THRS) THEN
       WRITE(LUPRI,'(/,1X,A,A16,/,A,2D12.7)')
     &  'CHECK_PCTO:  For property ', PROPNAME,
     &  'we have (occup,empty): ',DOCCUPB,DEMPTYB   
       WRITE(LUPRI,'(2X,A)')
     & 'WARNING: Non-zero sum of values in empty block!!!'
      ENDIF

      CALL QEXIT('CHECK_PCTO')
      RETURN
      END
