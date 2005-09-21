C  /* Deck make_bss */
       SUBROUTINE MAKE_BSS(VMAT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C On input : 
C =========
C       VMAT - Lowdin transformation matrix for 4-comp. Fock matrix
C
C On output :
C ===========
C  if wished, VMAT - Lowdin transformation matrix adapted for 2-comp. SCF calculations
C          
C
C    During this routine H2c AO matrix elements as well as
C picture change transformation matrixes are written to the file no LUBSS 'BSSMAT'
C
C Written by Miro Ilias & Hans Jorgen Aa Jensen, Odense, 2002, 
C         reprogrammed by Miro Ilias, Strasburg, 2005.
C
C****************************************************************************
C
C    '  Dobrorec dusa moja Hospodinu, a cele moje vnutro 
C     Jeho svatemu menu !' (Zalm 103,1)
C
C    '  Praise the Lord, my soul ! All my beeing,
C     praise his holy name!'  (Psalm 103,1)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbdhf.h"

        DIMENSION VMAT(*),WORK(LWORK)
        LOGICAL DK2SO1,ONLYSO1,SAVESO1,SF_BEG,FPB,FNDLAB,
     &         GETPCTM,ADD_AMFI,ONLYTRA,ELIM_SF

        CALL QENTER('MAKE_BSS')
#include "memint.h"


C     Print section
C
      CALL HEADER('*** Output from MAKE_BSS ***',-1)

C    ... reference
      WRITE(LUPRI,'(A)')
     & 'H. J. Aa. Jensen and M. Ilias,',
     & '  "Two-component relativistic methods based on the quaternion'//
     &   ' modified Dirac equation. I:',
     & '   from Douglas-Kroll-He{\ss} second order method to the infi'//
     &   'nite order two-component method.",',
     & '  J. Chem. Phys., submitted (2005).'

C extract information about the BSS Hamiltonian 
      CALL EXTR_BSS_INFO(IBSS,IPZ,SPINFR,DK2SO1,SAVESO1,ADD_AMFI,
     &       ONLYSO1,ONLYTRA,SF_BEG,ELIM_SF,NOPRTR,BLOCKD,
     &       DO4C2C,USEDF,YREQ1,IC_AMFI,GETPCTM,IPRHAM)
      FPB = .NOT.NOPRTR

C get transf. matrix (RKB)
      CALL MEMGET('REAL',KTBUF,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMAT,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NTBAS(0),WORK,KFREE,LFREE)

      IF (IPRHAM.GE.5) THEN
       CALL HEADER('MAKE_BSS: Entering Lowdin matrix',-1)
       CALL PRQMAT(VMAT,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),1,
     &            IPQTOQ(1,0),LUPRI)
      ENDIF

C ... get the RKB transformation matrix into WORK(KTMAT)
      CALL MODHAM(WORK(KTMAT),VMAT,WORK(KTBUF),WORK(KEIG),.NOT.FPB,
     &            WORK(KFREE),LFREE)
C      ... NTMO is set up !!!!
C .... this setting are needed (DFDIAG!!!) for the RKB option !!!
C .. NZT is in dcbbas.h, after MODHAM for RKB is equal to  NZ !
      N2TMT   = 0
      N2TMOTQ = 0
      NORBT  = 0
      N2ORBT = 0
      DO IFRP = 1,NFSYM
        NISHMF(IFRP) = NISH(IFRP)
        NOCCMF(IFRP) = NOCC(IFRP)

        I2TMT(IFRP) = N2TMT
C       N2TMT = N2TMT + NFBAS(IFRP,0)*NORB(IFRP)*NZT
        N2TMT = N2TMT + NFBAS(IFRP,0)*NORB(IFRP)*NZ

        I2TMOT(IFRP) = N2TMOTQ
        N2TMOTQ = N2TMOTQ + NTMO(IFRP)*NTMO(IFRP)*NZ

        IORB(IFRP)  = NORBT
C          ....  Count total number of orbitals...
        NORBT    = NORBT  + NORB(IFRP)
C         .....
        I2ORBT(IFRP)= N2ORBT*NZ
        N2ORBT   = N2ORBT + NORB(IFRP)*NORB(IFRP)
      ENDDO

C     N2TMT  = 0
C     NORBT  = 0
C     N2ORBT = 0
C     DO I = 1,NFSYM
C       I2TMT(I) = N2TMT
C       N2TMT    = N2TMT + NFORB(I,0)*NORB(I)*NZ
C       IORB(I)  = NORBT
C          ....  Count total number of orbitals...
C       NORBT    = NORBT  + NORB(I)
C       I2ORBT(I)= N2ORBT*NZ
C       N2ORBT   = N2ORBT + NORB(I)*NORB(I)
C     ENDDO


      LUBSS = 35
      CALL OPNFIL(LUBSS,'BSSMAT','UNKNOWN','GMOTRA')
C ... save Lowdin matrix into BSSMAT
      CALL NEWLAB('VMAT4C  ',LUBSS,LUPRI)
      CALL WRITT(LUBSS,N2BBASX,VMAT)
      WRITE(LUPRI,'(/2X,A)')
     &'MAKE_BSS: 4c Lowdin AO to MO transf. matrix is '//
     &' written to the BSSMAT file.'

C ... write RKB matrix into BSSMAT
      IF (IPRHAM.GE.5) THEN
        CALL HEADER('MAKE_BSS: RKB transformation matrix 4c AO->MO ',-1)
        IMAT = KTMAT
        DO I = 1,NFSYM                                              
        IF(NORB(I).GT.0) THEN                                    
          WRITE(LUPRI,'(A,I1,A,I1)')
     &    '* Fermion ircop no.',I,'/',NFSYM                          
          CALL PRQMAT(WORK(IMAT),NFBAS(I,0),NORB(I),              
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)        
          IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ                         
        ENDIF                                                   
        ENDDO                                                         
      ENDIF

C   .... utilized TMT !!!             
      CALL NEWLAB('TMAT4C  ',LUBSS,LUPRI)
      CALL WRITT(LUBSS,N2TMT,WORK(KTMAT))
      CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
      WRITE(LUPRI,'(/2X,A)')
     &'MAKE_BSS: RKB transf. matrix has been '//
     &' written to the BSSMAT file, EOFLABEL renewed.'

C      ... allocate memmory 
      CALL MEMGET('INTE',KIBEIG,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFACKIN,NTBAS(1),WORK,KFREE,LFREE)

C ===============================================================
C ===       Do prepare the free particle (FP) basis "k"
C ===============================================================
      IF (FPB) THEN
C generate free particle one-electron matrix in RKB ON basis "theta"
C    ... into WORK(KTBUF), reuse VMAT...
C   ... merge into 1 routine ! - get + diagonalize !
C diagonalize the FP matrix in the ON "theta" basis
        CALL MAKE_FPM(WORK(KTBUF),WORK(KTMAT),VMAT,WORK(KEIG),
     &             WORK(KIBEIG),WORK(KFREE),LFREE)

      ELSE
C ===============================================================
C ===   We have the RKB basis, prepare the boson irreps info
C ===============================================================
       IF (SPINFR) THEN
        CALL INI_SF_RKB(WORK(KIBEIG),IPRHAM)
       ENDIF
      ENDIF

C NOTE: If SPINFR=.true. we utilize later the WORK(KIBEIG) array with boson info !!!

C get A factors into WORK(KFACKIN)
C        VMAT - entering ON basis in "xhi"
CCC    IF (.NOT.FPB) CALL DCOPY(N2BBASXQ,WORK(KTMAT),1,VMAT,1)
       CALL GET_A_FAC(WORK(KFACKIN),WORK(KEIG),VMAT,WORK(KTMAT),
     &              WORK(KTBUF),FPB,IPRHAM,WORK(KFREE),LFREE)
C      IF (FPB) CALL DCOPY(N2BBASXQ,VMAT,1,WORK(KTMAT),1)
C      CALL DCOPY(N2BBASXQ,VMAT,1,WORK(KTMAT),1)

C ===============================================================
C    ... Do prepare the MO->AO transformation matrix ...
C ===============================================================
C        ...entering WORK(KTMAT) contains the transformation matrix in AO
       CALL MKSAOMO_Le(WORK(KTBUF),WORK(KTMAT),VMAT,
     &         WORK(KFACKIN),FPB,IPRHAM,WORK(KFREE),LFREE)

C get H1 (Dirac bare nucleus) Hamiltonian in the ON basis into WORK(KBUF)
C    ... that is  the RKB  or FP "k" basis; WORK(KTMAT) entering ... !
      CALL GET_H1_ONB(WORK(KTBUF),WORK(KTMAT),VMAT,WORK(KIBEIG),
     &                   SF_BEG,FPB,SAVESO1,WORK(KFREE),LFREE)


C get H2C Hamiltonian(inf,DKH2,DKH2SO1,DKH1)  in FP k basis 

C     ... by default do prepare the picture change transformation matrix...
C     GETPCTM=.FALSE.
      IF (IPZ.EQ.9) THEN
C       ... infinite irder ... get also picture change transformation matrix
C       GETPCTM=.TRUE.
       
        IF (BLOCKD) THEN
         CALL GET_HINF_BD(WORK(KTBUF),WORK(KTMAT),VMAT,WORK(KIBEIG),
     &                       GETPCTM,SF_BEG,WORK(KFREE),LFREE,IPRHAM)
        ELSE
C        ... we need additional integer array...
         ISIZE=NORB(1)
         IF (NFSYM.EQ.2) THEN
          IF (ISIZE.LT.NORB(2)) ISIZE=NORB(2)
         ENDIF
         CALL MEMGET('INTE',KIPIVOT,ISIZE,WORK,KFREE,LFREE)
C        ... both for FPM "k" & RKB 
         CALL GET_HINF_R(WORK(KTBUF),WORK(KTMAT),VMAT,WORK(KFACKIN),
     &                   WORK(KEIG),WORK(KIPIVOT),FPB,YREQ1,
     &               GETPCTM,SF_BEG,WORK(KFREE),LFREE,IPRHAM)
         CALL MEMREL('AFTER_INF_R',WORK,KWORK,KIPIVOT,KFREE,LFREE)
        ENDIF

      ELSE IF (IPZ.EQ.2) THEN
C       ... second order DKH - needed FREE PARTICLE BASIS !!!
        CALL GET_HDK2_FPM(WORK(KTBUF),WORK(KTMAT),VMAT,WORK(KEIG),
     &                    WORK(KIBEIG),DK2SO1,
     &                    IPRHAM,WORK(KFREE),LFREE)
      ELSE
C       ...and very rough first order, or SO1 integrals only... needed FREE PARTICLE BASIS !!!
        CALL GET_HDK1_FPM(WORK(KTBUF),VMAT, WORK(KIBEIG),ONLYSO1,
     &                            IPRHAM,WORK(KFREE),LFREE)
      ENDIF

C  Do  transform H2c RKB/"k" -> H2c "xhi" LL
      CALL H2CFINAL(WORK(KTBUF),WORK(KTMAT),VMAT,
     &              WORK(KFACKIN),WORK(KIBEIG),SAVESO1,
     &            SPINFR,SF_BEG,GETPCTM,FPB,IPRHAM,WORK(KFREE),LFREE)

C      ... after preparing final H2c Hamiltonian do renew the full relatvity !
      IF (SAVESO1.AND.SPINFR) THEN
       write(lupri,'(/,2X,A)') 'MAKE_BSS: SPINFR set to false !'
       SPINFR = .FALSE.
       SF_BEG = .FALSE.
      ENDIF

C      ... do prepare picture change transformation matrixes...
      IF (GETPCTM) THEN
        CALL GETPCTMAT(WORK(KTBUF),WORK(KTMAT),VMAT,FPB,
     &               SPINFR,SF_BEG,IPRHAM,WORK(KFREE),LFREE)

      ENDIF

C transform H2C Hamitonian from FP "k" basis to the SO-AO  basis,
C prepare the final picture change transformation matrix if wished,
C save integrals (and pctmtx) into the file
      IF (.NOT.ONLYSO1) THEN
C       ... compare eigenvalues of 4c and 2c 1-el.Hamiltonians
        CALL CMP_EIGVAL(VMAT,WORK(KTMAT),WORK(KTBUF),WORK(KEIG),
     &                  WORK(KIBEIG),SPINFR,SF_BEG,
     &                  ONLYSO1,SAVESO1,IPRHAM,WORK(KFREE),LFREE)
C  .... VMAT contains 2c adapted Lowdind matrix (VMAT2C label in BSSMAT)
      ENDIF

C     ...      before calling AMFI stuff do release ALL
C     ...    allocated memmory (except entering VMAT)
      CALL MEMREL('AMFI',WORK,KWORK,KFACKIN,KFREE,LFREE)
      CALL MEMREL('AMFI',WORK,KWORK,KIBEIG,KFREE,LFREE)
      CALL MEMREL('AMFI',WORK,KWORK,KEIG,KFREE,LFREE)
      CALL MEMREL('AMFI',WORK,KWORK,KTMAT,KFREE,LFREE)
      CALL MEMREL('AMFI',WORK,KWORK,KTBUF,KFREE,LFREE)

      IF (ADD_AMFI) THEN
C    ... transformation cartesian -> spherical, calling the AMFI stuff,
C        fitting AMFI integrals in the DIRAC order, saving to the BSSMAT
       IF (SPINFR.AND.ELIM_SF) SPINFR=.FALSE.
C       CALL AMFI_STUFF(WORK(KTBUF),ONLYSO1,ONLYTRA,IC_AMFI,
C      ... VMAT is reused as helping array...
        CALL AMFI_STUFF(VMAT,ONLYSO1,ONLYTRA,IC_AMFI,
     &                 WORK(KFREE),LFREE,IPRHAM,IMFCH)
C       ... do read back VMAT 2c from the file !
       REWIND LUBSS
       IF (FNDLAB('VMAT2C  ',LUBSS)) THEN
         CALL READT(LUBSS,N2BBASX,VMAT)
         IF (IPRHAM.GE.2) WRITE(LUPRI,'(/2X,A)')
     &  'MAKE_BSS: VMAT2C read into VMAT.'
       ELSE
        CALL QUIT(
     &  'MAKE_BSS: VMAT2C is not found in BSSMAT file!')
       END IF
      ENDIF

C ============================================================
C  Finally set up number of positronic shells equal to zero
C ============================================================
      NSYM = 4/NZ
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

C    .. close the BSSMAT file...
      CLOSE(LUBSS,STATUS='KEEP')
C **********************************************************************************
C    ... Number of transf. matrixes is set to 1 - we use only the Lowdin AO->MO
      NZT = 1
C **********************************************************************************
      WRITE(LUPRI,'(2X,A)') 
     &'MAKE_BSS: Two-component BSS mode - positronic shells deleted !'

C   ... release the rest of allocted memmory ...
      CALL MEMREL('MAKE_BSS',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('MAKE_BSS')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck ini_sf_rkb */
      SUBROUTINE INI_SF_RKB(IBEIG,IPRINT)
C****************************************************************************
C   Initialize the IBEIG(*) array with boson info for 2-component
C   calculations.
C****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION IBEIG(*)

      DO I = 1, NFSYM
        CALL INISUB
        CALL INIBOS(IBEIG(1+IORB(I)),I,.TRUE.,IPRINT)
      ENDDO

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck TEST_UNIT_MTX */
      SUBROUTINE TEST_UNIT_MTX(S,N,LR,LC,NZ)
C****************************************************************************
C
C   PURPOSE: Check the enterning matrix - is it unit mtx ?
C
C  Written by Miro ILIAS, Odense,2003
C             
C****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00, THRNULL = 1.00D-09 )

      DIMENSION S(LR,LC,NZ)

      SUM_DE  = SUM_DIAG_EL(S,N,LR,LC,NZ)
      SUM_DEN = SUM_DE/DFLOAT(N)
      SUM_ODE = SUM_OFFDIAG_EL(S,N,LR,LC,NZ)

      WRITE(LUPRI,'(/2X,A,1P,D12.5)') 
     &'Sum of diagonal elements /N  :',SUM_DEN
      WRITE(LUPRI,'(2X,A,1P,D12.5)') 
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
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
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
      SUBROUTINE EXTR_BSS_INFO(IBSS,IPZ,SPINFR,DK2SO1,SAVESO1,ADD_AMFI,
     &                    ONLYSO1,ONLYTRA,SF0,ELIM_SF,
     &                    NOFPB,BLOCKD,DO4C2C,USEDF,
     &                    YREQ1,IC_AMFI,GETPCTM,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C PURPOSE:
C ========
C     Set control variables based on IBSS coding
C     IBSS = axyz
C
C        a - 1 to 6, defines AMFI flags..
C
C        x - .gt. 0: define spinfree before (as 4-comp)
C            .eq. 0: define spinfree at end
C            ... x is active with the SPINFREE option & with the y=1 (means only SO1)
C           
C        y - order of spin-orbit term:
C            0 - spinfree, 1 - first order for DK2/BSS, 2+ - full.
C
C        z - order of Douglas-Kroll transformation
C            (0, 1, 2, or 9), where 9 is full.
C
C     Pick up the x,y,z variables from IBSS = xyz
C
C On input: IBSS, SPINFR, IPRINT, ....
C
C On output: IPZ - order of Hamiltonian
C 
C            DK2SO1,SF0
C
C            ONLYSO1 - flagfo BSS-DKH-SO1 integrals to be compared with AMFI DKH-SO1 integrals
C
C            SAVESO1  - flag to decide to save SO1 elements (from H1) into the BSSMAT file
C  to be added at the later stage to the final spin-free 2c Hamilt.(BSSsfB,BSSsfE,DKH2sfB,DKH2sfE)
C
C            ONLYTRA - flag for C->S->C transformation
C
C            ELIM_SF
C
C            IC_AMFI - flag for AMFI subroutine
C
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      LOGICAL SPINFR,DK2SO1,SF0,ONLYSO1,ONLYTRA,ADD_AMFI,ELIM_SF,
     &      NOFPB,BLOCKD,YREQ1,DO4C2C,USEDF,SAVESO1,GETPCTM

      CALL QENTER('EXTR_BSS_INFO')

C     ... do extract a,x,y,z digits of the entering IBSS=axyz value
      CALL EXTR_DIG(IBSS,IPA,IPX,IPY,IPZ,IPRINT)

      IF (DO4C2C) THEN
       WRITE(LUPRI,'(/2X,A/)')
     & 'After the Dirac-Fock procedure do continue'//
     & ' in the TWO-COMPONENT MODE !'
       IF (USEDF) THEN
         WRITE(LUPRI,'(/2X,A/)')
     & '...transforming the converged four-component'//
     & ' DIRAC-FOCK operator !'
       ENDIF
      ELSE
        WRITE(LUPRI,'(/2X,A/)')
     & 'Transforming the one-electron DIRAC bare nucleus operator.'
      ENDIF


      WRITE (LUPRI,'(A)')
     & ' * Order of the Douglas-Kroll-Hess transformation :'

      IF (IPZ.EQ.9) THEN
       write(lupri,'(5X,A/)')
     &   '...infinite order ! ( IOTC = BSSinf = DKHinf )'
      ELSE IF (IPZ.EQ.2) THEN
       write(lupri,'(5X,A/)') '...second order ! ( DKH2 )'
      ELSE IF (IPZ.EQ.1) THEN
       write(lupri,'(5X,A/)') '...first order ! ( DKH1 )'
      ELSE IF (IPZ.EQ.0) THEN
       write(lupri,'(5X,A/)')
     & '...only testing of the SO1 approach, then quit ...'
      ELSE
       write(lupri,*) 'IPZ=',IPZ
       write(lupri,*) 'Can be only 9,2,1 or 0!'
       CALL QUIT('MAKE_BSS: EXTR_BSS_INFO: Wrong value of'//
     &           ' IZ (should be 0,1,2,9) !')
      ENDIF

C      ... treat the NOFPB & BLOCKD
      IF (IPZ.EQ.9) THEN
       IF (NOFPB) THEN
        WRITE(LUPRI,*) 'Only the RKB basis ! No preliminary free-'//
     &  'particle basis transformation according to Barysz !'
       ELSE
        WRITE(LUPRI,*)
     & '...the preliminary free-particle (or "k") basis'//
     & ' transformation of the Dirac bare nuclei is applied !'
       ENDIF
       BLOCKD = .FALSE.
       IF (BLOCKD) THEN
        WRITE(LUPRI,*) 'The block diagonalization applied (QJACOBI) !'
       ELSE
        WRITE(LUPRI,*) 'R-equation solving applied !'
        IF (YREQ1) THEN
          WRITE(LUPRI,'(3X,A)') '...solving YR-equation 1 (HJAaJ) !'
        ELSE
          WRITE(LUPRI,'(3X,A)') '...solving YR-equation 2 (Trond) !'
        ENDIF
       ENDIF
      ELSE
C     ... for lower orders do always use the FP "k"- basis !
         NOFPB = .FALSE.
         BLOCKD = .FALSE.
      ENDIF

      ELIM_SF=.FALSE.
      IF (SPINFR.AND.IPY.GE.1) THEN 
       WRITE(LUPRI,'(A)') ' * "SPINFR" option, no spin-orbit terms!' 
       IPY=0
      ENDIF

      DK2SO1 = .FALSE.
         SF0 = .FALSE.
      IF (IPY .EQ. 0) THEN
         IF (.NOT.SPINFR.AND.IPA.EQ.0) CALL QUIT(
     & 'You forgot to specify .SPINFR in the Hamiltonian input!')
         DK2SO1 = .FALSE.
         IF (SPINFR.AND.IPA.GT.0) CALL QUIT(
     & 'Wrong .SPINFR in the Hamiltonian input with the AMFI on!')
         IF (.NOT.SPINFR.AND.IPA.GT.0) THEN
          SPINFR=.TRUE.
          ELIM_SF=.TRUE.
         ENDIF
      ELSE IF (IPY.EQ.1.AND.IPZ.EQ.2.AND.IPX.EQ.1
     &           .AND.(.NOT.SPINFR)) THEN
C        ... DKH2sfB + SO1
        DK2SO1 = .TRUE.
        WRITE (LUPRI,'(A)')
     &  ' * spin-orbit terms only to the first order -'//
     &  ' (within the DKH2 Hamilt. generation)'//
     &  ' - DKH2sfBSO1 Hamiltonian'
      ELSE IF (.NOT.SPINFR) THEN
         DK2SO1 = .FALSE.
         WRITE (LUPRI,'(A)')
     &   ' * all spin-orbit terms included (full relativity)'
      END IF

      IF (IPX .GE. 1 .AND. SPINFR) THEN  
         SF0 = .TRUE.
         WRITE (LUPRI,'(A)')
     & ' * Elimination of spin-orbit terms from the Dirac bare '//
     & 'nuclei in RKB/"k" basis (spinfree "BEFORE / AT THE BEGINNING")' 
      ELSE IF(IPX .EQ. 0 .AND. SPINFR) THEN
         SF0 = .FALSE.
         IF (SPINFR) WRITE (LUPRI,'(A)')
     & ' * Elimination of spin-orbit terms from the 2comp.Ham. '//
     & ' in RKB/"k" basis (spinfree "AFTER / AT THE END")'
      ENDIF

C     ... by default do not include SO1 addition to spin-free Hamiltonians DKH2,BSS
      SAVESO1 = .FALSE.
      IF (IPY.EQ.1.AND..NOT.SPINFR.AND..NOT.DK2SO1
     &            .AND.IPZ.GT.1.AND..NOT.NOFPB) THEN
        SAVESO1=.TRUE.
        SPINFR=.TRUE.
        SF0=.FALSE.
        IF (IPX.GE.1) SF0=.TRUE.
        WRITE(LUPRI,'(/2X,A)')
     & ' Addition of SO1 terms (extracted from H1_ONB)'//
     & ' to the desired spin-free Hamiltonian !'
        IF (IPZ.EQ.2.AND.SF0) THEN
          write(lupri,'(5X,A/)') 'DKH2sfB + SO1 Hamiltonian'
        ELSE IF (IPZ.EQ.2.AND..NOT.SF0) THEN
          write(lupri,'(5X,A/)') 'DKH2sfE + SO1 Hamiltonian'
        ELSE IF (IPZ.GT.2.AND.SF0) THEN
          write(lupri,'(5X,A/)') 'IOTCsfB + SO1 Hamiltonian'
        ELSE IF (IPZ.GT.2.AND..NOT.SF0) THEN
          write(lupri,'(5X,A/)') 'IOTCsfE + SO1 Hamiltonian'
        ENDIF
      ENDIF

      IF (ELIM_SF) WRITE (LUPRI,'(A)')                                
     & ' * SPIN-FREE flag is valid only for the generation of the BSS'//
     & ' Hamiltonian, later removed '//
     & ' due to the AMFI spin-orbit integrals.'
     
C      .... flag IBSS=010 - activate the ONLYSO1 flag
      ONLYSO1=.FALSE.
      ADD_AMFI=.FALSE.
      IF (IPX.EQ.0.AND.IPY.EQ.1.AND.IPZ.EQ.0) THEN
         WRITE (LUPRI,'(/A)')
     &  ' * only SO1 terms (for comparison with AMFI, then quit)'
        ONLYSO1=.TRUE.
C       ... AMFI gives only SO1 integrals
        IC_AMFI = 1
        ADD_AMFI=.TRUE.
      ENDIF

      ONLYTRA=.FALSE.
      IF (IPA.EQ.1) THEN
        ONLYTRA=.TRUE.
        ADD_AMFI=.TRUE.
C     ... AMFI gives SO1+MFSO2 integrals,
C        but these are added to the Hamiltionian with ZERO weight
        IC_AMFI = 3
        WRITE(LUPRI,'(/A)')
     &  ' * C->S->C transformation + AMFI run(test), but no adding'//
     &  ' of AMFI integrals.'  
      ELSE IF (IPA.EQ.2) THEN
        ADD_AMFI=.TRUE.
        IC_AMFI = 2
        WRITE(LUPRI,'(/A)')
     &  ' * AMFI: MFSSO2 terms (spin-same orbit only) added.'  
      ELSE IF (IPA.EQ.3) THEN
        ADD_AMFI=.TRUE.
        IC_AMFI = 3
        WRITE(LUPRI,'(/A)')
     &  ' * AMFI: MFSO2 terms (both spin-same&spin-other orbit) added.'
      ELSE IF (IPA.EQ.4) THEN
        ADD_AMFI=.TRUE.
        IC_AMFI = 4
        WRITE(LUPRI,'(/A)')
     &  ' * AMFI: SO1+MFSSO2 terms (only spin-same orbit) added.'
      ELSE IF (IPA.EQ.5) THEN
        ADD_AMFI=.TRUE.
        IC_AMFI = 5
        WRITE(LUPRI,'(/A)')
     &  ' * AMFI: SO1+MFSO2 terms (spin-same&spin-other orbit) added.'
      ELSE IF (IPA.EQ.6) THEN
        ADD_AMFI=.TRUE.
        IC_AMFI = 1
        WRITE(LUPRI,'(/A)')
     &  ' * AMFI: SO1 terms added.'  
      ENDIF

      IF (.NOT.ADD_AMFI) THEN
        WRITE (LUPRI,'(A/)')
     &  ' * AMFI contribution is NOT included.'
      ENDIF

C ... treat the GETPCTM flag
      GETPCTM = .FALSE.
      IF (.NOT.SPINFR.AND..NOT.SAVESO1.AND.IPZ.EQ.9) THEN
       GETPCTM = .TRUE. 
       WRITE(LUPRI,'(2X,A)')
     & '...do prepare the picture-change transformation matrix !'
      ELSE
       WRITE(LUPRI,'(2X,A)')
     & '...do NOT make the picture-change transformation matrix !'
      ENDIF

C     ... consistency checks - do not repeat 1 el SO interactions !!!
C =========================================================================
      IF (IPY.GT.0.AND.(IPA.EQ.6.OR.IPA.EQ.5.OR.IPA.EQ.4)) THEN
       write(lupri,*)
       WRITE(LUPRI,*) 'EXTR_BSS_INFO: IPY=',IPY,' IPA=',IPA
       WRITE(LUPRI,*) 'EXTR_BSS_INFO: Double inclusion of 1el'//
     & ' SO terms - AMFI & BSS Ham. '
       CALL  QUIT('EXTR_BSS_INFO: Double inclusion of 1el SO terms! ')
      ENDIF
      WRITE(LUPRI,*)

      CALL QEXIT('EXTR_BSS_INFO')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck EXTR_DIG */
      SUBROUTINE EXTR_DIG(IBSS,IPA,IPX,IPY,IPZ,IPRINT)
C****************************************************************************
C
C  PURPOSE: For entering integer number of type xyz,
C            extract a,x,y,z into IPA,IPX,IPY,IPZ
C             
C  Written: M Ilias, Odense, aug. 2002             
C
C****************************************************************************
#include "implicit.h"
#include "priunit.h"

      IF (IBSS.LT.1.OR.IBSS.GT.9999)
     & CALL QUIT('EXTR_DIG: The entering number not in range 1-9999! ')
  
C     ... save the IBSS
      IBSS_SAVE = IBSS
      IBSS0     = IBSS

      IF (IBSS.GE.1000) THEN
       IPA  = IBSS/1000
       IBSS0 = MOD(IBSS,1000)
       IBSS = MOD(IBSS0,100)
       IPX  = IBSS0/100
       IPY = IBSS/10
       IPZ = MOD(IBSS,10)
      ELSE IF (IBSS.GE.100) THEN 
       IBSS = MOD(IBSS,100)
       IPA = 0
       IPX  = IBSS0/100
       IPY = IBSS/10
       IPZ = MOD(IBSS,10)
      ELSE IF (IBSS.GE.10) THEN
       IPA = 0
       IPX = 0
       IPY = IBSS/10
       IPZ = MOD(IBSS,10)
      ELSE IF (IBSS.GE.0) THEN
       IPA = 0
       IPX = 0
       IPY = 0
       IPZ = IBSS
      ENDIF

      IF (IPRINT.GE.3) THEN
       WRITE (LUPRI,'(/A,I5,A,4I3)')
     &'EXTR_DIG: IBSS=axyz=',IBSS_SAVE,
     &' ...extracted digits  a,x,y,z:',IPA,IPX,IPY,IPZ
      ENDIF

C     ... restore the IBSS
      IBSS = IBSS_SAVE
C     ... check the extraction procedure
      IBSS1 = (1000*IPA)+(100*IPX)+(10*IPY)+IPZ
      IF (IBSS1.NE.IBSS) THEN 
        WRITE(LUPRI,*) 'EXTR_DIG: input IBSS=axyz',IBSS
        WRITE(LUPRI,*) 'EXTR_DIG: extracted digits:',IPA,IPX,IPY,IPZ
        CALL QUIT('EXTR_DIG: Error in extraction of a,x,y,z from axyz!')
      ENDIF

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck spf_hk */
      SUBROUTINE SPF_HK(HK,SPFR,IBS1,IBS2,
     &                 NROW,NCOL,LROW,LCOL,NZ,IPRINT)
C****************************************************************************
C
C  PURPOSE: Eliminate spin-orbit off-diagonal elements of entering matrix
C           in MO basis (SPFR=.true.) or eliminate spin-free elements(SPFR=.false.) 
C
C  On input: HK - quaternion hermitian matrix
C            IBS1 - boson symetries of first index
C            IBS2 - boson symetries of second index
C            SPFR - if true, elim.SO terms, otherwise elim.SF terms
C
C  On output: HK - the same matrix without spin-orbit/spin-free components
C             
C  Written: M Ilias, Odense, aug. 2002             
C           last update MI,2005,Strasbourg             
C
C****************************************************************************
#include "implicit.h"
      PARAMETER ( D0 = 0.00D00 )
#include "priunit.h"

      DIMENSION HK(LROW,LCOL,NZ),IBS1(*),IBS2(*)
      LOGICAL SPFR

      CALL QENTER('SPF_HK')

      IF (NROW.NE.NCOL) THEN
        WRITE(LUPRI,'(2X,A)') 'SPF_HK: Warning, NROW<>NCOL!'
        CALL QUIT('SPF_HK: NROW<>NCOL! How it is possible?')
      ENDIF


      IF (IPRINT. GE. 5) THEN
       write(lupri,*)
       IF (SPFR) THEN
        write(lupri,*)
     &  '..eliminating SPIN-ORBIT TERMS..'
       ELSE
        write(lupri,*)
     &  '..eliminating SPIN-FREE TERMS..'
       ENDIF
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

      IF (SPFR) THEN
C     ... eliminate spin-orbit components...
       IZ = 1
       DO I=1,NROW
       DO J=1,NCOL
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

      ELSE
C     ... eliminate spin-free components
       IZ = 1
       DO I=1,NROW
       DO J=1,NCOL
         IF (IBS1(I).EQ.IBS2(J)) HK(I,J,IZ)=D0
       ENDDO
       ENDDO

      ENDIF

      CALL QEXIT('SPF_HK')

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
      LOGICAL FNDLAB,DOTR
      CHARACTER  PROPNAME*16,TREV*1

      CALL QENTER('DOPCTRA')
#include "memint.h"

      DOTR=(INDXPR.EQ.-1)

C ... check if one needs to do the pic.change transf.
      IF (.NOT.DOPCT.AND..NOT.DOTR) GOTO 789

      IPRINT = IPRPRP

      IF (.NOT.DOTR) THEN
      IPRINT = IPRPRP
      IREP = IPRPSYM(INDXPR) - 1
C ... determine the time reversal symmetry
C     of the property operator
      ITIM = IPRPTIM(INDXPR)
      ELSE
       IREP=0
      ENDIF


C ,,, read (load) the picture change transformation matrix:
C       4c AO-> 2c AO(LL)
      CALL MEMCHK('in DOPCT',WORK,1)

      CALL OPNFIL(LUBSS,'BSSMAT','OLD','DOPCTRA')
C     ... set up at the beginning of the file
      REWIND LUBSS
      IF (FNDLAB('U_PICTRM',LUBSS)) THEN
        ISIZE = NFBAS(1,0)*NFBAS(1,1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NFBAS(2,1)*NZ
C ... array for the picture change transformation matrix (PCTM)...
       CALL MEMGET('REAL',KPCTM,ISIZE,WORK,KFREE,LFREE)
        CALL READT(LUBSS,ISIZE,WORK(KPCTM))
C       CALL READT(LUBSS,N2BBASXQ,WORK(KPCTM))
        IF (IPRPRP.GE.2) WRITE(LUPRI,'(/2X,A)')
     &   'DOPCTRA: Picture change transformation matrix was read.'
      ELSE
        IF (DOTR) THEN
          write(lupri,*) 'DOPCTRA: U_PICTRM not found in BSSMAT file!' 
          write(lupri,*)
     &      'deleting the SS block of entering matrix only...' 
          CLOSE (LUBSS,STATUS='KEEP')
          DO I = 1, NFSYM
             IS1 = 1 + (NTBAS(0)*NFBAS(I,1))+NFBAS(I,1) +
     &                  I2BASX(I,I)
            DO J = 1, NFBAS(I,2)
              ISS = IS1 + (J-1)*NTBAS(0)
              CALL DZERO(PRPAO(ISS),NFBAS(I,2))
            ENDDO
          ENDDO
          GOTO 789
        ELSE
          CALL QUIT('DOPCTRA:Error, U_PICTRM not found in BSSMAT file!')
        ENDIF
      END IF
      CLOSE (LUBSS,STATUS='KEEP')

      IF (IPRINT.GE.5) THEN
       CALL HEADER('Output from DOPCTRA',-1)
       IF (.NOT.DOTR) THEN
       WRITE(LUPRI,'(A,A16,/,10X,A,I2,A,I1)')
     & 'Total entering FOUR-COMPONENT AO-matrix of property ',
     & PRPNAM(INDXPR),' time rev.symm:',ITIM,
     & ' number of matrixes:',MZ
       ENDIF
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
        IF (.NOT.DOTR) THEN
        WRITE(LUPRI,'(/A,1X,A16,/,10X,A,I2,A,I1)')
     &'Symm. blocked FOUR-COMPONENT AO-matrix of property',
     &  PRPNAM(INDXPR),' time reversal symm:',ITIM,
     &  ' number of matrixes: ',MZ
        ELSE
        WRITE(LUPRI,'(/A,1X,A,I1)')
     &'Symm. blocked FOUR-COMPONENT AO-matrix ',
     &  ' number of matrixes: ',MZ
        ENDIF
       WRITE(LUPRI,'(5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I1,'/',NFSYM
       CALL PRQMAT(PRPAO(1+I2BASX(I1,I2)),NFBAS(I1,0),NFBAS(I2,0),
     &        NTBAS(0),NTBAS(0),MZ,IPQTOQ(1,IREP),LUPRI)

        IF (.NOT.DOTR) THEN
        WRITE(LUPRI,'(/A,1X,A16,/,10X,A,I2,A,I1)')
     &'Symm. blocked TRANSFORMED TWO-COMPONENT AO-matrix of property',
     &  PRPNAM(INDXPR),' time reversal symm:',ITIM,
     &  ' number of matrixes: ',MZ
        ELSE
        WRITE(LUPRI,'(/A)')
     &'Symm. blocked TRANSFORMED FOUR-COMPONENT AO-matrix '
        ENDIF
       WRITE(LUPRI,'(5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I1,'/',NFSYM
       CALL PRQMAT(WORK(KPCTO+I2BASX(I1,I2)),NFBAS(I1,1),NFBAS(I2,1),
     &        NTBAS(0),NTBAS(0),MZ,IPQTOQ(1,IREP),LUPRI)
        ENDIF      

       IPCTM1 = IPCTM1 + NFBAS(I1,0)*NFBAS(I1,1)*NZ
       IPCTM2 = IPCTM2 + NFBAS(I2,0)*NFBAS(I2,1)*NZ

      ENDDO 

      IF (IPRPRP.GE.5) THEN
       IF (.NOT.DOTR) THEN
       WRITE(LUPRI,'(/A,1X,A16)')
     &'Total (blocked) TRANSFORMED TWO-COMPONENT AO-matrix of property',
     &  PRPNAM(INDXPR)
        ELSE
       WRITE(LUPRI,'(/A)')
     &'Total (blocked) TRANSFORMED TWO-COMPONENT AO-matrix: '
        ENDIF
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
      CALL CHECK_PCTOP(IREP,IOPSY,PROPNAME,PRPAO,MZ,IPRPRP)

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
     & '*** CHECK_PCTOP: Total AO-matrix of 2comp. transf. property',
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


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck MAKE_FPM */      
       SUBROUTINE MAKE_FPM(TBUF,TMAT,VMAT,EIG,IBEIG,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:
C    ========
C      get FP (kinetic) matrix in the ON "theta" basis into TBUF
C
C   On input:
C  -----------
C
C   TMAT - transformation matrix
C
C  On output:   
C  -----------   
C
C   TMAT FP (kinetic) matrix in the SA-AO "xhi" basis 
C
C
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C                                                         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0,
     &           THRNULL = 1.0D-10, THRNULL1 = 1.0D-4)
C  ...  THRNULL is round off error with 15 digits relative to 2c^2 = 3.7D4
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbdhf.h"
      DIMENSION TBUF(*),VMAT(*),TMAT(*),EIG(*),IBEIG(*),WORK(LWORK)
      LOGICAL DBGDIAG, DOJACO_SAVE, FNDLAB

      CALL QENTER('MAKE_FPM')
#include "memint.h"


C ... get FP matrix in SA-AO into TBUF
      CALL FREEMT(TBUF,IPRINT,WORK(KFREE),LFREE) 

C         Transform the free particle matrix
C         TBUF in AO basis
C         in modified Dirac form (orhonormal "theta" basis); 
C
C      VMAT =  TMAT+ * TBUF * TMAT 

      IVMAT = 1
      ITMAT = 1
      DO I = 1,NFSYM
      IF(NORB(I).GT.0) THEN

        IF (IPRINT.GE.7) THEN
           CALL HEADER(
     &  'MAKE_FPM: Free particle matrix in SA-AO basis "xhi"'//
     &  ' (after call FREEMT)',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TBUF(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                 NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)

           CALL HEADER(
     &  'MAKE_FPM: ENTERING RKB transf. matrix "xhi"->"theta", TMAT',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NORB(I),
     &                 NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

        ENDIF

        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                             NORB(I),NORB(I),
     &      TBUF(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &           VMAT(IVMAT),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &           TMAT(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           TMAT(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRINT)
C
C       Zero the small elements in real part of WORK(KBUF)
C       which are non-zero because of round-off errors:
C       check that quaternion blocks are numerically zero
C
        DO J = 0,NORB(I)*NORB(I)-1
          IF (DABS(VMAT(IVMAT+J)) .LT. THRNULL) VMAT(IVMAT+J) = D0
        END DO

        IF (NZ .GT. 1) THEN
           IVMAT2 = IVMAT + NORB(I)*NORB(I)
           NVMAT2 = (NZ-1)*NORB(I)*NORB(I)
           TIMAGNORM =  DNRM2(NVMAT2,VMAT(IVMAT2),1)
        ELSE
           TIMAGNORM =  D0
        END IF

C  Write out TIMAGNORM

        IF (TIMAGNORM.GT.THRNULL1) THEN

          WRITE(LUPRI,'(2X,A,A,I1,A,I1,/,2X,A,D12.6)')
     & 'Free particle matrix in the ON basis "theta" ',
     & 'in fermion corep ',I,'/',NFSYM,
     & '*** Norm of the i,j,k imag. parts *** :',TIMAGNORM
          WRITE(LUPRI,'(/2X,A,D13.6)')
     & 'MAKE_FPM: WARNING, TIMAGNORM > THRNULL1 (=',THRNULL1

        ENDIF

        IF (IPRINT.GE.6) THEN
          CALL HEADER(
     &'MAKE_FPM: Free particle matrix in the ON'//
     & ' RKB basis "theta",VMAT',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &'*** Fermion corep ',I,'/',NFSYM
          WRITE(LUPRI,'(/A,1P,D12.6)')
     & ' Norm of the i,j,k imag. parts of this matrix:',TIMAGNORM
          CALL PRQMAT(VMAT(IVMAT),NORB(I),NORB(I),
     &                NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IVMAT = IVMAT + NORB(I)*NORB(I)*NZ
        ITMAT = ITMAT + NFBAS(I,0)*NORB(I)*NZ

      ENDIF
      ENDDO


C  ... Td_kin is set to be pure real matrix !
      NZEKIN = 1 
      NZSAVE = NZ
      NZ = NZEKIN 

C     note: DFDIAG zeroes all of TBUF, even if NZ is reset to 1
C     (it uses N2BBASXQ); so we can continue with the quaternion TBUF afterwards.

      DOJACO_SAVE = DOJACO
      DOJACO = .TRUE.
C     ... use JACOBI diagonalization because it does not mix degenerate eigenvectors
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

      IF (IPRHAM.GE.3) THEN
       WRITE(LUPRI,'(/,3X,A,A/)') 
     & 'MAKE_FPM: Diagonalizing (Jacobi) the real part',
     & ' of the full FP matrix in the ON "theta" basis !'
      ENDIF


C ... after diagonalization we get the IBEIG array which is used later !
      IF (IPRHAM.GE.5) THEN 
        DBGDIAG=.TRUE.
      ELSE
        DBGDIAG=.FALSE.
      ENDIF

C  DO copy FPM from VMAT to TBUF - do restore block !

      CALL DZERO(TBUF,N2BBASXQ)
     
C ... check I2ORBT and I2TMOT...
C    ... subroutine in dirscf.F for diagonalizing the Dirac-Fock matrix
C  FPK in TBUF
C   TMAT - eigenvectors in "k", EIG - eigenvalues, IBEIG - boson irreps
C     CALL DZERO(TMAT,N2BBASXQ)

C      ... do the diagonalization at once ! Block structure is needed !!!
C      ... needed NTMO, I2TMOT, IORB offsets ...
C     CALL DFDIAG (VMAT,EIG,IBEIG,TMAT,DBGDIAG,
      CALL DFDIAG (VMAT,EIG,IBEIG,TBUF,DBGDIAG,
     &             WORK,KFREE,LFREE)


C   ... read 4ctmat to the TBUF
C... check it !
C     REWIND(LUBSS)
C     IF (FNDLAB('TMAT4C  ',LUBSS)) THEN
C        CALL READT(LUBSS,N2TMT,TBUF)
C        IF (IPRHAM.GE.2) WRITE(LUPRI,'(/2X,A)')
C    &'DIAG_FPK: TMAT4C read into TBUF.'
C     ELSE
C       CALL QUIT('DIAG_FPK:Error, TMAT4C not found in BSSMAT file!')
C     END IF

C     CALL DZERO(VMAT,N2BBASXQ)
    
      IF (N2TMT.GT.N2BBASXQ) THEN
        CALL QUIT('DIAG-FPM: N2TMT.GT.N2BBASX!')
      ENDIF

      NZ = NZSAVE
      DOJACO = DOJACO_SAVE
      SUB_BL = .FALSE.

      CVAL2 = CVAL*CVAL

      IMAT = 1
      DO  I = 1,NFSYM
      IF(NORB(I).GT.0) THEN

C   ... do always check the RKB free-particle eigenvalues !

        IF (IPRHAM.GE.2) THEN
         CALL HEADER(
     &'MAKE_FPM: Free particle (in ON basis "theta") eigenvalues :',-1)
         WRITE(LUPRI,'(3X,A,I2,A,I2)') '*** Fermion corep ',I,'/',NFSYM
         WRITE(LUPRI,'(A,3X,A,2X,A,2X,A,3X,A)')
     &   '   Shell no ','FP eigval',
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
        WRITE(LUPRI,'(2X,A,I1,A,I1,A,D12.6)')
     & 'Average difference for e-p pairing of '//
     & 'free particle eigenvalues (fsym ',I,
     & '/',NFSYM,') : ',TOTDIFF

       IF (IPRHAM.GE.7) THEN
        CALL HEADER(
     & 'MAKE_FPM: |k> matrix (FPM eigenvectors)'//
     & ' in the ON basis "theta", TBUF ',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(TBUF(IMAT),NORB(I),NORB(I),
     &                NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C-------------------------------------------------------------------------------
C  Backtransform ELECTRONIC and POSITRONIC solutions from the ON basis "theta"
C  to the AO SO-basis "xhi".
C (we later do the no-pair projection to the electronic solutions after
C  we have done the DKH folding of the nuclear attraction into the
C  electronic part).
C
C  VMAT = TBUF*TMAT 
C-------------------------------------------------------------------------------
        CALL BCKTRA(VMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),
     &              TBUF(IMAT),NORB(I),NORB(I),
     &              NORB(I),NZ,NORB(I),1,NFBAS(I,0),
     &              TMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),NZ,
     &              IPRHAM)

        IF(IPRHAM.GE.7) THEN
          CALL HEADER(
     &    'MAKE_FPM: FP |k> matrix, VMAT (output-TMAT)'//
     &    ' backtransformed to SA-AO basis "xhi" ',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     &    '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(VMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

         IMAT = IMAT + NORB(I)*NORB(I)*NZ

      ENDIF    
      ENDDO 

C     ... do export "k" in "xchi" in TMAT
      CALL DCOPY(N2TMT,VMAT,1,TMAT,1)

C  ... do save VMAT - "k" in "xhi" into the file BSSMAT
C ... save Lowdin matrix into BSSMAT
      REWIND LUBSS
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
       BACKSPACE(LUBSS)
C        ... N2TMT is NFBAS(1,0)*NORB(1) + NFBAS(2,0)*NORB(2)
       CALL NEWLAB('FPMAT4C ',LUBSS,LUPRI)
       CALL WRITT(LUBSS,N2TMT,TMAT)

C ===============================================================
C        ... write help transformation matrix if needed ....
       ISIZE = NORB(1)*NORB(1)*NZ
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
       CALL NEWLAB('FPMAT_TH',LUBSS,LUPRI)
       CALL WRITT(LUBSS,ISIZE,TBUF)
       WRITE(LUPRI,'(/2X,A)')
     & 'MAKE_FPM: 4c "k"-> 4c "theta" transf. matrix FPMAT_TH, TBUF'//
     & ' written to the BSSMAT file.'
C ===============================================================

       CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
       WRITE(LUPRI,'(/2X,A)')
     & 'MAKE_FPM: 4c "k"-xhi transf. matrix, TMAT, is '//
     & ' written to the BSSMAT file, EOFLABEL renewed.'
      ELSE
       CALL QUIT(
     & 'MAKE_FPM:Error, EOFLABEL not found in BSSMAT file!')
      END IF

C     CALL QUIT('DIAG_FPK')
      CALL QEXIT('MAKE_FPM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_a_fac */      
      SUBROUTINE GET_A_FAC(FACKIN,EIG,VMAT,TMAT,TBUF,
     &                      FPB,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:  Construct (and check for PPB) the diagonal elements of the A factor 
C    ========
C                                                         
C    On input: TMAT - ON RKB/FP "k" in "xhi" SO-AO basis
C
C              EIG  - free particle eigenvalues  if FPB=.true.                                                   
C              FPB
C
C   Output:  FACKIN - normalization factors                                                         
C                                                         
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C                                                         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0, D2 = 2.0D0, THRNULL = 1.0D-10 ) 
#include "dcbgen.h"  
#include "dcbbas.h"  
#include "dcborb.h" 
#include "dgroup.h"   
      DIMENSION EIG(*),VMAT(*),TMAT(*),TBUF(*),FACKIN(*),WORK(LWORK)    
      LOGICAL FPB

      CALL QENTER('GET_A_FAC')
#include "memint.h"         


C       ... TBUF subst  WORK(KBUF)
      CALL GTOVLT(TBUF,D0,WORK(KFREE),LFREE)

      CVAL2 = CVAL*CVAL

      IMAT = 1
      IFC = 0
      IBUF = 1
      IKUN1 = 1
      DO I = 1, NFSYM
      IF (NORB(I).GT.0) THEN

         IF(IPRINT.GE.7) THEN
          CALL HEADER(
     &  'GET_A_FAC: Entering ON RKB/|k> matrix, TMAT,'//
     &  ' in SO basis "xhi" ',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')'* Fermion irp no.',I,'/',NFSYM
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          CALL HEADER(
     &    'GET_A_FAC: Obtained (GTOVLT) LL Sovlp matrix'//
     &    ' in SO basis "xhi"',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')'* Fermion irp no.',I,'/',NFSYM
          CALL PRQMAT(TBUF(IBUF),NFBAS(I,0),NFBAS(I,0),
     &                NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),LUPRI)
         ENDIF

** DO C+(e) S(LL) C(e) into TMAT(KKUN) **

        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &       NESH(I),NESH(I),
     &       TBUF(IBUF),NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),
     &       VMAT(IKUN1),NESH(I),NESH(I),1,IPQTOQ(1,0),
     &       TMAT(IMAT+NFBAS(I,0)*NPSH(I)),
     &       NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &       TMAT(IMAT+NFBAS(I,0)*NPSH(I)),
     &       NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRINT)

         IF(IPRINT.GE.2) THEN
          CALL HEADER(
     &    'GET_A_FAC: Diagonal (MUST BE!) RKB/<k_e|k_e>=A+.A'//
     &    ' overlap matrix, VMAT',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')'* Fermion irp no.',I,'/',NFSYM
          CALL TEST_DIAG_MTX(VMAT(IKUN1),NESH(I),NESH(I),NESH(I),1)
          IF(IPRINT.GE.7) THEN
           CALL PRQMAT(VMAT(IKUN1),NESH(I),NESH(I),
     &            NESH(I),NESH(I),1,IPQTOQ(1,0),LUPRI)
          ENDIF
         ENDIF

*************************************************************
*        Extract the A kinematic factors into FACKIN(1,*)   *
*************************************************************
         IF (IPRINT.GE.2) THEN 
           WRITE(LUPRI,'(/9X,A,I1,A,I1,//,12X,A/)')
     & 'GET_A_FAC: "A" kinem. electr. factors; ferm irp:',I,'/',NFSYM,
     & 'A^2_ovlp          A^2_calc          diff;'
         ENDIF

         SUMDIFF = D0
         DO J = 1, NESH(I)
          DAA = VMAT(IKUN1+(J-1)*NESH(I)+J-1)
          FACKIN(IFC+J) = DAA
          X1 = (EIG(IORB(I)+NESH(I)+J)/CVAL2) + D2
          X2 = D2*(EIG(IORB(I)+NESH(I)+J)/CVAL2) + D2
          XD = (X1/X2)
          DIFF = FACKIN(IFC+J)-XD
          SUMDIFF = SUMDIFF + DABS(DIFF)

          IF (IPRINT.GE.2) THEN
          WRITE(LUPRI,"(5x,2F18.10,4X,D11.5)")
     &         FACKIN(IFC+J),XD,DIFF
          ENDIF

         ENDDO 
         SUMDIFF = SUMDIFF/DFLOAT(NESH(I))         

         WRITE(LUPRI,'(/2X,A,I1,A,I1,A,D14.8)')
     & 'Average differences for A-factors (fsym ',
     & I,'/',NFSYM,') : ',SUMDIFF
         CALL FLSHFO(LUPRI)

         IF (SUMDIFF.GE.THRNULL) THEN
          WRITE(LUPRI,'(2X,A,D14.8)')
     & 'GET_A_FAC:  DIFFERENCES A^2_ovlp-A^2_calc/N DIFFER'//
     & ' TOO MUCH FROM "ZERO"=!',THRNULL
         ENDIF

         IFC   = IFC + NESH(I)
         IKUN1 = IKUN1 + NESH(I)*NESH(I)
         IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
         IBUF = IBUF + NFBAS(I,0)*NFBAS(I,0)

      ENDIF
      ENDDO

C     CALL QUIT('GET_A_FAC')
      CALL QEXIT('GET_A_FAC')

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck mksaomo_le  */
      SUBROUTINE MKSAOMO_Le(TBUF,TMAT,VMAT,FACKIN,FPB,
     &                    IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     PURPOSE:
C    ==========
C        Construct half transformed overlap matrix, 
C         Saomo = Sao(LL) * C(L,el) * Ad(ee)^(-1) (FPB=.true.)
C         Saomo = Sao(LL) * C(L,el)   (FPB=.false.)
C
C    On input:
C   -----------
C       TMAT - transformation matrix 
C       FACKIN  Ad(ee) _ if FPB=.true.
C
C       TBUF -  
C             
C
C    On output:
C   ------------
C
C      Written:  M.Ilias, Odense, 2003  (based on MKSAOMO by L. Visscher)
C                Last update by M.Ilias, Strasburg, 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00 )
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION TMAT(*),VMAT(*),TBUF(*),FACKIN(*),WORK(LWORK)
      LOGICAL  FNDLAB,FPB

      CALL QENTER('MKSAOMO_Le')
#include "memint.h"

C     Get AO_LL overlap matrix , S(LL), WITH SS part ZEROED, into VMAT
      CALL GTOVLT(VMAT,D0,WORK(KFREE),LFREE)
C     CALL GTOVLT(VMAT,D1,WORK(KFREE),LFREE)

C  ... read "k" in "xhi" from BSSMAT into TMAT
C     REWIND LUBSS
C     IF (FNDLAB('FPMAT4C ',LUBSS)) THEN
C        ISIZE = NFBAS(1,0)*NORB(1)*NZ
C        IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NORB(2)*NZ
C        CALL READT(LUBSS,ISIZE,TMAT)
C        IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
C    &   'FPMAT4C  read.'
C     ELSE
C       CALL QUIT('FPMAT4C  not found !')
C     END IF

C
C   S_LL(NFBAS(I,0),NFBAS(I,0)) * TMAT_e(NFBAS(I,0),NESH(I))
C

      IOFF = 1
      JOFFESH = 0
      ITMT = 1
      ISAO = 1
      DO I = 1,NFSYM
      IF ( NESH(I) .GT. 0 ) THEN

         IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &    'MKSAOMO_Le: ENTERING ON RKB/"k" matrix, TMAT,'//
     &    ' in 4c SA-AO basis "xhi" ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TMAT(ITMT+NFBAS(I,0)*NESH(I)),NFBAS(I,0),NESH(I),
     &               NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     &   'MKSAOMO_Le: Read SOVLP matrix(LL) in SO basis "xhi" ',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(ISAO),NFBAS(I,0),NFBAS(I,0),
     &                NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),LUPRI)
         ENDIF

C            SAOMO_Le = S_LL * TMAT
C             TBUF   =  VMAT * TMAT

         CALL QGEMM(NFBAS(I,0),NESH(I),NFBAS(I,0),D1,
     &         'N','N',IPQTOQ(1,0),VMAT(ISAO),
     &         NFBAS(I,0),NFBAS(I,0),1,
     &         'N','N',IPQTOQ(1,0),
     &         TMAT(ITMT+NFBAS(I,0)*NESH(I)),
     &         NFBAS(I,0),NESH(I),NZ,
     &         D0,IPQTOQ(1,0),TBUF(IOFF),
     &         NFBAS(I,0),NESH(I),NZ)
C

            IF ( IPRINT .GE. 6) THEN
             CALL HEADER('MKSAOMO_Le: Right index transf. overlap'//
     &         ' matrix Saomo_Le',-1)
            WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
             CALL PRQMAT(TBUF(IOFF),NFBAS(I,0),NESH(I),
     &            NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
             CALL FLSHFO(LUPRI)
            END IF

C            Now multiply by diagonal Ad^-1, FACKIN(j)  = Ad^2(j,j)
C               SAOMO_Le = SAOMO_Le * Ad^-1
          IF (FPB) THEN
           ICX = 0
           DO IZ = 1,NZ
           DO IESH = 1, NESH(I)
           DO IFBAS = 1, NFBAS(I,0) 
            TBUF(IOFF+ICX) = 
     &      TBUF(IOFF+ICX)/DSQRT(FACKIN(IESH+JOFFESH))
            ICX = ICX + 1
           ENDDO
           ENDDO 
           ENDDO 
          ENDIF

         IF ( IPRINT .GE. 5 ) THEN
         CALL HEADER('MKSAOMO_Le: FINAL transf. matrix '//
     &   'H2c SA_AO "xhi"  -> H2c ON RKB/FP "k", TBUF ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
            CALL PRQMAT(TBUF(IOFF),NFBAS(I,0),NESH(I),
     &            NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
            CALL FLSHFO(LUPRI)
          END IF

        JOFFESH = JOFFESH + NESH(I)
        IOFF = IOFF + NFBAS(I,0)*NESH(I)*NZ
        ITMT = ITMT + NFBAS(I,0)*NORB(I)*NZ
        ISAO = ISAO + NFBAS(I,0)*NFBAS(I,0)

      ENDIF
      ENDDO

C  ... write TBUF to the file...
      REWIND LUBSS
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
         BACKSPACE(LUBSS)
         CALL NEWLAB('TM2CAOFP',LUBSS,LUPRI)
         ISIZE = NFBAS(1,0)*NESH(1)*NZ
         IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NESH(2)*NZ
         CALL WRITT(LUBSS,ISIZE,TBUF)
         CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
         IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &   'MKSAOMO_Le: TM2CAOFP written into BSSMAT, EOFLABEL renewed.'
      ELSE
        CALL QUIT('MKSAOMO_Le: EOFLABEL not found !')
      END IF


C     CALL QUIT('MKSAOMO_Le')
      CALL QEXIT('MKSAOMO_Le')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_h1_fpm */      
       SUBROUTINE GET_H1_ONB(TBUF,TMAT,VMAT,IBEIG,
     &                  SF_BEG,FPB,SAVESO1,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: Get the H1 Hamiltonian in the FP ON "k" basis (FPB=.TRUE.)
C    ======== or in the ON RKB MO basis, 
C             if wished (SF_BEG=.true.),eliminate spin-orbit terms 
C             from the beginning using the IBEIG boson info
C
C  On input:
C -----------
C   IBEIG - boson info on MO shells (needed for the spin-free case)
C
C   TMAT  - "k" eigenvectors in "xhi" (when FPB=.true.), or RKB basis
C                                                         
C   FPB  - is .true. the "k" FP basis is used, otherwise RKB basis                                                         
C
C  On output:
C ------------
C     TBUF - H1 Hamiltonian in the FP ON "k" basis (spin-free,if SF_BEG .true.)                                            
C                           in RKB ON basis (when FPB=.false.) 
C 
C 
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0 ) 
#include "dcbgen.h"  
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"  
C ... needed due to the BSS,USEDF 
#include "dcbham.h"  
#include "dcbdhf.h"  
      DIMENSION TBUF(*),TMAT(*),VMAT(*),IBEIG(*),WORK(LWORK)    
      LOGICAL SF_BEG, SAVESO1, FPB, FNDLAB

      CALL QENTER('GET_H1_ONB')
#include "memint.h"         

C     ... it is in the dcbham.h block...
      IF (.NOT.USEDF) THEN

C First get the bare nucleus Dirac Hamiltonian (BNDH) in the SA-AO basis "xhi"
C into VMAT array

       BSS=.FALSE.
       CALL ONEFCK(VMAT,IPRHAM,WORK(KFREE),LFREE) 
       BSS=.TRUE.

C
C    transform the BNDH from its SO-basis "xhi" to the ON basis "|k>"
C    (reuse TBUF)
C  where H is the bare nucleus Hamiltonian in the SO basis "xhi"
C
C        TBUF = VMAT+ * TMAT  * VMAT; 
C  H1(in |k>) =   U1+ * H0 * U1 (notation in our paper)
C
C  if desired, eliminate spin-orbit components from H1 in C1 symmetry
C  - spin-free at the "beginning"
C

      IF (SF_BEG) THEN
         MZ = 1
      ELSE
         MZ = NZ
      END IF

      IOFT = 1
      IMAT = 1

      DO I = 1,NFSYM
      IF(NORB(I).GT.0) THEN

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &   'GET_H1_ONB: Entering the ON RKB/FP "k" matrix,'//
     &   ' TMAT, in 4c SO-AO basis "xhi" ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'GET_H1_ONB: Read Dirac bare nucleus matrix,'//
     &   ' TMAT, in 4c SO-AO basis "xhi" ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(1+I2BASX(I,I)),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &       NORB(I),NORB(I),
     &       VMAT(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &       TBUF(IOFT),NORB(I),NORB(I),MZ,IPQTOQ(1,0),
     &       TMAT(IMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &       TMAT(IMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'GET_H1_ONB: Obtained Dirac H1 matrix in '//
     &  '4c "k"/ON basis, TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C       ... update the offsets
        IOFT = IOFT + NORB(I)*NORB(I)*NZ
        IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ

      ENDIF
      ENDDO

      ELSE
C ...  read converged Fock-Dirac MO matrix from the DFFOCK file
        OPEN (LUFOCK, FILE='DFFOCK', FORM='UNFORMATTED',
     +        ACCESS='DIRECT',RECL=8*N2TMOTQ,STATUS='OLD')
C  ... HEEELP ! IREC=??
C ... RESORT NEEDED !!!!
       IREC = -1
       IF(MXDIIS.EQ.0) THEN
          IREC = 1
        ELSE
          IREC = MOD(NITER-2,MXDIIS-1) + 1
        ENDIF
         write(lupri,*) 'N2TMOTQ=',N2TMOTQ,' IREC=',IREC
        CALL READAC(LUFOCK,N2TMOTQ,TBUF,IREC)
        IOFT = 1
        DO I = 1, NFSYM
        IF (NORB(I).GT.0) THEN
         IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'GET_H1_ONB: Read(DFFOCK) Fock-Dirac matrix in '//
     &  '4c ON basis, TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
        
         IOFT = IOFT + NORB(I)*NORB(I)*NZ
        ENDIF
        ENDDO
C      call quit('GET_H1_ONB: USEDF ... !')
      ENDIF

C     ... do extract (ee) elements from TBUF into VMAT 
C         and eliminate SF terms to get pure SO1 matrix !
      IF (SAVESO1) THEN
  
       IOFT = 1
       IEE  = 1
       DO I=1,NFSYM
       IF (NORB(I).GT.0) THEN

         N2ESHI = NESH(I)*NESH(I)
         N2T    = NORB(I)*NORB(I)
         DO IZ = 1, NZ
           IBUF = IEE + (IZ-1)*N2ESHI
           JH1EE = IOFT + (NPSH(I)*NORB(I))+NPSH(I)
     &          + (IZ-1)*N2T
           DO J=1,NESH(I)
           CALL DCOPY(NESH(I),TBUF(JH1EE),1,VMAT(IBUF),1)
           JH1EE = JH1EE + NORB(I)
           IBUF  = IBUF  + NESH(I)
           ENDDO
         ENDDO

         IF (IPRHAM.GE.4) THEN
           CALL HEADER(
     & 'GET_H1_ONB: (ee) part of the H1 matrix in ON RB/"k" basis '//
     & ' before eliminating spin-free terms (SAVESO1=.true.),VMAT',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(VMAT(IEE),NESH(I),NESH(I),
     &         NESH(I),NESH(I), NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

         CALL SPF_HK(VMAT(IEE),.FALSE.,IBEIG(IORB(I)+1+NPSH(I)),
     &            IBEIG(IORB(I)+1+NPSH(I)),
     &            NESH(I),NESH(I),NESH(I),NESH(I),NZ,IPRHAM)

         IF (IPRHAM.GE.4) THEN
           CALL HEADER(
     & 'GET_H1_ONB: (ee) part of the H1 matrix in ON RB/"k" basis '//
     & ' after eliminating spin-free terms-pure SO1 Ham., VMAT',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(VMAT(IEE),NESH(I),NESH(I),
     &         NESH(I),NESH(I), NZ,IPQTOQ(1,0),LUPRI)
         ENDIF
         IEE  = IEE + NESH(I)*NESH(I)*NZ
         IOFT = IOFT + NORB(I)*NORB(I)*NZ
       ENDIF
       ENDDO
C      ... do save H_SO1(ee) term into the file BSSMAT ...
      REWIND LUBSS
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
       BACKSPACE(LUBSS)
       ISIZE = NESH(1)*NESH(1)*NZ
       IF (NFSYM.EQ.2) ISIZE = ISIZE + (NESH(2)*NESH(2)*NZ)
       CALL NEWLAB('HSO1_ONB',LUBSS,LUPRI)
       CALL WRITT(LUBSS,ISIZE,VMAT)
       CALL NEWLAB('H1_ONB  ',LUBSS,LUPRI)
       ISIZE=NORB(1)*NORB(1)*NZ
       IF (NFSYM.EQ.2) ISIZE=ISIZE+NORB(2)*NORB(2)*NZ
       CALL WRITT(LUBSS,ISIZE,TBUF)
       CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
       WRITE(LUPRI,'(/2X,A)')
     & 'GET_H1_ONB: H1_SO1 & H1_ONB in ON RKB/"k" basis have been '//
     & 'written to the BSSMAT file, EOFLABEL is renewed.'
      ENDIF
C      call quit('GET_H1_ONB: H1_SO1 ...')
      ENDIF

C .... if spin-free, eliminate spin-orbit dependent terms
      IF (SPINFR.AND.SF_BEG) THEN
        IOFT = 1
        DO I = 1, NFSYM
        IF (NORB(I).GT.0) THEN
         CALL SPF_HK( TBUF(IOFT),.TRUE.,IBEIG(IORB(I)+1),
     &                IBEIG(IORB(I)+1),
     &             NORB(I),NORB(I),NORB(I),NORB(I),NZ,IPRHAM )
         IF(IPRHAM.GE.5) THEN
          CALL HEADER(
     &  'GET_H1_ONB: H1 Dirac matrix in 4c ON RKB/"k" basis '//
     &  ' after elim. spin-orbit dependent terms "BEFORE", TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          write(lupri,*) 'boson irreps: ifsym=',I,':',
     &                   ( (IBEIG(IORB(I)+J)),J=1,NORB(I))
         ENDIF
         IOFT = IOFT + NORB(I)*NORB(I)*NZ
        ENDIF
        ENDDO
      ENDIF

C ... write the TBUF - H1 ON RKB/"k" to the file...
      REWIND LUBSS
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
       BACKSPACE(LUBSS)
       CALL NEWLAB('H1_ONBAS',LUBSS,LUPRI)
       ISIZE=NORB(1)*NORB(1)*NZ
       IF (NFSYM.EQ.2) ISIZE=ISIZE+NORB(2)*NORB(2)*NZ
       CALL WRITT(LUBSS,ISIZE,TBUF)
       CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
       WRITE(LUPRI,'(/2X,A)')
     &'GET_H1_ONB:  H1_ONBAS (in ON RKB/"k") matrix is '//
     &' written to the BSSMAT file, EOFLABEL renewed.'
       CALL FLSHFO(LUPRI)
      ELSE
        CALL QUIT(
     & 'GET_H1_ONB:Error, EOFLABEL not found in BSSMAT file!')
      END IF


CMI   CALL QUIT('GET_H1_ONB')
      CALL QEXIT('GET_H1_ONB')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_hinf_fpm */      
      SUBROUTINE GET_HINF_BD(TBUF,TMAT,VMAT,IBEIG,
     &                   GETPCTM,SF_BEG,WORK,LWORK,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: Get the H_inf 2c Hamiltonian in the FP ON "k" basis
C    ========  
C             - save H "k"_block diagonal (full, both pp and ee blocks)
C             - save 4c "k"->2c "k" (full) picture change transf.matrix
C             - save 4c "xhi" -> 2c "k" (only ee) picture change transf.matrix 
C             - save Hinf(2c) in "k" basis to the file
C
C  On input:
C --------------
C    TBUF - BNDH in "k" basis
C    VMAT "k" in SA-AO "xhi" basis
C
C  On output:
C --------------
C
C Written by:  M.Ilias, HJ Aa Jensen, Odense, 2002,
C            last revison: M.Ilias, Strasbourg, february 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER( D1 = 1.0D0, DM1= -1.0D0, D0 = 0.0D0,THRNULL = 1.0D-10) 
#include "dcbgen.h"  
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"   
      DIMENSION TBUF(*),TMAT(*),VMAT(*),IBEIG(*),
     &          WORK(LWORK)    
      DIMENSION IBL(1,4)
      LOGICAL   SPINFR_SAVE,DBGDIAG,GETPCTM,SF_BEG,FNDLAB

      CALL QENTER('GET_HINF_BD')
#include "memint.h"         

      NBL = 1

      IOFT = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN

        IF(IPRINT.GE.5) THEN
          CALL HEADER(
     &  'GET_HINF_FPM: Entering H1 Dirac matrix in 4c "k" basis '//
     &  'TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
       
C     .. starting point
        IBL(1,1)=1
        IBL(1,2)=NPSH(I)+1
C     ... ending point
        IBL(1,3)=NPSH(I)
        IBL(1,4)=NORB(I)

C     ... problem ! it takes too long !
        CALL QJACOBI(TBUF(IOFT),TMAT(IOFT),NORB(I),NZ,NBL,IBL,.TRUE.,1)

        IF(IPRINT.GE.5) THEN
          CALL HEADER(
     &  'GET_HINF_FPM: H1 FULL BLOCK DIAGONAL Dirac matrix '//
     &  'in 4c "k" basis, TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     &  'GET_HINF_FPM: Picture change tranf.matrix,  '//
     &  'full 4c "k" -> 2c "k", TMAT',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TMAT(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)

        ENDIF
   
        IOFT = IOFT + NORB(I)*NORB(I)*NZ

      ENDIF
      ENDDO

      ISIZE1=NFBAS(1,0)*NORB(1)*NZ
      ISIZE2=NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2) THEN
        ISIZE1=ISIZE1+NFBAS(2,0)*NORB(2)*NZ
        ISIZE2=ISIZE2+NORB(2)*NORB(2)*NZ
      ENDIF

C ... save H1 block diagonal as well as picture change matrix
      CALL NEWLAB('H1FPM_BD',LUBSS,LUPRI)
      CALL WRITT(LUBSS,ISIZE2,TBUF)
      WRITE(LUPRI,'(/2X,A)')
     &'H1 "k" block-diagonal  matrix is '//
     &' written to the BSSMAT file.'

      CALL NEWLAB('PCTM_FPM',LUBSS,LUPRI)
      CALL WRITT(LUBSS,ISIZE1,TMAT)
      WRITE(LUPRI,'(/2X,A)')
     &'4c "k" -> 2c "k" picture change  '//
     &'transformation matrix written to the BSSMAT file.'

C  ... get the picture change transformation matric 4c "xhi" -> 2c "k" K.U 
C  ... into the TBUF = VMAT*TMAT

C     CALL DZERO(TBUF,N2BBASXQ)

      IMAT=1
      IPC =1 
      IOFT = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN

        IF(IPRINT.GE.5) THEN
          CALL HEADER(
     & 'GET_HINF_FPM: Entering |k> matrix, VMAT, in SO basis "xhi" ',-1)
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',I
          CALL PRQMAT(VMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        CALL QGEMM(NFBAS(I,0),NESH(I),NORB(I),D1,
     &         'N','N',IPQTOQ(1,0),VMAT(IMAT),
     &         NFBAS(I,0),NORB(I),NZ,
     &         'N','N',IPQTOQ(1,0),
     &         TMAT(IPC+NORB(I)*NPSH(I)),
     &         NORB(I),NESH(I),NZ,
     &         D0,IPQTOQ(1,0),TBUF(IOFT),
     &         NFBAS(I,0),NESH(I),NZ)

        IF(IPRINT.GE.5) THEN
          CALL HEADER(
     & 'GET_HINF_FPM: 4c "xhi" - 2c "k" picture change'//
     & ' transformation matrix ',-1)
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',I
          CALL PRQMAT(TBUF(IOFT),NFBAS(I,0),NESH(I),
     &             NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
        IPC  = IPC  + NORB(I)*NORB(I)*NZ
        IOFT = IOFT + NFBAS(I,0)*NESH(I)*NZ

      ENDIF
      ENDDO
 
C  ... save the picture change transformation matrix, TBUF...
      ISIZE1=NFBAS(1,0)*NESH(1)*NZ
      IF (NFSYM.EQ.2) ISIZE1=ISIZE1+NFBAS(2,0)*NESH(2)*NZ
      CALL NEWLAB('PCTM4C2K',LUBSS,LUPRI)
      CALL WRITT(LUBSS,ISIZE1,TBUF)
C     ... write only label for marking end of the file...
      CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
      WRITE(LUPRI,'(/2X,A)')
     &'4c SA-AO"xhi" -> 2c "k" picture change '//
     &'transformation matrix (PCTM4C2K) is written to the BSSMAT file.'
       WRITE(LUPRI,'(2X,A)') 'EOFLABEL is added.'

C  ... read the block diagonal matrix to extract the (ee) block only !
      REWIND LUBSS
      IF (FNDLAB('H1FPM_BD',LUBSS)) THEN
         CALL READT(LUBSS,ISIZE2,TBUF)
         IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &   'GET_HINF_FPM: H1FPM_BD read.'
      ELSE
        CALL QUIT(
     &  'GET_HINF_FPM: H1FPM_BD not found in BSSMAT file!')
      END IF

      IBUF = 1
      IMAT = 1
      DO I = 1, NFSYM

       IF(IPRINT.GE.5) THEN
          CALL HEADER(
     & 'GET_HINF_FPM: H1 "k" BLOCK DIAGONAL '//
     & ' Hamiltonian, TBUF ',-1)
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',I
          CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),
     &          NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C      ... copy TBUF -> VMAT
       DO IZ=1,NZ
       DO J=1,NESH(I)
        CALL DCOPY(NESH(I),
     &    TBUF(IBUF+(NORB(I)*NPSH(I))+NESH(I)+
     &        (IZ-1)*NORB(I)*NORB(I)+(J-1)*NORB(I)),1,
     &    VMAT(IMAT+(IZ-1)*NESH(I)*NESH(I)+(J-1)*NESH(I)),1)
       ENDDO
       ENDDO

       IF(IPRINT.GE.5) THEN
          CALL HEADER(
     & 'GET_HINF_FPM: H2c "k" (ee)'//
     & ' Hamiltonian, VMAT ',-1)
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',I
          CALL PRQMAT(VMAT(IMAT),NESH(I),NESH(I),
     &          NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

       IBUF = IBUF + NORB(I)*NORB(I)*NZ
       IMAT = IMAT + NESH(I)*NESH(I)*NZ
      ENDDO

      REWIND LUBSS
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
       BACKSPACE(LUBSS)
       CALL NEWLAB('H2C_FPM ',LUBSS,LUPRI)
       ISIZE3 = NESH(1)*NESH(1)*NZ
       IF (NFSYM.EQ.2) ISIZE3 = ISIZE3 + NESH(2)*NESH(2)*NZ
       CALL WRITT(LUBSS,ISIZE3,VMAT)
       CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
       WRITE(LUPRI,'(/2X,A)')
     &'H2c "k" (ee) written to BSSMAT, EOFLABEL renewed '
      ELSE
        CALL QUIT(
     &  'GET_HINF_FPM:Error, EOFLABEL not found in BSSMAT file!')
      END IF


C     CALL QUIT('GET_HINF_FPM')
      CALL QEXIT('GET_HINF_BD')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_hinf_fpm */      
      SUBROUTINE GET_HINF_R(TBUF,TMAT,VMAT,FACKIN,EIG,IBEIG,FPB,
     &                   YREQ1,GETPCTM,SF_BEG,WORK,LWORK,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  PURPOSE: Get the H_inf 2c Hamiltonian in the  ON RKB/"k" basis
C  ========  
C  Use the method of "R" solving to get picture change (in ON basis)
C  transformation matrix U and after applying U^+ H1 U get the H_inf 2c
C  Hamiltonian in a given ON basis.
C
C    On input:
C --------------
C    TBUF - BNDH in ON RKB / FP "k" basis
C
C   On output:
C --------------
C    TBUF - 2 comp. infinite order Hamiltonian in ON RKB / FP "k" basis (ee block)
C
C    Written by:  M.Ilias, HJ Aa Jensen, Odense, 2002,
C  last revison: M.Ilias, T.Saue, Strasbourg, february 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER( D1 = 1.0D0, DM1= -1.0D0, D0 = 0.0D0,THRNULL = 1.0D-10) 
#include "dcbgen.h"  
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"   
      DIMENSION TBUF(*),TMAT(*),VMAT(*),FACKIN(*),EIG(*),IBEIG(*),
     &          WORK(LWORK)    
      LOGICAL  SPINFR_SAVE,DBGDIAG,GETPCTM,SF_BEG, YREQ1,
     &          FNDLAB,FPB,FP2TH,FP2TH1,FP2TH2,USEJACO

      CALL QENTER('GET_HINF_R')
#include "memint.h"         


C     USEJACO=.TRUE.
      USEJACO=.FALSE.

      IBUF = 1
      IEIG = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN     

        IF(IPRINT.GE.5) THEN
         CALL HEADER(
     &  'GET_HINF_R (0): ENTERING Dirac H1 bare nucl. matrix in '//
     &  ' ON RKB/FP "k" basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),NORB(I),
     &                NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C  ... do diagonalize the H1 in ON RKB/"k" basis
        IF (USEJACO.AND.NZ.EQ.1) THEN
C       ... use Jacobi diagonalization...
         CALL RSJACO(NORB(I),NORB(I),NORB(I),TBUF(IBUF),
     &                  EIG(IEIG),1,1,0,TMAT(IBUF),
     &                  WORK(KFREE),LFREE)
        ELSE
         CALL QDIAG(NZ,NORB(I),TBUF(IBUF),NORB(I),NORB(I),
     &        EIG(IEIG),1,TMAT(IBUF),NORB(I),NORB(I),
     &        WORK(KFREE),LFREE,IERR)
        ENDIF

        IF(IPRINT.GE.5) THEN
          CALL HEADER(
     &   'GET_HINF_R (0): Eigenvalues of H1 bare nucl. matrix in '//
     &   ' ON "k"/RKB basis, EIG',-1)
         IX=0
         DO J=1,NORB(I)
          WRITE(LUPRI,'(2X,I4,2X,F16.6)') J-1-NPSH(I)+IX,EIG(IEIG+J-1)
          IF (J.EQ.NESH(I)) IX=1
         ENDDO

         CALL HEADER(
     &   'GET_HINF_R (0): Eigenvectors of H1 bare nucl. matrix in '//
     &   'ON "k"/RKB basis, TMAT',-1)
         CALL PRQMAT(TMAT(IBUF),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IBUF = IBUF + NORB(I)*NORB(I)*NZ
       IEIG = IEIG + NORB(I)

      ENDIF
      ENDDO

C  .... next (control) cycle ... do transform HDBN "k" eigenvectors 
C       to the HDBN "theta" eigenvectors
C     FP2TH=.TRUE.
      FP2TH=.FALSE.
C     FP2TH1=.FALSE.
      IF (FP2TH.AND.FPB) THEN
C ... read the FPMAT_TH 4c "k"-> 4c "theta" transf. matrix into TBUF
        REWIND(LUBSS)
        IF (FNDLAB('FPMAT_TH',LUBSS)) THEN
        ISIZE = NORB(1)*NORB(1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
        CALL READT(LUBSS,ISIZE,TBUF)
        WRITE(LUPRI,'(/2X,A)')
     & 'GET_HINF_R (TEST): 4c "k"-> 4c "theta" transf. '//
     & '  matrix FPMAT_TH was read from the BSSMAT file into TBUF.'
        ELSE
         CALL QUIT(
     &  'GET_HINF_R: Error, EOFLABEL not found in BSSMAT file!')
        END IF

        IBUF = 1
        
        DO I=1,NFSYM
        IF (NORB(I).GT.0) THEN

          IF(IPRINT.GE.5) THEN

          CALL HEADER(
     &   'GET_HINF_R (T1): Eigenvectors of H1 bare nucl. matrix in '//
     &   'ON "k"/RKB basis, TMAT',-1)
         CALL PRQMAT(TMAT(IBUF),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     &  'GET_HINF_R (T1): 4c "k"-> 4c "theta" transf.'//
     & ' matrix FPMAT_TH  in TBUF',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),NORB(I),
     &                NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

C       ... VMAT: =  TBUF*TMAT
        CALL QGEMM(NORB(I),NORB(I),NORB(I),D1,
     &     'N','N',IPQTOQ(1,0),
     &     TBUF(IBUF),NORB(I),NORB(I),NZ,
     &     'N','N',IPQTOQ(1,0),
     &      TMAT(IBUF),NORB(I),NORB(I),NZ,
     &      D0,IPQTOQ(1,0),VMAT(IBUF),NORB(I),NORB(I),NZ)

          IF(IPRINT.GE.5) THEN
            CALL HEADER(
     &  'GET_HINF_R (T1): Eigenvectors of H1("k") transformed'//
     &  ' to the "theta" basis, VMAT',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(VMAT(IBUF),NORB(I),NORB(I),NORB(I),
     &                NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
            IBUF = IBUF + NORB(I)*NORB(I)*NZ
        ENDIF
        ENDDO
C      CALL QUIT('FP2TH ')
      ENDIF


      IBUF = 1
      IBUFE = 1
      IEIG = 1
      IBMAT = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN     
       NREF   = NESH(I)
       N2REF  = NREF*NREF
       N2EQ   = NREF*NREF*NZ
       N2REFQ = N2REF*NZ

       IF (YREQ1) THEN
C ... make A and B matrixes as in the original (HJAaJ) approach...
       CALL DUNIT2(TBUF(IBUFE),NREF,NREF,NREF,NZ)
C ... do A := A - Yp*(Yp)+ into TBUF(IBUFE), A = 1 - Y-+.Y-+^+
       CALL QGEMM(NREF,NREF,NREF,DM1, 
     &    'N','N',IPQTOQ(1,0),
     &     TMAT(IBUF+(NPSH(I)*NORB(I))),NORB(I),NORB(I),NZ,
     &    'H','N',IPQTOQ(1,0),
     &     TMAT(IBUF+(NPSH(I)*NORB(I))),NORB(I),NORB(I),NZ,
     &     D1,IPQTOQ(1,0),TBUF(IBUFE),NREF,NREF,NZ)

C  -- make b matrix (store in VMAT); b := Yp * (Ye)+ = Y-+.Y++^+
       CALL QGEMM(NREF,NREF,NREF,D1,
     &     'N','N',IPQTOQ(1,0),
     &     TMAT(IBUF+(NPSH(I)*NORB(I))),NORB(I),NORB(I),NZ,
     &     'H','N',IPQTOQ(1,0),
     &      TMAT(IBUF+(NPSH(I)*NORB(I))+NPSH(I)),NORB(I),NORB(I),NZ,
     &      D0,IPQTOQ(1,0),VMAT(IBUFE),NREF,NREF,NZ)

      ELSE
C       .....  Trond's suggestion (YR1=.false.)  .....
C        A = Y--.Y--^+
       CALL QGEMM(NREF,NREF,NREF,D1, 
     &    'N','N',IPQTOQ(1,0),
     &     TMAT(IBUF),NORB(I),NORB(I),NZ,
     &    'H','N',IPQTOQ(1,0),
     &     TMAT(IBUF),NORB(I),NORB(I),NZ,
     &     D0,IPQTOQ(1,0),TBUF(IBUFE),NREF,NREF,NZ)

C  -- make b matrix (store in VMAT); B := - Y--.Y+-^+
       CALL QGEMM(NREF,NREF,NREF,DM1,
     &     'N','N',IPQTOQ(1,0),
     &     TMAT(IBUF),NORB(I),NORB(I),NZ,
     &     'H','N',IPQTOQ(1,0),
     &      TMAT(IBUF+NPSH(I)),NORB(I),NORB(I),NZ,
     &      D0,IPQTOQ(1,0),VMAT(IBUFE),NREF,NREF,NZ)

       ENDIF

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &   'GET_HINF_R (1): A matrix for getting R in ON basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'GET_HINF_R (1): B matrix for getting R in ON basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

       IBUF = IBUF + NORB(I)*NORB(I)*NZ
       IBUFE = IBUFE + N2EQ
       IEIG = IEIG + NORB(I)
      ENDIF
      ENDDO

C ... now we have eigenvectors: do prepare for solving R equation !
C        (1-Y-+.Y-+^+) R = Y-+.Y++^+ or
C        (Y--.Y--^+) R = -Y--.Y+^+

      IEIG = 1
      IBUFE =1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN
       NREF   = NESH(I)
       N2REF = NREF*NREF
       N2EQ = NREF*NREF*NZ
       N2REFQ = N2REF*NZ

       CALL MEMGET('REAL',KW,NREF,WORK,KFREE,LFREE)  

C ... DOES NOT INFLUENCE...
C      THRNULL1 = 1.0D-20
C      THRNULL1 = 0.0D0
       JOB = 1
       CALL QCHOLD(TBUF(IBUFE),NREF,NZ,NREF,NREF,EIG(IEIG),
C    &                  THRNULL1,NEFF,JOB,IBEIG(IEIG))
     &                  THRNULL,NEFF,JOB,IBEIG(IEIG))

       IF(NEFF.LT.NREF) THEN
         WRITE(LUPRI,'(/2A)') '* GET_HINF_R WARNING: ',
     &   'linear dependencies detected in Cholesky decomposition !'
         WRITE(LUPRI,'(11X,A,I5,A,I5)')
     &   'Reduced dimensionality: ',NREF,' --> ',NEFF
         CALL QUIT('GET_HINF_R: Cholesky decompos.problem !')
C         NEFF = NREF ! Do not do that !
       ENDIF


       CALL QCHOLS(TBUF(IBUFE),NREF,NEFF,NREF,NZ,NREF,NREF,
     &           EIG(IEIG),VMAT(IBUFE),NREF,NREF,
     &                     TMAT(IBUFE),NREF,NREF,
     &           JOB,IBEIG(IEIG),WORK(KW))

C   ... Do print out the R
       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &   'GET_HINF_R (1): Obtained R after CHoleski, TMAT',-1)
         IF (.NOT.YREQ1) THEN
           write(lupri,*) "Trond's form of the YR equation."
         ELSE
           write(lupri,*) "HJAaJ's form of the YR equation."
         ENDIF
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL MEMREL('W    ',WORK,KWORK,KW,KFREE,LFREE)

       IEIG = IEIG + NREF
       IBUFE = IBUFE + N2EQ

      ENDIF
      ENDDO

C ... save the full R, TMAT, into the file BSSMAT !
      REWIND LUBSS
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
       BACKSPACE(LUBSS)
       CALL NEWLAB('R_MTX   ',LUBSS,LUPRI)
       ISIZE = NESH(1)*NESH(1)*NZ
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NESH(2)*NESH(2)*NZ
       CALL WRITT(LUBSS,ISIZE,TMAT(1))
       CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
       WRITE(LUPRI,'(/2X,A)')
     &'R matrix (ON basis), TMAT, written to BSSMAT, EOFLABEL renewed'
      ELSE
        CALL QUIT(
     & 'GET_HINF_FPM:Error, EOFLABEL not found in BSSMAT file!')
      END IF

C =========== TESTING "k" AND RKB approaches ======================
C      ... does not work for R !!!

      FP2TH1=.FALSE.
C     FP2TH1=.TRUE.
      IF (FP2TH1.AND.FPB) THEN
C ... read the FPMAT_TH 4c "k"-> 4c "theta" transf. matrix into TBUF
        REWIND(LUBSS)
        IF (FNDLAB('FPMAT_TH',LUBSS)) THEN
        ISIZE = NORB(1)*NORB(1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
        CALL READT(LUBSS,ISIZE,TBUF)
        WRITE(LUPRI,'(/2X,A)')
     & 'GET_HINF_R (T2): 4c "k"-> 4c "theta" transf. matrix FPMAT_TH '//
     & ' was read from the BSSMAT file into TBUF.'
        ELSE
         CALL QUIT(
     &  'GET_HINF_R: Error, EOFLABEL not found in BSSMAT file!')
        END IF

        IMAT = 1
        IBUF = 1
        DO I=1,NFSYM
        IF (NORB(I).GT.0) THEN
C      ... do TBUF. TMAT. TBUF^+ = VMAT (MOAO transf.)

C        IBUFA = NPSH(I)
         IBUFA = IBUF + NPSH(I)*NORB(I) 

         IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'GET_H_INF_R (T2): R in ON "k" basis ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,' /',NFSYM
           CALL PRQMAT(TMAT(IMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)

           CALL HEADER(
     &     'GET_H_INF_R (T2): Transf.mtx FPM'//
     &     ' "k"->RKB "theta" basis ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,' /',NFSYM
           CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),NORB(I),NORB(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)

           CALL HEADER(
     &     'GET_H_INF_R (T2): Transf.mtx FPM'//
     &     ' "k"->RKB "theta" basis, ep-block ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,' /',NFSYM
           CALL PRQMAT(TBUF(IBUFA),NPSH(I),NESH(I),NORB(I),NORB(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

         CALL QTRANS('MOAO','S',D0,NESH(I),NESH(I),
     &               NESH(I),NESH(I),
     &            VMAT(IMAT),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &            TMAT(IMAT),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &            TBUF(IBUFA),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &            TBUF(IBUFA),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)

         IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'GET_H_INF_R (T2): R in ON "theta" (transformed'//
     &     ' from "k") basis ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,' /',NFSYM
           CALL PRQMAT(VMAT(IMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                   NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

         IMAT = IMAT + NESH(I)*NESH(I)*NZ
         IBUF = IBUF + NORB(I)*NORB(I)*NZ

        ENDIF
        ENDDO
        CALL QUIT('R(SL): "k"->"theta" transform')
      ENDIF


C ... do generate the Omega+ matrix  
C ... first get TBUF = 1 + R^+.R matrix
      IBUFE = 1
      IEIG =  1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN
       NREF   = NESH(I)
       N2REF = NREF*NREF
       N2EQ = N2REF*NZ

C ... first make the unit matrix TBUF(IBUFE)
C .... perhaps try DUNIT2 ...
C      CALL DZERO(TBUF(IBUFE),N2REFQ)
C      CALL DUNIT(TBUF(IBUFE),NREF)
       CALL DUNIT2(TBUF(IBUFE),NREF,NREF,NREF,NZ)

C ... then do TBUF = TBUF + R^+.R
            CALL QGEMM(NREF,NREF,NREF,D1, 
     &          'H','N',IPQTOQ(1,0),TMAT(IBUFE),NREF,NREF,NZ,
     &          'N','N',IPQTOQ(1,0),TMAT(IBUFE),NREF,NREF,NZ,
     &          D1,IPQTOQ(1,0),TBUF(IBUFE),NREF,NREF,NZ)

C ... do print out 1 + R^+.R matrix
       IF (IPRINT.GE.7) THEN
         CALL HEADER( 'GET_HINF_R (2): 1 + R^+.R matrix, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C ... do diagonalize  1 + R^+.R, TBUF
       CALL QDIAG(NZ,NREF,TBUF(IBUFE),NREF,NREF,
     &                 EIG(IEIG),1,VMAT(IBUFE),NREF,
     &                 NREF,WORK(KFREE),LFREE,IERR)
       IF (IERR.NE.0) CALL QUIT('GET_HINF_R: QDIAG failed!')

       IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &  'GET_HINF_R (2): Eigenvalues of (1 + R+ * R)',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         DO J=1,NREF
           WRITE(LUPRI,'(I5,F25.10)') J,EIG(IEIG-1+J)
         ENDDO
       ENDIF

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &'GET_HINF_R (2): Eigenvectors of (1 + R+ * R)',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IBUFE),NREF,NREF,NREF,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C ... do copy eigenvectors into TBUF which is no longer needed
       CALL DCOPY(N2EQ,VMAT(IBUFE),1,TBUF(IBUFE),1)

C ... do rescale eigenvectors
       JEVEC = IBUFE
       DO IZ = 1, NZ
       DO J = 1,NREF
         FACTOR = D1 / DSQRT(EIG(IEIG-1+J))
C              ....  d^-1/2 * U(IZ)
         CALL DSCAL(NREF,FACTOR,TBUF(JEVEC),1)
         JEVEC = JEVEC + NREF 
       ENDDO
       ENDDO

C     now generate Omega+ = (X+ * X)^(-1/2) into TMAT (is reused)
C       .... U * d^-1/2 * U+ 
       CALL QGEMM(NREF,NREF,NREF,D1,  
     &   'N','N',IPQTOQ(1,0),TBUF(IBUFE),NREF,NREF,NZ,
     &   'H','N',IPQTOQ(1,0),VMAT(IBUFE),NREF,NREF,NZ,
     &    D0,IPQTOQ(1,0),TMAT(IBUFE),NREF,NREF,NZ)

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &'GET_HINF_R (2): Omega+ =  (1 + R+ * R)^(-1/2), TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IBUFE),NREF,NREF,NREF,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IBUFE = IBUFE + N2EQ
       IEIG  = IEIG + NREF

      ENDIF
      ENDDO

C ... read R into VMAT
      REWIND LUBSS
      IF (FNDLAB('R_MTX   ',LUBSS)) THEN
       ISIZE = NESH(1)*NESH(1)*NZ
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NESH(2)*NESH(2)*NZ
       CALL READT(LUBSS,ISIZE,VMAT)
       WRITE(LUPRI,'(/2X,A)')
     &'R matrix (in ON basis) was read from BSSMAT into VMAT.'
      ELSE
        CALL QUIT(
     & 'GET_HINF_FPM:Error, R_MTX not found in BSSMAT file!')
      ENDIF


C ... do make VMAT -> TBUF = (1 R)  !
      IBUFE  = 1
      IBUF2E = 1
      DO I=1, NFSYM
      IF (NORB(I).GT.0) THEN

       NREF  = NESH(I)
       N2REF = NESH(I)*NESH(I)
       NREF2 = 2*NESH(I)
       N2EQ = NREF*NREF*NZ

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &   'GET_HINF_R (3): R matrix read from BSSMAT, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C       now insert unit matrix in bottom half of X:
C    ... do make (1 R) into TBUF (R is in VMAT)
       JOFF =  IBUFE + NREF
       DO IZ = 1,NZ
        DO J = 1,NREF
C      ... do copy NREF into TBUF
         IVMAT = IBUFE + (IZ-1)*N2REF + (J-1)*NREF
         IBUF  = IBUF2E + (IZ-1)*2*N2REF + (J-1)*NREF2
         CALL DCOPY(NREF,VMAT(IVMAT),1,TBUF(IBUF),1)
C      ... add 1 and O
         IBUF2 = IBUF + NREF
         CALL DZERO(TBUF(IBUF2),NREF)
         IF (IZ.EQ.1) TBUF(IBUF2+J-1) = D1
        END DO
       END DO

       IBUFE = IBUFE + N2EQ
       IBUF2E = IBUF2E + 2*N2EQ

      ENDIF
      ENDDO


      IBUFE  = 1
      IBUF2E = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN

       NREF  = NESH(I)
       N2REF = NESH(I)*NESH(I)
       NREF2 = 2*NESH(I)
       N2EQ = NREF*NREF*NZ

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &'GET_HINF_R (4): (1 + R) matrix, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF2E),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C      and multiply it to form X Omega+ = (1  R)t * Omega+:
       CALL QGEMM(NREF2,NREF,NREF,D1,
     &          'N','N',IPQTOQ(1,0),TBUF(IBUF2E),NREF2,NREF,NZ,
     &          'N','N',IPQTOQ(1,0),TMAT(IBUFE),NREF ,NREF,NZ,
     &          D0,IPQTOQ(1,0),VMAT(IBUF2E),NREF2,NREF,NZ)


       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &   'GET_HINF_R (4): Final U1= (1 R)t*Omega+ picture change '//
     &   'transformation matrix in ON RKB/"k" basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IBUF2E),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IBUFE = IBUFE + N2EQ
       IBUF2E = IBUF2E + 2*N2EQ

      ENDIF
      ENDDO

C ... TEST: Do U+.U into TMAT
      IBUF2E=1
      IBUFE=1
      IBUF4E=1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN
        NREF = NESH(I)
        NREF2 = 2*NESH(I)

       CALL QGEMM(NREF,NREF,2*NREF,D1, 
     &    'H','N',IPQTOQ(1,0),
     &     VMAT(IBUF2E),2*NREF,NREF,NZ,
     &    'N','N',IPQTOQ(1,0),
     &     VMAT(IBUF2E),2*NREF,NREF,NZ,
     &     D0,IPQTOQ(1,0),TMAT(IBUFE),NREF,NREF,NZ)

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &   'GET_HINF_R : U1+ * U1 multiplied '//
     &   'transformation matrix in ON RKB/"k" basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IBUFE),NREF,NREF,NREF,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL QGEMM(NREF2,NREF2,NREF,D1, 
     &    'N','N',IPQTOQ(1,0),
     &     VMAT(IBUF2E),NREF2,NREF,NZ,
     &    'H','N',IPQTOQ(1,0),
     &     VMAT(IBUF2E),NREF2,NREF,NZ,
     &     D0,IPQTOQ(1,0),TMAT(IBUF4E),NREF2,NREF2,NZ)

C    ... this does not work !!! 
       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &   'GET_HINF_R : U1 * U1+ multiplied '//
     &   'transformation matrix in ON RKB/"k" basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IBUF4E),NREF2,NREF2,NREF2,NREF2,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

          IBUFE =  IBUFE + NESH(I)*NESH(I)*NZ
         IBUF2E = IBUF2E + 2*NESH(I)*NESH(I)*NZ
         IBUF4E = IBUF4E + NREF2*NREF2*NZ

      ENDIF
      ENDDO
C     CALL QUIT('U1+.U1')

C =================================================================================
C   ... do save 4c ON-> 2c ON picture change transformation matrix into the file
C =================================================================================
C  ... do save the picture change transformation matrix into BSSMAT
      REWIND LUBSS
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
       BACKSPACE(LUBSS)
       CALL NEWLAB('U1_ONBAS',LUBSS,LUPRI)
       ISIZE = 2*NESH(1)*NESH(1)*NZ
       IF (NFSYM.EQ.2) ISIZE = ISIZE + 2*NESH(2)*NESH(2)*NZ
       CALL WRITT(LUBSS,ISIZE,VMAT)
       CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
       WRITE(LUPRI,'(/2X,A)')
     & 'GET_HINF_R: U1_ONBAS (picture change transf.matrix '//
     & '4c ONB->2c ONB) is written to the BSSMAT'//
     & ' file, EOFLABEL is renewed.'
      ELSE
        CALL QUIT(
     & 'GET_HINF_R:Error, EOFLABEL not found in BSSMAT file!')
      END IF

C =============================================================================
C ... do the  V U1("k") (Ad-1 V+) transformation to get the U1 in "theta" basis
C
C    V-left is full transf.matrix, U1 is ep.e matrix and V-right is only (ee) block !
C
C        We do not get U1' matrix identical with its RKB counterpart !
C  (checked with Octave for 1s(identical) and 2s basis (not identical) )
C =============================================================================
      FP2TH2=.FALSE.
C     FP2TH2=.TRUE.
      IF (FP2TH2.AND.FPB) THEN
C ... read the FPMAT_TH 4c "k"-> 4c "theta" transf. matrix into TBUF
        REWIND(LUBSS)
        IF (FNDLAB('FPMAT_TH',LUBSS)) THEN
        ISIZE = NORB(1)*NORB(1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
        CALL READT(LUBSS,ISIZE,TBUF)
        WRITE(LUPRI,'(/2X,A)')
     & 'GET_HINF_R (T3): 4c "k"-> 4c "theta" transf. matrix FPMAT_TH '//
     & ' was read from the BSSMAT file into TBUF.'
        ELSE
         CALL QUIT(
     &  'GET_HINF_R (T3): Error, EOFLABEL not found in BSSMAT file!')
        END IF

        IBUF = 1
        IBUF2E = 1
        JOFFESH = 0
        DO I=1,NFSYM
        IF (NORB(I).GT.0) THEN
            
         NREF = NESH(I)
         NREF2 = NORB(I)

         IF(IPRINT.GE.5) THEN
          CALL HEADER(
     &   'GET_HINF_R (T3): Final U1= (1 R)t*Omega+ picture change '//
     &   'transformation matrix in ON RKB/"k" basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IBUF2E),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     &  'GET_HINF_R (T3): 4c "k"-> 4c "theta" transf.'//
     & ' matrix FPMAT_TH  in TBUF',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),NORB(I),
     &                NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

             write(lupri,*) 'FACKIN:',
     &     ( FACKIN(JESH+JOFFESH),JESH=1,NESH(I))

C         .... renormalize only ee block ????
            DO IZ = 1,NZ
C           write(lupri,*) 'IZ=',IZ
            DO IESH = 1, NREF ! run over columns...
CC           IBUFX = IBUF+(NREF*NREF2)+NREF+
             IBUFX = IBUF+(NREF*NREF2)+
     &         ((IESH-1)*NREF2)+(IZ-1)*NORB(I)*NORB(I)
            DO JESH = 1, NREF2
              TBUFOLD = TBUF(IBUFX+JESH-1)
C             TBUF(IBUFX+JESH-1) = 
C    &        TBUFOLD/DSQRT(FACKIN(IESH+JOFFESH))
             write(lupri,*) 'elem:',IBUFX+JESH-1,
     &               TBUFOLD,TBUF(IBUFX+JESH-1),FACKIN(IESH+JOFFESH)
            ENDDO
            ENDDO 
            ENDDO 

        IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &  'GET_HINF_R (T3): 4c "k"-> 4c "theta" transf.'//
     & ' matrix FPMAT_TH in TBUF after RENORMALIZATION',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),NORB(I),
     &                NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C   DO the V U1 V+ into TMAT

        CALL QTRANS('MOAO','S',D0,NREF2,NREF,
     &                            NREF2,NREF,
     &            TMAT(IBUF2E),NREF2,NREF,NZ,IPQTOQ(1,0),
     &            VMAT(IBUF2E),NREF2,NREF,NZ,IPQTOQ(1,0),
     &            TBUF(IBUF),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &            TBUF(IBUF),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)

        IF(IPRINT.GE.5) THEN
         CALL HEADER(
     &   'GET_HINF_R (T3): Final U1 picture change transformation '//
     &   'matrix TRANSFORMED BACK from ON "k" to RKB basis, TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IBUF2E),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)

        ENDIF

         JOFFESH = JOFFESH + NESH(I)
         IBUF = IBUF + NORB(I)*NORB(I)*NZ
         IBUF2E = IBUF2E + NESH(I)*NORB(I)*NZ

        ENDIF
        ENDDO
        call quit('GET_HINF_R: Tra. U1 transformation to "k"')
      ENDIF


C ... do prepare the resulting Hamiltonian in ON RKB/"k" basis
C   after the picture change transformation

C ... read the Dirac Hamiltonian in the ON RKB/"k" basis into  TMAT
      REWIND(LUBSS)
      IF (FNDLAB('H1_ONBAS',LUBSS)) THEN
        ISIZE = NORB(1)*NORB(1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
        CALL READT(LUBSS,ISIZE,TMAT)
        IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &  'GET_HINF_R: H1_ONBAS read into TMAT.'
      ELSE
        CALL QUIT(
     & 'GET_HINF_R:Error, H1_ONBAS not found in BSSMAT file!')
      END IF

      IBUF = 1
      IE   = 1
      ITMAT = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN
       NREF  = NESH(I)
       N2REF = NESH(I)*NESH(I)
       NREF2 = 2*NESH(I)
       N2EQ = NREF*NREF*NZ

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &   'GET_HINF_R (5): Final U1= (1 R)t*Omega+ picture change '//
     &   'transformation matrix in ON RKB/"k" basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(ITMAT),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'GET_HINF_R (5): H1 Dirac bare nucleus matrix '//
     &   'in ON RKB/"k" basis, TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IBUF),NREF2,NREF2,NREF2,NREF2,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL QTRANS('AOMO','S',D0,NREF2,NREF2,NREF,NREF,
     &        TMAT(IBUF),NREF2,NREF2,NZ ,IPQTOQ(1,0),
     &        TBUF(IE),  NREF ,NREF ,NZ ,IPQTOQ(1,0),
     &        VMAT(ITMAT),NREF2,NREF ,NZ, IPQTOQ(1,0),
     &        VMAT(ITMAT),NREF2,NREF ,NZ, IPQTOQ(1,0),
     &        WORK(KFREE),LFREE,IPRINT)

       IF (IPRINT.GE.4) THEN
           CALL HEADER(
     &   'GET_HINF_R (5): Obtained U1+ H1 U1 = H2c_infinite'//
     &   ' in ON RKB/"k" basis, TBUF (output)',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TBUF(IE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
        END IF

        IE  = IE + N2EQ
        ITMAT = ITMAT + 2*N2EQ 
        IBUF = IBUF + NORB(I)*NORB(I)*NZ

      ENDIF
      ENDDO

C ... save H2c_ON basis into BSSMAT
      REWIND LUBSS
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
       BACKSPACE(LUBSS)
       CALL NEWLAB('H2C_ONBA',LUBSS,LUPRI)
       ISIZE = NESH(1)*NESH(1)*NZ
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NESH(2)*NESH(2)*NZ
       CALL WRITT(LUBSS,ISIZE,TBUF)
       CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
       WRITE(LUPRI,'(/2X,A)')
     &'H2c_infinite (in ON basis) written to BSSMAT, EOFLABEL renewed '
      ELSE
        CALL QUIT(
     & 'GET_HINF_R:Error, EOFLABEL not found in BSSMAT file!')
      END IF

C     CALL QUIT('GET_HINF_R')
      CALL QEXIT('GET_HINF_R')
      RETURN
      END


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck h2cfinal */
      SUBROUTINE H2CFINAL(TBUF,TMAT,VMAT,FACKIN,IBEIG,SAVESO1,SPFR,
     &                    SPFR_BEG,GETPCTM,FPB,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     PURPOSE:
C    ==========
C
C     Do the transformation of H2c(ee) from ON MO basis to SA-AO(LL) basis.
C
C    On input:
C   -----------
C       TBUF -  H2c(ee) in ON MO basis RKB/"k"
C
C    On output: H2c in SA-AO is written to the file...
C   ------------
C
C    Written:  M.Ilias, Strasbourg 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00 )
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION TMAT(*),VMAT(*),TBUF(*),FACKIN(*),IBEIG(*),WORK(LWORK)
      LOGICAL FNDLAB, SAVESO1, SPFR,SPFR_BEG, GETPCTM, TRAFPRKB, FPB

      CALL QENTER('H2CFINAL')
#include "memint.h"


C     ... read H2C in ON RKB/FPM "k" basis into VMAT
C     REWIND LUBSS
C     IF (FNDLAB('H2C_ONBA',LUBSS)) THEN
C        ISIZE1 = NESH(1)*NESH(1)*NZ
C        IF (NFSYM.EQ.2) ISIZE1 = ISIZE1 + NESH(2)*NESH(2)*NZ
C        CALL READT(LUBSS,ISIZE1,VMAT)
C        IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
C    &   'H2C_ONBA read into VMAT.'
C     ELSE
C       CALL QUIT('H2CFINAL: H2C_ONBA not found !')
C     END IF


C ... assumes that H2C_ONBA is in FP "k" basis !
C     TRAFPRKB = .TRUE.
      TRAFPRKB = .FALSE.
C .... NEED TO INTERCHANGE  VARIABLES !!!
      IF (TRAFPRKB.AND.FPB) THEN
C  ... read the transformation matrix "k" -> "theta"
      REWIND LUBSS
      IF (FNDLAB('FPMAT_TH',LUBSS)) THEN
         ISIZE = NORB(1)*NORB(1)*NZ
         IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
         CALL READT(LUBSS,ISIZE,TBUF)
         IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &   'FPMAT_TH read into TBUF.'
      ELSE
        CALL QUIT('H2CFINAL: FPMAT_TH not found !')
      END IF

      IVMAT = 1
      IBUF = 1
      ITMAT = 1
      JOFFESH = 0
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN

         IBUFA = IBUF + (NORB(I)*NESH(I)) + NESH(I)

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'H2CFINAL: H2c in ON "k" basis ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',I,' /',NFSYM
            CALL PRQMAT(VMAT(IVMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)

           CALL HEADER(
     &     'H2CFINAL: Transf.mtx FPM "k"->RKB "theta" basis ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',I,' /',NFSYM
            CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),NORB(I),NORB(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

           write(lupri,*) 'FACKIN:',
     &     ( FACKIN(JESH+JOFFESH),JESH=1,NESH(I) )

C         .... renormalize only the "ee" block !
            DO IZ = 1,NZ
C           write(lupri,*) 'IZ=',IZ
            DO IESH = 1, NESH(I)
             IBUFX = IBUFA+((IESH-1)*NORB(I))+(IZ-1)*NORB(I)*NORB(I)
            DO JESH = 1, NESH(I) 
             TBUF(IBUFX+JESH-1) = 
     &       TBUF(IBUFX+JESH-1)/DSQRT(FACKIN(IESH+JOFFESH))
C            write(lupri,*) 'elem:',ICX,IBUFX+ICX,TBUF(IBUFX+ICX)
            ENDDO
            ENDDO 
            ENDDO 

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'H2CFINAL: RENORMALIZED FPM "k"->RKB "theta" basis'//
     &     ' ee block',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',I,' /',NFSYM
            CALL PRQMAT(TBUF(IBUFA),NESH(I),NESH(I),NORB(I),NORB(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

C         write(lupri,*) 'FACKIN:',(FACKIN(J),J=1,NESH(I))

        CALL QTRANS('MOAO','S',D0,NESH(I),NESH(I),
     &               NESH(I),NESH(I),
     &            TMAT(IVMAT),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &            VMAT(IVMAT),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &            TBUF(IBUFA),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &            TBUF(IBUFA),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'H2CFINAL: H2c in ON "theta" (transformed'//
     &     ' from "k") basis ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',I,' /',NFSYM
            CALL PRQMAT(TMAT(IVMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                 NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

        JOFFESH = JOFFESH  + NESH(I)
        IVMAT = IVMAT + NESH(I)*NESH(I)*NZ
        IBUF = IBUF + NORB(I)*NORB(I)*NZ

      ENDIF
      ENDDO 
C      CALL QUIT('H2CFINAL: See H2c "theta" froml H2c "k"')
      ENDIF


C      ... do add SO1 in MO into TBUF !!!
      IF (SAVESO1) THEN
        REWIND LUBSS
        IF (FNDLAB('HSO1_ONB',LUBSS)) THEN
          ISIZE = NESH(1)*NESH(1)*NZ
          IF (NFSYM.EQ.2) ISIZE = ISIZE + NESH(2)*NESH(2)*NZ
          CALL READT(LUBSS,ISIZE,TMAT)
          IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &    'H2CFINAL: HSO1_ONB read into TMAT.'
        ELSE
          CALL QUIT('H2CFINAL: HSO1_ONB not found in BSSMAT !')
        END IF

        IMAT = 1
        DO I=1,NFSYM
        IF (NESH(I).GT.0) THEN

         N2EQ = NESH(I)*NESH(I)*NZ

         IF (IPRINT.GE.5) THEN
           CALL HEADER(
     & 'H2CFINAL: H1_SO1 in ON  RKB/"k" basis, TMAT',-1)
           WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     & '*** Fermion corep ',I,' /',NFSYM
           CALL PRQMAT(TMAT(IMAT),NESH(I),NESH(I),
     &     NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF
C      ... do eliminate spinèorbit terms !
         IF (SPFR.AND..NOT.SPFR_BEG) THEN

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'H2CFINAL: H2c in RKB/|k> - before elimination '//
     &     ' of spin_orbit terms "AFTER"',-1)
           WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     &      '*** Fermion corep ',I,' /',NFSYM
            CALL PRQMAT(TBUF(IMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
         ENDIF

          CALL SPF_HK(TBUF(IMAT),.TRUE.,IBEIG(IORB(I)+1+NPSH(I)),
     &            IBEIG(IORB(I)+1+NPSH(I)),
     &           NESH(I),NESH(I),NESH(I),NESH(I),NZ,IPRINT)

          IF (IPRINT.GE.5) THEN
          CALL HEADER(
     & 'H2CFINAL: H2c (SPIN-FREE!) in the ON RKB/"k"basis, TBUF',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     & '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(TBUF(IMAT),NESH(I),NESH(I),
     &          NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

C        ... do add HSO1 into TBUF !!!
         CALL DAXPY(N2EQ,D1,TMAT(IMAT),1,TBUF(IMAT),1)

         IMAT = IMAT + N2EQ

        ENDIF
        ENDDO

      ENDIF


C  ... read transformation matrix 4c AO_LL -> 4c ON(ee)
      REWIND LUBSS
      IF (FNDLAB('TM2CAOFP',LUBSS)) THEN
         ISIZE = NFBAS(1,0)*NESH(1)*NZ
         IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NESH(2)*NZ
         CALL READT(LUBSS,ISIZE,TMAT)
         IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &   'H2CFINAL: TM2CAOFP read into TMAT.'
      ELSE
        CALL QUIT('H2CFINAL: TM2CAOFP not found !')
      END IF

C ....  do the transformation into VMAT
      IVMAT = 1
      IBUF = 1
      ITMAT = 1
      DO I=1,NFSYM
      IF (NESH(I).GT.0) THEN
C  ********************************************************************
C    ... if wished, eliminate SO terms from H2FP(IH2FP) "at the end"
C  ********************************************************************
        IF (SPFR.AND..NOT.SPFR_BEG.AND..NOT.SAVESO1) THEN

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'H2CFINAL: H2c in RKB/|k> - before elimination '//
     &     ' of spin_orbit terms "AFTER"',-1)
           WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     &      '*** Fermion corep ',I,' /',NFSYM
            CALL PRQMAT(TBUF(IVMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

          CALL SPF_HK(TBUF(IVMAT),.TRUE.,IBEIG(IORB(I)+1+NPSH(I)),
     &            IBEIG(IORB(I)+1+NPSH(I)),
     &           NESH(I),NESH(I),NESH(I),NESH(I),NZ,IPRINT)

          IF (IPRINT.GE.5) THEN
            CALL HEADER(
     &     'H2CFINAL: H2c in ON RKB/"k" - eliminated'//
     &     ' spin_orbit terms "AFTER"',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     &     '*** Fermion corep ',I,' /',NFSYM
           CALL PRQMAT(TBUF(IVMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

        ENDIF

        IF (IPRINT.GE.4) THEN
          CALL HEADER(
     & 'H2CFINAL: H2c in ON basis',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     & '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(TBUF(IVMAT),NESH(I),NESH(I),
     &          NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     & 'H2CFINAL: Transf.matrix 2c ON MO RKB/"k" -> 2c SA-AO, TMAT',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     & '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(TMAT(IBUF),NFBAS(I,0),NESH(I),
     &          NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C       Row dimension of H2AO: only L [NFBAS(I,1)] or L+S [NFBAS(I,0)]?
C       currently we use L+S, to use only L we should delete small Xhi
C       after MAKE_DKH /5-Aug-2002 hjaaj

        CALL QTRANS('MOAO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &               NESH(I),NESH(I),
     &            VMAT(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            TBUF(IVMAT),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &            TMAT(IBUF),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &            TMAT(IBUF),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)

        IF (IPRINT.GE.4) THEN
          CALL HEADER(
     & 'H2CFINAL: FINAL H2c in SA-AO basis "xhi_LL", VMAT',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     & '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(VMAT(1+I2BASX(I,I)),NFBAS(I,0),NFBAS(I,0),
     &          NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IVMAT = IVMAT + NESH(I)*NESH(I)*NZ
        IBUF = IBUF + NFBAS(I,0)*NESH(I)*NZ

      ENDIF
      ENDDO

C  ... write H2c SA-AO to the file...
      REWIND LUBSS
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
        BACKSPACE(LUBSS)
        CALL NEWLAB('H1AO_DK ',LUBSS,LUPRI)
        CALL WRITT(LUBSS,N2BBASXQ,VMAT)
        CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
        IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &  'H2CFINAL:  H1AO_DK written, EOFLABEL renewed.'
      ELSE
        CALL QUIT('H2CFINAL: EOFLABEL  not found !')
      END IF

C     CALL QUIT('H2CFINAL')
      CALL QEXIT('H2CFINAL')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck getpctmat */      
      SUBROUTINE GETPCTMAT(TBUF,TMAT,VMAT,FPB,
     &               SPFR,SPFR_BEG,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:  
C    ========  
C
C    On input/output - none, all quantities are read/written from/into BSSMAT.
C                                                         
C  ...  construct the picture change transformation matrixes depending on the mode
C ===================================================================================
C
C.  1. full: H4c AO ->  H2c AO
C
C   2.  sfB: H4c AO ->  H4c "k"  &   H4c"k" -> H2c AO
C
C   3.  sfE: H4c AO ->  H2c "k"  &   H2c"k" -> H2c AO
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00 )
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION TMAT(*),VMAT(*),TBUF(*),WORK(LWORK)
      LOGICAL   FNDLAB, FPB, SPFR,SPFR_BEG

      CALL QENTER('GETPCTMAT')
#include "memint.h"


      IF (.NOT.SPFR) THEN
C     ... do construct full picture change transf.matrix hd 4c AO -> h2c AO (LL)


       REWIND LUBSS
       IF (FNDLAB('U1_ONBAS',LUBSS)) THEN
          ISIZE = NORB(1)*NESH(1)*NZ
          IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NESH(2)*NZ
          CALL READT(LUBSS,ISIZE,TBUF)
        IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &  'GETPCTMAT: U1_ONBAS (4c ON->2c ON) read into TBUF.'
       ELSE
        CALL QUIT(
     &  'GETPCTMAT: U1_ONBAS not found in the BSSMAT file!')
       END IF

      REWIND LUBSS
      IF (FPB) THEN
       IF (FNDLAB('FPMAT4C ',LUBSS)) THEN
         ISIZE2 = NFBAS(1,0)*NORB(1)*NZ
         IF (NFSYM.EQ.2) ISIZE2 = ISIZE2 + NFBAS(2,0)*NORB(2)*NZ
C           ...  ISIZE2 is N2TMT !!!! ...
         CALL READT(LUBSS,ISIZE2,VMAT)
         IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &  'GETPCTMAT: FPMAT4C read into VMAT.'
       ELSE
        CALL QUIT(
     &  'GETPCTMAT: FPMAT4C not found in the BSSMAT file!')
       END IF

      ELSE
C        ... resticted kinetic balance only ...
       IF (FNDLAB('TMAT4C  ',LUBSS)) THEN
C        ISIZE2 = NFBAS(1,0)*NORB(1)*NZ
C        IF (NFSYM.EQ.2) ISIZE2 = ISIZE2 + NFBAS(2,0)*NORB(2)*NZ
C           ...  ISIZE2 is N2TMT !!!! ...
C        CALL READT(LUBSS,ISIZE2,VMAT)
         CALL READT(LUBSS,N2TMT,VMAT)
         IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &  'GETPCTMAT: TMAT4C read into VMAT.'
       ELSE
        CALL QUIT(
     &  'GETPCTMAT: TMAT4C is not found in the BSSMAT file!')
       END IF
      ENDIF

      IBUF = 1
      IVMAT = 1
      ITMAT = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN
        IF (IPRINT.GE.6) THEN
          CALL HEADER(
     & 'GETPCTMAT: Picture change transf.matrix '//
     & ' H4c in ON -> H2c in ON',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IBUF),NORB(I),NESH(I),
     &                    NORB(I),NESH(I),
     &                     NZ,IPQTOQ(1,0),LUPRI)
          CALL HEADER(
     &'GETPCTMAT: 4c ON basis (RKB/"k") in SA-AO basis "xhi"',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(IVMAT),NFBAS(I,0),NORB(I),
     &        NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C ...   construct picture change transf. matrix,  VMAT*TBUF = TMAT
      CALL QGEMM(NFBAS(I,0),NESH(I),NORB(I),D1, 
     &       'N','N',IPQTOQ(1,0),
     &        VMAT(IVMAT),NFBAS(I,0),NORB(I),NZ,
     &       'N','N',IPQTOQ(1,0),
     &        TBUF(IBUF),NORB(I),NESH(I),NZ,
     &        D0,IPQTOQ(1,0),
     &        TMAT(ITMAT),NFBAS(I,0),NESH(I),NZ)

        IF (IPRINT.GE.6) THEN
          CALL HEADER(
     & 'GETPCTMAT: Picture change transf.matrix '//
     & 'H4c SA-AO "xhi" -> H2c ON RKB/"k", TMAT',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NESH(I),
     &    NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)

        ENDIF

        IBUF = IBUF + NORB(I)*NESH(I)*NZ
        IVMAT = IVMAT + NFBAS(I,0)*NORB(I)*NZ
        ITMAT = ITMAT + NFBAS(I,0)*NESH(I)*NZ
 
      ENDIF
      ENDDO

      REWIND LUBSS
      IF (FNDLAB('TM2CAOFP',LUBSS)) THEN
         ISIZE2 = NFBAS(1,0)*NESH(1)*NZ
         IF (NFSYM.EQ.2) ISIZE2 = ISIZE2 + NFBAS(2,0)*NESH(2)*NZ
         CALL READT(LUBSS,ISIZE2,VMAT)
         IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &  'GETPCTMAT: TM2CAOFP read into VMAT.'
      ELSE
        CALL QUIT(
     &  'GETPCTMAT: TM2CAOFP not found in the BSSMAT file!')
      END IF

      IBUF = 1
      IVMAT = 1
      ITMAT = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN

        IF (IPRINT.GE.6) THEN
         CALL HEADER(
     &   'GETPCTMAT: Transf.matrix H2c ON(e) -> H2c "xhi_L" '//
     &   ' (TM2CAOFP), VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IVMAT),NFBAS(I,0),NESH(I),
     &         NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C  ...  construct full picture change transfprmatin matrix,  TMAT*VMAT+ = TBUF
      CALL QGEMM(NFBAS(I,0),NFBAS(I,1),NESH(I),D1, 
     &          'N','N',IPQTOQ(1,0),
     &           TMAT(ITMAT) ,NFBAS(I,0),NESH(I),NZ,
     &          'H','N',IPQTOQ(1,0),
     &           VMAT(IVMAT),NFBAS(I,0),NESH(I),NZ,
     &          D0,IPQTOQ(1,0),
     &          TBUF(IBUF),NFBAS(I,0),NFBAS(I,1),NZ)

        IF (IPRINT.GE.6) THEN
          CALL HEADER(
     & 'GETPCTMAT: Picture change transf.matrix '//
     & '  H4c SA-AO "xhi" -> H2c SA-AO "xhi"',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IBUF),NFBAS(I,0),NFBAS(I,1),
     &              NFBAS(I,0),NFBAS(I,1),
     &             NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IBUF =  IBUF + NFBAS(I,0)*NFBAS(I,1)*NZ
       ITMAT = ITMAT + NFBAS(I,0)*NESH(I)*NZ
       IVMAT = IVMAT + NFBAS(I,0)*NESH(I)*NZ

      ENDIF
      ENDDO

C .. do save final picture change transformation matrix into the file
      REWIND(LUBSS)
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
        BACKSPACE(LUBSS)
        CALL NEWLAB('U_PICTRM',LUBSS,LUPRI)
        ISIZE = NFBAS(1,0)*NFBAS(1,1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NFBAS(2,1)*NZ
        CALL WRITT(LUBSS,ISIZE,TBUF)
        CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
        WRITE(LUPRI,'(/2X,A)')
     &  'GETPCTMAT: Final picture change transformation matrix '//
     &  ' - U_PICTRM - is written to the BSSMAT file, EOFLABEL renewed.'
      ELSE
        CALL QUIT(
     & 'GETPCTMAT: Error, EOFLABEL not found in the BSSMAT file!')
      END IF

C ... prepare picture change matrixes for SPINFREE "BEFORE"
C   2.  sfB: H4c AO ->  H4c "k"  &   H4c"k" -> H2c AO
      ELSE IF (SPFR.AND.SPFR_BEG) THEN

C ... prepare picture change matrix for SPINFREE "AFTER"
      ELSE IF (SPFR.AND..NOT.SPFR_BEG) THEN

      ELSE
        CALL QUIT('GETPCTMAT: wrong end up within the "if" branch !')
      ENDIF


C     CALL QUIT('GETPCTMAT')
      CALL QEXIT('GETPCTMAT')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck cmp_eigval */      
       SUBROUTINE CMP_EIGVAL(VMAT,TMAT,TBUF,EIG,
     &                       IBEIG,SPFR,SPFR_BEG,ONLYSO1,SAVESO1,
     &                       IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:
C    ========
C  Do comparison of the eigenvalues of the Dirac and H2c Hamiltonians
C   trasnform each to the ON basis and do the diagonalization
C
C  On input: only arrays
C
C  On output: VMAT reorganized for H2c                                                    
C                                                         
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0 ) 
#include "dcbgen.h"  
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"
#include "dummy.h"
   
      DIMENSION VMAT(*),TMAT(*),
     &          EIG(*),IBEIG(*),TBUF(*),WORK(LWORK)    
      LOGICAL SPFR,SPFR_BEG,ONLYSO1,SAVESO1,FNDLAB,DOJACOBI
      CHARACTER*8 LABELH1

      CALL QENTER('CMP_EIGVAL')
#include "memint.h"         

C      ... this option does not matter....
C        DOJACOBI = .TRUE.
         DOJACOBI = .FALSE.

C   ... get AND reorganize the Lowdin transformation matrix VMAT
      REWIND LUBSS
      IF (FNDLAB('VMAT4C  ',LUBSS)) THEN
         CALL READT(LUBSS,N2BBASX,VMAT)
         IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &  'CMP_EIGVAL: VMAT4C read.'
      ELSE
        CALL QUIT(
     &  'CMP_EIGVAL: VMAT4C is not found in BSSMAT file!')
      END IF


      IF (IPRINT.GE.6) THEN
        CALL HEADER(
     &'CMP_EIGVAL: The total Lowdin matrix before BSS reordering',-1)
        IVMAT=1
        DO I=1, NFSYM
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(IVMAT),
     &          NFBAS(I,0),NFORB(I,0),NFBAS(I,0),NFORB(I,0),
     &          1,IPQTOQ(1,0),LUPRI)
          IVMAT = IVMAT + NFBAS(I,0)*NFORB(I,0)
        ENDDO
      ENDIF
      IF (NFSYM.EQ.2) THEN
C    HJJ:For TMAT: copy the VMAT for the L block ifsym=2 to
C        follow immediately after the L block ifsym=1.
C       MI: what a simple solution...
C          NDIM2 = NFBAS(2,0)*NORB(2)
           NDIM2 = NFBAS(2,0)*NESH(2)
           JV2   = 1 + NFBAS(1,0)*NFORB(1,0)
C          JT2   = 1 + NFBAS(1,0)*NORB(1)
           JT2   = 1 + NFBAS(1,0)*NESH(1)
           CALL DCOPY(NDIM2,VMAT(JV2),1,VMAT(JT2),1)
      ENDIF

C .... write the 2c-adapted Lowdin matrix into the file
      REWIND LUBSS
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
        BACKSPACE(LUBSS)
        CALL NEWLAB('VMAT2C  ',LUBSS,LUPRI)
        CALL WRITT(LUBSS,N2BBASX,VMAT)
        CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
        IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &  'CMP_EIGVAL: Adapted Lowdin matrix, VMAT2C, was written'//
     &  ' into the BSSMAT file, EOFLABEL was renewed.'
      ELSE
        CALL QUIT('CMP_EIGVAL: EOFLABEL not found !')
      END IF


C ... read H2c AO into TBUF
      REWIND LUBSS
      IF (FNDLAB('H1AO_DK ',LUBSS)) THEN
         CALL READT(LUBSS,N2BBASXQ,TBUF)
         IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &  'CMP_EIGVAL: H1AO_DK read.'
      ELSE
        CALL QUIT(
     &  'CMP_EIGVAL: H1AO_DK is not found in BSSMAT file!')
      END IF


C   ... transform H2c AO through Lowdin matrix into ON MO ee basis
C   ... transform H4c AO through  matrix into ON MO  basis

      IH2MO = 1
      IVMAT = 1
      IH1FP = 1
      DO I = 1, NFSYM
      IF (NORB(I).GT.0) THEN

       IF (IPRINT.GE.4) THEN

C    ... print out H2c matrix in AO "xhi" basis
        CALL HEADER(
     &'CMP_EIGVAL: Final H(2c)_LL in "xhi" SA-AO basis ' ,-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TBUF(1+I2BASX(I,I)),NFBAS(I,1),NFBAS(I,1),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)


C    ... print out Lowdin matrix
        CALL HEADER(
     &  'CMP_EIGVAL: Adapted Lowd. transf. matrix (Le)'//
     &  ' for the 2c SCF procedure:',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(VMAT(IVMAT),NFBAS(I,0),NESH(I),NFBAS(I,0),
     &              NESH(I),1,IPQTOQ(1,0),LUPRI)

       ENDIF

C   ... transform H2c AO through Lowdin matrix into ON MO ee basis
        CALL QTRANS('AOMO','S',D0,NFBAS(I,1),NFBAS(I,1),
     &               NESH(I),NESH(I),
     &    TBUF(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &    TMAT(IH2MO),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &    VMAT(IVMAT),NFBAS(I,0),NESH(I),1,IPQTOQ(1,0),
     &    VMAT(IVMAT),NFBAS(I,0),NESH(I),1,IPQTOQ(1,0),
     &        WORK(KFREE),LFREE,IPRINT)

       IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &   'CMP_EIGVAL: H2c in (MO_ee) ON Lowdin basis "Theta" to'//
     &   ' be diagonalized',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TMAT(IH2MO),NESH(I),NESH(I),
     &             NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

          IH2MO = IH2MO + NESH(I)*NESH(I)*NZ 
C       ... offset for Lowdin matrix
          IVMAT = IVMAT + NFBAS(I,0)*NESH(I)
          IH1FP = IH1FP + NORB(I)*NORB(I)*NZ

      ENDIF
      ENDDO

C ... WE WANT THE EIGENVECTORS AS WELL .... ADD ONE MORE CYCLE ONLY 
C      FOR THE DIAGONALIZATION...REUSE TBUF for eigenvectors

      IH2MO = 1
      IEIG = 1
      DO I = 1, NFSYM
      IF (NORB(I).GT.0) THEN

C      ... do the diagonalization of the H2c in MO(ee) basis !
       IF (DOJACOBI.AND.NZ.EQ.1) THEN
         write(lupri,*) 'CMP_EIGVAL: RSJACO used for diagonalization' 
         CALL RSJACO(NESH(I),NESH(I),NESH(I),TMAT(IH2MO),
     &                  EIG(IEIG),1,1,0,TBUF(IH2MO),
     &                  WORK(KFREE),LFREE)
       ELSE
        CALL QDIAG(NZ,NESH(I),TMAT(IH2MO),NESH(I),NESH(I), 
     &             EIG(IEIG),1,          
     &             TBUF(IH2MO),NESH(I),       
     &             NESH(I),WORK(KFREE),LFREE,IERR)
        IF (IERR.NE.0) CALL QUIT('CMP_EIGVAL: FAILED QDIAG 1 !...')
       ENDIF

       IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &   'CMP_EIGVAL: Eigenvectors of the diagonalized H2c'//
     &   ' in ON MO_ee Lowdin basis "Theta": ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IH2MO),NESH(I),NESH(I),
     &             NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

        IEIG  = IEIG  + NESH(I)
        IH2MO = IH2MO + NESH(I)*NESH(I)*NZ 
      ENDIF
      ENDDO

C .... read and diagonalize H1 in ON RKB/FP "k" basis
C ... write the TBUF - H1 "k" to the file...
      ISIZE=NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2) ISIZE=ISIZE+NORB(2)*NORB(2)*NZ
C ... get H1_ONB into TBUF
      REWIND LUBSS
      LABELH1 = 'H1_ONBAS'
      IF (SAVESO1) LABELH1 = 'H1_ONB  '
      IF (FNDLAB(LABELH1,LUBSS)) THEN
         CALL READT(LUBSS,ISIZE,TBUF)
         IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &  'CMP_EIGVAL: H1_ONBAS or H1_ONB read into TBUF.'
      ELSE
        CALL QUIT(
     &  'CMP_EIGVAL: H1_ONBAS or H1_ONB is not found in BSSMAT file!')
      END IF

       CALL HEADER(
     & 'CMP_EIGVAL: Electronic eigenv. of H(2c)/H(4c) in ON basis',-1)

      IH1FP = 1
      IH2MO = 1
      IEIG4  = 1
      IEIG2  = 1
      DO I = 1, NFSYM
      IF (NORB(I).GT.0) THEN
C        ... eliminate spin-free terms utilizing IBEIG(*)...
       IF (SPFR.AND..NOT.SPFR_BEG) THEN

         IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &    'CMP_EIGVAL: H4c in ON MO basis RKB/FP "k" before '//
     &    'elimination of SO terms (we have 2c spin-free "AFTER") ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IH1FP),NORB(I),NORB(I),
     &             NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

         CALL SPF_HK(TBUF(IH1FP),.TRUE.,IBEIG(IORB(I)+1),
     &    IBEIG(IORB(I)+1),
     &    NORB(I),NORB(I),NORB(I),NORB(I),NZ,IPRINT )

       ENDIF

        IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &    'CMP_EIGVAL: H4c in ON MO basis RKB/FP "k" before '//
     &    ' the diagonalization ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IH1FP),NORB(I),NORB(I),
     &             NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF


        CALL QDIAG(NZ,NORB(I),TBUF(IH1FP),NORB(I),NORB(I), 
     &             TMAT(IEIG4),0,          
     &             DUMMY,NORB(I),       
     &             NORB(I),WORK(KFREE),LFREE,IERR)
        IF (IERR.NE.0) CALL QUIT('CMP_EIGVAL: FAILED QDIAG 2!')

         WRITE(LUPRI,'(8X,A,I1,A,I1/)')
     &  '*** Fermion corep ',I,'/',NFSYM

C ... first positronic states first
        WRITE(LUPRI,'(3X,I5,F25.10)')
     &  ( J-NPSH(I)-1,TMAT(IEIG4+J-1),J=1,NPSH(I))

C ... then continue with electronic states
        DO J = 1, NESH(I)
          VAL2C = EIG(IEIG2+J-1)
          VAL4C = TMAT(IEIG4+NPSH(I)+J-1)
        DIFF =  VAL4C - VAL2C
        WRITE(LUPRI,'(3X,I5,2F25.10,3X,D12.5)') 
     &    J, VAL4C,VAL2C,DIFF
        ENDDO
        write(lupri,*)

        IEIG4  = IEIG4  + NORB(I)
        IEIG2  = IEIG2  + NESH(I)
        IH2MO = IH2MO + NESH(I)*NESH(I)*NZ 
        IH1FP = IH1FP + NORB(I)*NORB(I)*NZ

      ENDIF
      ENDDO
      write(lupri,*)

      CALL QEXIT('CMP_EIGVAL')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_hdk2_fpm */      
      SUBROUTINE GET_HDK2_FPM(TBUF,TMAT,VMAT,EIG,IBEIG,
     &                     DK2SO1,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:
C    ========
C
C  Do comparison of the eigenvalues of the Dirac and H2c Hamiltonians
C   trasnform each to the ON basis and do the diagonalization
C
C  On input:  TBUF H1 in "k" basis
C             EIG - free particle eigenvalues needed for denominators
C             IBEIG - boson irreps if spin-free, or DK2SO1
C                                                         
C  On output:  TBUF 
C             the desired DK2 Hamiltonian in "k" basis is saved into BSSMAT
C                                                         
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJAaJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0, D2 = 2.0D0, DM1 = -1.0D0 ) 
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dcbgen.h"   
#include "dgroup.h"
   
      DIMENSION TBUF(*),TMAT(*),VMAT(*),EIG(*),IBEIG(*),WORK(LWORK)
      LOGICAL SPFR_BEG, DK2SO1, FNDLAB

      CALL QENTER('GET_HDK2_FPM')
#include "memint.h"         

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
C     (taken from Barysz&Sadlej, JMS(Theochem) 573(2001) 181, eq. 46 )    
C

      IF (IPRINT.GE.2) THEN
        CALL HEADER(
     &'GET_HDK2_FPM: DKH2 Hamiltonian to be generated',-1)
      ENDIF

      CVAL2 = CVAL*CVAL
      JOFT1 = 1
      JBUF  = 1
      DO  I = 1,NFSYM  
      IF(NORB(I).GT.0) THEN

        N2T    = NORB(I)*NORB(I)
        N2ESHI = NESH(I)*NESH(I)
C
C         1) generate R_DK2(pe) in TBUF(pe) block:
C
         IOFT  = JOFT1 + NORB(I)*NPSH(I) - 1
         JR2PE = IOFT  + 1

         DO IE = NPSH(I) + 1,  NORB(I)
           DO IP = 1, NPSH(I)
C ... eig_p and eig_e are separated by 2c^2; use rather average
C ... eig(p) = -eig(e) - 2c^2
C         .... Calculate the average value
             EIG_P = 
     & ((-EIG(NORB(I)-IP+1+IORB(I))-(D2*CVAL2)) + 
     &   EIG(IP+IORB(I)))/D2

             EIG_E =
     & ( EIG(IE+IORB(I)) + 
     &   (-EIG(NORB(I)-IE+1+IORB(I))-(D2*CVAL2)) )/D2

             FAC = DM1/ (EIG_P - EIG_E)
             DO IZ = 0,NZ-1
               TBUF(IOFT+IP+IZ*N2T) =
     &         TBUF(IOFT+IP+IZ*N2T) * FAC
             END DO
           END DO
           IOFT = IOFT + NORB(I)
         END DO

          JH1EP = JOFT1 + NPSH(I)
 
C ----------------------------------------------------------------------------y
C  If desired, eliminate SO1 terms from both H1_ep and R_pe - DK2SO1 approach
C -----------------------------------------------------------------------------
          IF (DK2SO1) THEN
            IF (IPRINT.GE.2) THEN
              WRITE(LUPRI,'(3X,A)')
     & 'GET_HDK2_FPM: DKH2sfBSO1 option - elimination of spin-orbit'//
     & ' terms from H1(ep) and R(pe) blocks'
            ENDIF
           CALL SPF_HK( TBUF(JH1EP),.TRUE.,IBEIG(IORB(I)+NPSH(I)+1),
     &                 IBEIG(IORB(I)+1),
     &                 NESH(I),NPSH(I),NORB(I),NORB(I),NZ,IPRINT )
           CALL SPF_HK( TBUF(JR2PE),.TRUE.,IBEIG(IORB(I)+1),
     &                 IBEIG(IORB(I)+1+NPSH(I)),
     &                 NPSH(I),NESH(I),NORB(I),NORB(I),NZ,IPRINT )
          ENDIF

          IF (IPRINT.GE.6) THEN
           CALL HEADER(
     &    'GET_HDK2_FPM: 1. order wave function R_DK2(pe)'//
     &    ' in "|k>" basis',-1)
           WRITE(LUPRI,'(5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM

           CALL PRQMAT(TBUF(JR2PE),NPSH(I),NESH(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
           CALL HEADER(
     &'GET_HDK2_FPM: H(ep) part of H1 in "|k>" basis',-1)
           WRITE(LUPRI,'(5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TBUF(JH1EP),NESH(I),NPSH(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

C
C         2) Generate  "X" = H_1(ep) R_DK2(pe) into VMAT(JBUF)
C

          CALL QGEMM(NESH(I),NESH(I),NPSH(I),D1,
     &          'N','N',IPQTOQ(1,0),TBUF(JH1EP),NORB(I),NORB(I),NZ,
     &          'N','N',IPQTOQ(1,0),TBUF(JR2PE),NORB(I),NORB(I),NZ,
     &          D0,IPQTOQ(1,0),VMAT(JBUF), NESH(I),NESH(I),NZ)
C
C         3) Generate  X = 0.5 ( X + X(dagger)) 
C
          CALL QGETHE(VMAT(JBUF),NESH(I),NZ,NESH(I),NESH(I))

          IF (IPRINT.GE.7) THEN
            CALL HEADER(
     &'GET_HDK2_FPM: 2. order correction to H1 in "|k_e>" basis',-1)
           WRITE(LUPRI,'(5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
            CALL PRQMAT(VMAT(JBUF),NESH(I),NESH(I),NESH(I),
     &                  NESH(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
*  Add to first order (ee-blocks!): H_1(ee) =  H_1(ee) + X

          DO IZ = 1,NZ
            IBUF  = JBUF + (IZ-1)*N2ESHI
            JH1EE = JOFT1 + NPSH(I)*NORB(I) + NPSH(I)
     &            + (IZ-1)*NORB(I)*NORB(I)
            DO J = 1, NESH(I)
              CALL DAXPY(NESH(I),D1,TBUF(JH1EE),1,VMAT(IBUF),1)
              JH1EE = JH1EE + NORB(I)
              IBUF  = IBUF  + NESH(I)
            END DO
          END DO

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &    'GET_HDK2_FPM: Full H2c_DKH2(or DKH2sfBSO1) Hamiltonian'//
     &    ' in ON FP "k" basis, VMAT',-1)
          WRITE(LUPRI,'(5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(JBUF),NESH(I),NESH(I),NESH(I),
     &                  NESH(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
C         ... update the offsets...
          JBUF  = JBUF  + NESH(I)*NESH(I)*NZ
          JOFT1 = JOFT1 + N2T*NZ
      ENDIF
      ENDDO

C ... do save Hamiltonian matrix (VMAT) into the file BSSMAT
      REWIND(LUBSS)
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
        BACKSPACE(LUBSS)
        CALL NEWLAB('H2C_ONBA',LUBSS,LUPRI)
        ISIZE = NESH(1)*NESH(1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NESH(2)*NESH(2)*NZ
        CALL DCOPY(ISIZE,VMAT,1,TBUF,1)
        CALL WRITT(LUBSS,ISIZE,TBUF)
        CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
        WRITE(LUPRI,'(/2X,A)')
     &  'GET_HDK2_FPM: H2c_DKH2 Hamiltonian matrix (in "k" basis) is '//
     &  ' written to the BSSMAT file, EOFLABEL renewed.'
      ELSE
        CALL QUIT(
     & 'GET_HDK2_FPM: Error, EOFLABEL not found in the BSSMAT file!')
      END IF

      CALL QEXIT('GET_HDK2_FPM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_hdk1_fpm */      
      SUBROUTINE GET_HDK1_FPM(TBUF,VMAT,IBEIG,ONLYSO1,
     &                         IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: Do construct the simplest H1(ee) Hamiltonian in "k" basis
C    ========
C
C  On input:  TBUF -  H1 in "k" basis
C             IBEIG - boson irreps if pure spin-orbit integrals (ONLYSO1=.true.)
C             ONLYSO1 - flag for removing all scalar terms from Hamiltonian                                                         
C                                                         
C  On output: TBUF - H2c in "k" basis
C           ...  the desired DK1 Hamiltonian in "k" basis is written to the file BSSMAT
C                                                         
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER( DM1= -1.0D0 )
#include "dcbgen.h"  
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"   

      DIMENSION TBUF(*),VMAT(*),IBEIG(*),WORK(LWORK)
      LOGICAL ONLYSO1, FNDLAB

      CALL QENTER('GET_HDK1_FPM')
#include "memint.h"         
  
C  ... do copy the ee-part of the full H1 Hamiltonian into VMAT
      JBUF  = 1
      JOFT1 = 1
      DO I = 1, NFSYM
      IF (NESH(I).GT.0) THEN

        N2ESHI = NESH(I)*NESH(I)
        N2T    = NORB(I)*NORB(I)

        DO IZ = 1, NZ
         IBUF  = JBUF + (IZ-1)*N2ESHI
         JH1EE = JOFT1 + (NPSH(I)*NORB(I)) + NPSH(I)
     &          + (IZ-1)*N2T
         DO J = 1, NESH(I)
           CALL DCOPY(NESH(I),TBUF(JH1EE),1,VMAT(IBUF),1)
           JH1EE = JH1EE + NORB(I)
           IBUF  = IBUF  + NESH(I)
         ENDDO
        ENDDO

C      ... if desired, remove all scalar terms to have only SO1 integrals...
C      ======================================================================
        IF (ONLYSO1) THEN

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &  'GET_HDK1_FPM: Full H_DK1 matrix in ON FP "k" basis '//
     &  ' before eliminating spin-free terms (ONLYSO1=.true.)',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(VMAT(JBUF),NESH(I),NESH(I),
     &         NESH(I),NESH(I), NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

          CALL SPF_HK(VMAT(JBUF),.TRUE.,IBEIG(IORB(I)+1+NPSH(I)),
     &               IBEIG(IORB(I)+1+NPSH(I)),
     &              NESH(I),NESH(I),NESH(I),NESH(I),NZ,IPRINT)

        ENDIF

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &  'GET_HDK1_FPM: Final H_DK1-SO1 matrix in "k" basis ',-1)
           IF (ONLYSO1) WRITE(LUPRI,'(3X,A)')
     &  '...spin-free terms are deleted!'
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(VMAT(JBUF),NESH(I),NESH(I),
     &         NESH(I),NESH(I), NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

C       ... update the offsets...
        JBUF  = JBUF  + N2ESHI*NZ
        JOFT1 = JOFT1 + N2T*NZ
      ENDIF 
      ENDDO

C ... do write the H1ee Hamiltonian into the file...
      REWIND(LUBSS)
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
        BACKSPACE(LUBSS)
        CALL NEWLAB('H2C_ONBA',LUBSS,LUPRI)
        ISIZE = NESH(1)*NESH(1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NESH(2)*NESH(2)*NZ
        CALL DCOPY(ISIZE,VMAT,1,TBUF,1)
        CALL WRITT(LUBSS,ISIZE,TBUF)
        CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
        WRITE(LUPRI,'(/2X,A)')
     &  'GET_HDK1_FPM: H2c_DKH1 Hamiltonian matrix (in "k" basis) is '//
     &  ' written to the BSSMAT file, EOFLABEL renewed.'
      ELSE
        CALL QUIT(
     & 'GET_HDK2_FPM: Error, EOFLABEL not found in the BSSMAT file!')
      END IF

      CALL QEXIT('GET_HDK1_FPM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck amfiin */
      SUBROUTINE AMFIIN(H2SAO,NB,NTNUC,IPTYPAO,INDXLM,IPNUC,IZTNUCS,
     &                  IPTA,INDXLMA,IAMFU,IANUCS,ISTNUCTYP,CHPNUC,
     &                   NSPH,LABINT,ONLYSO1,ONLYTRA,
     &             IMFCHA,IC_AMFI,IQPOS,CVAL_AMFI,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
C  Interface to B.Schimmelpfennig's atomic mean-field integral code, AMFI   
C
C  On input:  H2SAO
C             NB - total number of basis functions(L+S)
C             NTNUC - total number (all) of nucleis 
C              
C          ... IZTNUCS - array for charges (int) of all nuclei of the system
C             
C          ... for each sperical AO in H2CAO matrix: 
C             IPTYPAO - typ of spherical AO
C             INDXLM - index of typ
C             IPNUC - index nuclei (1-NTNUC)
C             ISTNUCTYP - integer array for storing ICENT, INDXLM 
C 
C 
C          ... for each sperical AO in AMFI triangular matrixes
C             IPTA - typ of spherical AO for AMFI
C             INDXLMA - index of AMFI type
C             IAMFU - mapping function AMFI ints -> H2SAO positions
C
C             LABINT - spherical AO labels
C             ONLYSO1 - is .true. do only comparison of SO1 BSS/AMFI integrals 
C
C          ... for AMFI program (subroutine)
C              IMFCHA
C              IC_AMFI
C              CVAL (from dcbgen.h) - actual speed of light
C
C
C On output: H2SAO with mean-field contributions (only if wished)
C            - only comparison made
C
C     Written by Miro Ilias, Strasbourg, 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER(THRNULL=1.0D-7,D0=0.00D00)
C#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dummy.h"
#include "nuclei.h"
#include "shells.h"
C#include "dcbsph.h"
#include "chrnos.h"

C
      LOGICAL BREIT,PROPRI,OLDNRM,ISONLIST,ONLYSO1,ONLYTRA,EXAMF
      DIMENSION WORK(LWORK),H2SAO(NB,NB,4),IQPOS(1:8)  
      DIMENSION IPTYPAO(NSPH),INDXLM(NSPH),IPNUC(NSPH),
     &       IPTA(*),INDXLMA(*),IAMFU(NSPH),ISTNUCTYP(NSPH,2),
     &       IZTNUCS(NTNUC),IANUCS(NTNUC),CHPNUC(NTNUC)
      DIMENSION ncontrac(0:4)

C     CHARACTER LABINT(MXQN**2)*4
      CHARACTER*(*) LABINT(*)
      character*8 xa,ya,za
      DIMENSION xa(4),ya(4),za(5)
      CHARACTER*16 AFN

      CALL QENTER('AMFIIN')
#include "memint.h"

      CALL HEADER('*** Output from AMFIIN ***',-1)

      write(lupri,*) 'NUCIND=',NUCIND
      write(lupri,*) 'NLRGSH=',NLRGSH
      write(lupri,*) 'charge:',(charge(i),i=1,nucind)
      write(lupri,*) 'nuc_deg:',(nucdeg(i),i=1,nucind)
      write(lupri,*) 'Total number of centers:',NTNUC


C   ... fill IZTNUCS - LIST OF ALL NUCLEIS
      WRITE(LUPRI,'(2X,A)') 'Total list of nucleis:'
      IX = 1
      DO I=1,NUCIND
      DO J=1,NUCDEG(I)
       IZTNUCS(IX) = NINT(CHARGE(I))
       write(lupri,*) IX,'. nuclei of charge:',IZTNUCS(IX) 
       IX = IX + 1
      ENDDO
      ENDDO

C     CALL QUIT('amfiin')

C    .. make ready the IPTYPAO, INDXLM,IPNUC
      IF (IPRINT.GE.1) THEN
       WRITE(LUPRI,'(/2X,A)') 'List of shells of the system:' 
      ENDIF
      IOFF = 1
      JOFF = 1
      IX = 1
      DO ISHELL = 1, NLRGSH
        L     = NHKT(ISHELL)
        ICENT = NCENT(ISHELL)
        NDEG  = NUCDEG(ICENT)
        INDX=NUMCF(ISHELL)
        DO IDEG = 1, NDEG
C        ... find out the nuclei
         INUCPOS = 0
         DO IC=1,ICENT
           IF (IC.LT.ICENT) THEN
             IDEG1=NUCDEG(IC)
           ELSE
             IDEG1=IDEG
           ENDIF
           DO IDG = 1, IDEG1
             INUCPOS = INUCPOS + 1
           ENDDO
         ENDDO

         IF (IPRINT.GE.1) THEN
          write(lupri,*) 'SHELL L=',NHKT(ISHELL),
     &  ' ICENT=',NCENT(ISHELL),' charge:',CHARGE(ICENT),
     &  ' IDEG/NDEG=',IDEG,'/',NUCDEG(ICENT),' centre=',INUCPOS,
     &  ' INDXLM=',INDX
         ENDIF
        
C      ... specify IPTYPAO,  INDX_L_ML and center
        IF (L.EQ.1) THEN
          IPTYPAO(IX)=1
C          write(lupri,*) 'IX=',IX,' >',WORK(KIPTYPAO+IX)
          INDXLM(IX)=INDX
          IPNUC(IX)=INUCPOS
          IX = IX + 1 
C          write(lupri,*) 'IX set to:',IX
        ELSE IF (L.EQ.2) THEN
          IPTYPAO(IX)=2 ! 2px
C          write(lupri,*) 'IX=',IX,' >',WORK(KIPTYPAO+IX)
          INDXLM(IX)=INDX
          IPNUC(IX)=INUCPOS
          IPTYPAO(IX+1)=3  ! 2py
C          write(lupri,*) 'IX=',IX+1,' >',WORK(KIPTYPAO+IX+1)
          INDXLM(IX+1)=INDX
          IPNUC (IX+1)=INUCPOS
          IPTYPAO(IX+2)=4  ! 2pz
C          write(lupri,*) 'IX=',IX+2,' >',WORK(KIPTYPAO+IX+2)
          INDXLM(IX+2)=INDX
          IPNUC (IX+2)=INUCPOS
          IX = IX + 3
C          write(lupri,*) 'IX set to:',IX
        ELSE IF (L.EQ.3) THEN
          IPTYPAO(IX  )=5 ! 3d2-
          INDXLM(IX+0)=INDX
          IPNUC(IX)=INUCPOS

          IPTYPAO(IX+1)=6 ! 3d1-
          INDXLM(IX+1)=INDX
          IPNUC (IX+1)=INUCPOS

          IPTYPAO(IX+2)=7 ! 3d0
          INDXLM(IX+2)=INDX
          IPNUC (IX+2)=INUCPOS

          IPTYPAO(IX+3)=8 ! 3d1+
          INDXLM(IX+3)=INDX
          IPNUC (IX+3)=INUCPOS

          IPTYPAO(IX+4)=9 ! 3d2+
          INDXLM(IX+4)=INDX
          IPNUC (IX+4)=INUCPOS

          IX = IX + 5
        ELSE IF (L.EQ.4) THEN

          IPTYPAO(IX  )=10 ! 4f3-
          INDXLM(IX+0)=INDX
          IPNUC (IX+0)=INUCPOS

          IPTYPAO(IX+1)=11 ! 4f2-
          INDXLM(IX+1)=INDX
          IPNUC (IX+1)=INUCPOS

          IPTYPAO(IX+2)=12 ! 4f1-
          INDXLM(IX+2)=INDX
          IPNUC (IX+2)=INUCPOS

          IPTYPAO(IX+3)=13 ! 4f0 
          INDXLM(IX+3)=INDX
          IPNUC (IX+3)=INUCPOS

          IPTYPAO(IX+4)=14 ! 4f1+
          INDXLM(IX+4)=INDX
          IPNUC (IX+4)=INUCPOS

          IPTYPAO(IX+5)=15 ! 4f2+
          INDXLM(IX+5)=INDX
          IPNUC (IX+5)=INUCPOS

          IPTYPAO(IX+6)=16 ! 4f3+
          INDXLM(IX+6)=INDX
          IPNUC (IX+6)=INUCPOS

          IX = IX + 7
        ELSE IF (L.EQ.5) THEN

          IPTYPAO(IX+0)=17 ! 5g4-
          INDXLM(IX+0)=INDX
          IPNUC (IX+0)=INUCPOS

          IPTYPAO(IX+1)=18 ! 5g3-
          INDXLM(IX+1)=INDX
          IPNUC (IX+1)=INUCPOS

          IPTYPAO(IX+2)=19 ! 5g2-
          INDXLM(IX+2)=INDX
          IPNUC (IX+2)=INUCPOS

          IPTYPAO(IX+3)=20 ! 5g1-
          INDXLM(IX+3)=INDX
          IPNUC (IX+3)=INUCPOS

          IPTYPAO(IX+4)=21 ! 5g0
          INDXLM(IX+4)=INDX
          IPNUC (IX+4)=INUCPOS

          IPTYPAO(IX+5)=22 ! 5g1+
          INDXLM(IX+5)=INDX
          IPNUC (IX+5)=INUCPOS

          IPTYPAO(IX+6)=23 ! 5g2+
          INDXLM(IX+6)=INDX
          IPNUC (IX+6)=INUCPOS

          IPTYPAO(IX+7)=24 ! 5g3+
          INDXLM(IX+7)=INDX
          IPNUC (IX+7)=INUCPOS

          IPTYPAO(IX+8)=25 ! 5g4+
          INDXLM(IX+8)=INDX
          IPNUC (IX+8)=INUCPOS

          IX = IX + 9
        ELSE
          CALL QUIT('AMFI MAX UP TO g FUNCTIONS!')
        ENDIF
      ENDDO
      ENDDO

      IF ((IX-1).NE.NSPH) THEN
        write(lupri,*) 'AMFIIN error: IX-1=',IX-1, 'NSPH=',NSPH
        CALL QUIT('AMFIIN: IX-1<>NSPH ! CHECK IT !')
      ENDIF

C ... rearange INDXLM for splitted basis sets...
      CALL UNITE_INDXS(NSPH,LABINT,IPNUC,IPTYPAO,
     &                   INDXLM,ISTNUCTYP,IPRINT)

C  ... cancel control print out as IQDEF is def. in dgroup.h
C     IF (IPRINT.GE.6) THEN
C       CALL HEADER(
C    &'AMFIIN: TOTAL H(2c)_LL in "xhi" S AO basis ' ,-1)
C       CALL PRQMAT(H2SAO,NB,NB,NB,NB,4,IQDEF,LUPRI)    
C     ENDIF

C     Open Mean-field input file generated in READIN      
      LUTMP2 = 0
      MNFPRP = 0
      CALL GPOPEN(LUTMP2,'MNF.INP','OLD',' ','FORMATTED',
     &            IDUMMY,.FALSE.)

      OLDNRM=.FALSE.
      BREIT = .FALSE.

C  ++++ Find out individual charges
      IF (IMFCHA.NE.0) THEN

C     ... find out degeneracy of that nuclei IATOM with charge ICHARGE
        IMINCHA=IZTNUCS(1)
        INDXMINCH=1
        ISUMCH=0
        DO I = 1, NTNUC
C        .... NOTE: The minimal charge MUST NOT be zero !
         IF ( (IMINCHA.GT.IZTNUCS(I).AND.IZTNUCS(I).NE.0)
     &        .OR.(IMINCHA.EQ.0)) THEN
             IMINCHA=IZTNUCS(I)
             INDXMINCH = I
         ENDIF
         ISUMCH = ISUMCH + IZTNUCS(I)
        ENDDO

        IF (ISUMCH.GE.ABS(IMFCHA)) THEN
           WRITE(LUPRI,*) 'Total sum of nucl.charges:',ISUMCH,
     &   ' IMFCHA is',IMFCHA,
     &  ' WARNING, charge too high!'
        ENDIF

C     .... find out the individual charges on atoms
          WRITE(LUPRI,*)
     &     'The total nonzero charge of the system:',IMFCHA
           CHSUM = D0
          DO I=1,NTNUC
            CHPNUC(I)=DFLOAT(IZTNUCS(I))/DFLOAT(IZTNUCS(INDXMINCH))
            CHSUM=CHSUM+CHPNUC(I)
C           write(lupri,*) I,' -',CHPNUC(I)
          ENDDO
          FACCHA=DFLOAT(IMFCHA)/CHSUM
          write(lupri,*) ' factor is :',FACCHA
          CHSUM = D0
          DO I=1,NTNUC
            CHPNUC(I)=FACCHA*CHPNUC(I)
            CHSUM = CHSUM + CHPNUC(I)
            write(lupri,*) I,'.atom-nucleus charge:',IZTNUCS(I),
     &      '  partial charge:',CHPNUC(I)
          ENDDO
          write(lupri,*) '     Sum of all charges (real):',CHSUM
          write(lupri,*) 'Total charge of the system is :',IMFCHA
          IF (NINT(CHSUM).NE.IMFCHA) CALL 
     &      QUIT('AMFIIN: error in partial charges !')
      ENDIF

C ===================================
C     We loop over unique centers
C ===================================

      REWIND (LUTMP2)
C     ... skip text
      READ(LUTMP2,*) 
C     ... read number of unique nuclei
      READ(LUTMP2,*) NUNUC
      write(lupri,*) 'MNF.INP - read number of unique nuclei:',NUNUC

      TOTDIFX=D0
      TOTDIFY=D0
      TOTDIFZ=D0
      ITD=0

      DO IATOM = 1, NUNUC

       READ(LUTMP2,*) ICHARGE

       WRITE(LUPRI,'(2X,A,I2,A,I5)')
     &   '>>> AMFI is about to calculate ',
     &   IATOM,'. unique nucleus of charge ',ICHARGE

C     ... file name for AMFI integrals of center
       IF (ICHARGE.LT.10) THEN
        AFN='AOPROPER_MNF.'//CHRNOS(ICHARGE)
       ELSE IF (ICHARGE.LT.100) THEN
        AFN='AOPROPER_MNF.'//CHRNOS(ICHARGE/10)//CHRNOS(MOD(ICHARGE,10))
       ELSE IF (ICHARGE.LT.1000) THEN
        AFN='AOPROPER_MNF.'//CHRNOS(ICHARGE/100)//
     &                       CHRNOS(MOD(ICHARGE,100)/10)//
     &                       CHRNOS(MOD( MOD(ICHARGE,100), 10) )
       ELSE
        CALL QUIT('AMFIIN: nucl.charge > 1000! This is not supported!')
       ENDIF

       WRITE(LUPRI,'(2X,A,A)')
     &  'AMFIIN: The file name with AMFI integrals: ',AFN


C      ... decide whether to calculate or not the AOPROPER_MNF.XY file
       INQUIRE(FILE=AFN,EXIST=EXAMF)

        IF (.NOT.EXAMF) THEN
          WRITE(LUPRI,'(4X,A)')
     &  '...this file does not exist, AMFI is utilized ! '
          CALL GPOPEN(MNFPRP,AFN,'UNKNOWN',' ','UNFORMATTED',
     &             IDUMMY,.FALSE.)
        ELSE
          WRITE(LUPRI,'(4X,A)')
     &  '...this file does exist, AMFI is not utilized ! '
          CALL GPOPEN(MNFPRP,AFN,'OLD',' ','UNFORMATTED',
     &             IDUMMY,.FALSE.)
        ENDIF


C      ... find out the correct index of nuclear exponent
       IF (GAUNUC) THEN
        IATOM1 = 0
        DO J = 1, NUCIND
          IF (NINT(CHARGE(J)).EQ.ICHARGE) IATOM1=J
        ENDDO
        IF (IATOM1.EQ.0) THEN
         CALL QUIT('AMFIIN: IATOM1 not found')
        ELSE
         write(lupri,*) ' indx=',IATOM1,' CHARGE=',
     &        CHARGE(IATOM1),ICHARGE,
     &   ' nucl. exponent:',GNUEXP(IATOM1)
        ENDIF

C     ... gauss.nucleus (for SO1) does not work in AMFI !!!!
C     ... Bernd Schim., help us !
        IF (IC_AMFI.EQ.1.OR.IC_AMFI.EQ.4.OR.IC_AMFI.EQ.5) THEN
          CALL
     &  QUIT('AMFIIN: GNUC is not working in AMFI!')
        ENDIF
       ENDIF

C     ... find out degeneracy of that nuclei IATOM with charge ICHARGE
        NICHARGE = 0
        INDXMINCH=0
        DO I = 1, NTNUC
         IF (IZTNUCS(I).EQ.ICHARGE) THEN
          NICHARGE=NICHARGE+1
          INDXMINCH=I
         ENDIF
        ENDDO
        IF (NICHARGE.EQ.0) THEN
         CALL QUIT('AMFIIN: The nucleus does not exist is system!')
        ENDIF  
        write(lupri,*) 'The AMFI unique nuclei Z=',ICHARGE,
     &              ' index=',INDXMINCH,
     &              ' is degenerate: ',NICHARGE
        write(lupri,*) ' partial charge is (real):',CHPNUC(INDXMINCH)
        IACHNUC=NINT(CHPNUC(INDXMINCH))
        write(lupri,*) ' partial charge is (round int):',IACHNUC

C     .... call the AMFI program of B.Schimmelpfennig 
C     ... if file with integrals is present, return back without calcul.
         CALL AMFI(LUTMP2,MNFPRP,OLDNRM,BREIT,GAUNUC,GNUEXP(IATOM1),
     &           CVAL_AMFI,EXAMF,IC_AMFI,IACHNUC,WORK(KFREE),LFREE)

C      ... read AOPROPER.MNF file...
C    numballcart... total no of functions per centre...
      REWIND (MNFPRP)

      READ(MNFPRP) xa
      READ(MNFPRP) numballcart

C     write(lupri,*) 'numballcart=',numballcart

C  ... allocate memmory for this AMFI center
      length3=((numballcart*numballcart)+numballcart)/2
      CALL MEMGET('REAL',KSOX,length3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSOY,length3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSOZ,length3,WORK,KFREE,LFREE)

      CALL MEMCHK('after AMFI 2',WORK,1)

      WRITE(LUPRI,'(/,2X,A)')
     & 'AMFI order of (centre specific) basis function: '
      DO I=1,numballcart
        READ(MNFPRP) IPTA(I),INDXLMA(I)
        write(lupri,*) ' amfi fu no:',I,' type(indx)',LABINT(IPTA(I)),
     &                '(',INDXLMA(I),')'
      ENDDO

      READ(MNFPRP)  Lhigh,(ncontrac(I),I=0,Lhigh) 

       write(lupri,*) 'Lhigh=',Lhigh,' ncontrac:',
     &     (ncontrac(I),I=0,Lhigh)
       iscontr=0
       DO I=0,Lhigh
        iscontr=iscontr+(ncontrac(I)*(I+I+1))
       ENDDO
       IF (iscontr.NE.numballcart) THEN
         write(lupri,*) 'iscontr,numballcart:',iscontr,numballcart
         CALL QUIT('AMFIIN: iscontr <> numballcart!!')
       ENDIF

C        ... find NICHARGE centers on IZTNUCS of H2CAO with ICHARGE !
C          use IANUCS(1-NTNUC)
        CALL IZERO(IANUCS,NTNUC)
        IX = 0
        DO I=1,NTNUC
         IF (ICHARGE.EQ.IZTNUCS(I)) THEN
           IX = IX + 1
           IANUCS(IX)=I
         ENDIF
        ENDDO
        IF (IX.NE.NICHARGE) THEN
          CALL QUIT('AMFIIN: nuclei do not correspond !')
        ENDIF

        WRITE(LUPRI,'(/2X,A,$)') 'AMFIIN: corresponding H2CAO nuclei:'
        WRITE(LUPRI,*) (IANUCS(I),I=1,NICHARGE)


C       ... prepare the mapping array AMFI-H2C basis functions...
        CALL IZERO(IAMFU,NSPH)
        DO I=1, NICHARGE
        DO J=1, numballcart
           K=0
 23        K=K+1
            ISONLIST=.FALSE.
            ISONLIST=(IPTA(J).EQ.IPTYPAO(K)).AND.
     &             (INDXLMA(J).EQ.INDXLM(K)).AND.
     &              (IANUCS(I).EQ.IPNUC(K))
C           write(lupri,*) K,IPTA(J),IPTYPAO(K),INDXLMA(J),
C    &                    INDXLM(K),IANUCS(I),IPNUC(K),ISONLIST
C          IF (ISONLIST) THEN
C           write(lupri,*) 'MATCH! K=',K,' I=',I,' J=',J
C          ENDIF
          IF (K.LT.NSPH.AND..NOT.ISONLIST) GOTO 23 
          IF (K.EQ.NSPH.AND..NOT.ISONLIST) THEN
           CALL QUIT('AMFIIN: isonlist 1')
          ELSE IF (ISONLIST.AND.K.LE.NSPH) THEN
            JI=J+(numballcart*(I-1))
            IAMFU(JI) = K 
C           write(lupri,*) 'stored k=',K,' at pos:',JI
          ELSE
           CALL QUIT('AMFIIN: isonlist 2')
          ENDIF

        ENDDO
        ENDDO

C   ... write out mapping function correspondance....
      IF (IPRINT.GE.3) THEN
        write(lupri,*) '*** Mapping AMFI - H2CAO *** '
       IX = 0
       DO I=1,NICHARGE
       DO J=1,numballcart
         IX = IX + 1
        write(lupri,*) 'amfi nuc, charge:',IANUCS(I),ICHARGE,
     &    ' amfi fu:',LABINT(IPTA(J)),'(',
     &                INDXLMA(J),')',
     &   '  H2CAO: pos= ',IAMFU(IX), ' nucl. charge:',
     &      IPNUC(IAMFU(IX)),IZTNUCS(IPNUC(IAMFU(IX))),
     &   ' fu:',LABINT(IPTYPAO(IAMFU(IX))),
     &      '(',INDXLM(IAMFU(IX)),')'
        IF (ICHARGE.NE.IZTNUCS(IPNUC(IAMFU(IX)))) THEN
         CALL QUIT('AMFIIN: Different charges !')
        ENDIF
       ENDDO
       ENDDO
        write(lupri,*)
      ENDIF

C     ... read AMFI integrals...
C     READ(MNFPRP) xa
      READ(MNFPRP) (WORK(KSOX+irun),irun=0,length3-1)
C      write(lupri,*) ' AMFI X integrals:'
C      write(lupri,*) (WORK(KSOX+irun),irun=0,length3-1) 
      READ(MNFPRP) ya
      READ(MNFPRP) (WORK(KSOY+irun),irun=0,length3-1)
C      write(lupri,*) ' AMFI Y integrals:'
C      write(lupri,*) (WORK(KSOY+irun),irun=0,length3-1) 
C     READ(MNFPRP) za
      READ(MNFPRP) (WORK(KSOZ+irun),irun=0,length3-1)
C      write(lupri,*) ' AMFI Z integrals:'
C      write(lupri,*) (WORK(KSOZ+irun),irun=0,length3-1) 

      CALL MEMCHK('after AMFI2',WORK,1)


C ************************************************************************
C     ... do comparison of SO1 integrals if wished (ONLYSO1=.true.)
C       or add AMFI contributions into H2CAO
C ************************************************************************

      IF (ONLYSO1.AND.ICHARGE.GT.0) THEN
       WRITE(LUPRI,'(/2X,A,I5/)') 
     & 'AMFIIN:  Doing comparison of SO1 integrals for'//
     & '  nucleus of charge ',ICHARGE
      ENDIF

      IXPOS=IQPOS(4)
      IYPOS=IQPOS(3)
      IZPOS=IQPOS(2)

      write(lupri,*) '...quaternions x,y,z:',IQPOS(4),IQPOS(3),IQPOS(2)
      IF (IXPOS.EQ.IYPOS.OR.IXPOS.EQ.IZPOS.OR.IYPOS.EQ.IZPOS.OR.
     & (IXPOS+IYPOS+IZPOS).NE.9) THEN
        CALL QUIT('AMFIIN: Wrong x,y,z quaternions !')
      ENDIF   

      DO IANUC=1,NICHARGE
      IX = 0
      DO  I = 1, numballcart
      DO  J = 1, I

C       ... get SO integrals from AMFI  
         SOXINT=-WORK(KSOX+IX)
         SOYINT=-WORK(KSOY+IX)
         SOZINT=-WORK(KSOZ+IX)

C      ... indexes of H2CAO(I1,J1,2-4)
       I1 = IAMFU(I+(numballcart*(IANUC-1)))
       J1 = IAMFU(J+(numballcart*(IANUC-1)))

       IF (IPRINT.GE.4) THEN
       write(lupri,*) 'amfi: <',LABINT(IPTA(I)),'(',INDXLMA(I),')',
     &  '[',IANUCS(IANUC),']','|O|',LABINT(IPTA(J)),'(',INDXLMA(J),')',
     &   '[',IANUCS(IANUC),']> H2CAO: <',LABINT(IPTYPAO(I1)),
     &     '(',INDXLM(I1),
     &    ')','[',IPNUC(I1),']','|O|',LABINT(IPTYPAO(J1)),
     &   '(',INDXLM(J1),
     &   ')[',IPNUC(J1),']>'
       ENDIF


C    ... do only comparison of between AMFI SO1 - BSS SO1 integrals
C   ... only for nucleis with charge > 0 !!!
      IF (ONLYSO1.AND.ICHARGE.GT.0) THEN
C     IF (ONLYSO1) THEN

C      DIFSOX = SOXINT-H2SAO(I1,J1,4)
       DIFSOX = SOXINT-H2SAO(I1,J1,IXPOS)
C      DIFSOY = SOYINT-H2SAO(I1,J1,3)
       DIFSOY = SOYINT-H2SAO(I1,J1,IYPOS)
C      DIFSOZ = SOZINT-H2SAO(I1,J1,2)
       DIFSOZ = SOZINT-H2SAO(I1,J1,IZPOS)
     
       TOTDIFX=TOTDIFX + DABS(DIFSOX)
       TOTDIFY=TOTDIFY + DABS(DIFSOY)
       TOTDIFZ=TOTDIFZ + DABS(DIFSOZ)

       ITD = ITD + 1

       IF (DABS(DIFSOX).GE.THRNULL) THEN
           ratio=SOXINT/H2SAO(I1,J1,4)
           write(lupri,*) '<',LABINT(IPTYPAO(I1)),'(',
     &        INDXLM(I1),')','|SOX|',
     &         LABINT(IPTYPAO(J1)),'(',INDXLM(J1),')','>'
           write(lupri,*) SOXINT,H2SAO(I1,J1,4),DIFSOX,ratio
         ENDIF

         IF (DABS(DIFSOY).GE.THRNULL) THEN
           ratio=SOYINT/H2SAO(I1,J1,3)
           write(lupri,*) '<',LABINT(IPTYPAO(I1)),'(',
     &        INDXLM(I1),')','|SOY|',
     &         LABINT(IPTYPAO(J1)),'(',INDXLM(J1),')','>'
           write(lupri,*) SOYINT,H2SAO(I1,J1,3),DIFSOY,ratio
         ENDIF

         IF (DABS(DIFSOZ).GE.THRNULL) THEN
           ratio=SOZINT/H2SAO(I1,J1,2)
           write(lupri,*) '<',LABINT(IPTYPAO(I1)),'(',
     &        INDXLM(I1),')','|SOZ|',
     &         LABINT(IPTYPAO(J1)),'(',INDXLM(J1),')','>'
           write(lupri,*) SOZINT,H2SAO(I1,J1,2),DIFSOZ,ratio
         ENDIF

         IF (I1.NE.J1) THEN

C           DIFSOX = SOXINT+H2SAO(J1,I1,4)
            DIFSOX = SOXINT+H2SAO(J1,I1,IXPOS)
C           DIFSOY = SOYINT+H2SAO(J1,I1,3)
            DIFSOY = SOYINT+H2SAO(J1,I1,IYPOS)
C           DIFSOZ = SOZINT+H2SAO(J1,I1,2)
            DIFSOZ = SOZINT+H2SAO(J1,I1,IZPOS)

            TOTDIFX=TOTDIFX + DABS(DIFSOX)
            TOTDIFY=TOTDIFY + DABS(DIFSOY)
            TOTDIFZ=TOTDIFZ + DABS(DIFSOZ)

            ITD = ITD + 1

C            write(lupri,*) ITD,'<<<',TOTDIFX,TOTDIFY,TOTDIFZ

          IF (DABS(DIFSOX).GE.THRNULL) THEN
           ratio=SOXINT/H2SAO(J1,I1,4)
           write(lupri,*) '<',LABINT(IPTYPAO(J1)),'(',
     &        INDXLM(J1),')','|SOX|',
     &         LABINT(IPTYPAO(I1)),'(',INDXLM(I1),')','>'
           write(lupri,*) SOXINT,-H2SAO(J1,I1,4),DIFSOX,ratio
          ENDIF

          IF (DABS(DIFSOY).GE.THRNULL) THEN
           ratio=SOYINT/H2SAO(J1,I1,3)
           write(lupri,*) '<',LABINT(IPTYPAO(J1)),'(',
     &        INDXLM(J1),')','|SOY|',
     &         LABINT(IPTYPAO(I1)),'(',INDXLM(I1),')','>'
           write(lupri,*) SOYINT,-H2SAO(J1,I1,3),DIFSOY,ratio
           ENDIF

           IF (DABS(DIFSOZ).GE.THRNULL) THEN
            ratio=SOZINT/H2SAO(J1,I1,2)
           write(lupri,*) '<',LABINT(IPTYPAO(J1)),'(',
     &        INDXLM(J1),')','|SOZ|',
     &         LABINT(IPTYPAO(I1)),'(',INDXLM(I1),')','>'
           write(lupri,*) SOZINT,-H2SAO(J1,I1,2),DIFSOZ,ratio
           ENDIF

         ENDIF
       ELSE IF (.NOT.ONLYTRA) THEN
C  ==================================================
C      ... add AMFI contribution to the H2CAO
C        when onlytra=.false.
C  ==================================================
C       H2SAO(I1,J1,2) = H2SAO(I1,J1,2) + SOZINT 
        H2SAO(I1,J1,IZPOS) = H2SAO(I1,J1,IZPOS) + SOZINT 
C       H2SAO(I1,J1,3) = H2SAO(I1,J1,3) + SOYINT 
        H2SAO(I1,J1,IYPOS) = H2SAO(I1,J1,IYPOS) + SOYINT 
C       H2SAO(I1,J1,4) = H2SAO(I1,J1,4) + SOXINT 
        H2SAO(I1,J1,IXPOS) = H2SAO(I1,J1,IXPOS) + SOXINT 

C       H2SAO(J1,I1,2) = H2SAO(J1,I1,2) - SOZINT 
        H2SAO(J1,I1,IZPOS) = H2SAO(J1,I1,IZPOS) - SOZINT 
C       H2SAO(J1,I1,3) = H2SAO(J1,I1,3) - SOYINT 
        H2SAO(J1,I1,IYPOS) = H2SAO(J1,I1,IYPOS) - SOYINT 
C       H2SAO(J1,I1,4) = H2SAO(J1,I1,4) - SOXINT 
        H2SAO(J1,I1,IXPOS) = H2SAO(J1,I1,IXPOS) - SOXINT 
       ENDIF
       IX = IX + 1
      ENDDO
      ENDDO
        write(lupri,*) '*******************************************'
      ENDDO

C     ... when you have integrals, compare them..
C        ... release memmory...
      CALL MEMREL('KSOZ',WORK,KWORK,KSOZ,KFREE,LFREE)
      CALL MEMREL('KSOY',WORK,KWORK,KSOY,KFREE,LFREE)
      CALL MEMREL('KSOX',WORK,KWORK,KSOX,KFREE,LFREE)
C     CALL MEMREL('AMFI',WORK,KWORK,KAMFI,KFREE,LFREE)

      CALL GPCLOSE(MNFPRP,'KEEP')

      END DO

C    .. close the MNF.INP file
      CALL GPCLOSE(LUTMP2,'KEEP')

      IF (ONLYSO1) THEN
C       write(lupri,*) ITD,'>>>',TOTDIFX,TOTDIFY,TOTDIFZ
        TOTDIFX=TOTDIFX/DFLOAT(ITD)
        TOTDIFY=TOTDIFY/DFLOAT(ITD)
        TOTDIFZ=TOTDIFZ/DFLOAT(ITD)

        WRITE(LUPRI,'(/3X,A/)')
     &  '**** AMFIIN: average differences between AMFI and '//
     &  'BSS:SO1 spin-orbit integrals (point nucleus) ****'
        WRITE(LUPRI,'(2X,A,D20.10)') '....for SO1X integrals:',TOTDIFX
        WRITE(LUPRI,'(2X,A,D20.10)') '....for SO1Y integrals:',TOTDIFY
        WRITE(LUPRI,'(2X,A,D20.10)') '....for SO1Z integrals:',TOTDIFZ

        CALL QUIT('AMFIIN: Only BSS-AMFI SO1 integral comparisons !')
      ENDIF

      CALL QEXIT('AMFIIN')
      RETURN
      END

      INTEGER FUNCTION INDX_LML(IL,IML,I)
CMI  ... called from AMFI/contandmult
#include "implicit.h"
#include "priunit.h"
      DIMENSION IL(*),IML(*)

       IX = I
        J = 0
      LX = IL(I)
      MLX=IML(I)

 20   CONTINUE
      IF (LX.EQ.IL(IX) .AND. MLX.EQ.IML(IX) .AND. IX.GT.0) THEN
         IX = IX - 1
         J = J + 1
         GOTO 20
      ENDIF
     
      INDX_LML = J 

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck UNITE_INDXS */
      SUBROUTINE UNITE_INDXS(NSPH,LABINT,IPNUC,
     &                   IPTYPAO,INDXLM,ISTNUCTYP,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  PURPOSE: If the basis set is splitted, do unite the INDXLM so that their fit
C           for entering of the AMFI integrals.
C
C  Written by Miro Ilias, Strasbourg, 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"

      DIMENSION IPNUC(NSPH),IPTYPAO(NSPH),INDXLM(NSPH),ISTNUCTYP(NSPH,2)
      CHARACTER LABINT(MXQN**2)*4
      LOGICAL ITWAS

      CALL QENTER('UNITE_INDXS')

      IF (IPRINT.GE.2) THEN
       write(lupri,'(/2X,A)')
     & 'Entering order of basis functions in spher. AO basis:'
       write(lupri,*) ' '
       DO I=1, NSPH
        write(lupri,*) I,'.fu center:',IPNUC(I),' type:',
     &  LABINT(IPTYPAO(I)),'(',INDXLM(I),')'
       ENDDO
       write(lupri,*) ' '
      ENDIF

C    ... list of stored INUC&IPTYP
      IT = 0 
      IL = 0
      I = 0
 234  I = I + 1

      INUC=IPNUC(I)
      IPTYP=IPTYPAO(I)

C     ...find out if INUC&IPTYP are not in the list
      ITWAS=.FALSE.
      DO K=1,IL
       IF (INUC.EQ.ISTNUCTYP(K,1).AND.IPTYP.EQ.ISTNUCTYP(K,2))
     &    ITWAS=.TRUE.
      ENDDO
      IF (ITWAS) GOTO 234

C     ... store INUC&IPTYP into the store_list_array
      IL=IL+1
      ISTNUCTYP(IL,1)=INUC 
      ISTNUCTYP(IL,2)=IPTYP 
C     ... check out indexes from J
      JJ = 0
      J = I-1
 235  J = J + 1

      IF (INUC.EQ.IPNUC(J).AND.IPTYP.EQ.IPTYPAO(J).AND.J.LE.NSPH) THEN
      JJ = JJ + 1
      IT = IT + 1
        IF (IPRINT.GE.2) THEN
        write(lupri,*) 'I=',I,' J=',J,
     &   ' for INUC=',INUC,' IPTYP=',IPTYP,
     &   ' indexes:', JJ,INDXLM(J)
        ENDIF

        IF (INDXLM(J).NE.JJ) THEN
         IF (IPRINT.GE.2) THEN
          write(lupri,*) 'difference in indexes !'
         ENDIF
C          ... change the index !
           INDXLM(J)=JJ
        ENDIF
      ENDIF
      IF (J.LT.NSPH) GOTO 235 
      IF (I.LT.NSPH) GOTO 234

      IF (IT.NE.NSPH) THEN
        write(lupri,*) 'IT=',IT,' NSPH=',NSPH
        CALL QUIT('UNITE_INDXS: IT<>NSPH !')
      ENDIF

C     ... when indexes are changed, print them out....
      IF (IPRINT.GE.1) THEN
       write(lupri,'(2X,A)')
     &  'UNITE_INDXS: Order of spher.atomic basis functions:'
       DO I=1, NSPH
         write(lupri,*) I,'.fu center:',IPNUC(I),' at.fu.type(indx):',
     &   LABINT(IPTYPAO(I)),'(',INDXLM(I),')'
       ENDDO
      ENDIF


      CALL QEXIT('UNITE_INDXS')
      RETURN
      END


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck amfi_stuff */
      SUBROUTINE AMFI_STUFF(H2CAO,ONLYSO1,
     &                      ONLYTRA,IC_AMFI,WORK,LWORK,IPRINT,IMFCHA)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: 
C   ==========
C Routine to include AMFI integrals, or simply do comparison 
C or only the transformation 
C
C On input: H2CAO - entering 2c integrals in SA-AO basis...
C                  
C           ONLYSO1 - flag for making AMFI-SO1 and BSS-SO1 integrals comparison
C           ONLYTRA - flag announcing only CSC transformation (testing purposes)
C           IC_AMFI - flag for the AMFI stuff
C           IMFCHA  - total charge of the system - important for the mean-field summation
C
C On output: H2CAO - CSC trasformed (+ AMFI, if wished) 2c SA-AO integrals.
C
C Written by Miro Ilias, Strasbourg, 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbqal.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"

      DIMENSION H2CAO(*),WORK(LWORK)
      CHARACTER LABINT(MXQN**2)*4
      LOGICAL ONLYSO1,ONLYTRA,FNDLAB,TESTCSC

      CALL QENTER('AMFI_STUFF')
#include "memint.h"

C        ... no dcbham.h common block !
       IPRHAM=IPRINT

C ========================================================
C       ... read H2C_AO into the H2CAO file
C ========================================================
      REWIND LUBSS
      IF (FNDLAB('H1AO_DK ',LUBSS)) THEN
         CALL READT(LUBSS,N2BBASXQ,H2CAO)
         IF (IPRINT.GE.2) WRITE(LUPRI,'(/2X,A)')
     &  'AMFI_STUFF: H1AO_DK read into H2CAO.'
      ELSE
        CALL QUIT(
     &  'AMFI_STUFF: H1AO_DK is not found in BSSMAT file!')
      END IF


      IF (IPRHAM.GE.5) THEN
C    ... print out H2c matrix in AO "xhi" basis
        CALL HEADER(
     &'AMFI_STUFF: READ TOTAL H(2c)_LL in "xhi" SA-AO basis ' ,-1)
        CALL PRQMAT(H2CAO(1),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C ... insert phases
      IREP=0
      IF (NZ.LT.4) THEN
      DO IZ = 1,NZ                                                 
         IQ = IPQTOQ(IZ,IREP)                                       
         CALL Q2BPHASE('F',IQ,1,H2CAO(1+N2BBASX*(IZ-1)))
      ENDDO            
      ENDIF           

C   ... transform to BU basis  
      CALL BSTOBU(H2CAO,NZ,WORK(KFREE),LFREE)

      CALL MEMGET('REAL',KBUF1,N2BBASX*4,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF2,N2BBASX*4,WORK,KFREE,LFREE)

C  ... transform to the AO basis ... as in dirone.F/PR1CEX2
      IPAR = JBTOF(IREP,1)
      NB = NTBAS(0)
      CALL DZERO(WORK(KBUF1),4*N2BBASX)
      DO IZ =1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL MTSOAO(H2CAO(1+N2BBASX*(IPQ-1)),
     &             WORK(KBUF1+N2BBASX*(IZ-1)),
     &             NB,IREPD,IPRHAM)
      ENDDO

C     CALL DCOPY(N2BBASXQ,WORK(KH2AOS1),1,WORK(KH2AOS),1)

      IF (IPRHAM.GE.4) THEN
        CALL HEADER(
     &'AMFI_STUFF: TOTAL H(2c)_LL in "xhi" C AO basis ' ,-1)
C        ... use IQDEF instead of IPQTOQ only when beeing in AO basis!
        CALL PRQMAT(WORK(KBUF1),NB,NB,NB,NB,4,IQDEF,LUPRI)    
      ENDIF

C      ... generate transf.coefficients...
      CALL SPHGEN(1,2,.FALSE.,WORK(KFREE),LFREE,IPRHAM)

C      ... transform C AO, WORK(KBUF1) -> S AO, WORK(KBUF2) 
      CALL DZERO(WORK(KBUF2),4*N2BBASX)
      LCOMP=0 
      DO IZ = 1, 4
       CALL MTAOCS(WORK(KBUF1+N2BBASX*(IZ-1)),NB,
     &         WORK(KBUF2+N2BBASX*(IZ-1)),NB,LCOMP,
     &         WORK(KFREE),LFREE,IPRHAM)
      ENDDO

      IF (IPRHAM.GE.3) THEN
        CALL HEADER(
     &'AMFI_STUFF: TOTAL ENTERING H(2c)_LL in "xhi" S AO basis ' ,-1)
        CALL PRQMAT(WORK(KBUF2),NB,NB,NB,NB,4,IQDEF,LUPRI)    
      ENDIF

C   ... find out total number of L-component spherical functions
      NSPH = NSPHCM(1)

      CALL MEMGET('INTE',KHINLM,NSPH,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPTYPAO,NSPH,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPNUC,NSPH,WORK,KFREE,LFREE)
      CALL MEMCHK('XXX',WORK,1)
      CALL SPHLAB(MXQN-1,LABINT)
       write(lupri,*) 'SPHLAB done, NSPH=',NSPH
C      write(lupri,*) 'NLRGSH=',NLRGSH


C     ...allocations for AMFI integrals: typ and index
       CALL MEMGET('INTE',KAIPT,NSPH,WORK,KFREE,LFREE)
       CALL MEMGET('INTE',KAILM,NSPH,WORK,KFREE,LFREE)
       CALL MEMGET('INTE',KAMFU,NSPH,WORK,KFREE,LFREE)
       CALL MEMGET('INTE',KSTNUCTYP,2*NSPH,WORK,KFREE,LFREE)


C      ... find out total number of atomic centers of the system
       NTNUC = 0
       DO I=1,NUCIND
         NTNUC = NTNUC + NUCDEG(I)
       ENDDO
       write(lupri,*) 'Total number of centers:',NTNUC

       CALL MEMGET('INTE',KZTNUCS,NTNUC,WORK,KFREE,LFREE)
       CALL MEMGET('INTE',KANUCS,NTNUC,WORK,KFREE,LFREE)
C      .... charges per nucleis
       CALL MEMGET('REAL',KCHPNUC,NTNUC,WORK,KFREE,LFREE)
       CALL DZERO(WORK(KCHPNUC),NTNUC)

       CALL MEMCHK('before AMFI2',WORK,1)

       CVAL_AMFI = CVAL

       CALL AMFIIN(WORK(KBUF2),NB,NTNUC,WORK(KIPTYPAO),WORK(KHINLM),
     &  WORK(KIPNUC),WORK(KZTNUCS),WORK(KAIPT),WORK(KAILM),WORK(KAMFU),
     &  WORK(KANUCS),WORK(KSTNUCTYP),WORK(KCHPNUC),
     &  NSPH,LABINT,ONLYSO1,ONLYTRA,IMFCHA,IC_AMFI,JM4POS,
     &  CVAL_AMFI,IPRHAM,WORK(KFREE),LFREE)

      CALL MEMCHK('before S->C',WORK,1)

      IF (IPRHAM.GE.3) THEN
        CALL HEADER(
     &  'AMFI_STUFF: TOTAL H(2c)_LL (after AMFIIN)'//
     &  ' in "xhi" S AO basis' ,-1)
        CALL PRQMAT(WORK(KBUF2),NB,NB,NB,NB,4,IQDEF,LUPRI)    
      ENDIF

C     write(lupri,*) '>>>>>  IQDEF:',(IQDEF(i),i=1,4)

C ===================================================================
C     now transform back: S AO, WORK(KBUF2) -> C AO, WORK(KBUF1) 
C ===================================================================
      CALL DZERO(WORK(KBUF1),N2BBASX*4)
C     ... generate with appropriare factors !!
      CALL SPHGEN(1,2,.TRUE.,WORK(KFREE),LFREE,IPRHAM)
      DO IZ = 1, 4
       CALL MTAOSC(WORK(KBUF1+N2BBASX*(IZ-1)),NB,
     &          WORK(KBUF2+N2BBASX*(IZ-1)),NB,LCOMP,
     &          WORK(KFREE),LFREE,IPRHAM)
      ENDDO

C   ... point to check H2c in C AO basis again...
      IF (IPRHAM.GE.4) THEN
        CALL HEADER(
     &'AMFI_STUFF: TOTAL H(2c)_LL BACK in "xhi" C AO basis ' ,-1)
        CALL PRQMAT(WORK(KBUF1),NB,NB,NB,NB,4,IQDEF,LUPRI) 
      ENDIF

      TESTCSC=.FALSE.
      IF (TESTCSC) THEN
       CALL DZERO(WORK(KBUF2),N2BBASX*4)
C       ... needed to go back to S
       CALL SPHGEN(1,2,.FALSE.,WORK(KFREE),LFREE,IPRHAM)
       DO IZ = 1, 4
         CALL MTAOCS(WORK(KBUF1+N2BBASX*(IZ-1)),NB,
     &         WORK(KBUF2+N2BBASX*(IZ-1)),NB,LCOMP,
     &         WORK(KFREE),LFREE,IPRHAM)
       ENDDO

        IF (IPRHAM.GE.6) THEN
          CALL HEADER(
     &'AMFI_STUFF: TOTAL CSC_H(2c)_LL AGAIN in "xhi" S AO basis ' ,-1)
          CALL PRQMAT(WORK(KBUF2),NB,NB,NB,NB,4,IQDEF,LUPRI)    
        ENDIF

       CALL DZERO(WORK(KBUF1),N2BBASX*4)
      CALL SPHGEN(1,2,.TRUE.,WORK(KFREE),LFREE,IPRHAM)

C        ... and from S to C
        DO IZ = 1, 4
         CALL MTAOSC(WORK(KBUF1+N2BBASX*(IZ-1)),NB,
     &         WORK(KBUF2+N2BBASX*(IZ-1)),NB,LCOMP,
     &         WORK(KFREE),LFREE,IPRHAM)
        ENDDO

        IF (IPRHAM.GE.6) THEN
          CALL HEADER(
     &'AMFI_STUFF: TOTAL CSC_H(2c)_LL AGAIN in "xhi" C AO basis ' ,-1)
          CALL PRQMAT(WORK(KBUF1),NB,NB,NB,NB,4,IQDEF,LUPRI)    
        ENDIF

        CALL QUIT('CSC  again')
      ENDIF

C ============================================================
C   ... FINALLY: transform BACK from C AO(BU) to the SA-AO (BS) basis
C            WORK(KBUF1) -> H2CAO(1)
C ============================================================
      CALL DZERO(H2CAO,N2BBASXQ)
      IREP=0
      IPAR = JBTOF(IREP,1)
      NB = NTBAS(0)
      DO IZ = 1, 4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
       CALL MTAOSO(WORK(KBUF1+N2BBASX*(IZ-1)),
     &             H2CAO(1+N2BBASX*(IPQ-1)),
     &             NB,IREPD,IPRHAM)
      ENDDO

      CALL BUTOBS(H2CAO,NZ,WORK(KFREE),LFREE) 

      IF(NZ.LT.4) THEN                              
        DO IZ = 1,NZ                                
          IQ = IPQTOQ(IZ,IREP)                    
          CALL Q2BPHASE('F',IQ,1,H2CAO(1+N2BBASX*(IZ-1)))  
        ENDDO                                    
      ENDIF
C        ... print statement...
       IF (IPRHAM.GE.6) THEN
C    ... print out H2c matrix in AO "xhi" basis
        CALL HEADER(
     &'AMFI_STUFF: AGAIN TOTAL H(2c)_LL in "xhi" SA-AO basis ' ,-1)
        CALL PRQMAT(H2CAO,NTBAS(0),NTBAS(0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

      CALL MEMREL('AMFI_STUFF',WORK,KWORK,KWORK,KFREE,LFREE)

C  SAVE THE CSC TRANSFOMED INTEGRALS 
C     ... include/dcbgen.h..LUBSS
C ...save hamiltonian matrix elements
      REWIND(LUBSS)
      IF (FNDLAB('EOFLABEL',LUBSS)) THEN
       BACKSPACE(LUBSS)
       CALL NEWLAB('H1AO_DK2',LUBSS,LUPRI)
       CALL WRITT(LUBSS,N2BBASXQ,H2CAO)
       WRITE(LUPRI,'(/2X,A)')
     &'AMFI_STUFF: New BSS2 relativistic atomic integrals'//
     &' - H1AO_DK2 - are written to the BSSMAT file.'
      ELSE
       CALL QUIT(
     &  'AMFI_STUFF: Error, EOFLABEL not found in BSSMAT file!')
      END IF

      CALL QEXIT('AMFI_STUFF')
      RETURN
      END

      subroutine matrix_inversion ( A, NP, IPIV, INDXR, INDXC )
      !-------------------------------------------------------------------------
      !
      !	      Taken from "Numeric recipes".  The original program was
      !       GAUSSJ which solves linear equations by the Gauss_Jordon
      !       elimination method.  Only the parts required to invert
      !	      matrices have been retained.
      !
      !	      J.P. Griffith  6/88
      !
CMI       A - input/output....
      !-------------------------------------------------------------------------
#include "implicit.h" 
#include "priunit.h"  

C     PARAMETER (NMAX=50)

      DIMENSION A(NP,NP), IPIV(*), INDXR(*), INDXC(*)

      n = np

      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
CMI             PAUSE 'Singular matrix'
                CALL QUIT('matrix_inversion: Singular matrix!')
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE

        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
CMI     IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
        IF (A(ICOL,ICOL).EQ.0.) CALL QUIT('Singular matrix.')
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
