C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck make_bss */
       SUBROUTINE MAKE_BSS(VMAT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C On input : VMAT - Lowdin transformation matrix for 4-comp. Fock matrix
C
C On output: VMAT - Lowdin transformation matrix modified for 2-comp. Fock matrix
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
        LOGICAL DK2SO1,ONLYSO1,SF_BEG,GETPCTM

        CALL QENTER('MAKE_BSS')
#include "memint.h"


C     Print section
C
      CALL HEADER('*** Output from MAKE_BSS ***',-1)

C    ... refrence
      WRITE(LUPRI,'(2X,A/,A/,A/)')
     & 'Jensen H.-J.Aa., Ilias M., "Two-component'//
     & ' relativistic methods based on the quaternion modified',
     & '  dirac equation: from the Douglas--Kroll to the'//
     & ' Barysz--Sadlej--Snijders infinite order";',
     & '  in preparation (2005)'

C extract information about the BSS Hamiltonian 
      CALL EXTR_BSS_INFO(IBSS,IPZ,SPINFR,DK2SO1,ONLYSO1,SF_BEG,IPRHAM)

C generate free particle one-electron matrix in Lowdin-basis "Theta"
      CALL MEMGET('REAL',KTBUF,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMAT,N2BBASXQ,WORK,KFREE,LFREE)
C    ... into WORK(KTMAT)
      CALL GET_FP_LOWD(WORK(KTBUF),WORK(KTMAT),
     &              VMAT,IPRHAM,WORK(KFREE),LFREE)

C generate the transformation matrix Lowdin ON basis "Theta" -> ON "theta" basis
C        ... into the WORK(KTBUF)...
      CALL MEMGET('REAL',KEIG ,NTBAS(0),WORK,KFREE,LFREE)
      CALL MAKE_RKB(WORK(KTBUF),WORK(KTMAT),WORK(KEIG),
     &                 IPRHAM,WORK(KFREE),LFREE)

C get FP (kinetic) matrix in the ON "theta" basis into WORK(KTBUF)
C get the transformation matrix AO "xhi" -> ON "theta" basis into WORK(KTMAT)
C     ....
      CALL GET_FP_RKB(WORK(KTBUF),WORK(KTMAT),VMAT,
     &                    IPRHAM,WORK(KFREE),LFREE)

C diagonalize the FP matrix in the ON "theta" basis
C .. get WORK(KIBEIG)=boson irreps
      CALL MEMGET('INTE',KIBEIG,NTBAS(0),WORK,KFREE,LFREE)
      CALL DIAG_FP(WORK(KTBUF),WORK(KTMAT),WORK(KEIG),
     &             WORK(KIBEIG),WORK(KFREE),LFREE)

C get A factors 
      CALL MEMGET('REAL',KFACKIN,NTBAS(1),WORK,KFREE,LFREE)
      CALL GET_A_FAC(WORK(KFACKIN),WORK(KEIG),
     &              WORK(KTMAT),IPRHAM,WORK(KFREE),LFREE)

C get H1 (Dirac bare nucleus) Hamiltonian in the FP ON "k" basis 
      CALL MEMGET('REAL',KHDAO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL GET_H1_FPM(WORK(KTBUF),WORK(KTMAT),WORK(KIBEIG),
     &              WORK(KHDAO),SF_BEG,WORK(KFREE),LFREE)
C    ... save H1 in "k" basis for the last comparison of eigenvalues...
      CALL MEMGET('REAL',KH1FP,N2BBASXQ,WORK,KFREE,LFREE)
      CALL DCOPY(N2BBASXQ,WORK(KTBUF),1,WORK(KH1FP),1)

C ... allocate accurate memmory space for H2c in FP "k" basis
      NH2FP = NESH(1)*NESH(1)*NZ
      IF (NFSYM.EQ.2) NH2FP=NH2FP+NESH(2)*NESH(2)*NZ
      CALL MEMGET('REAL',KH2FP,NH2FP,WORK,KFREE,LFREE)

C get H2C Hamiltonian(inf,DK2,DK2SO1,DK1)  in FP k basis 

C     ... by default do prepare the picture change transformation matrix...
      GETPCTM=.FALSE.
      IF (IPZ.EQ.9) THEN
C       ... infinite irder ... get also picture change transformation matrix
        CALL MEMGET('REAL',KPCTM,N2BBASXQ,WORK,KFREE,LFREE)

C       ... for the infinite order, do always prepare pict.change transf.matrix
        GETPCTM=.TRUE.
        CALL GET_HINF_FPM(WORK(KTMAT),WORK(KTBUF),WORK(KIBEIG),
     &             WORK(KHDAO),WORK(KH2FP),WORK(KPCTM),
     &             GETPCTM,SF_BEG,WORK(KFREE),LFREE)
      ELSE IF (IPZ.EQ.2) THEN
C       ... second order DKH
        CALL GET_HDK2_FPM(WORK(KH2FP),WORK(KTBUF),WORK(KEIG),
     &                    WORK(KIBEIG),DK2SO1,
     &                    IPRHAM,WORK(KFREE),LFREE)
      ELSE
C       ...and very rough first order, or SO1 integrals only...
        CALL GET_HDK1_FPM(WORK(KH2FP),WORK(KTBUF),WORK(KIBEIG),
     &                    ONLYSO1,IPRHAM,WORK(KFREE),LFREE)
      ENDIF

C  Make contravariant transformation matrix
C   (C+el)-1 = Sao(LL) C(el) Ad^(-1) to transform H2c "k"-> H2c "xhi" LL
      NCKLQ = NFBAS(1,0)*NESH(1)*NZ
      IF (NFSYM.EQ.2) NCKLQ = NCKLQ + NFBAS(2,0)*NESH(2)*NZ
      CALL MEMGET('REAL',KCTM,NCKLQ,WORK,KFREE,LFREE)

      CALL MKSAOMO_Le(WORK(KTBUF),WORK(KCTM),WORK(KTMAT),WORK(KFACKIN),
     &            SPINFR,SF_BEG,GETPCTM,IPRHAM,WORK(KFREE),LFREE)

C transform H2C Hamitonian from FP "k" basis to the SO-AO  basis,
C prepare the final picture change transformation matrix if wished
      CALL MAKE_H2C_AO(WORK(KCTM),WORK(KTMAT),WORK(KTBUF),
     &                 WORK(KH2FP),WORK(KHDAO),WORK(KPCTM),
     &                 WORK(KIBEIG),
     &                 GETPCTM,SPINFR,SF_BEG,
     &                 IPRHAM,WORK(KFREE),LFREE)

C       ... compare eigenvalues of 4c and 2c 1-el.Hamiltonians
      CALL CMP_EIGVAL(VMAT,WORK(KTMAT),WORK(KH1FP),
     &                  WORK(KEIG),WORK(KTBUF),WORK(KIBEIG),
     &                  SPINFR,SF_BEG,IPRHAM,WORK(KFREE),LFREE)

C  Finally set up number of positronic shells equal to zero
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

C **********************************************************************************
C    ... Number of transf. matrixes is set to 1 - it is only the Lowdin AO->MO
      NZT = 1
C **********************************************************************************

C     WRITE(LUPRI,'(/2X,A/)') 
C    &'MAKE_DKH:Two-component BSS mode - positronic shells deleted !'

      CALL MEMREL('MAKE_BSS',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('MAKE_BSS')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck get_fp_lowd */
       SUBROUTINE GET_FP_LOWD(TBUF,TMAT,VMAT,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    PURPOSE:
C  ===========
C  Get the free particle one-electron matrix in the ON Lowdin-basis "Theta"
C
C    On input: VMAT ... IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)
C           
C              TBUF ... buffer
C
C    On output: TMAT ... IOFT = IOFT + NFORB(I,0)*NFORB(I,0)*NZ
C
C Written by Miro Ilias & Hans Jorgen Aa Jensen, Odense, 2002, 
C         reprogrammed by Miro Ilias, Strasburg, 2005.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0 )
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION TBUF(*),VMAT(*),TMAT(*),WORK(LWORK)

      CALL QENTER('GET_FP_LOWD')     
#include "memint.h"

C=====================================================================
C First construct quaternion one-electron free particle in SO-basis 'xhi'
C=====================================================================
      CALL FREEMT(TBUF,IPRINT,WORK(KFREE),LFREE) 

C=====================================================================
C Then transform to Lowdin basis "Theta" (from TBUF to TMAT):
C
C                TMAT = VMATt *  TBUF  * VMAT
C=====================================================================
      IOFT = 1
      IOFV = 1
      DO I = 1,NFSYM
        IF(NFORB(I,0).GT.0) THEN

        IF (IPRINT.GE.7) THEN
           CALL HEADER(
     &  'GET_FP_LOWD: Free particle matrix in SA-AO basis "xhi"',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TBUF(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                 NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)

           CALL HEADER(
     &  'GET_FP_LOWD: Lowdin transf. matrix "xhi"->"Theta"',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(VMAT(IOFV),NFBAS(I,0),NFORB(I,0),
     &                 NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),LUPRI)

        ENDIF
C=====================================================================
C                TMAT  =  VMATt *  TBUF  * VMAT
C=====================================================================
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0), 
     &                         NFORB(I,0),NFORB(I,0),
     &            TBUF(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            TMAT(IOFT),NFORB(I,0),NFORB(I,0),NZ,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRHAM)

        IF (IPRINT.GE.7) THEN
          CALL HEADER(
     &'GET_FP_LOWD: Free particle matrix in ON Lowdin basis "Theta"',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
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
C  /* Deck TEST_UNIT_MTX */
      SUBROUTINE TEST_UNIT_MTX(S,N,LR,LC,NZ)
C****************************************************************************
C
C  PURPOSE: Check the enterning matrix - is it unit mtx ?
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
      SUBROUTINE EXTR_BSS_INFO(IBSS,IPZ,SPINFR,DK2SO1,
     &                         ONLYSO1,SF0,IPRINT)
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
C
C On output: IPZ - order of Hamiltonian
C            DK2SO1,SF0
C            ONLYSO1
C
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"

      LOGICAL SPINFR,DK2SO1,SF0,ONLYSO1

      CALL EXTR_DIG(IBSS,IPX,IPY,IPZ,IPRINT)

      WRITE (LUPRI,'(A,I2)')
     &  '* Order of Barysz-Sadlej_Snijders transformation :',IPZ
      IF (SPINFR.AND.IPY.GE.1) THEN 
       WRITE(LUPRI,'(A)') '* "SPINFR" option, no spin-orbit terms!' 
       IPY=0
      ENDIF
      IF (IPY .EQ. 0) THEN
         IF (.NOT.SPINFR) CALL QUIT(
     & 'You forgot to specify .SPINFR in the Hamiltonian input!')
         DK2SO1 = .FALSE.
      ELSE IF (IPY .EQ. 1 .AND. IPZ .EQ. 2 .AND.(.NOT.SPINFR)) THEN
         DK2SO1 = .TRUE.
         WRITE (LUPRI,'(A)')
     &   '* spin-orbit terms only to the first order -'//
     &   ' DKH2sfB+SO1 Hamiltonian'
      ELSE IF (.NOT.SPINFR) THEN
         DK2SO1 = .FALSE.
         WRITE (LUPRI,'(A)')
     &   '* all spin-orbit terms included (full relativity)'
      END IF
      IF (IPX .GE. 1 .AND. SPINFR) THEN  
         SF0 = .TRUE.
         WRITE (LUPRI,'(A)')
     & '* Elimination of spin-orbit terms from the Dirac bare '//
     & 'nuclei in "k" basis (spinfree "BEFORE / AT THE BEGINNING")' 
      ELSE
         SF0 = .FALSE.
         IF (SPINFR) WRITE (LUPRI,'(A)')
     & '* Elimination of spin-orbit terms from the 2comp.Ham. in "k"'//
     & ' basis (spinfree "AFTER / AT THE END")'
      ENDIF
     
      ONLYSO1=.FALSE.
      IF (IPX.EQ.0.AND.IPY.EQ.1.AND.IPZ.EQ.0) THEN
         WRITE (LUPRI,'(/A)')
     &  '* only SO1 terms (for comparison with AMFI, then quit)'
        ONLYSO1=.TRUE.
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
#include "implicit.h"
#include "priunit.h"

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
      ENDIF

      IF (IPRINT. GE. 5) THEN
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

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck */      
       SUBROUTINE MAKE_RKB(TBUF,TMAT,EIG,IPRINT,WORK,LWORK)
C---------------------------------------------------------------------
C    PURPOSE:
C    ========
C   Generate the transformation matrix:
C                 Lowdin ON basis "Theta" -> ON "theta" basis
C                                                         
C On input: TMAT ... Free particle matrix in "Lowdin" basis "Theta"
C
C On output: TBUF - desired transformation matrix
C               ...  IBUF = IBUF + NFORB(I,0)*NORB(I)*NZ
C                                                         
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C                                                         
C----------------------------------------------------------------------
#include "implicit.h"   
#include "priunit.h" 
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0 ) 
#include "dcbbas.h"  
#include "dcbgen.h"  
#include "dcborb.h"
#include "dgroup.h"   
#include "pgroup.h"   
#include "dcbdhf.h"
      DIMENSION TBUF(*),TMAT(*),EIG(*),WORK(LWORK)    

      CALL QENTER('MAKE_RKB')
#include "memint.h"                              

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

        NSYM = 4/NZ
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

C
C           Make (sigma.p)(sigma.p) = p2
C

            IREP = JFSYM(ISYM,I) - 1
            JSL = IMAT + IBOL * NFORB(I,0) + NFORB(I,1)
            JLS = IMAT + NFORB(I,1) * NFORB(I,0) + IBOL

            CALL QGEMM(NBORB(ISYM,I,1),NBORB(ISYM,I,1),NFORB(I,2),D1,
     &          'N','N',IPQTOQ(1,0),TMAT(JLS),NFORB(I,0),NFORB(I,0),NZ,
     &          'N','N',IPQTOQ(1,0),TMAT(JSL),NFORB(I,0),NFORB(I,0),NZ,
     &          D0,IPQTOQ(1,0),WORK(KSMT),
     &          NBORB(ISYM,I,1),NBORB(ISYM,I,1),1)

            IF (IPRINT.GE.7) THEN
              WRITE(LUPRI,'(A,A3)')
     &        'MAKE_RKB: * Boson irrep: ',REP(IREP)
              CALL HEADER('Kinetic energy matrix in Lowdin L basis',-1)
              CALL PRQMAT(WORK(KSMT),NBORB(ISYM,I,1),NBORB(ISYM,I,1),
     &                    NBORB(ISYM,I,1),NBORB(ISYM,I,1),
     &                    1,IPQTOQ(1,0),LUPRI)
            ENDIF

C Do canonical orthonormalization of small components
            CALL LOWGEN(WORK(KSMT),NBORB(ISYM,I,1),WORK(IVMT),NEFFI,
     &                  NBORB(ISYM,I,1),STOL(2),IPRINT,EIG,
     &                  WORK(KFREE),LFREE)
            IF (NEFFI.LT.NBORB(ISYM,I,1)) WRITE(LUPRI,'(A)') 
     &      'MAKE_RKB : *** WARNING *** Linear dependence in RKB basis'
            NEFF = NEFF + NEFFI
            NBORB(ISYM,I,2) = NEFFI


             IF (IPRINT.GE.7) THEN
              CALL HEADER(
     &'MAKE_RKB: Lowdin transf. matrix to ON modif. Dirac small basis',
     &        -1)
              CALL PRQMAT(WORK(IVMT),NBORB(ISYM,I,1),NBORB(ISYM,I,2),
     &                    NBORB(ISYM,I,1),NBORB(ISYM,I,2),
     &                    NZ,IPQTOQ(1,0),LUPRI)
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

C .... We have it !
          IF (IPRINT.GE.8) THEN
           CALL HEADER(
     &'MAKE_RKB:Transf. matrix from Lowdin "Theta" to ON basis "theta"',
     &     -1)
           CALL PRQMAT(TBUF(IBUF),NFORB(I,0),NORB(I),
     &                 NFORB(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

          IMAT = IMAT + NFORB(I,0)*NFORB(I,0)*NZ
          IBUF = IBUF + NFORB(I,0)*NORB(I)*NZ

C   ... release the memmory allocated in the cycle...
       CALL MEMREL('MAKE_RKB',WORK,KWORK,KVMT,KFREE,LFREE)
       CALL MEMREL('MAKE_RKB',WORK,KWORK,KSMT,KFREE,LFREE)

   40   CONTINUE

       CALL QEXIT('MAKE_RKB')
       RETURN
       END


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck GET_FP_RKB */      
       SUBROUTINE GET_FP_RKB(TBUF,TMAT,VMAT,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:
C    ========
C      get FP (kinetic) matrix in the ON "theta" basis into TBUF
C      get the transformation matrix AO "xhi" -> ON "theta" basis  TMAT
C
C On input:
C-----------
C  TBUF contains the transformation from the Lowdin ON basis "Theta"
C          to the ON "theta" basis.
C
C  VMAT still contains the transformation from SO basis "xhi" to
C       the Lowdin ON "Theta"
C
C  TMAT contains the free particle matrix in Lowdin basis "Theta".
C                                                         
C On output:   TMAT transformation matrix AO "xhi" -> ON RKB "theta" 
C-----------   TBUF FP (kinetic) matrix in the ON "theta" basis 
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
#include "dcborb.h"   
#include "dgroup.h"
      DIMENSION TBUF(*),VMAT(*),TMAT(*),WORK(LWORK)    

      CALL QENTER('GET_FP_RKB')
#include "memint.h"

      N2TMT  = 0
      NORBT  = 0
      N2ORBT = 0
      DO I = 1,NFSYM
        I2TMT(I) = N2TMT
        N2TMT    = N2TMT + NFORB(I,0)*NORB(I)*NZ
        IORB(I)  = NORBT
        NORBT    = NORBT  + NORB(I) ! Count total number of orbitals...
        I2ORBT(I)= N2ORBT*NZ
        N2ORBT   = N2ORBT + NORB(I)*NORB(I)
      ENDDO

C ... allocate array for the FP in "theta"
      CALL MEMGET('REAL',KBUF,N2BBASXQ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KBUF),N2BBASXQ)
C
C         Transform the free particle matrix
C         from TMAT in Lowdin "Theta" basis to WORK(IBUF)
C         in modified Dirac form (orhonormal "theta" basis); 
C
C                WORK(IBUF) =  TBUF+ * TMAT * TBUF 
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
     &           TBUF(1+I2TMT(I)),NFORB(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           TBUF(1+I2TMT(I)),NFORB(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRINT)
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

        IF (TIMAGNORM.GT.THRNULL1) THEN

          WRITE(LUPRI,'(2X,A,A,I1,A,I1,/,2X,A,D12.6)')
     & 'Free particle matrix in the ON basis "theta" ',
     & 'in fermion corep ',I,'/',NFSYM,
     & '*** Norm of the i,j,k imag. parts *** :',TIMAGNORM
          WRITE(LUPRI,'(/2X,A,D13.6)')
     & 'GET_FP_RKB: WARNING, TIMAGNORM > THRNULL1 (=',THRNULL1

        ENDIF

        IF (IPRINT.GE.6) THEN
          CALL HEADER(
     &'GET_FP_RKB:Free particle matrix in the ON basis "theta"',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &'*** Fermion corep ',I,'/',NFSYM
          WRITE(LUPRI,'(A,1P,D12.6)')
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
     &              VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1, IPRINT)

        IF (IPRINT.GE.7) THEN
          CALL HEADER(
     &  'GET_FP_RKB:Transf. mtx TMAT H4c "xhi"->final'//
     & ' ON basis "theta" (restr.kin.bal.)',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IOFT = IOFT + NFORB(I,0)*NFORB(I,0)*NZ
        IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
        IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)

 50   CONTINUE

C  *** Store the WORK(KBUF) into the TBUF(1) ***
      CALL DCOPY(N2BBASXQ,WORK(KBUF),1,TBUF(1),1)
C  *** Dealocate WORK(KBUF) into the TBUF(1) ***
      CALL MEMREL('MAKE_DKH 1',WORK,KWORK,KBUF,KFREE,LFREE)


      CALL QEXIT('GET_FP_RKB')
      RETURN
      END
 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck DIAG_FP */      
       SUBROUTINE DIAG_FP(TBUF,TMAT,EIG,IBEIG,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:
C    ========
C
C      Diagonalize TBUF with kinetic energy matrix (free particle, FP)
C      in the orthonormal "theta" basis as a pure real metrix
C
C     HJ:
C       Diagonalize this TMAT to get the Douglas-Kroll 4c |k> basis:
C       (the "exact" instead of Hess's so-called "SQR" approximation;
C       Hess states that the difference is small, and we h tested that
C       with a 25 s L basis set where we also got a small difference;
C       our model corresponds to Hess's "SQRA".)
C                                                         
C                                                         
C On input:
C 
C      TBUF reused for |k> eigenvectors in ON "theta" basis
C      TMAT - transf. matrix ON RKB "theta" -> SA-AO "xhi"
C
C On output:
C
C      EIG - eigenvalues 
C        - used also for checking the A prefactors
C        - used for 'denominator' in DKH2 method R_DK2;  eps(el)-eps(pos)
C
C      IBEIG - (approximate) boson symmetries of shells
C
C      TMAT - eigenvectors of FP ("k" basis) in SO-AO basis "xhi"
C
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0 ) 
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
      DIMENSION TBUF(*),TMAT(*),EIG(*),IBEIG(*),WORK(LWORK)
      LOGICAL DBGDIAG, DOJACO_SAVE

      CALL QENTER('DIAG_FP')
#include "memint.h"    


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

      IF (IPRHAM.GE.3) THEN
       WRITE(LUPRI,'(/,3X,A,A/)') 
     & 'DIAG_FP: Diagonalizing (Jacobi) the real part',
     & ' of the FP matrix in the ON "theta" basis !'
      ENDIF

C ... after diagonalization we get the IBEIG array which is used later !
      IF (IPRHAM.GE.5) THEN 
        DBGDIAG=.TRUE.
      ELSE
        DBGDIAG=.FALSE.
      ENDIF

      CALL MEMGET('REAL',KBUF,N2BBASXQ,WORK,KFREE,LFREE)
      CALL DCOPY(N2BBASXQ,TBUF(1),1,WORK(KBUF),1)

      IF (IPRHAM.GE.6) THEN
       DO IFSYM=1,NFSYM
          IBUF=KBUF+ I2ORBT(IFSYM)
          CALL HEADER(
     &'DIAG_FP:Free particle matrix in the ON basis "theta"',-1)
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',IFSYM
          CALL PRQMAT(WORK(IBUF),NORB(IFSYM),NORB(IFSYM),
     &              NORB(IFSYM),NORB(IFSYM),NZ,IPQTOQ(1,0),LUPRI)

       ENDDO
      ENDIF

C    ... subroutine in dirscf.F for diagonalizing the Dirac-Fock matrix
C   TBUF - eigenvectors in "k", EIG - eigenvalues, IBEIG - boson irreps
      CALL DFDIAG (WORK(KBUF),EIG,IBEIG,TBUF,DBGDIAG,
     &             WORK,KFREE,LFREE)

      NZ = NZSAVE
      DOJACO = DOJACO_SAVE
      SUB_BL = .FALSE.

      CVAL2 = CVAL*CVAL

      N2TMT = 0
      IMAT = 1
      IBUF = KBUF
      DO 60 I = 1,NFSYM
       IF(NORB(I).EQ.0) GOTO 60
       I2TMT(I) = N2TMT
       N2TMT = NFBAS(I,0)*NORB(I)*NZ

C   ... do always check the RKB free-particle eigenvalues !

        IF (IPRHAM.GE.2) THEN
         CALL HEADER(
     &'DIAG_FP: Free particle (in ON basis "theta") eigenvalues :',-1)
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
     &  'DIAG_FP:DKH |k> matrix TBUF in the ON basis "theta" ',-1)
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
C  TMAT = TMAT*TBUF (WORK(IBUF) reused as scratch area)
C-------------------------------------------------------------------------------
        CALL BCKTRA(WORK(IBUF),NFBAS(I,0),NORB(I),
     &              TBUF(IMAT),NORB(I),NORB(I),
     &              NORB(I),NZ,NORB(I),1,NFBAS(I,0),
     &              TMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),NZ,
     &              IPRHAM)
C   ... copy the result into the TMAT
        CALL DCOPY(N2TMT,WORK(IBUF),1,TMAT(I2TMT(I)+1),1)


        IMAT = IMAT + NORB(I)*NORB(I)*NZ ! Set up offset !!!
        IBUF = IBUF + IMAT  ! Upgrade the IBUF !

 60   CONTINUE  


      IF (IPRHAM.GE.2) THEN

C   Obtain overlap integrals in SO basis 'xhi' into WORK(KBUF)
      CALL MEMGET('REAL',KKUN,N2BBASXQ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KKUN),N2BBASXQ) ! Store here Ad^2 diagonal matrix

      CALL GTOVLT(WORK(KBUF),D1,WORK(KFREE),LFREE)
      IMAT = 1
      IBUF = KBUF
      IKUN1 = KKUN
      DO 80 I = 1,NFSYM
        IF (NORB(I).EQ.0) GOTO 80

C ... control print out:
         IF(IPRHAM.GE.7) THEN
          CALL HEADER(
     &    'DIAG_FP:DKH |k> matrix TMAT in SO basis "xhi" ',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     &    '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     &    'DIAG_FP: Full SOVLP matrix in SO basis "xhi" ',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     &    '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(WORK(IBUF),NFBAS(I,0),NFBAS(I,0),
     &                NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),LUPRI)
         ENDIF


**  DO C+ S C into WORK(KKUN)
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &       NORB(I),NORB(I),
     &       WORK(IBUF),NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),
     &       WORK(IKUN1),NORB(I),NORB(I),1,IPQTOQ(1,0),
     &       TMAT(IMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &       TMAT(IMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)

         CALL HEADER(
     &  'DIAG_FP: Full unit <k|k>=C+.S.C matrix  ',-1)
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

C      ... release all used memmory...
      CALL MEMREL('DIAG_FP',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('DIAG_FP')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_a_fac */      
       SUBROUTINE GET_A_FAC(FACKIN,EIG,TMAT,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: Check the diagonal elements of the A factor 
C    ========
C                                                         
C    On input: TMAT - FP eigenvector "k" in "xhi" SO-AO basis
C              EIG  - free particle eigenvalues                                                       
C                                                         
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
      DIMENSION EIG(*),TMAT(*),FACKIN(*),WORK(LWORK)    

      CALL QENTER('GET_A_FAC')
#include "memint.h"         


        CALL MEMGET('REAL',KBUF,N2BBASXQ,WORK,KFREE,LFREE)
        CALL GTOVLT(WORK(KBUF),D0,WORK(KFREE),LFREE)

        CALL MEMGET('REAL',KKUN,N2BBASXQ,WORK,KFREE,LFREE)
C       CALL DZERO(WORK(KKUN),N2BBASXQ) ! Store here Ad^2 diagonal matrix

      CVAL2 = CVAL*CVAL

      IMAT = 1
      IFC = 0
      IBUF = KBUF 
      IKUN1 = KKUN
      DO 82 I = 1,NFSYM
       IF (NORB(I).EQ.0) GOTO 82

         IF(IPRINT.GE.7) THEN
          CALL HEADER(
     &  'GET_A_FAC: Entering |k> matrix TMAT in SO basis "xhi" ',-1)
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',I
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          CALL HEADER(
     &        'GET_A_FAC: Obtained (GTOVLT) LL Sovlp matrix'//
     &        ' in SO basis "xhi"',-1)
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',I
          CALL PRQMAT(WORK(IBUF),NFBAS(I,0),NFBAS(I,0),
     &                NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),LUPRI)
         ENDIF

** DO C+(e) S(LL) C(e) into WORK(KKUN) **

        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &       NESH(I),NESH(I),
     &       WORK(IBUF),NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),
     &       WORK(IKUN1),NESH(I),NESH(I),1,IPQTOQ(1,0),
     &       TMAT(IMAT+NFBAS(I,0)*NPSH(I)),
     &       NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &       TMAT(IMAT+NFBAS(I,0)*NPSH(I)),
     &       NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRINT)

         IF(IPRINT.GE.2) THEN
          CALL HEADER(
     &    'MAKE_DKH: Diagonal <k_e|k_e>=A+.A matrix',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')'* Fermion irp no.',I,'/',NFSYM
          CALL TEST_DIAG_MTX(WORK(IKUN1),NESH(I),NESH(I),NESH(I),1)
          IF(IPRINT.GE.7) THEN
           CALL PRQMAT(WORK(IKUN1),NESH(I),NESH(I),
     &            NESH(I),NESH(I),1,IPQTOQ(1,0),LUPRI)
          ENDIF
         ENDIF

*************************************************************
*        Extract the A kinematic factors into FACKIN(1,*)   *
*************************************************************
         IF (IPRINT.GE.2) THEN 
           WRITE(LUPRI,'(/9X,A,I1,A,I1,//,12X,A/)')
     & 'MAKE_DKH: "A" kinem. electr. factors; ferm irp:',I,'/',NFSYM,
     & 'A^2_ovlp          A^2_calc          diff;'
         ENDIF

         SUMDIFF = D0
         DO J = 1, NESH(I)
          DAA = WORK(IKUN1+(J-1)*NESH(I)+J-1)
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

         WRITE(LUPRI,'(2X,A,I1,A,I1,A,D14.8)')
     & 'Average differences for A-factors (fsym ',
     & I,'/',NFSYM,') : ',SUMDIFF

         IF (SUMDIFF.GE.THRNULL) THEN
          WRITE(LUPRI,'(2X,A,D14.8)')
     & 'GET_A_FAC:  DIFFERENCES A^2_ovlp-A^2_calc/N DIFFER'//
     & ' TOO MUCH FROM "ZERO"=!',THRNULL
         ENDIF

         IFC   = IFC + NESH(I)
         IKUN1 = IKUN1 + NESH(I)*NESH(I)
         IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
         IBUF = IBUF + NFBAS(I,0)*NFBAS(I,0)

 82     CONTINUE

C   .... release all the allocated memmory at once...
      CALL MEMREL('GET_A_FAC',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('GET_A_FAC')

      RETURN
      END


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_h1_fpm */      
       SUBROUTINE GET_H1_FPM(TBUF,TMAT,IBEIG,HDAO,SF_BEG,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: Get the H1 Hamiltonian in the FP ON "k" basis
C    ========
C
C  On input:
C    TMAT - transformation matrix "xhi" -> "k" 
C   IBEIG - boson info on MO shells (needed for the spin-free case)
C                                                         
C  On output:
C    TBUF - H1 Hamiltonian in the FP ON "k" basis (spin-free,if SF_BEG .true.)                                            
C     HDAO - bare nucleus Dirac Ham; in SAO-basis "xhi"                                                         
C
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0 ) 
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"   
#include "dcbham.h"   
      DIMENSION TBUF(*),TMAT(*),HDAO(*),IBEIG(*),WORK(LWORK)    
      LOGICAL SF_BEG

      CALL QENTER('GET_H1_FPM')
#include "memint.h"         


C First get the bare nucleus Dirac Hamiltonian (BNDH) in the SO-basis "xhi"
C in WORK(KBUF)

      CALL MEMGET('REAL',KBUF,N2BBASXQ,WORK,KFREE,LFREE) 
      BSS=.FALSE.
      CALL ONEFCK(WORK(KBUF),IPRHAM,WORK(KFREE),LFREE) 
      BSS=.TRUE.

C   ... store the  bare nucleus Dirac Hamiltonian
      CALL DCOPY(N2BBASXQ,WORK(KBUF),1,HDAO,1)

C
C    transform the BNDH from its SO-basis "xhi" to the ON basis "|k>"
C    (reuse TBUF)
C  where H is the bare nucleus Hamiltonian in the SO basis "xhi"
C
C                            TBUF = TMAT+ * H  * TMAT; 
C                      H1(in |k>) =   U1+ * H0 * U1
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
      DO 85 I = 1,NFSYM
      IF(NORB(I).EQ.0) GOTO 85
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &       NORB(I),NORB(I),
     &       WORK(KBUF+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &       TBUF(IOFT),NORB(I),NORB(I),MZ,IPQTOQ(1,0),
     &       TMAT(1+I2TMT(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &       TMAT(1+I2TMT(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'GET_H1_FPM: Obtained Dirac H1 matrix in '//
     &  '4c "k" basis, TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IF (SF_BEG) THEN
         CALL SPF_HK( TBUF(IOFT),.TRUE.,IBEIG(IORB(I)+1),
     &                 IBEIG(IORB(I)+1),
     &                NORB(I),NORB(I),NORB(I),NORB(I),NZ,IPRHAM )
         IF(IPRHAM.GE.5) THEN
          CALL HEADER(
     &  'GET_H1_FPM: H1 Dirac matrix in 4c "k" basis '//
     &  ' after elim. spin-orbit terms "BEFORE", TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF
        ENDIF

C       ... update the offset
        IOFT = IOFT + NORB(I)*NORB(I)*NZ

   85 CONTINUE

      CALL MEMREL('DKH',WORK,KWORK,KBUF,KFREE,LFREE)

      CALL QEXIT('GET_H1_FPM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_hinf_fpm */      
      SUBROUTINE GET_HINF_FPM(TMAT,TBUF,IBEIG,HDAO,H2FP,
     &             PCTM,GETPCTM,SF_BEG,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: Get the H_inf 2c Hamiltonian in the FP ON "k" basis
C    ========  Check the YR-condition and R-condition as well
C
C  On input:
C --------------
C TBUF - BNDH in "k" basis
C TMAT - 4c "k" -> 4c AO "xhi"
C HDAO - BNDH in SA-AO "xhi" basis
C IBEIG _ boson irreps if spin-free from beginning
C
C  On output:
C --------------
C H2FP H2c Hamiltonian in the "k" ON FP basis (no elimination of spin-free terms!)
C
C PCTM  picture change transformation matrix 4c SO-AO "xhi" -> 2c "k", if desired (GETPCTM=.true.)
C
C TBUF  picture change transformation matrix 4c "k" -> 2c "k"
C
C Written by:  M.Ilias, HJ Aa Jensen, Odense, 2002,
C            last revison: M.Ilias, Strasbourg, february 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER( D1 = 1.0D0, DM1= -1.0D0, D0 = 0.0D0,THRNULL = 1.0D-10) 
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"   
#include "dcbham.h"   
      DIMENSION TBUF(*),TMAT(*),HDAO(*),IBEIG(*),
     &          H2FP(*),PCTM(*),WORK(LWORK)    
      LOGICAL   SPINFR_SAVE,DBGDIAG,GETPCTM,SF_BEG

      CALL QENTER('GET_HINF_FPM')
#include "memint.h"         

C diagonalize the BNDH in "k" basis

      CALL MEMGET('INTE',KIBEIG,NORBT,WORK,KFREE,LFREE)
C        ... the size must be NTBAS(0) !
      CALL MEMGET('REAL',KEIG,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH1VEC,N2BBASXQ,WORK,KFREE,LFREE)

      CALL MEMGET('REAL',KHD,N2BBASXQ,WORK,KFREE,LFREE)
C     ... save the H1 in the "k" basis...
      CALL DCOPY(N2BBASXQ,TBUF(1),1,WORK(KHD),1)

      IF (IPRHAM.GE.6) THEN 
         DBGDIAG = .TRUE.
      ELSE
         DBGDIAG = .FALSE.
      ENDIF

      SPINFR_SAVE = SPINFR
      SPINFR=.FALSE.
C     ... diagonalize the H1 BNDH in the  FP ON "k" basis
      CALL DFDIAG (TBUF(1),WORK(KEIG),WORK(KIBEIG),WORK(KH1VEC),
     &                     DBGDIAG,WORK,KFREE,LFREE)
      SPINFR      = SPINFR_SAVE

      JOFT1 = 1
      IH2FP=1
      IPCTM = 1
      IBUF = 1
C    Run big cycle over the fermion symmetries
      DO 95 I = 1,NFSYM  
      IF(NORB(I).EQ.0) GOTO 95

C      ... off-set for the eigenvectors...
       JH1VEC = KH1VEC + I2ORBT(I)

       IF(IPRHAM.GE.4) THEN
C          ... Write out the eigenvalues in the |k> basis
        CALL HEADER(
     &    'MAKE_DKH: Bare nucleus eigenvalues:',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1,4X,A)')
     &    '*** Fermion corep ',I,'/',NFSYM,'boson irrep'
        WRITE(LUPRI,'(I5,F25.10,I5)')
     &   (J,WORK(KEIG+IORB(I)+J-1),IBEIG(IORB(I)+J), J=1,NORB(I))
       ENDIF

       IF(IPRHAM.GE.7) THEN
C          ... Write out the eigenvectors in the |k> basis
       CALL HEADER(
     &'MAKE_DKH: Bare nucleus eigenvectors |v> in the |k> basis',-1)
       WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
       CALL PRQMAT(WORK(JH1VEC),NORB(I),NORB(I),
     &        NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C*******************************************************************
C ... prepare for the solving of "YR" linear matrix equation
C*******************************************************************

       JVECP = JH1VEC + NPSH(I)*NORB(I)
       JVECE = JVECP + NPSH(I)
C--------------------------------------------------------------------
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
C--------------------------------------------------------------------
       NREF = NESH(I)
       NREF2 = 2*NREF  
       N2REF = NREF*NREF
       N2REFQ = N2REF*NZ
       N2EQ = NREF*NREF*NZ

       CALL MEMGET('REAL',KAMAT,2*N2REFQ,WORK,KFREE,LFREE)
C       ... factor 2 for Omega+
       CALL MEMGET('REAL',KXMAT,2*N2REFQ,WORK,KFREE,LFREE)
       CALL MEMGET('INTE',KIPVT,NREF,WORK,KFREE,LFREE)

C           -- make A matrix
       CALL DZERO(WORK(KAMAT),N2REFQ)
       CALL DUNIT(WORK(KAMAT),NREF)   ! A := 1
       CALL QGEMM(NREF,NREF,NREF,DM1, ! A := A - Yp*(Yp)+
     &      'N','N',IPQTOQ(1,0),WORK(JVECP),NORB(I),NORB(I),NZ,
     &      'H','N',IPQTOQ(1,0),WORK(JVECP),NORB(I),NORB(I),NZ,
     &       D1,IPQTOQ(1,0),WORK(KAMAT),NREF,NREF,NZ)

       IF (IPRHAM.GE.7) THEN
         CALL HEADER( 'MAKE_DKH: A matrix for getting R',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(WORK(KAMAT),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL MEMGET('REAL',KA,N2REFQ,WORK,KFREE,LFREE) ! For A mtx
       CALL DCOPY(N2REFQ,WORK(KAMAT),1,WORK(KA),1)

C         -- make b matrix (store in X); b := Yp * (Ye)+
       CALL QGEMM(NREF,NREF,NREF,D1,
     &     'N','N',IPQTOQ(1,0),WORK(JVECP),NORB(I),NORB(I),NZ,
     &     'H','N',IPQTOQ(1,0),WORK(JVECE),NORB(I),NORB(I),NZ,
     &      D0,IPQTOQ(1,0),WORK(KXMAT),NREF2,NREF,NZ)

       IF (IPRHAM.GE.7) THEN
         CALL HEADER( 'MAKE_DKH: b matrix for getting R',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(WORK(KXMAT),NREF,NREF,NREF2,NREF,
     &              NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL MEMGET('REAL',KB,N2REFQ,WORK,KFREE,LFREE) ! For B mtx
       DO IZ = 1, NZ
       DO IE = 1, NESH(I)
         CALL DCOPY(NREF,
     &           WORK(KXMAT+(IZ-1)*NREF2*NREF+(IE-1)*NREF2),1,
     &           WORK(KB+(IZ-1)*N2REF+(IE-1)*NREF),1)
       ENDDO
       ENDDO

       JDVEC = JH1VEC

C        ... reuse memory with positron solutions for
C              D(*) vector in Choleskyl decomp.
       JOB = 1
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

C=============================================================================
C  that was R, check the YR-equation and the R-condition (for ferm.symm. I)
C=============================================================================
         CALL R_CHECK(I,WORK(KHD+JOFT1-1),WORK(KXMAT),
     &               WORK(KA),WORK(KB),IPRHAM,WORK(KFREE),LFREE)

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

C=====================================================================
C ...find Omega+ matrix and multiply it into WORK(KXMAT)
C=====================================================================
            CALL MEMGET('REAL',KH1AO,N2EQ,WORK,KFREE,LFREE)

            CALL QGEMM(NREF,NREF,NREF2,D1,  ! R+.R = H1AO
     &          'H','N',IPQTOQ(1,0),WORK(KXMAT),NREF2,NREF,NZ,
     &          'N','N',IPQTOQ(1,0),WORK(KXMAT),NREF2,NREF,NZ,
     &          D0,IPQTOQ(1,0),WORK(KH1AO),NREF,NREF,NZ)

            IF (IPRHAM.GE.7) THEN
              CALL HEADER(
     &  '  MAKE_DKH: X+ = ( 1 R+)',-1)
              WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
              CALL PRQMAT(WORK(KXMAT),NREF2,NREF,NREF2,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
              CALL HEADER(
     &  '  MAKE_DKH: X+ * X = 1 + R+ * R',-1)
              WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
              CALL PRQMAT(WORK(KH1AO),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
            ENDIF

            CALL QDIAG(NZ,NREF,WORK(KH1AO),NREF,NREF,
     &                 WORK(JH1VEC),1,WORK(KAMAT),NREF,
     &                 NREF,WORK(KFREE),LFREE,IERR)

            CALL MEMGET('REAL',KAMAT1,2*N2REFQ,WORK,KFREE,LFREE)
            CALL DCOPY(N2REFQ,WORK(KAMAT),1,WORK(KAMAT1),1)

            IF (IPRHAM.GE.5) THEN
             CALL HEADER(
     &    'MAKE_DKH: Eigenvalues of (1 + R+ * R)',-1)
              WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
             DO J=1,NREF
              WRITE(LUPRI,'(I5,F25.10)') J,WORK(JH1VEC-1+J)
             ENDDO
            ENDIF

            IF (IPRHAM.GE.7) THEN
              CALL HEADER(
     &'MAKE_DKH: Eigenvectors of (1 + R+ * R)',-1)
              WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
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

C           now generate Omega+ = (X+ * X)^(-1/2) in WORK(KH1AO)

            CALL QGEMM(NREF,NREF,NREF,D1,  ! U * d^-1/2 * U+
     &          'N','N',IPQTOQ(1,0),WORK(KAMAT1),NREF,NREF,NZ,
     &          'H','N',IPQTOQ(1,0),WORK(KAMAT),NREF,NREF,NZ,
     &          D0,IPQTOQ(1,0),WORK(KH1AO),NREF,NREF,NZ)

            IF (IPRHAM.GE.6) THEN
              CALL HEADER(
     &'MAKE_DKH: Omega+ = (1+ R+ * R)^(-1/2) ',-1)
              WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
              CALL PRQMAT(WORK(KH1AO),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
            ENDIF


            IF (IPRHAM.GE.6) THEN
              CALL HEADER(
     &  '  MAKE_DKH: Once again ( 1 R)',-1)
              WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
              CALL PRQMAT(WORK(KXMAT),NREF2,NREF,NREF2,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
            ENDIF

C           and multiply it to form X Omega+ = (1  R)t * Omega+:
            CALL QGEMM(NREF2,NREF,NREF,D1,
     &          'N','N',IPQTOQ(1,0),WORK(KXMAT),NREF2,NREF,NZ,
     &          'N','N',IPQTOQ(1,0),WORK(KH1AO),NREF ,NREF,NZ,
     &          D0,IPQTOQ(1,0),WORK(KAMAT),NREF2,NREF,NZ)

C      ... save the WORK(KAMAT), 4 "k" -> 2 "k" pict.ch.tr.mat, into TBUF !

           CALL DCOPY(2*N2REFQ,WORK(KAMAT),1,TBUF(IBUF),1)

C ... we have the final transformation matrix in "k" basis !
         IF (IPRHAM.GE.6) THEN
           CALL HEADER(
     &'GET_HINF_FPM: Final U1= X Omega+ transformation'//
     &' matrix in ON FP |k> basis',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(WORK(KAMAT),NREF2,NREF,NREF2,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
           CALL HEADER(
     &'GET_HINF_FPM: Dirac bare nucleus in |k> basis',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(WORK(KHD+JOFT1-1),NREF2,NREF2,NREF2,NREF2,
     &                 NZ,IPQTOQ(1,0),LUPRI)
         END IF

C ======================================================================
C ... do the trasformation - get the H2c Hamiltonian in the "k" basis
C               U1+ H1(in "k") U1
C                WORK(KAMAT)+ *   WORK(KHD)*WORK(KAMAT)
C        (remember NREF = NESH(I) and NREF2 = NORB(I) = 2*NREF)
C=====================================================================
            CALL QTRANS('AOMO','S',D0,NREF2,NREF2,NREF,NREF,
     &           WORK(KHD+JOFT1-1),NREF2,NREF2,NZ ,IPQTOQ(1,0),
     &           H2FP(IH2FP),NREF ,NREF ,NZ ,IPQTOQ(1,0),
     &           WORK(KAMAT),NREF2,NREF ,NZ, IPQTOQ(1,0),
     &           WORK(KAMAT),NREF2,NREF ,NZ, IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRHAM)

         IF (IPRHAM.GE.4) THEN
           CALL HEADER(
     &   'GET_HINF_FPM: Obtained U1+ H1 U1 = H2c_infinite'//
     &   ' in FP electronic |k> basis',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(H2FP(IH2FP),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
         END IF

C**************************************************************************
C
C     Do the stuff with the picture change transformation matrix
C
C    special case: spin-sfree "from beginning"
C
C**************************************************************************
         IF (GETPCTM) THEN
C   ... get the picture change transf. mtx 4c H1DBN "xhi"-> 2c H |k>
          IF (IPRHAM.GE.6) THEN
           CALL HEADER(
     &     'GET_HINF_FPM: Transf. mtx TMAT 4c "xhi"- 4c |k>',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &       '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TMAT(I2TMT(I)+1),NFBAS(I,0),
     &               NORB(I),NFBAS(I,0),NORB(I),
     &               NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

          CALL QGEMM(NFBAS(I,0),NESH(I),NORB(I),D1,
     &       'N','N',IPQTOQ(1,0),TMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),NZ,
     &       'N','N',IPQTOQ(1,0),WORK(KAMAT),NREF2,NREF,NZ,
     &       D0,IPQTOQ(1,0),PCTM(IPCTM),NFBAS(I,0),NESH(I),NZ)

           IF (IPRHAM.GE.5) THEN
             CALL HEADER(
     &    'GET_HINF_FPM: Transf. matrix PCTM, 4c "xhi"-> 2c |k_e>',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &       '*** Fermion corep ',I,'/',NFSYM
             CALL PRQMAT(PCTM(IPCTM),NFBAS(I,0),NESH(I),
     &                  NFBAS(I,0),NESH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
           ENDIF

C ========================================================================
C ...  now check the picture change from 4c BNDH "xhi" to the H2c in "k"
C ========================================================================
           IF (IPRHAM.GE.5) THEN
             CALL HEADER(
     & 'GET_HINF_FPM: Dirac bare nucl.Ham. in SA-AO basis "xhi"',-1)
             WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &       '*** Fermion corep ',I,'/',NFSYM
             CALL PRQMAT(HDAO(1+I2BASX(I,I)),NFBAS(I,0),NFBAS(I,0),
     &                  NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
           ENDIF

C  ... do the picture change transformation 4c Dirac Ham. "xhi" -> H2c in "k"
       CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),NREF,NREF,
     &           HDAO(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &           WORK(KA)   ,NREF ,NREF ,NZ ,IPQTOQ(1,0),
     &           PCTM(IPCTM),NFBAS(I,0),NREF ,NZ,IPQTOQ(1,0),
     &           PCTM(IPCTM),NFBAS(I,0),NREF ,NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRHAM)

      IF (IPRHAM.GE.5) THEN
        CALL HEADER(
     &  'GET_HINF_FPM: H2c obtained again in electronic |k>'//
     &  '  FP basis',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(WORK(KA),NREF,NREF,NREF,NREF,
     &               NZ,IPQTOQ(1,0),LUPRI)
           ENDIF

       IF (SF_BEG) THEN
        CALL SPF_HK(WORK(KA),.TRUE.,IBEIG(IORB(I)+1+NPSH(I)),
     &            IBEIG(IORB(I)+1+NPSH(I)),
     &             NREF,NREF,NREF,NREF,NZ,IPRHAM)

        IF (IPRHAM.GE.5) THEN
          CALL HEADER(
     &  'GET_HINF_FPM: H2c obtained again in electronic |k>'//
     &  ' basis after elimin. spin-orbit terms',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(WORK(KA),NREF,NREF,NREF,NREF,
     &               NZ,IPQTOQ(1,0),LUPRI)
        ENDIF


       ENDIF

C     ... calculate the difference for the irrep....
       CALL DAXPY(N2EQ,DM1,H2FP(IH2FP),1,WORK(KA),1) 
       DN =  DNRM2(N2EQ,WORK(KA),1)

C      ... do always print differences...
       WRITE(LUPRI,'(2X,A,I1,A,I1,A,D13.6)')
     &'Difference between 2 times constructed H2c in ON FP |k> basis'//
     &'(fsym ',I,'/',NFSYM,'):',DN

        ENDIF

C    ... offsets
        JOFT1 = JOFT1 + NORB(I)*NORB(I)*NZ
C    ... for the picture change transf.matrix
        IPCTM = IPCTM + NESH(I)*NFBAS(I,0)*NZ
C    ... update the offset for H2C Ham. in "k" basis
        IH2FP = IH2FP + NESH(I)*NESH(I)*NZ
C    ... offset for 4c "k" -> 2c "k" picture change transf.matrix
        IBUF = IBUF + 2*N2REFQ

C ... release memmory allocated in this cycle...
       CALL MEMREL('BSS 7',WORK,KWORK,KAMAT1,KFREE,LFREE)
       CALL MEMREL('BSS 6',WORK,KWORK,KH1AO,KFREE,LFREE)
       CALL MEMREL('BSS 5',WORK,KWORK,KB,KFREE,LFREE)
       CALL MEMREL('BSS 4',WORK,KWORK,KA,KFREE,LFREE)
       CALL MEMREL('BSS 3',WORK,KWORK,KIPVT,KFREE,LFREE)
       CALL MEMREL('BSS 2',WORK,KWORK,KXMAT,KFREE,LFREE)
       CALL MEMREL('BSS 1',WORK,KWORK,KAMAT,KFREE,LFREE)


   95 CONTINUE
C     end do 95 I = 1,NFSYM

C     ... release the rest of allocated memmory in this routine...
      CALL MEMREL('GET_HINF',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('GET_HINF_FPM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck R_CHECK */
      SUBROUTINE R_CHECK(ISYMM,HK,R,A,B,IPRHAM,WORK,LWORK)
C****************************************************************************
C
C     PURPOSE:  
C    ==========
C     Check the: 
C                a) A.R = B condition for getting R
C
C                b)  R-condition for BSS infinite order:
C
C                    H(pp)*R = -H(pe) + R*H(ee) + R*H(pp)*R
C
C    On input:  HK  - bare nucleus Hamiltonian matrix in the ON |k> basis
C                R  - the R matrix in |k> basis
C                A, B -  matrixes  A.R = B
C                ISYMM - the fermion symmetry (1 or 2) 
C
C   Output:  just check and print the result of the R-condition
C
C   Written:  M.Ilias, Odense, 2002,
C            last revison MI,february 2005,Strasburg
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00, DM1 = -1.00D00 )
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION HK(*),R(*),A(*),B(*),WORK(*)

      CALL QENTER('R_CHECK')
#include "memint.h"

C .... save here the fermion symmetry
      I = ISYMM

      IF (NPSH(I).NE.NESH(I)) CALL QUIT("CHECK_R:Ne<>Np!")
      IF (NORB(I).NE.(NESH(I)+NPSH(I))) CALL QUIT("CHECK_R:Ne<>Np!")


      IF ( IPRHAM .GE. 5 ) THEN
        CALL HEADER(
     & 'CHECK_R: Bare nucleus matrix in |k> basis',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(HK(1),NORB(I),NORB(I),NORB(I),NORB(I),
     &            NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     & 'CHECK_R: R matrix in |k> basis',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(R(1),NESH(I),NESH(I),2*NESH(I),NESH(I),
     &            NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     & 'CHECK_R: A matrix in |k> basis',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(A(1),NESH(I),NESH(I),NESH(I),NESH(I),
     &            NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     & 'CHECK_R: B matrix in |k> basis',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(B(1),NESH(I),NESH(I),NESH(I),NESH(I),
     &            NZ,IPQTOQ(1,0),LUPRI)

      END IF

C     ... this is the dimension...
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

      WRITE(LUPRI,'(2X,A,I1,A,I1,A,D13.6)')
     &'Norm of |A.R - B| (ifsym ',I,'/',NFSYM,') :',ARMBNULLNORM

      IF ( IPRHAM .GE. 5 ) THEN
        CALL HEADER(
     & 'CHECK_R: A * R - B zero matrix',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')'* Fermion irp ',I,'/',NFSYM
        CALL PRQMAT(WORK(K2),NESH(I),NESH(I),NESH(I),NESH(I),
     &               NZ,IPQTOQ(1,0),LUPRI)
      ENDIF


C .... check the Barysz "R-condition" equation

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
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(R(1),NESH(I),NESH(I),2*NESH(I),NESH(I),
     &              NZ,IPQTOQ(1,0),LUPRI)
        CALL HEADER(
     & 'CHECK_R: H(pp) part in |k> basis',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(HK(1),NESH(I),NESH(I),NORB(I),NORB(I),
     &         NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     & 'CHECK_R: Bare nucleus matrix - ep part in |k> basis',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
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
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
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
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
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
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
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
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
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

      WRITE(LUPRI,'(2X,A,I1,A,I1,A,D13.6)')
     & 'Norm of the R-condition, H(pp)*R+H(ep)-R*H(ee)-R*H(pe)*R'//
     & ' (ifsym ',I,'/',NFSYM,'):',RNULLNORM

      IF ( IPRHAM .GE. 5 ) THEN
       CALL HEADER(
     &   'CHECK_R: Final R-condition matrix '//
     &  ' H(pp)*R+H(ep)-R*H(ee)-R*H(pe)*R ',-1)
      WRITE(LUPRI,'(7X,A,I1,A,I1,3X,A,D17.12)')
     &' Fermion corep ',I,'/',NFSYM,
     & '* norm: ',RNULLNORM
        CALL PRQMAT(WORK(K1),NESH(I),NESH(I),NESH(I),NESH(I),
     &               NZ,IPQTOQ(1,0),LUPRI)
      ENDIF

C    ... release all used memmory in R_CHECK....
      CALL MEMREL('R_CHECK',WORK,KWORK,KWORK,KFREE,LFREE)
 
      CALL QEXIT('R_CHECK')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck mksaomo_le */
      SUBROUTINE MKSAOMO_Le(TBUF,SAOMO,TMAT,FACKIN,SPFR,
     &                      SPFR_BEG,GETPCTM,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     PURPOSE:
C    ==========
C        Construct half transformed overlap matrix, 
C
C         Saomo = Sao(LL) * C(L,el) * Ad(ee)^(-1)
C
C    On input:
C   -----------
C       TMAT - FP eigenvector "k" in "xhi" SA-AO basis
C       FACKIN  Ad(ee)
C
C       TBUF - 4c "k" to 2c "k" picture change transf.matrix, 
C              to be overwritten if desired...
C
C    On output:
C   ------------
C      SAOMO - transf.matrix form  2c SO-AO "xhi" -> 2c "k" NB: type AOMO !!!
C
C      TBUF - if spin-free from the beginning, it is the picture change transf.matrix
C                H4c "k" -> H2c "xhi_LL"
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
      DIMENSION SAOMO(*),TMAT(*),TBUF(*),WORK(LWORK),FACKIN(*)
      LOGICAL   SPFR, SPFR_BEG, GETPCTM

      CALL QENTER('MKSAOMO_Le')
#include "memint.h"
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

      IF (SPFR.AND.SPFR_BEG.AND.GETPCTM) THEN
C  ... array for the picture change transf. matrix 4c "k" -> 2c "xhi"
        ISI=NORB(1)*NFBAS(1,0)*NZ
        IF (NFSYM.EQ.2) ISI=ISI+(NORB(2)*NFBAS(2,0)*NZ)
        CALL MEMGET('REAL',KH4KH2AO,ISI,WORK,KFREE,LFREE)
        IBUF = 1
        IH4KH2AO = KH4KH2AO
      ENDIF

C
C   S_LL(NFBAS(I,0),NFBAS(I,0)) * TMAT_e(NFBAS(I,0),NESH(I))
C

      IOFF = 1
      JOFFESH = 0
      ITMT = 1
      ISAO = KSAO
      DO I = 1,NFSYM
      IF ( NESH(I) .GT. 0 ) THEN

         IF ( IPRINT .GE. 5 ) THEN
            CALL HEADER(
     &   'MKSAOMO_Le: Read SOVLP matrix(LL) in SO basis "xhi" ',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(WORK(ISAO),NFBAS(I,0),NFBAS(I,0),
     &                NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     &  'MKSAOMO_Le: Entering |k>_el matrix TMAT in SO basis "xhi" ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TMAT(ITMT+NFBAS(I,0)*NESH(I)),NFBAS(I,0),NESH(I),
     &                NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

C            SAOMO_Le = S_LL * TMAT

         CALL QGEMM(NFBAS(I,0),NESH(I),NFBAS(I,0),D1,
     &         'N','N',IPQTOQ(1,0),WORK(ISAO),
     &         NFBAS(I,0),NFBAS(I,0),1,
     &         'N','N',IPQTOQ(1,0),
     &         TMAT(ITMT+NFBAS(I,0)*NESH(I)),
     &         NFBAS(I,0),NESH(I),NZ,
     &         D0,IPQTOQ(1,0),SAOMO(IOFF),
     &         NFBAS(I,0),NESH(I),NZ)
C

            IF ( IPRINT .GE. 12) THEN
             CALL HEADER('MKSAOMO_Le: Right index transf. overlap'//
     &         ' matrix Saomo_Le',-1)
            WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
             CALL PRQMAT(SAOMO(IOFF),NFBAS(I,0),NESH(I),
     &            NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
             CALL FLSHFO(LUPRI)
            END IF

C            Now multiply by diagonal Ad^-1, FACKIN(j)  = Ad^2(j,j)
C               SAOMO_Le = SAOMO_Le * Ad^-1

             ICX = 0
             DO IZ = 1,NZ
             DO IESH = 1, NESH(I)
             DO IFBAS = 1, NFBAS(I,0) 
              SAOMO(IOFF+ICX) = 
     &        SAOMO(IOFF+ICX)/DSQRT(FACKIN(IESH+JOFFESH))
              ICX = ICX + 1
             ENDDO
             ENDDO 
             ENDDO 

         IF ( IPRINT .GE. 5 ) THEN
         CALL HEADER('MKSAOMO_Le: Final transf. matrix '//
     &   'H2c SA_AO "xhi"  -> H2c ON FP "k", SAOMO ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
            CALL PRQMAT(SAOMO(IOFF),NFBAS(I,0),NESH(I),
     &            NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
            CALL FLSHFO(LUPRI)
          END IF

C     ... take care of pict.change transf.matrixes if spinfree from "beginning"
        IF (SPFR.AND.SPFR_BEG.AND.GETPCTM) THEN

         IF ( IPRINT .GE. 5 ) THEN
          CALL HEADER('MKSAOMO_Le: Picture change transf. matrix'//
     &   ' H4c "k" -> H2c "k", TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IBUF),NORB(I),NESH(I),
     &           NORB(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
          CALL FLSHFO(LUPRI)
         END IF

C TBUF[NORB(I),NESH(I)] *  SAOMO [NFBAS(I,0),NESH(I)]+  =  WORK(IH4KH2AO)[NORB(I),NFBAS(I,0)]

         CALL QGEMM(NORB(I),NFBAS(I,0),NESH(I),D1,
     &         'N','N',IPQTOQ(1,0),TBUF(IBUF),
     &         NORB(I),NESH(I),NZ,
     &         'H','N',IPQTOQ(1,0),
     &         SAOMO(IOFF),NFBAS(I,0),NESH(I),NZ,
     &         D0,IPQTOQ(1,0),WORK(IH4KH2AO),
     &         NORB(I),NFBAS(I,0),NZ)
C

         IF (IPRINT.GE.5) THEN
          CALL HEADER('MKSAOMO_Le: Picture change transf. matrix'//
     &   ' H4c "k" -> H2c "xhi"',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(WORK(IH4KH2AO),NORB(I),NFBAS(I,0),
     &           NORB(I),NFBAS(I,0),NZ,IPQTOQ(1,0),LUPRI)
          CALL FLSHFO(LUPRI)
         END IF

C    ... offset for 4c "k" -> 2c "xhi" picture change transf.matrix
         IH4KH2AO = IH4KH2AO + NORB(I)*NFBAS(I,0)*NZ
C    ... offset for 4c "k" -> 2c "k" picture change transf.matrix, TBUF
        IBUF = IBUF + NORB(I)*NESH(I)*NZ

        ENDIF

      END IF

      JOFFESH = JOFFESH + NESH(I)
      IOFF = IOFF + NZ*NFBAS(I,0)*NESH(I)
      ITMT = ITMT + NFBAS(I,0)*NORB(I)
      ISAO = ISAO + NFBAS(I,0)*NFBAS(I,0)

      END DO

      IF (SPFR.AND.SPFR_BEG.AND.GETPCTM) THEN
C      ... save the 4c "k" -> 2c "xhi" picture change transform.matrix
        CALL DCOPY(ISI,WORK(KH4KH2AO),1,TBUF(1),1)
      ENDIF

      CALL MEMREL('SAOMOLe',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('MKSAOMO_Le')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck */      
       SUBROUTINE MAKE_H2C_AO(SAOMO,TMAT,TBUF,H2FP,H4AO,PCTM,IBEIG,
     &                        GETPCTM,SPFR,SPFR_BEG,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: Make the Hamiltonian in 2c AO basis
C    ========
C                                                         
C   On input: 
C      SAOMO - transf. matrix H2c "k" <- H2c "xhi_LL"
C      H2FP - H2c Ham. in the electronic FP "k" basis 
C            to be transformed to the SA-AO "xhi_LL" basis
C
C      PCTM - transf. matrix H4c "xhi" -> H2c "k"
C      H4AO - Dirac Ham. in "xhi"       
C                                                 
C      IBEIG - info on boson irreps
C
C
C     TBUF - temporary array, or H4c "k" -> H2c "xhi_LL" pict.change transformation matrix
C                               when spin-free from beginning
C
C     TMAT - temporary array, or H4c "xhi" -> H4c "k" transformation matrix, SF beginning
C                                                
C   On output:  
C     TMAT - H2c(BSSinf,DK2,DK1) in SA-AO "xhi" basis,
C
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C                                                         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, DM1 = -1.0D0, D0 = 0.0D0 ) 
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"
#include "dcbgen.h"
   
      DIMENSION SAOMO(*),TBUF(*),TMAT(*),PCTM(*),
     &          H4AO(*),H2FP(*),IBEIG(*),WORK(LWORK)    
      LOGICAL GETPCTM,SPFR,SPFR_BEG


      CALL QENTER('MAKE_H2C_AO')
#include "memint.h"         

C    Generate DKn  n'th order  H1(xhi_LL) matrix in H1AO ....
C         H1(xhi_L) = (C+)-1      H_DKn(k) C-1
C                   = Sao C Ad^-1 H_DKn(k) Ad^-1 C+ Sao
C       H_DKn = WORK(JBUF) ... H_ee in |k_e> basis

      CALL MEMGET('REAL',KH2AO,N2BBASXQ,WORK,KFREE,LFREE)

C     CALL DZERO(TMAT(1),N2BBASXQ)
      CALL DZERO(WORK(KH2AO),N2BBASXQ)

      IF (SPFR.AND.SPFR_BEG.AND.GETPCTM) THEN

       ISIZE1=NFBAS(1,0)*NORB(1)*NZ
       IF(NFSYM.EQ.2) ISIZE1=ISIZE1+(NFBAS(2,0)*NORB(2)*NZ)
       CALL MEMGET('REAL',KTMAT,ISIZE1,WORK,KFREE,LFREE)
C     ... save TMAT
       CALL DCOPY(ISIZE1,TMAT(1),1,WORK(KTMAT),1)
       ITMAT = KTMAT

       ISIZE2=NORB(1)*NORB(1)*NZ
       IF(NFSYM.EQ.2) ISIZE2=ISIZE2+(NORB(2)*NORB(2)*NZ)
       CALL MEMGET('REAL',KH4K,ISIZE2,WORK,KFREE,LFREE)
       IH4K = KH4K

C      ... starting offset for TBUF ....
       IBUF1 = 1
      ENDIF
C      ... needs cleaning for both SF_beg/SF_end !...
      CALL DZERO(TMAT(1),N2BBASXQ)

      IH2FP = 1
      ISAO  = 1
      IPCTM = 1
      IBUF  = 1
      IMAT  = 1
      DO 105 I = 1, NFSYM
      IF(NFBAS(I,0).EQ.0) GOTO 105

C    ... if wished, eliminate SO terms from H2FP(IH2FP) "at the end"
C  ********************************************************************
        IF (SPFR.AND..NOT.SPFR_BEG) THEN

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'MAKE_H2C_AO: H2c in |k> - before elimination '//
     &     ' of spin_orbit terms "AFTER"',-1)
           WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     &      '*** Fermion corep ',I,' /',NFSYM
            CALL PRQMAT(H2FP(IH2FP),NESH(I),NESH(I),NESH(I),NESH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

          CALL SPF_HK(H2FP(IH2FP),.TRUE.,IBEIG(IORB(I)+1+NPSH(I)),
     &            IBEIG(IORB(I)+1+NPSH(I)),
     &           NESH(I),NESH(I),NESH(I),NESH(I),NZ,IPRINT)

          IF (IPRINT.GE.5) THEN
            CALL HEADER(
     &     'MAKE_H2C_AO: H2c in |k> - eliminated'//
     &     ' spin_orbit terms "AFTER"',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     &     '*** Fermion corep ',I,' /',NFSYM
           CALL PRQMAT(H2FP(IH2FP),NESH(I),NESH(I),NESH(I),NESH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

        ENDIF

C       Row dimension of H2AO: only L [NFBAS(I,1)] or L+S [NFBAS(I,0)]?
C       currently we use L+S, to use only L we should delete small Xhi
C       after MAKE_DKH /5-Aug-2002 hjaaj

        IH2AO  = KH2AO+I2BASX(I,I)

        CALL QTRANS('MOAO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &               NESH(I),NESH(I),
     &            WORK(IH2AO),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            H2FP(IH2FP),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &            SAOMO(ISAO),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &            SAOMO(ISAO),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)

        IF (IPRINT.GE.4) THEN
          CALL HEADER(
     & 'MAKE_H2C_AO: FINAL H2c in SA-AO basis "xhi_LL"',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     & '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(WORK(IH2AO),NFBAS(I,0),NFBAS(I,0),
     &          NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C =================================================================================
C
C     ... prepare the total picture change transformation matrix, if wished,
C         and compare results of transformations....
C
C =================================================================================
        IF (GETPCTM) THEN

         IF (IPRINT.GE.6) THEN
           CALL HEADER(
     &'MAKE_H2C_AO: 4c bare nucleus in "xhi" AO basis',-1)
            WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(H4AO(1+I2BASX(I,I)),NFBAS(I,0),NFBAS(I,0),
     &                        NTBAS(0),NTBAS(0),
     &                        NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

         IF (IPRINT.GE.6) THEN
C  ... these two pict.change transf.matrixes are to be saved,
C      when given spin-free "from the end" !
          CALL HEADER(
     & 'MAKE_H2C_AO: Picture change transf.matrix '//
     & ' Dirac Ham. 4c "xhi" -> H 2c |k_e>',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(PCTM(IPCTM),NFBAS(I,0),NESH(I),
     &                   NFBAS(I,0),NESH(I),
     &                     NZ,IPQTOQ(1,0),LUPRI)
          CALL HEADER(
     &'MAKE_H2C_AO: Transf.matrix H2c |k_e> -> H2c "xhi_L"',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(SAOMO(ISAO),NFBAS(I,0),NESH(I),
     &        NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)

       ENDIF

       IF (.NOT.SPFR) THEN
C **** Full BSS ****
C      construct full picture change transfprmatin matrix,  PCTM*SAOMO+ = TBUF
      CALL QGEMM(NFBAS(I,0),NFBAS(I,1),NESH(I),D1, 
     &          'N','N',IPQTOQ(1,0),
     &           PCTM(IPCTM),NFBAS(I,0),NESH(I),NZ,
     &          'H','N',IPQTOQ(1,0),
     &           SAOMO(ISAO),NFBAS(I,0),NESH(I),NZ,
     &          D0,IPQTOQ(1,0),
     &          TBUF(IBUF),NFBAS(I,0),NFBAS(I,1),NZ)

C That we have: transformation matrix 4c AO "xhi" -> 2c AO "xhi"
       IF (IPRINT.GE.6) THEN
        CALL HEADER(
     &'  MAKE_H2C_AO: TBUF, final full picture change '//
     &'  transf. matrix H4c "xhi" -> H2c "xhi"',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TBUF(IBUF),NFBAS(I,0),NFBAS(I,1),
     &                         NFBAS(I,0),NFBAS(I,1),
     &                        NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C -----------------------------------------------------------------------------------
C  ..... when we have the picture change transf.matrix, do check the H2c in AO..
C              TBUF+ * WORK(IH2AO) * TBUF = TMAT
C         Now do U+ H4c(AO)  U = H2c(AO) and check it !
C .... do trasformation 4c bare nucleus AO "xhi" -> 2c AO "xhi"
C -----------------------------------------------------------------------------------
       CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                        NFBAS(I,1),NFBAS(I,1),
     &  H4AO(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &  TMAT(IMAT),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &  TBUF(IBUF),NFBAS(I,0),NFBAS(I,1),NZ,IPQTOQ(1,0),
     &  TBUF(IBUF),NFBAS(I,0),NFBAS(I,1),NZ,IPQTOQ(1,0),
     &    WORK(KFREE),LFREE,IPRINT)

        IF (IPRINT.GE.7) THEN
          CALL HEADER(
     &'MAKE_H2C_AO: H2c constructed (FULL picture change tranf.)'//
     & ' again in SA-AO basis "xhi_L", TMAT',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,1),NFBAS(I,1),
     &                        NTBAS(0),NTBAS(0),
     &                     NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

       ELSE IF (SPFR.AND..NOT.SPFR_BEG) THEN
C ****  spin-free at the end ****
C        ... PCTM+ . H4AO . PCTM = TBUF
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                        NESH(I),NESH(I),
     &  H4AO(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &  TBUF(IH2FP),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &  PCTM(IPCTM),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &  PCTM(IPCTM),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &    WORK(KFREE),LFREE,IPRINT)


C        ... control print out...
        IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &'MAKE_H2C_AO: H2c constructed (part. picture change tranf.)'//
     & ' again in FP ON "k" basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IH2FP),NESH(I),NESH(I),
     &                        NESH(I),NESH(I),
     &                     NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
C        ... eliminate SO terms from TBUF, H2c in "k"
        CALL SPF_HK(TBUF(IH2FP),.TRUE.,IBEIG(IORB(I)+1+NPSH(I)),
     &             IBEIG(IORB(I)+1+NPSH(I)),NESH(I),NESH(I),
     &             NESH(I),NESH(I),NZ,IPRINT)
C        ... control print out...
       IF (IPRINT.GE.7) THEN
        CALL HEADER(
     &'MAKE_H2C_AO: H2c again in FP ON "k" basis, TBUF,'//
     & ' after elimination of SO terms',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TBUF(IH2FP),NESH(I),NESH(I),
     &                       NESH(I),NESH(I),
     &                     NZ,IPQTOQ(1,0),LUPRI)
       ENDIF
C        ... SAOMO+ . TBUF . SAOMO = TMAT
        CALL QTRANS('MOAO','S',D0,NFBAS(I,1),NFBAS(I,1),
     &                        NESH(I),NESH(I),
     &  TMAT(IMAT),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &  TBUF(IH2FP),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &  SAOMO(ISAO),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &  SAOMO(ISAO),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &  WORK(KFREE),LFREE,IPRINT)

        IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &'MAKE_H2C_AO: H2c constructed again '//
     & ' - SPIN-FREE AT THE END -in SA-AO basis "xhi_L", TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IMAT),NFBAS(I,1),NFBAS(I,1),
     &                        NTBAS(0),NTBAS(0),
     &                     NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

       ELSE
C ****  spin-free at the beginning, utilize TMAT, TBUF ****
        IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &   'MAKE_H2C_AO: transformation matrix H4c "xhi"'//
     &   ' -> H4c "k"',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(WORK(ITMAT),NFBAS(I,0),NORB(I),
     &        NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
C     ... do the transformation WORK(ITMAT)+.H4AO.WORK(ITMAT)=WORK(IH4K)
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                        NORB(I),NORB(I),
     &  H4AO(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &  WORK(IH4K),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &  WORK(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &  WORK(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &   WORK(KFREE),LFREE,IPRINT)

        IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &   'MAKE_H2C_AO: Dirac Hamiltonian (before elim.SO terms) in'//
     &   ' FP ON basis "k"',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(WORK(IH4K),NORB(I),NORB(I),
     &          NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
C       ... eliminate SO components..
         CALL SPF_HK(WORK(IH4K),.TRUE.,IBEIG(IORB(I)+1),
     &                 IBEIG(IORB(I)+1),
     &                NORB(I),NORB(I),NORB(I),NORB(I),NZ,IPRINT )

        IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &   'MAKE_H2C_AO: Dirac Hamiltonian (AFTER elim.SO terms) in'//
     &   ' FP ON basis "k"',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(WORK(IH4K),NORB(I),NORB(I),
     &          NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'MAKE_H2C_AO: TBUF, picture change transf. matrix H4c "k"'//
     &   ' -> H2c "xhi"',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF1),NORB(I),NFBAS(I,0),
     &          NORB(I),NFBAS(I,0),NZ,IPQTOQ(1,0),LUPRI)

        ENDIF
        CALL DZERO(NFBAS(I,0)*NFBAS(I,0)*NZ,TMAT(IMAT) ) 
C       ... do the pict.change transformation, TBUF+.WORK(IH4K).TBUF = TMAT(IMAT)

        CALL QTRANS('AOMO','S',D0,
     &                    NORB(I),NORB(I),
     &                    NFBAS(I,0),NFBAS(I,0),
     &  WORK(IH4K),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &  TMAT(IMAT),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &  TBUF(IBUF1),NORB(I),NFBAS(I,0),NZ,IPQTOQ(1,0),
     &  TBUF(IBUF1),NORB(I),NFBAS(I,0),NZ,IPQTOQ(1,0),
     &  WORK(KFREE),LFREE,IPRINT)

        IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &'MAKE_H2C_AO: H2c constructed (part. picture change tranf.)'//
     & ' again in SA-AO basis "xhi_L",TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NFBAS(I,0),
     &                        NTBAS(0),NTBAS(0),
     &                     NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C      ... offset for 4c "xhi" -> 4c "k" picture change transf.matrix
        ITMAT = ITMAT + NFBAS(I,0)*NORB(I)*NZ
C      ... offset for 4c "k" -> 2c "xhi" picture change transf.matrix
        IBUF1 = IBUF1 + NFBAS(I,0)*NORB(I)*NZ
C      ... offset for H4c in "k" basis
        IH4K  = IH4K + NORB(I)*NORB(I)*NZ
       ENDIF

C .........................................................
C ... calculate the difference between H(2c) in AO "xhi"
C .........................................................
       ISIZE = NFBAS(I,0)*NFBAS(I,0)*NZ
       CALL DAXPY(ISIZE,DM1,WORK(IH2AO),1,TMAT(IMAT),1) 
       DN =  DNRM2(ISIZE,TMAT(IMAT),1)

C========================================================================
        WRITE(LUPRI,'(2X,A,I1,A,I1,A,D13.6)')
     & 'Difference between 2 times constructed'//
     & ' H2c in SA-AO basis (fsym ',I,'/',NFSYM,'): ',DN
C========================================================================

       ENDIF

C       ... update the offset for H2C Ham. in "k" basis
        IH2FP = IH2FP + NESH(I)*NESH(I)*NZ
C       ... update the offset for Saomo_Le * Ad(e)^-1 matrix
        ISAO  = ISAO + NZ*NFBAS(I,0)*NESH(I)
C       ... update the offset for the matrix 4c AO "xhi" -> 2c H2 "k"
        IPCTM = IPCTM  + NFBAS(I,0)*NESH(I)*NZ
C       ... offset for pict.change transf.matrix 4c AO "xhi" -> 2c  AO "xhi"
        IBUF = IBUF + NFBAS(I,0)*NFBAS(I,1)*NZ
C       ... offset for H2c in  AO "xhi", TMAT
        IMAT = IMAT + NFBAS(I,0)*NFBAS(I,0)*NZ

  105 CONTINUE

C    ... save the integrals for exporting them in TMAT
      CALL DCOPY(N2BBASXQ,WORK(KH2AO),1,TMAT(1),1)


C  SAVE THE INTEGRALS AND PICTURE CHANGE TRANSFORMATION MATRIX/MATRIXES, IF WISHED
C     ... include/dcbgen.h
      LUBSS = 35
C ...save hamiltonian matrix elements
      CALL OPNFIL(LUBSS,'BSSMAT','UNKNOWN','GMOTRA')
      CALL NEWLAB('H1AO_DK ',LUBSS,LUPRI)
      CALL WRITT(LUBSS,N2BBASXQ,WORK(KH2AO))
      WRITE(LUPRI,'(/2X,A)')
     &'BSS relativistic atomic integrals'//
     &' are written to the BSSMAT file.'
     
      IF (GETPCTM) THEN
C ...save the picture change transformation matrixes
       IF (.NOT.SPFR) THEN
        CALL NEWLAB('U_PICTRM',LUBSS,LUPRI)
        CALL WRITT(LUBSS,N2BBASXQ,TBUF(1))
        WRITE(LUPRI,'(2X,A)')
     &'BSS full picture change transformation matrix'//
     &' is written to the BSSMAT file.'
       ELSE IF (SPFR.AND..NOT.SPFR_BEG) THEN

        ISIZE = NFBAS(1,0)*NESH(1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NESH(2)*NZ

C         .... X4c-AO -> X2c-"k"
        CALL NEWLAB('PCTM_AOK',LUBSS,LUPRI)
        CALL WRITT(LUBSS,ISIZE,PCTM(1))

C         .... X2c-"k spin-free" -> X2c-AO
        CALL NEWLAB('PCTM_KAO',LUBSS,LUPRI)
        CALL WRITT(LUBSS,ISIZE,SAOMO(1))

C        ... we need to save also boson irreps of "k"
        CALL NEWLAB('PCTM_IBE',LUBSS,LUPRI)
        CALL WRITT(LUBSS,NTBAS(0),IBEIG(1))

        WRITE(LUPRI,'(2X,A)')
     &'BSS picture change transformation matrixes'//
     &' "spin-free AFTER" are written to the BSSMAT file.'

       ENDIF
      ENDIF

      CLOSE(LUBSS,STATUS='KEEP')

C ....   memmory deallocation
      CALL MEMREL('MAKE_BSS',WORK,KWORK,KWORK,KFREE,LFREE)


      CALL QEXIT('MAKE_H2C_AO')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck cmp_eigval */      
       SUBROUTINE CMP_EIGVAL(VMAT,H2AO,H1FP,EIG,TBUF,
     &                       IBEIG,SPFR,SPFR_BEG,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:
C    ========
C  Do comparison of the eigenvalues of the Dirac and H2c Hamiltonians
C   trasnform each to the ON basis and do the diagonalization
C
C  On input: H2AO - 2c Ham. in SA-AO "xhi" basis
C            H1FP - 4c Ham. in FP ON "k" basis with boson info in IBEIG
C                                                         
C            VMAT - Lowdin matrix - to be reorganized for H2c AO  
C                                                         
C  On output: VMAT reorganized for H2c                                                    
C                                                         
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0 ) 
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"
#include "dummy.h"
   
      DIMENSION H2AO(*),VMAT(*),H1FP(*),
     &          EIG(*),IBEIG(*),TBUF(*),WORK(LWORK)    
      LOGICAL SPFR,SPFR_BEG

      CALL QENTER('CMP_EIGVAL')
#include "memint.h"         

C   ... reorganize the Lowdin transformation matrix VMAT

      IF (IPRINT.GE.6) THEN
        CALL HEADER(
     &'GMOTRA: The total Lowdin matrix before BSS reordering',-1)
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

C   ... transform H2c AO through Lowdin matrix into ON MO ee basis
C   ... transform H4c AO through  matrix into ON MO  basis
      N1MO = NESH(1)*NESH(1)*NZ
      N2MO = NORB(1)*NORB(1)*NZ
 
      IF (NFSYM.EQ.2) THEN
       N1MO = N1MO + NESH(2)*NESH(2)*NZ
       N2MO = N2MO + NORB(2)*NORB(2)*NZ
      ENDIF

      CALL MEMGET('REAL',KH2MO,N1MO,WORK,KFREE,LFREE)

      IH2MO = KH2MO
      IVMAT = 1
      IH1FP = 1
      DO I = 1, NFSYM

       IF (IPRINT.GE.4) THEN
C    ... print out Lowdin matrix
        CALL HEADER(
     &'CMP_EIGVAL: Lowd. transf. matrix (Le) for 2c SCF procedure:',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(VMAT(IVMAT),NFBAS(I,0),NESH(I),NFBAS(I,0),
     &              NESH(I),1,IPQTOQ(1,0),LUPRI)
       ENDIF

C   ... transform H2c AO through Lowdin matrix into ON MO ee basis
        CALL QTRANS('AOMO','S',D0,NFBAS(I,1),NFBAS(I,1),
     &               NESH(I),NESH(I),
     &    H2AO(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &    WORK(IH2MO),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &    VMAT(IVMAT),NFBAS(I,0),NESH(I),1,IPQTOQ(1,0),
     &    VMAT(IVMAT),NFBAS(I,0),NESH(I),1,IPQTOQ(1,0),
     &        WORK(KFREE),LFREE,IPRINT)

       IF (IPRINT.GE.4) THEN
C    ... print out H2c matrix in AO "xhi" basis
        CALL HEADER(
     &'CMP_EIGVAL: Final H(2c)_LL in "xhi" SA-AO basis ' ,-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(H2AO(1+I2BASX(I,I)),NFBAS(I,1),NFBAS(I,1),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF


C =============================================================
C         ... remove spin-orbit dependent terms...
C  ... if spin-free "from the end"...
C =============================================================
       IF (SPFR.AND..NOT.SPFR_BEG) THEN

         IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &'CMP_EIGVAL: H4c in MO basis (before elim.SO terms) ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(H1FP(IH1FP),NORB(I),NORB(I),
     &             NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

         CALL SPF_HK(H1FP(IH1FP),.TRUE.,IBEIG(IORB(I)+1),
     &    IBEIG(IORB(I)+1),
     &    NORB(I),NORB(I),NORB(I),NORB(I),NZ,IPRINT )

       ENDIF

       IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &'CMP_EIGVAL: H2c in (MO_ee) ON Lowdin basis "Theta"',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(WORK(IH2MO),NESH(I),NESH(I),
     &             NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     &'CMP_EIGVAL: H4c in (MO) FP basis "k"',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(H1FP(IH1FP),NORB(I),NORB(I),
     &             NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF
        
C      ... do the diagonalization and compare eigenvalues ....

        CALL QDIAG(NZ,NORB(I),H1FP(IH1FP),NORB(I),NORB(I),   ! Diagonalize
     &             EIG(1),0,          ! Eigenvalues
     &             DUMMY,NORB(I),       ! Eigenvectors not necessary
     &             NORB(I),WORK(KFREE),LFREE,IERR)
        IF (IERR.NE.0) CALL QUIT('CMP_EIGVAL: FAILED QDIAG 1...')

        CALL QDIAG(NZ,NESH(I),WORK(IH2MO),NESH(I),NESH(I),   ! Diagonalize
     &             TBUF(1),0,          ! Eigenvalues
     &             DUMMY,NESH(I),       ! Eigenvectors not necessary
     &             NESH(I),WORK(KFREE),LFREE,IERR)
        IF (IERR.NE.0) CALL QUIT('CMP_EIGVAL: FAILED QDIAG 2...')


        CALL HEADER(
     &  'CMP_EIGVAL: Electronic eigenv. of H(2c)/H(4c) in ON basis',-1)
        WRITE(LUPRI,'(8X,A,I1,A,I1/)')
     &  '*** Fermion corep ',I,'/',NFSYM

       WRITE(LUPRI,'(3X,I5,F25.10)')
     &  ( J-NPSH(I)-1,EIG(J),J=1,NPSH(I))

       DO J =1, NESH(I)
        DIFF = TBUF(J) - EIG(NPSH(I)+J) 
        WRITE(LUPRI,'(3X,I5,2F25.10,3X,D12.5)') 
     &    J,
     &    EIG(NPSH(I)+J),
     &    TBUF(J),
     &    DIFF
       ENDDO
          IH2MO = IH2MO + NESH(I)*NESH(I)*NZ 
C       ... offset for Lowdin matrix
          IVMAT = IVMAT + NFBAS(I,0)*NESH(I)
          IH1FP = IH1FP + NORB(I)*NORB(I)*NZ
      ENDDO

      CALL MEMREL('CMP_EIGVAL',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('CMP_EIGVAL')

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_hdk2_fpm */      
      SUBROUTINE GET_HDK2_FPM(H2FP,TBUF,EIG,IBEIG,
     &                     DK2SO1,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:
C    ========
C  Do comparison of the eigenvalues of the Dirac and H2c Hamiltonians
C   trasnform each to the ON basis and do the diagonalization
C
C  On input:  TBUF H1 in "k" basis
C             EIG - free particle eigenvalues needed for denominators
C             IBEIG - boson irreps if spin-free, or DK2SO1
C                                                         
C  On output: H2FP(*) - the desired DK2 Hamiltonian in "k" basis
C                                                         
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0, D2 = 2.0D0, DM1 = -1.0D0 ) 
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dcbgen.h"   
#include "dgroup.h"
   
      DIMENSION H2FP(*),TBUF(*),EIG(*),IBEIG(*),WORK(LWORK)

      LOGICAL SPFR_BEG, DK2SO1

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
     &  'MAKE_DKH: DK2 Hamiltonian to be generated',-1)
      ENDIF

      CVAL2 = CVAL*CVAL
      JOFT1 = 1
      JBUF  = 1
      DO  I = 1,NFSYM  ! Run big cycle over the fermion symmetries
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
             EIG_P = ! Calculate the average value
     & ((-EIG(NORB(I)-IP+1+IORB(I))-(D2*CVAL2)) + 
     &   EIG(IP+IORB(I)))/D2

             EIG_E = ! Calculate the average value
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
     & 'DK2SO1 - elim. spin-orbit effects from H1(ep) and R(pe)'
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
C         2) Generate  "X" = H_1(ep) R_DK2(pe) in H2FP(JBUF)
C

          CALL QGEMM(NESH(I),NESH(I),NPSH(I),D1,
     &          'N','N',IPQTOQ(1,0),TBUF(JH1EP),NORB(I),NORB(I),NZ,
     &          'N','N',IPQTOQ(1,0),TBUF(JR2PE),NORB(I),NORB(I),NZ,
     &          D0,IPQTOQ(1,0),H2FP(JBUF), NESH(I),NESH(I),NZ)
C
C         3) Generate  X = 0.5 ( X + X(dagger)) 
C
          CALL QGETHE(H2FP(JBUF),NESH(I),NZ,NESH(I),NESH(I))

          IF (IPRINT.GE.7) THEN
            CALL HEADER(
     &'MAKE_DKH: 2. order correction to H1 in "|k_e>" basis',-1)
            WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
            CALL PRQMAT(H2FP(JBUF),NESH(I),NESH(I),NESH(I),
     &                  NESH(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
*  Add to first order (ee-blocks!): H_1(ee) =  H_1(ee) + X

          DO IZ = 1,NZ
            IBUF  = JBUF + (IZ-1)*N2ESHI
            JH1EE = JOFT1 + NPSH(I)*NORB(I) + NPSH(I)
     &            + (IZ-1)*NORB(I)*NORB(I)
            DO J = 1, NESH(I)
              CALL DAXPY(NESH(I),D1,TBUF(JH1EE),1,H2FP(IBUF),1)
              JH1EE = JH1EE + NORB(I)
              IBUF  = IBUF  + NESH(I)
            END DO
          END DO
C         ... update the offsets...
          JBUF  = JBUF  + NESH(I)*NESH(I)*NZ
          JOFT1 = JOFT1 + N2T*NZ
      ENDIF
      ENDDO


      CALL QEXIT('GET_HDK2_FPM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_hdk1_fpm */      
      SUBROUTINE GET_HDK1_FPM(H2FP,TBUF,IBEIG,ONLYSO1,
     &                         IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:
C    ========
C  Do comparison of the eigenvalues of the Dirac and H2c Hamiltonians
C   trasnform each to the ON basis and do the diagonalization
C
C  On input:  TBUF H1 in "k" basis
C             IBEIG - boson irreps if pure spin-orbit integrals (ONLYSO1=.true.)
C                                                         
C  On output: H2FP(*) - the desired DK2 Hamiltonian in "k" basis
C                                                         
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"   

      DIMENSION H2FP(*),TBUF(*),IBEIG(*),WORK(LWORK)
      LOGICAL ONLYSO1

      CALL QENTER('GET_HDK1_FPM')
#include "memint.h"         
  
C  ... store into WORK(IBUF) the ee-part of the full H1 Hamiltonian
      JBUF  = 1
      JOFT1 = 1
      DO I = 1, NFSYM
      IF (NESH(I).GT.0) THEN
        N2ESHI = NESH(I)*NESH(I)
        N2T    = NORB(I)*NORB(I)
        DO IZ = 1, NZ
          IBUF  = JBUF + (IZ-1)*N2ESHI
          JH1EE = JOFT1 + NPSH(I)*NORB(I) + NPSH(I)
     &           + (IZ-1)*N2T
           DO J = 1, NESH(I)
             CALL DCOPY(NESH(I),TBUF(JH1EE),1,H2FP(IBUF),1)
             JH1EE = JH1EE + NORB(I)
             IBUF  = IBUF  + NESH(I)
           ENDDO
        ENDDO
C      ... if desired, remove scalar terms to have only SO1 integrals...
        IF (ONLYSO1) THEN
          CALL SPF_HK(H2FP(JBUF),.FALSE.,IBEIG(IORB(I)+1+NPSH(I)),
     &               IBEIG(IORB(I)+1+NPSH(I)),
     &              NESH(I),NESH(I),NESH(I),NESH(I),NZ,IPRINT)
        ENDIF
C       ... update the offsets...
        JBUF  = JBUF  + N2ESHI*NZ
        JOFT1 = JOFT1 + N2T*NZ
      ENDIF 
      ENDDO

      CALL QEXIT('GET_HDK1_FPM')
      RETURN
      END
