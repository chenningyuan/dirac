!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck onegen */
      SUBROUTINE ONEGEN(WORK,LWORK)
C*****************************************************************************
C
C     This routine goes through the list of one-electron integrals specified
C     in the COMMON block XCBCLS. 
C
C     Call other routines for integrals generation.
C
C     Called from: MINPRP
C                  PAMSET
C                  LSQMAT
C
C     Written by T.Saue - May 1996
C
C*****************************************************************************
! checking
       use fde_types
       use fde_mod
!!     use FDE_ImportPotentials

#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
      DIMENSION WORK(LWORK)
      CHARACTER COMB*4
      LOGICAL TEST
      type(fde_import) :: itmp
C
#include "dcbgen.h"
#include "dcbprl.h"
#include "dcbham.h"
#include "dcbcls.h"
#include "cbihr1.h"
#include "cbiher.h"
#include "nuclei.h"
      CALL QENTER('ONEGEN')
C.....With modified Hamiltonian this routine is skipped for now
      IF(MDIRAC) GOTO 20
      !mi: always print nuclear contribution to the dipole moment (costs nothing)
      CALL DIPNUC(1,.FALSE.)
C
C     Loop through list of class labels; calculate only the active ones(CLSCAL)
C     =========================================================================
C
      DO 10 I = 1,NPRPCLS
      IF(.NOT.CLSCAL(I)) GOTO 10
        IF (IPRONE.GE.2) WRITE(LUPRI,'(/2A,2X,A)')
     &       '* ONEGEN: class operator and class comb ',
     &       CLSINT(I), CLSCMB(I)
Chj     CLSCMB and PDOINT is now '+', '-', or '0',
Chj     corresponding to a factor of +1, -1, or 0 on that block.
C       The four blocks (AO basis): LL, SL, LS, SS
        DO J = 1,4
          IF (CLSCMB(I)(J:J) .EQ. '0') THEN
             COMB(J:J) = 'F'
          ELSE
             COMB(J:J) = 'T'
          END IF
        END DO  
C
C     Nuclear attraction integrals (molecular field)
C
      IF    (CLSINT(I).EQ.'MOLFLD ') THEN
        IF (NOSMLV) THEN
C          ... zero SS integrals
           CALL NOSSNUCATT
           COMB(4:4) = 'F'
           IF(IPRONE.GE.3) THEN
               WRITE(LUPRI,'(2X,A,A4)')
     &      'FYI: no SS-nuclear attraction integrals'// 
     &      '(MOLFLD), adapted COMB=',COMB
           ENDIF
        END IF
        CALL NUCATT(WORK,LWORK,COMB,IPRONE)
C
C     Effective electronic density on nuclei 
C
      ELSEIF (CLSINT(I).EQ.'EFFDEN ') THEN
        CALL EFFECTIVE_DENS_NUC(WORK,LWORK,COMB)
C
C     Contribution from (static) embedding potential
C
      ELSEIF (CLSINT(i).EQ.'FDEVEMB') THEN 
         IF (NOSMLV) THEN
C          ... zero SS integrals
           CALL NOSSNUCATT
           COMB(4:4) = 'F'
         END IF
        call fde_get_import_info(itmp)
        if (itmp%im_vemb) then
            CALL FDE_SaveAOPROPER_StaticEmbPot(COMB,IPRONE)
         endif
C
C     Contribution from beta-matrix
C
      ELSEIF(CLSINT(I).EQ.'BETAMAT') THEN
        CALL BETAMT(COMB,WORK,LWORK,IPRONE)
C
C     Contribution from CAP integrals 
C
      ELSEIF(CLSINT(I).EQ.'CAP_RE '.OR.CLSINT(I).EQ.'CAP_IM ' .OR.
     &       CLSINT(I).EQ.'CAP_OVL'.OR.CLSINT(I).EQ.'CAPD1R ') THEN
CMI     ... this serves only for 'catching' CAP... labels because
CMI     these integrals are not computed in the abacus part but are from external source
       CONTINUE

CMI  Contribution from Heff_EDM integrals
!     ELSEIF (CLSINT(I).EQ.'EDM') THEN
!       CALL EDM_EFF_OPER(WORK,LWORK,IPRONE)

      ELSE
C
C Process all other integrals, which are coded in the abacus part !
C
        TEST = ALLATM.AND.
     &         ((CLSINT(I).EQ.'NEFIELD').OR.
     &          (CLSINT(I).EQ.'ELFGRDC').OR.
     &          (CLSINT(I).EQ.'RM1N1H ').OR.
     &          (CLSINT(I).EQ.'NUCSNLO').OR.
     &          (CLSINT(I).EQ.'FERMI C').OR.
     &          (CLSINT(I).EQ.'PVC    ').OR.
     &          (CLSINT(I).EQ.'NSTCGO '))
!       write(lupri,*) 'Test output of CLSINT, TEST, NUCIND',
!                      CLSINT(I),TEST,NUCIND
        IF(TEST) THEN
          ALLATM = .FALSE.
          NPATOM = 1
          DO K = 1,NUCIND
            IPATOM(1)=K
            CALL OP1GEN(CLSINT(I),COMB,IORDCL(I),WORK,LWORK,IPRONE)
          ENDDO
          ALLATM = .TRUE.
          IPATOM(1) = 0
          NPATOM    = 0
        ELSE
          CALL OP1GEN(CLSINT(I),COMB,IORDCL(I),WORK,LWORK,IPRONE)
        ENDIF
      ENDIF
      CLSCAL(I) = .FALSE.
   10 CONTINUE
C
      CLOSE(LU1INT,STATUS='KEEP')
 20   CONTINUE
      CALL QEXIT('ONEGEN')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EFFECTIVE_DENS_NUC(WORK,LWORK,COMB)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Generates integrals used for the calculation of the effective electronic
C density at nuclei.  
C
C Called from: ONEGEN
C 
C Written by Robert van Meer, Sept. 2010, last edited Mar. 2011 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "cbihr1.h"
#include "cbiher.h"
#include "dcbham.h"
#include "nuclei.h"
#include "dcbbas.h"
#include "dcbprl.h"
#include "dcbcls.h"
#include "dcbxpr.h"
C     GEXPEFFDE2 contains nuclear gaussian exponents related to size change
C     EXPPLC = Exponent placeholder
      REAL*8 GEXPEFFDE2, EXPPLC, PREFACTOR 
      CHARACTER RTNLBL(2)*8,COMB*4,OMITVNUC(2)*4, 
     &          LABINT(15*MXCENT)*8
      LOGICAl DOINT(4), DOATOM(NUCIND)
      INTEGER ITYP
      real*8 prefac2(4)
      DIMENSION WORK(LWORK), GEXPEFFDE2(4), INTREP(9*MXCENT)
#include "memint.h"
      CALL QENTER('EFFECTIVE_DENS_NUC')

C     Memory allocation
      CALL MEMGET('REAL',KONEIN,4*NNBBASX,WORK,KFREE,LFREE)

      IF (BSS.OR.x2C.OR.ZORA.OR.LEVYLE) THEN
        WRITE(LUPRI,*) 'EFFECTIVE_DENS_NUC: BSS, x2C, ZORA and LEVY'//
     &  '-LEBLOND modes not supported'
        CALL QUIT('EFFECTIVE_DENS_NUC: not supported for current'//
     &  ' Hamiltonian')
      END IF

      READ(COMB,'(4L1)',ERR=1000) (DOINT(I),I=1,4)
      OMITVNUC(1) = 'FFFF'
      OMITVNUC(2) = 'FFFF'
C     Fetch and set labels/variables (-99 is used to turn on .SELECT)
      CALL SETATM(DOATOM,NATOM,-99)
      CALL EFFDE2TYP(NOPTYP,INTREP,LABINT,DOATOM,NATOM)
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = 'SY'
      WRITE(RTNLBL(2)(3:4),'(I2)') 1
      RTNLBL(2)(5:8) = COMB
      ITYP = 0

C     Get the integrals, exponent is changed for a single nucleus 
C     and reset to the original value after each calculation  
      DO K = 1,NUCIND
        IF (DOATOM(K)) THEN
C         Generate nuclear exponents for nuclear size displacments and prefactor 
          CALL EFF_DEN_NUC_EXP(GEXPEFFDE2, K, PREFACTOR)
          EXPPLC = GNUEXP(K)
C         Calculate property matrices for nuclear displacements
          DO L = 1, 4 
            GNUEXP(K) = GEXPEFFDE2(L)
C           Calculate the integrals  
            CALL NUCAT1(WORK(KONEIN+(L-1)*NNBBASX),WORK(KFREE),LFREE,
     &                  DOINT,OMITVNUC,IPRHAM)
!                                      iprint is not defined - stefan april 2011
!    &                  DOINT,OMITVNUC,IPRINT)
          ENDDO
          GNUEXP(K) = EXPPLC
C         Create total property matrix with correct prefactor
          CALL MEMGET('REAL',KHELP,NNBBASX,WORK,KFREE,LFREE)
          prefac2(1) = +1.0d0
          prefac2(2) = -8.0d0
          prefac2(3) = +8.0d0
          prefac2(4) = -1.0d0
          call dscal(NNBBASX,PREFACTOR*prefac2(1),WORK(KONEIN),1)
          do l = 2, 4
          call dcopy(NNBBASX,WORK(KONEIN+(L-1)*NNBBASX),1,work(KHELP),1)
          call daxpy(NNBBASX,PREFACTOR*prefac2(l),work(KHELP),1,
     &               WORK(KONEIN),1)
          end do
!         WORK(KONEIN:KONEIN+NNBBASX-1) = PREFACTOR *
!    &    (
!    &            WORK(KONEIN          :KONEIN+  NNBBASX-1) 
!    &      - 8 * WORK(KONEIN+  NNBBASX:KONEIN+2*NNBBASX-1)
!    &      + 8 * WORK(KONEIN+2*NNBBASX:KONEIN+3*NNBBASX-1)
!    &      -     WORK(KONEIN+3*NNBBASX:KONEIN+4*NNBBASX-1)
!    &    )
          CALL MEMREL('EFFECTIVE_DENS_NUC',WORK,KWORK,KHELP,KFREE,LFREE)
C         Store matrix in AOproper file 
          ITYP = ITYP + 1
          CALL WRTPRO(WORK(KONEIN),NNBBASX,LABINT(ITYP),RTNLBL,IPRHAM)
        ENDIF
      ENDDO

C     Memory deallocation
      CALL MEMREL('EFFECTIVE_DENS_NUC',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('EFFECTIVE_DENS_NUC')
      RETURN
      
 1000 CONTINUE
C     Not able to read DOINT information from COMB
      WRITE(LUPRI,'(A,A)') 'EFFECTIVE_DENS_NUC: Not able to read COMB ='
     & ,COMB
      CALL QUIT('EFFECTIVE_DENS_NUC: Not able to read COMB')
      RETURN
      END 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EFF_DEN_NUC_EXP(GEXPEFFDE2, NUCLEUS, PREFACTOR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Fills Mossbauer exponent array with exponent values
C  
C Called from: MOSSBAUER_ISOMER
C 
C Written by Robert van Meer, Sept. 2010 last edited Mar. 2011
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "pi.h"
      REAL*8 GEXPEFFDE2, RADIUS
      DIMENSION GEXPEFFDE2(4)
C     Generate nuclear radius from the gaussian exponent    
      RADIUS = SQRT(3D0 / (2D0 * GNUEXP(NUCLEUS)) )
C     Generate new gaussian parameters dR=0.005R
      GEXPEFFDE2(1) = 3D0 / (2D0 * 0.990D0 * 0.990D0 * RADIUS * RADIUS) 
      GEXPEFFDE2(2) = 3D0 / (2D0 * 0.995D0 * 0.995D0 * RADIUS * RADIUS)
      GEXPEFFDE2(3) = 3D0 / (2D0 * 1.005D0 * 1.005D0 * RADIUS * RADIUS)
      GEXPEFFDE2(4) = 3D0 / (2D0 * 1.010D0 * 1.010D0 * RADIUS * RADIUS)  
C     Prefactor = 3 / (48*Z*Pi*R*dR) 
      PREFACTOR = 12.5D0 / (RADIUS * RADIUS * CHARGE(NUCLEUS) * PI)
      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE NOSSNUCATT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Delete marking SS-block of in nucl.attraction integrals as empty
C
C Called from: PAMSET (when restart from AOPROPER)
C              ONEGEN
C when NOSMLV=.true.
C
C Written by Miro Ilias, Dec. 2008, Tel Aviv
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      LOGICAL MOLFLD_FOUND
#include "dcbprl.h"
#include "dcbcls.h"
#include "dcbxpr.h"

      MOLFLD_FOUND=.FALSE.

      DO I=1,NPRPCLS
        IF (CLSINT(I).EQ.'MOLFLD ') THEN
          MOLFLD_FOUND=.TRUE.
          CLSCMB(I)(4:4) = '0'
          ILBL = IPRPLBL(1,I)
          PDOINT(ILBL)(4:4) = '0'
        ELSE IF (CLSINT(I).EQ.'FDEVEMB') THEN
          CLSCMB(I)(4:4) = '0'
          ILBL = IPRPLBL(1,I)
          PDOINT(ILBL)(4:4) = '0'
        ENDIF
      ENDDO

      IF (MOLFLD_FOUND) THEN
        WRITE(LUPRI,'(2X,A,A4)')
     & 'The SS-block of the potential energy matrix and (if present)'//
     & ' embedding potential are neglected (Levy-Leblond Hamiltonian)'
      ELSE
        WRITE(LUPRI,'(2X,A)') 'Nuclear attraction integrals not'//
     &   ' found in operator list!'
        WRITE(LUPRI,*) 'NOSSNUCATT routine was supposed to zero'//
     &   'SS-block of them !'
        CALL QUIT('NOSSNUCATT: end here')
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck op1gen */
      SUBROUTINE OP1GEN(WORD,COMB,IORDER,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Driver for the generation of one-electron integrals
C    to be placed on the file AOPROPER.
C
C     On input:  WORD     ... 
C                COMB     ... integral characteritics
C                IORDER   ...
C
C     Called from ONEGENdirone.F
C
C     Written by Trond Saue - May 1996
C     Last revisions:  May 19 1996
C                  MI(+TEC): minor changes for LAO based RPA, dec.2002
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "dcbgen.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbiher.h"
      LOGICAL DOINT(4),TOFILE
      CHARACTER WORD*7,COMB*4,LABINT(5*MXCOOR)*8,MTFORM*6
      DIMENSION WORK(LWORK)
C
      CALL QENTER('OP1GEN')
#include "memint.h"
C     Write integrals to file
      TOFILE   = .TRUE.
      NCOMP    = 0
      MAXTYP = 5*MXCOOR
      IF ((WORD.EQ.'ELFGRDC').OR.(WORD.EQ.'ELFGRDS')) THEN
         NINTAD = 9*NUCIND*(MAXREP+1)
      ELSE IF (WORD.EQ.'DSO    ') THEN
CMI ... to be fixed for the future LAO implementation...
         NINTAD = 9*NUCDEP*NUCDEP
      ELSE
         NINTAD = MAXTYP
      END IF
      READ(COMB,'(4L1)',ERR=1000) (DOINT(I),I=1,4)
C
C     Memory allocation
C
      CALL MEMGET('INTE',KINTRP,MAXTYP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINTAD,NINTAD,WORK,KFREE,LFREE)
C
C  extra assignement of characteristics for very narrow group of integrals
C
      NCOMP = 0
      IF (WORD.EQ.'S1MAGR ') THEN 
         TRIANG = .FALSE.
         MTFORM = 'SQUARE'
      ELSE
         TRIANG = .TRUE.
         MTFORM = 'TRIANG'
      END IF

C
C     Generate integrals in abacus
C
C     ... calling of ../abacus/her1pro.F
      CALL PR1INT_1(WORK(KFREE),LFREE,WORK(KINTRP),WORK(KINTAD),
     &            LABINT,WORD,IORDER,NPQUAD,TRIANG,
     &            PROPRI,IPRINT,DUMMY,NCOMP,TOFILE,MTFORM, 
     &            DOINT)
C
C     Memory deallocation
C
      CALL MEMREL('OP1GEN',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('OP1GEN')
C
      RETURN
 1000 CONTINUE
C
C    Not able to read DOINT information from COMB
C
      WRITE(LUPRI,'(A,A)') 'OP1GEN: Class label:',WORD
      WRITE(LUPRI,'(A,A)') 'OP1GEN ERROR: Not able to read COMB =',COMB
      CALL QUIT('OP1GEN: Not able to read COMB')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck betamt */
      SUBROUTINE BETAMT(COMB,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Driver for the generation of beta-matrix integrals
C     to be placed on file AOPROPER.
C
C     Written by Trond Saue - May 1996
C     Last revision May 19 1996
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dummy.h"
      LOGICAL DOINT(4)
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4
      DIMENSION WORK(LWORK)
#include "dcbgen.h"
#include "dcbbas.h"
C
      CALL QENTER('BETAMT')
#include "memint.h"
C
C     Initialization
C
      NCOMP    = 1
      READ(COMB,'(4L1)',ERR=1000) (DOINT(I),I=1,4)
C
C     Memory allocation
C
      CALL MEMGET('REAL',KONEIN,NNBBASX,WORK,KFREE,LFREE)
C
C     Generate integrals
C
      CALL RGETINT('OVERLAP',WORK(KONEIN),DOINT,NCOMP,IDUMMY,IDUMMY,
     &            'TRIANG',.FALSE.,IPRINT,WORK(KFREE),LFREE)
C
      IF(IPRINT.GE.6) THEN
        CALL HEADER('Small component overlap integrals',-1)
        CALL OUTPAK(WORK(KONEIN),NTBAS(0),1,LUPRI)
      ENDIF
C
C     Generate integral labels
C
      LABEL = 'BETAMAT '
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = 'SY'
      WRITE(RTNLBL(2)(3:4),'(I2)') 1
      RTNLBL(2)(5:8) = COMB
C
C     Write integrals to file
C
      CALL WRTPRO(WORK(KONEIN),NNBBASX,LABEL,RTNLBL,IPRINT)
C
C     Memory deallocation
C
      CALL MEMREL('BETAMT',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('BETAMT')
C
      RETURN
 1000 CONTINUE
C
C     Not able to read DOINT information from COMB
C
      WRITE(LUPRI,'(A,A)') 'BETAMT ERROR: Not able to read COMB =',COMB
      CALL QUIT('BETAMT: Not able to read COMB')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck nucatt */
      SUBROUTINE NUCATT(WORK,LWORK,COMB,IPRINT)
C*****************************************************************************
C
C     Generate nuclear attraction integrals
C
C     IPRINT is IPRONE
C
C     Called from: ONEGEN
C
C     Written by T.Saue May 1995
C     Last revision May 19 1996 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      LOGICAl DOINT(4)
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4,OMITVNUC(2)*4
      DIMENSION WORK(LWORK)
#include "dcbbas.h"
C
      CALL QENTER('NUCATT')
#include "memint.h"
C     Memory allocation
      CALL MEMGET('REAL',KONEIN,NNBBASX,WORK,KFREE,LFREE)
C
C     Generate integrals
C
      READ(COMB,'(4L1)',ERR=1000) (DOINT(I),I=1,4)
      OMITVNUC(1) = 'FFFF'
      OMITVNUC(2) = 'FFFF'
      CALL NUCAT1(WORK(KONEIN),WORK(KFREE),LFREE,
     &            DOINT,OMITVNUC,IPRINT)
      IF(IPRINT.GE.4) THEN
        CALL HEADER('Nuclear attraction Integrals',-1)
        CALL OUTPAK(WORK(KONEIN),NTBAS(0),1,LUPRI)
      ENDIF
C
C     Generate integral labels
C
      LABEL = 'MOLFIELD'
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = 'SY'
      WRITE(RTNLBL(2)(3:4),'(I2)') 1
      RTNLBL(2)(5:8) = COMB
C
C     Write integrals to file
C
      CALL WRTPRO(WORK(KONEIN),NNBBASX,LABEL,RTNLBL,IPRINT)
C
C     Memory deallocation
      CALL MEMREL('NUCATT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('NUCATT')
      RETURN
 1000 CONTINUE
C
C     Not able to read DOINT information from COMB
C
      WRITE(LUPRI,'(A,A)') 'NUCATT: Not able to read COMB =',COMB
      CALL QUIT('NUCATT: Not able to read COMB')
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck nucat1 */
      SUBROUTINE NUCAT1(OP1INT,WORK,LWORK,DOINT,OMITVNUC,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
C
      LOGICAL DOINT(2,2)
      CHARACTER*4 OMITVNUC(2)
      DIMENSION OP1INT(NNBBASX),WORK(LWORK)
#include "dcbbas.h"
#include "nuclei.h"
#include "cbiher.h"
#include "cbihr1.h"
#include "cbisol.h"
#include "dcbgen.h" 
#include "dcbham.h"
C
      CALL QENTER('NUCAT1')
#include "memint.h"
C
      CALL MEMGET('REAL',KSTHMA,3*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFACIN,2*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      LWRK   = LWORK - KFREE + 1
C
      CALL ONEDR1(WORK(KSTHMA),DUMMY,DUMMY,WORK(KFACIN),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KNCENT),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),DUMMY,WORK(KFREE),LWRK,
     &            IPRINT,.FALSE.,0,.TRUE.,.TRUE.,.FALSE.,NTBAS(0),
     &            NNBBASX,NNBBAST,DOINT,OMITVNUC,.FALSE.)
C
C     Symmetry unpack nuclear attraction integrals
C
      KNUCS = KSTHMA + NNBBASX

#ifdef MOD_ECP
c
C      *********************************************
C      **                                         **
C      **      this part adds ARGOS ARPP          **
C      **    nuc-elec attraction integrals.       **
C      **                                         **
C      *********************************************
      IF (ECPCALC) then
         OPEN(59,FORM='UNFORMATTED',FILE='RECP_INT_C')
         CALL MEMGET('REAL',KARPP,NNBBAST,WORK,KFREE,LFREE)
         CALL READT(59,NNBBAST,WORK(KARPP))
         close(59)
         DO I=0,NNBBAST-1
            WORK(KNUCS+I)=WORK(KNUCS+I)+WORK(KARPP+I)
         ENDDO
         CALL MEMREL('NUCAT1.ECPCALC',WORK,KWORK,KARPP,KFREE,LFREE)
      ENDIF
#endif
      CALL SYMUPK(WORK(KNUCS),OP1INT,1,NNBBASX)

C
C     Memory deallocation
C
      CALL MEMREL('NUCAT1.ONEDR1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('NUCAT1')
      RETURN
C
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EDM_EFF_OPER(WORK,LWORK,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  prepare imaginary, antisymmetric KINENERG integrals into LS/SL blocks,
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
      DIMENSION WORK(LWORK)
!     ... local variables
      LOGICAl DOINT(4)
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4,OMITVNUC(2)*4
#include "memint.h"

C     Memory allocation
      CALL MEMGET('REAL',K_EDM,NNBBASX,WORK,KFREE,LFREE)
 
      write(LUPRI,'(/,2X,A,/)') '....integral routine for EDM '
       
C     Memory deallocation
      CALL MEMREL('EDM',WORK,KWORK,KWORK,KFREE,LFREE)
      !CALL QUIT('EDM end')        
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck onefck */
      SUBROUTINE ONEFCK(FOCK,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C   PURPOSE:
C ===========
C     Generate one-electron Fock matrix, using the one-electron operators
C     specified in COMMON block CBIHAM.
C
C     It can be (in SA-AO basis): 
C         i) Dirac 1-el. operator
C        ii) BSS/IOTC two-component operator (either in the 4comp. or in the reduced 2comp. framework)
C       iii) Non-relativistic operator (in reduced two-component framework)
C
C     Routine is called from many places.
C
C     NOTE: If BSS=.true., get BSS integrals (TWOCOMP does not matter later)
C              BSS=.false., TWOCOMP does matter (Dirac or nonrelativistic Ham.) 
C
C     Written by Trond Saue May 1996
C     Last revision: May 19 1996
C
C MI&HJAaJ (2003): If BSS=.TRUE., read the TWO-COMPONENT
C                   relativistic bare nuclei (in AO basis) from the file
C                   BSSMAT (label='H1AO_DK')
C
C   MI/Strasbg,2005 added new labels containing
C        various modified atomic integrals of H2c AO
C   MI/Strasbg,Febr. 2006 H2cAO Hamiltonian reading placed into separate routine
C
C*****************************************************************************
      use x2cmod_cfg
      use x2c_utils, only:
     &    read_1fock_x2c
#ifdef MOD_INTEREST
      use module_interest_interface
#endif
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "cbihr1.h"
      DIMENSION FOCK(*),WORK(LWORK)
C
      CALL QENTER('ONEFCK')
#include "memint.h"

CMI   ... control output ...
      IF (IPRHAM.GE.4) THEN
        IF (BSS .or. x2c) THEN
         if(x2c)then
          CALL HEADER(
     &    'ONEFCK: One-electron Fock matrix from - '//
     &    ' X2C 2c Hamiltonian ',-1)
         end if
         IF (TWOCOMPBSS) THEN
          CALL HEADER(
     &    'ONEFCK: One-electron Fock matrix from - '//
     &    ' BSS 2c Hamiltonian (in pure two-comp mode, '//
     &    ' saving memory)',-1)
         ELSE
          CALL HEADER(
     &    'ONEFCK: One-electron Fock matrix from - '//
     &    ' BSS 2c Hamiltonian (in four-component mode, '//
     &    'wasting memory)',-1)
         ENDIF
        ELSE
          IF (TWOCOMP) THEN
         CALL HEADER(
     &   'ONEFCK: One-electron Fock matrix - nonrelativistic'//
     &   '  Hamiltonian (BSS=.false.,TWOCOMP=.true.)',-1)
          ELSEIF(MDIRAC) THEN
         CALL HEADER(
     &   'ONEFCK: One-electron Fock matrix from the'//
     &   ' modified Dirac Hamiltonian (BSS=.false.,TWOCOMP=.false.)',-1)
          ELSE
         CALL HEADER(
     &   'ONEFCK: One-electron Fock matrix from the'//
     &   ' Dirac Hamiltonian (BSS=.false.,TWOCOMP=.false.)',-1)
          ENDIF
       ENDIF
      ENDIF

      if(bss)then
!MI If BSS is active read one-electron DKn(BSS) AO integrals from the
!   BSSMAT file (lu=LUBSS) according to the parameter I2COFK
       CALL GETH2CAO(I2COFK,FOCK,WORK(KFREE),LFREE,IPRHAM)
      else if(x2c.and..not.x2c_4c_fock_mtx_defh1)then
        print *, 'blubb blubb x2c read'
        call read_1fock_x2c(fock,n2bbasx,nz,i2cofk)
      else
C
C     Memory allocation
C
cisl-b    getting NZ*N2BBASX instead of NZ*NNBBASX
C
C     ************************************************
C     ****** Construct one-electron Fock matrix ******
C     ************************************************
C
       IF(MDIRAC) THEN
#ifdef MOD_INTEREST
         fock(1:n2bbasxq) = 0.0d0
         call process_kinetic_energy_integrals( Fock )
         call process_nuclear_attraction_integrals( Fock ) 
#endif
       ELSE
         CALL MEMGET('REAL',KFTRI ,NZ*N2BBASX,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KONEIN,   N2BBASX,WORK,KFREE,LFREE)
         CALL MEMGET('LOGI',KFIRST,NZ        ,WORK,KFREE,LFREE)
         CALL ONEFC1(FOCK,WORK(KFTRI),WORK(KONEIN),WORK(KFIRST),
     &              WORK(KFREE),LFREE)
C
C        Memory deallocation
         CALL MEMREL('ONEFCK.ONEFC1',WORK,KWORK,KWORK,KFREE,LFREE)
       ENDIF

      end if
C
C     Output section
C     ==============
C
#define blubb_blubb_so
#ifdef blubb_blubb_so
      iprham_save = iprham
      iprham      = 4
      print *, 'bss, x2c, twocomp',bss,x2c,twocomp
#endif
      IF(IPRHAM.GE.4) THEN
        IF (BSS) THEN
         IF (TWOCOMPBSS) THEN
          CALL HEADER(
     &    'ONEFCK: One-electron Fock (symm.blocked SA-AO)'//
     &    ' matrix from - '//
     &    ' BSS 2c Hamiltonian (in pure two-comp mode, '//
     &    ' saving memory)',-1)
         ELSE
          CALL HEADER(
     &    'ONEFCK: One-electron Fock (symm.blocked SA-AO)'// 
     &    'matrix from - '//
     &    ' BSS 2c Hamiltonian (LL blocks in four-component mode, '//
     &    'wasting memory)',-1)
         ENDIF
        ELSE IF(x2c)THEN
          CALL HEADER(
     &    'ONEFCK: One-electron Fock (symm.blocked SA-AO)'//
     &    ' matrix from - X2C 2c Hamiltonian)',-1)
        ELSE
          IF (TWOCOMP) THEN
         CALL HEADER(
     &   'ONEFCK: One-electron Fock matrix - nonrelativistic'//
     &   '  Hamiltonian (BSS=.false.,TWOCOMP=.true.)',-1)
          ELSEIF(MDIRAC) THEN
         CALL HEADER(
     &   'ONEFCK: One-electron Fock matrix from the'//
     &   ' modified Dirac Hamiltonian (BSS=.false.,TWOCOMP=.false.)',-1)
          ELSE
         CALL HEADER(
     &   'ONEFCK: One-electron Fock matrix from the'//
     &   ' Dirac Hamiltonian (BSS=.false.,TWOCOMP=.false.)',-1)
          ENDIF
        ENDIF
        DO 10 I = 1,NFSYM
          IF(NFBAS(I,0).EQ.0) GOTO 10
          WRITE(LUPRI,'(A,I1,A,I1)')
     &     '* Fermion ircop no.',I,'/',NFSYM   
          CALL PRQMAT(FOCK(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   10   CONTINUE
      ENDIF
#ifdef blubb_blubb_so
      iprham = iprham_save
#undef blubb_blubb_so
#endif
C
      CALL QEXIT('ONEFCK')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck onefc1 */
      SUBROUTINE ONEFC1(FOCK,FTRI,OP1INT,FIRST,WORK,LWORK)
C*****************************************************************************
C
C     Generate one-electron Fock matrix, using the one-electron operators
C     specified in COMMON block CBIHAM.
C
C     Called from: ONEFCK
C
C     Written by Trond Saue May 1996
C     Last revision: May 19 1996
C                  MI / Febr.2006
C
C*****************************************************************************
      use dirac_cfg
      use fde_types
      use fde_mod

#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "dummy.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dcbxpr.h"
#include "dcbprl.h"
C
      type(fde_import) :: itmp
      LOGICAL FIRST(NZ),TOBE
      DIMENSION FOCK(N2BBASX,NZ),FTRI(NNBBASX,NZ),OP1INT(NNBBASX),
     &          WORK(LWORK)
C
      CALL QENTER('ONEFCK1')
#include "memint.h"
C
C     Initialize
C
#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
      CALL LSET(NZ,.TRUE.,FIRST)
C
C     Nuclear attraction integrals
C     ============================
C
      IOP = IPMOLFLD
      CALL PRPMAO(LU1INT,IOP,NNBBASX,FTRI,
     &            OP1INT,FIRST,IPRHAM)
CMI     ... dcbham.h is needed
      IF (SMLV1C) THEN
         CALL CATCORR(FTRI,DUMMY,DUMMY,WORK,LFREE,IPRHAM)
         CALL PR1CEX1 (FTRI,FIRST,FOCK,WORK,KFREE,LFREE,IPRHAM)
      ENDIF

C     Embedding potential integrals, if FDE is requested
C     ==================================================
C
C aspg: i wonder whether or not there should be more control to allow 
C       for the possibility of not including the static potential to
C       the one-electron matrix prior to the scf...

      call fde_get_import_info(itmp)
      if (dirac_cfg_fde .and. itmp%im_vemb) then
         IOP = IPVEMB0
         CALL PRPMAO(LU1INT,IOP,NNBBASX,FTRI,
     &          OP1INT,FIRST,IPRHAM)
      ENDIF

C
C     Kinetic energy
C     ==============
C
CMI/Febr.2006 For TWOCOMP=.true. this is the exact nonrelativistic kinetic energy
C
      IOP = IPKINERG
      CALL PRPMAO(LU1INT,IOP,NNBBASX,FTRI,
     &       OP1INT,FIRST,IPRHAM)

C
C     Beta matrix
C     ===========
C
CMI/Febr.2006 Beta matrix is not initialized in OP1INI for the two-component mode !

      IF (.NOT.TWOCOMP) THEN

        IOP = IPBETAMT
        CALL PRPMAO(LU1INT,IOP,NNBBASX,FTRI,
     &           OP1INT,FIRST,IPRHAM)

      ENDIF

C
C     Finite field: 
C     Loop through list of one-electron operators in CBIHAM;
C     ======================================================
C     ...skip Fock-matrices for now ...
CMI   ... only ordinary (symmetric) operators added to the Hamiltonian
      TOBE = .FALSE.
      DO I = 1,N1OPER
        IOP = IPR1OP(I)
        IOPTYP = IPRPTYP(IOP)
        IF(IOPTYP.EQ.0) THEN
C.........Fock matrix
          TOBE = .TRUE.
        ELSE
          CALL PRPMAO(LU1INT,IOP,NNBBASX,FTRI,OP1INT,FIRST,
     &                IPRHAM)
        ENDIF
      ENDDO
C
C
C     Make full matrix
C     ================
C
      DO IZ = 1,NZ
        IF(FIRST(IZ)) THEN
          CALL DZERO(FOCK(1,IZ),N2BBASX)
        ELSE
          IQ = IPQTOQ(IZ,0)
          IH = IHQMAT(IQ,1)
          IF(IH.EQ.1) THEN        
            CALL DSPTSI(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DAPTGE(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ENDIF
        ENDIF
      ENDDO
c
c    now add SO integrals
c     SOPP integrals
c     ==============
c     
      if (ECPCALC) then
         IOP = IPSOPP
c      making it square by using N2BBASX
         CALL PRPMAOso(LU1INT,IOP,N2BBASX,FOCK,
     &          OP1INT,FIRST,IPRHAM)
      endif
C
      CLOSE(LU1INT,STATUS = 'KEEP')
C
C
C     Reindex to sorted basis
C     =======================
C
      CALL BUTOBS(FOCK,NZ,WORK(KFREE),LFREE)
C   
C     Add Fock matrices when requested
C
      IF(TOBE) THEN
        CALL MEMGET('REAL',KFBUF,N2BBASXQ,WORK,KFREE,LFREE)
        DO I = 1,N1OPER
          IOP = IPR1OP(I)
          IOPTYP = IPRPTYP(IOP)
          IF(IOPTYP.EQ.0) THEN
            KLBL  = IPRPLBL(1,IOP)
            INQUIRE(FILE=PRPLBL(KLBL),EXIST=FIRST(1))
            IF(FIRST(1)) THEN
              CALL OPNFIL(2,PRPLBL(KLBL),'OLD','ONEFC1')
              CALL REAFCK(2,WORK(KFBUF),.FALSE.,1)
              CLOSE(2,STATUS = 'KEEP')
              PFAC = FACPRP(1,IOP)
              CALL DAXPY(N2BBASXQ,PFAC,WORK(KFBUF),1,FOCK,1)
              CALL LSET(NZ,.FALSE.,FIRST)
            ELSE
              WRITE(LUPRI,'(A,A8,A,A16,A)')
     &        'Property file ',PRPLBL(KLBL),' for operator ',
     &        PRPNAM(IOP),' not found !'
              CALL QUIT('ONEFC1: Fock matrix file not found !')
            ENDIF
          ENDIF
        ENDDO
        CALL MEMREL('ONEFC1.fock',WORK,KWORK,KFBUF,KFREE,LFREE)
      ENDIF
      CALL QEXIT('ONEFCK1')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpmao */
CMI ... TODO: to be later renamed to PRPMBU ...
Chjaaj: NO! routine names should not be too cryptic, BU is too special to me.
C       Current name is made from: PRPMAO - PRoPerty Matrix AO basis
C       If we also want to tell that the order of the AO basis is BU, then
C       e.g. as new nanme: PRPMAO_BU /nov.08-hjaaj
      SUBROUTINE PRPMAO(LU1INT,IOP,IDIMENS,PMAT,
     &                  OP1INT,FIRST,IPRINT)
C*****************************************************************************
C
C     Generates  lower triangular/square  matrix of one-electron integrals
C     corresponding to property IOP of COMMON block XCBPRP
C
C     On input:  LU1INT  - lu of integral file
C                IOP - index of the property
C                IDIMENS - dimension of array (either triangular or square)
C                FIRST - array, must be set to true
C
C     On output: PMAT - accumulated property matrix (Hermit sorting)
C
C
C     Written by Trond Saue May 1996
C     Last revision: May 19 1996 - tsaue
C     Dec.2005, MI - generalized for square properties (IDIMENS = NNBBASX or N2BBASX)
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,PTOL = 1.0D-12)
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbxpr.h"
#include "dcbprl.h"
C
      CHARACTER MXFORM*6, FMT*6, RTNLBL(2)*8,JCHR*1
      LOGICAL FIRST(*), DOINT(4), SQUARE
      DIMENSION PMAT(IDIMENS,*),OP1INT(IDIMENS),IDOINT(4)
C

CMI   ... find out type of AO matrix - triangular or square
      IF (IDIMENS.NE.N2BBASX.AND.IDIMENS.NE.NNBBASX) THEN
       write(lupri,*) 'PRPMAO: NNBBASX,N2BBASX,IDIMENS=',
     &                NNBBASX,N2BBASX,IDIMENS
       call QUIT('PRPMAO: wrong dimension !!!')
      ENDIF

CMI   ... this initialization is needed
      SQUARE = .FALSE.
      SQUARE = (IDIMENS.EQ.N2BBASX)
C
C     Extract information about property
C     ==================================
C
      IOPTYP = IPRPTYP(IOP)
      IF(IPRINT.GE.3) THEN
        CALL HEADER('Output from PRPMAO',-1)
        WRITE(LUPRI,*)
     & 'Generating property matrix for the property: ',PRPNAM(IOP)
        IF (SQUARE) THEN
         write(lupri,*)
     &   '.... it is the SQUARE matrix ' 
        ELSE
         write(lupri,*)
     &   '.... it is the LOWER TRIANGULAR matrix ' 
        ENDIF
      ENDIF
C
C     Loop over components
C     ====================
C
      IBRPOP = IPRPSYM(IOP)-1
      NCMP   = MCMP(IOPTYP)
      DO 10 I = 1,NCMP
        PFAC  = FACPRP(I,IOP)
        IF(ABS(PFAC).LT.PTOL) GOTO 10
C
C       Extract information about component
C
C       SCALAR OPERATOR:
C       ----------------
C       KLBL - pointer to list of labels
C       ISYM - scalar operator is symmetric(1)/antisymmetric(2) about diagonal
C       IBREP - boson irrep of scalar operator
C
C       MATRIX OPERATOR:
C       ----------------
C
C       IM4  - pointer to matrix operator
C       IFAC - sign of coefficient in linear combination
C       IFAS - sign of matrix operator
C         Anti-symmetric scalar operators are provided with an imaginary phase 
C         in order to retain Hermiticity
C       IMAT - quaternion position of matrix operator
C
        KLBL  = IPRPLBL(I,IOP)
C
C       Determine what matrix the operator enters
C       =========================================
C
C       Is this a LL+SS operator (IC=1) or a LS+SL operator (IC=2) ?
        IF(PDOINT(KLBL)(2:3) .EQ. '00') THEN
          IC = 1
        ELSE IF((PDOINT(KLBL)(1:1) .EQ. '0') .AND.
     &          (PDOINT(KLBL)(4:4) .EQ. '0')) THEN 
          IC = 2
        ELSE
          write(lupri,'(3A)') 'PRPMAO: illegal PDOINT = "',PDOINT,'"'
          CALL QUIT('PRPMAO: illegal PDOINT!')
        ENDIF
CMI     ... assignements
        IBRP = IPRLREP(KLBL)
        IH   = IPRLTYP(KLBL)
        IM4  = JM4(I,IOPTYP)
        IZ   = JM4POS(IM4)
        IQ   = IQMULT(1,JQBAS(IBRP,IC),IZ)
        IQP  = IQTOPQ(IQ,IBRPOP)
CMI fix
!       IF (IOPTYP.EQ.21) THEN 
!          IQP=1
!          WRITE(LUPRI,'(/,2X,A)')
!     &    'Operator of type 21, IQP set to 1.'
!       ENDIF
! check
        IF(IQP.LE.0) THEN
          WRITE(LUPRI,'(/A)') '***** PRPMAO ERROR *****'
          WRITE(LUPRI,'(A,I3,A,I2,3X,A,A8)')
     &     'Component nr. ',I,'/',NCMP,' label:',PRPLBL(KLBL)
          WRITE(LUPRI,'(A,I2,A,I2,A,I2,A)')
     &    'IQP(',IQP,')=IQTOPQ(',IQ,',',IBRPOP,').LE.0 !'
          CALL WRIXPR(0,IOP)
          WRITE(LUPRI,'(A,I5)') 'Q: Scalar :',JQBAS(IBRP,IC)
          WRITE(LUPRI,'(A,I5)') 'Q: Matrix :',IZ
          WRITE(LUPRI,'(A,9I3)')
     &        'IM4,IBRP,IH,IZ,IQ,IBRPOP,IQP,IC,IOPTYP :',
     &         IM4,IBRP,IH,IZ,IQ,IBRPOP,IQP,IC,IOPTYP
          WRITE(LUPRI,'(A)') 'IQTOPQ(1:4,0:NBSYM-1):'
          DO IZ = 1,4
             WRITE(LUPRI,'(8I5)') (IQTOPQ(IZ,J),J=0,NBSYM-1)
          END DO
          WRITE(LUPRI,'(A)') 'IPQTOQ(1:NZ,0:NBSYM-1):'
          DO IZ = 1,NZ
             WRITE(LUPRI,'(8I5)') (IPQTOQ(IZ,J),J=0,NBSYM-1)
          END DO
          CALL QUIT('PRPMAO: IQP.LE.0 !!!')
        END IF
        IF(IQP.GT.NZ) THEN
          WRITE(LUPRI,'(/A)') '***** PRPMAO ERROR *****'
          WRITE(LUPRI,'(A,I5,3X,A8)') 'Component nr.',I,PRPLBL(KLBL)
          WRITE(LUPRI,'(A,I2,A,I2,A)') 'IQP(',IQP,').GT.NZ(',
     &          NZ,')'
          CALL WRIXPR(0,IOP)
          WRITE(LUPRI,'(A,I5)') 'Q: Scalar :',JQBAS(IBRP,IC)
          WRITE(LUPRI,'(A,I5)') 'Q: Matrix :',IZ
          CALL QUIT('PRPMAO: IQP.GT.NZ !!!')
        ENDIF
        IF(IPRINT.GE.4) THEN
          WRITE(LUPRI,'(A,I5,3X,A8)') 'Component nr.',I,PRPLBL(KLBL)
          FMT = MXFORM(PFAC,20)
          WRITE(LUPRI,'(A,4X,'//FMT//',A,A4)')
     &      'Factor: ',PFAC,' and PDOINT:',PDOINT(KLBL)
          WRITE(LUPRI,'(A,I5)') 'Integrals enter matrix :',IQP
          WRITE(LUPRI,'(A,L1,A,I2)')
     &    'FIRST = ',FIRST(IQP),' and IC(1 or 2):',IC
        ENDIF
C =============================================================
C               Get integrals from the AOPROPER file
C =============================================================
        CALL PRPREA(LU1INT,PRPLBL(KLBL),RTNLBL,PDOINT(KLBL),
     &              OP1INT,IDIMENS,IPRINT)
        IF(IPRINT.GE.8) THEN
          WRITE(LUPRI,'(A,I3)')
     &    '* Integrals from AOPROPER, size of matrix ',NTBAS(0)
          IF (NTBAS(0).EQ.NTBAS(1)) THEN
           WRITE(LUPRI,'(2X,A)')
     &   '...two-component framework (only LL block is present) '
          ELSE
           WRITE(LUPRI,'(2X,A)')
     &   '...four-component framework  '
          ENDIF
          IF (SQUARE) THEN
           CALL OUTPUT(OP1INT,1,NTBAS(0),1,NTBAS(0),
     &                NTBAS(0),NTBAS(0),-1,LUPRI)
          ELSE
           CALL OUTPAK(OP1INT,NTBAS(0),-1,LUPRI)
          ENDIF
        ENDIF        
C
C       Determine factor to scale matrix with
C       =====================================
C
        IFAC  = JCOM(I,IOPTYP)
        PFAC  = FACPRP(I,IOP)*IFAC
        IF(PFAC.NE.D1) THEN
          CALL DSCAL(IDIMENS,PFAC,OP1INT,1)
          IF(IPRINT.GE.8) THEN
            WRITE(LUPRI,'(/A,D13.6)')
     &    '* scaled integrals with the factor=',PFAC
           IF (SQUARE) THEN
            CALL OUTPUT(OP1INT,1,NTBAS(0),1,NTBAS(0),
     &              NTBAS(0),NTBAS(0),-1,LUPRI)
           ELSE
             CALL OUTPAK(OP1INT,NTBAS(0),-1,LUPRI)
           ENDIF
          ENDIF
        ENDIF        
C
C       Phase transformation/PDOINT factor if necessary.
C       The four blocks (AO basis): LL, SL, LS, SS
C
Chj     PDOINT is now '+', '-', or '0',
Chj     corresponding to a factor of +1, -1, or 0 on that block.
        NM1 = 0
        DO J = 1,4
           JCHR = PDOINT(KLBL)(J:J)
           IF (JCHR .EQ. '0') THEN
              IDOINT(J) = 0
              DOINT(J) = .FALSE.
           ELSE IF (JCHR .EQ. '+') THEN
              IDOINT(J) = 1
              DOINT(J) = .TRUE.
           ELSE IF (JCHR .EQ. '-') THEN
              NM1 = 1
              IDOINT(J) = -1
              DOINT(J) = .TRUE.
           ELSE
              CALL QUIT('PRPMAO: Illegal factor in PDOINT')
           END IF
        END DO
C
        IF (NZ.LT.4 .OR. NM1.GT.0) THEN
          IF (.NOT.SQUARE) THEN
           CALL QNPHASE(IZ,IBRP,OP1INT,IDOINT)
          ELSE
           CALL QMPHASE('S',IZ,IBRP,OP1INT,DOINT)
          ENDIF

          IF(IPRINT.GE.8) THEN
             WRITE(LUPRI,'(A/3A)')
     &          '* Phase and/or PDOINT transformed integrals:',
     &          '  ( PDOINT = ',PDOINT(KLBL),' )'
             IF (SQUARE) THEN
              CALL OUTPUT(OP1INT,1,NTBAS(0),1,NTBAS(0),
     &             NTBAS(0),NTBAS(0),-1,LUPRI)
             ELSE
              CALL OUTPAK(OP1INT,NTBAS(0),-1,LUPRI)
             ENDIF
          ENDIF
        ENDIF
C
C       Determine where integrals go and add them on
C
        IF(FIRST(IQP)) THEN
          CALL DCOPY(IDIMENS,OP1INT,1,PMAT(1,IQP),1)
          FIRST(IQP) = .FALSE.
        ELSE
C         ... accumulate the property ...
          CALL DAXPY(IDIMENS,D1,OP1INT,1,PMAT(1,IQP),1)
        ENDIF
 10   CONTINUE

      IF(IPRINT.GE.8) THEN
        DO IZ = 1, NZ
          WRITE(LUPRI,'(A,I2,A)') '* Property matrix (IZ=',IZ,')'
          IF (FIRST(IZ)) THEN
             WRITE(LUPRI,'(/A)') '  --- zero matrix ---'
          ELSE
           IF (SQUARE) THEN
             CALL OUTPUT(PMAT(1,IZ),1,NTBAS(0),1,NTBAS(0),
     &                   NTBAS(0),NTBAS(0),-1,LUPRI)
           ELSE
             CALL OUTPAK(PMAT(1,IZ),NTBAS(0),-1,LUPRI)
           ENDIF
          END IF
        ENDDO
      ENDIF

      RETURN
C
!radovan: the rest of the routine is unreachable
 1000 CONTINUE
C
C     Not able to read IDOINT information from PDOINT
C
      WRITE(LUPRI,'(A,A)') 'PRPMAO: Class label:',PRPLBL(KLBL)
      WRITE(LUPRI,'(A,A)') 'PRPMAO ERROR: Not able to read COMB =',
     &                     PDOINT(KLBL)
      CALL QUIT('PRPMAO: Not able to read COMB')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
C  /* Deck prpmaoso */
      SUBROUTINE PRPMAOso(LU1INT,IOP,IDIMENS,PMAT,
     &                  OP1INT,FIRST,IPRINT)
C*************************************************************************
****    
C     
C     Generates lower triangular/square  matrix of one-electron integrals
C     corresponding to property IOP of COMMON block XCBPRP
C
C     On input:  LU1INT
C                   IOP
C     
C     On output: PMAT - property matrix (Hermit sorting)
C
C        
C     Written by Trond Saue May 1996
C     Last revision: May 19 1996 - tsaue
C     Dec.2005, MI - adapted for square properties (IDIMENS = NNBBASX or N2BBASX)
C
C*************************************************************************
****
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,PTOL = 1.0D-12)
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbxpr.h"
#include "dcbprl.h"
C
      CHARACTER MXFORM*6, FMT*6, RTNLBL(2)*8,JCHR*1
      LOGICAL FIRST(*), DOINT(4), SQUARE
      DIMENSION PMAT(IDIMENS,4),OP1INT(IDIMENS),IDOINT(4)
C


CMI   ... find out type of AO matrix - triangular or square
      IF (IDIMENS.NE.N2BBASX.AND.IDIMENS.NE.NNBBASX) THEN
       write(lupri,*) 'PRPMAO: NNBBASX,N2BBASX,IDIMENS=',
     &                NNBBASX,N2BBASX,IDIMENS
       call QUIT('PRPMAO: wrong dimension !!!')
      ENDIF

CMI   ... this initialization is needed
      SQUARE = .true.
C
C     Extract information about property
C     ==================================
C
      IOPTYP = IPRPTYP(IOP)
      IF(IPRINT.GE.3) THEN
        CALL HEADER('Output from PRPMAO',-1)
        WRITE(LUPRI,*)
     & 'Generating property matrix for property ',PRPNAM(IOP)
        IF (SQUARE) THEN
         write(lupri,*)
     &   '.... it is the SQUARE matrix '
        ELSE
         write(lupri,*)
     &   '.... it is the LOWER TRIANGULAR matrix '
        ENDIF
      ENDIF
C
C     Loop over components
C     ====================
C
      IBRPOP = IPRPSYM(IOP)-1
      NCMP   = MCMP(IOPTYP)
      DO 11 I = 1,NCMP
        PFAC  = FACPRP(I,IOP)
        IF(ABS(PFAC).LT.PTOL) GOTO 11
C
C       Extract information about component
C
C       SCALAR OPERATOR:
C       ----------------
C       KLBL - pointer to list of labels
C       ISYM - scalar operator is symmetric(1)/antisymmetric(2) about diagonal
C       IBREP - boson irrep of scalar operator
C
C       MATRIX OPERATOR:
C       ----------------
C
C       IM4  - pointer to matrix operator
C       IFAC - sign of coefficient in linear combination
C       IFAS - sign of matrix operator
C         Anti-symmetric scalar operators are provided with an imaginary phase
C         in order to retain Hermiticity
C       IMAT - quaternion position of matrix operator
C
        KLBL  = IPRPLBL(I,IOP)
C
C       Determine what matrix the operator enters
C       =========================================
C
C       Is this a LL+SS operator (IC=1) or a LS+SL operator (IC=2) ?
        IF(PDOINT(KLBL)(2:3) .EQ. '00') THEN
          IC = 1
        ELSE IF((PDOINT(KLBL)(1:1) .EQ. '0') .AND.
     &          (PDOINT(KLBL)(4:4) .EQ. '0')) THEN
          IC = 2
        ELSE
          write(lupri,'(3A)') 'PRPMAOso: illegal PDOINT = "',PDOINT,'"'
          CALL QUIT('PRPMAOso: illegal PDOINT!')
        ENDIF
        IBRP = IPRLREP(KLBL)
        IH   = IPRLTYP(KLBL)
        IM4  = JM4(I,IOPTYP)
        IZ   = JM4POS(IM4)
        IQ   = IQMULT(1,JQBAS(IBRP,IC),IZ)
        IQP  = IQTOPQ(IQ,IBRPOP)
        IF(IQP.LE.0) THEN
          WRITE(LUPRI,'(/A)') '***** PRPMAO ERROR *****'
          WRITE(LUPRI,'(A,I5,3X,A8)') 'Component nr.',I,PRPLBL(KLBL)
          WRITE(LUPRI,'(A,I2,A)') 'IQP(',IQP,').LE.0'
          CALL WRIXPR(0,IOP)
          WRITE(LUPRI,'(A,I5)') 'Q: Scalar :',JQBAS(IBRP,IC)
          WRITE(LUPRI,'(A,I5)') 'Q: Matrix :',IZ
          WRITE(LUPRI,'(A,7I5)')
     &        'IM4,IBRP,IH,IZ,IQ,IBRPOP,IQP :',
     &         IM4,IBRP,IH,IZ,IQ,IBRPOP,IQP
          WRITE(LUPRI,'(A)') 'IQTOPQ(1:4,0:NBSYM-1):'
          DO IZ = 1,4
             WRITE(LUPRI,'(8I5)') (IQTOPQ(IZ,J),J=0,NBSYM-1)
          END DO
          WRITE(LUPRI,'(A)') 'IPQTOQ(1:NZ,0:NBSYM-1):'
          DO IZ = 1,NZ
             WRITE(LUPRI,'(8I5)') (IPQTOQ(IZ,J),J=0,NBSYM-1)
          END DO
          CALL QUIT('PRPMAO: IQP.LE.0 !!!')
        END IF
        IF(IQP.GT.NZ) THEN
          WRITE(LUPRI,'(/A)') '***** PRPMAO ERROR *****'
          WRITE(LUPRI,'(A,I5,3X,A8)') 'Component nr.',I,PRPLBL(KLBL)
          WRITE(LUPRI,'(A,I2,A,I2,A)') 'IQP(',IQP,').GT.NZ(',
     &          NZ,')'
          CALL WRIXPR(0,IOP)
          WRITE(LUPRI,'(A,I5)') 'Q: Scalar :',JQBAS(IBRP,IC)
          WRITE(LUPRI,'(A,I5)') 'Q: Matrix :',IZ
          CALL QUIT('PRPMAO: IQP.GT.NZ !!!')
        ENDIF
        IF(IPRINT.GE.4) THEN
          WRITE(LUPRI,'(A,I5,3X,A8)') 'Component nr.',I,PRPLBL(KLBL)
          FMT = MXFORM(PFAC,20)
          WRITE(LUPRI,'(A,4X,'//FMT//',A,A4)')
     &      'Factor: ',PFAC,' and PDOINT:',PDOINT(KLBL)
          WRITE(LUPRI,'(A,I5)') 'Integrals enter matrix :',IQP
          WRITE(LUPRI,'(A,L1,A,I2)')
     &    'FIRST = ',FIRST(IQP),' and IC(1 or 2):',IC

        ENDIF
C
C       Get integral from file
C
        CALL PRPREA(LU1INT,PRPLBL(KLBL),RTNLBL,PDOINT(KLBL),
     &              OP1INT,IDIMENS,IPRINT)
        IF(IPRINT.GE.8) THEN
          WRITE(LUPRI,'(A)') '* Integrals from AOPROPER'
          IF (SQUARE) THEN
           CALL OUTPUT(OP1INT,1,NTBAS(0),1,NTBAS(0),
     &                NTBAS(0),NTBAS(0),-1,LUPRI)
          ELSE
           CALL OUTPAK(OP1INT,NTBAS(0),-1,LUPRI)
          ENDIF
        ENDIF
C
C       Determine factor to scale matrix with
C       =====================================
C
        IFAC  = JCOM(I,IOPTYP)
        PFAC  = FACPRP(I,IOP)*IFAC
        IF(PFAC.NE.D1) THEN
          CALL DSCAL(IDIMENS,PFAC,OP1INT,1)
          IF(IPRINT.GE.8) THEN
            WRITE(LUPRI,'(/A,D13.6)')
     &    '* Scaled integrals with factor=',PFAC
           IF (SQUARE) THEN
            CALL OUTPUT(OP1INT,1,NTBAS(0),1,NTBAS(0),
     &              NTBAS(0),NTBAS(0),-1,LUPRI)
           ELSE
             CALL OUTPAK(OP1INT,NTBAS(0),-1,LUPRI)
           ENDIF
          ENDIF
        ENDIF
C
C       Phase transformation/PDOINT factor if necessary.
C       The four blocks (AO basis): LL, SL, LS, SS
C
c        this is changed so that
c        PDOINT is now '+', '-', or '0',
c        corresponding to a factor of -1, +1, or 0 on that block.
        NM1 = 0
        DO J = 1,4
           JCHR = PDOINT(KLBL)(J:J)
           IF (JCHR .EQ. '+') THEN
              IDOINT(J) = 0
              DOINT(J) = .FALSE.
           ELSE IF (JCHR .EQ. '0') THEN
              IDOINT(J) = 1
              DOINT(J) = .TRUE.
           ELSE IF (JCHR .EQ. '-') THEN
              NM1 = 1
              IDOINT(J) = -1
              DOINT(J) = .TRUE.
           ELSE
              CALL QUIT('PRPMAO: Illegal factor in PDOINT')
           END IF
        END DO
C
        IF (NZ.LT.4 .OR. NM1.GT.0) THEN
          IF (.NOT.SQUARE) THEN
           CALL QNPHASE(IZ,IBRP,OP1INT,IDOINT)
          ELSE
           CALL QMPHASE('S',IZ,IBRP,OP1INT,DOINT)
          ENDIF

          IF(IPRINT.GE.8) THEN
             WRITE(LUPRI,'(A/3A)')
     &          '* Phase and/or PDOINT transformed integrals:',
     &          '  ( PDOINT = ',PDOINT(KLBL),' )'
             IF (SQUARE) THEN
              CALL OUTPUT(OP1INT,1,NTBAS(0),1,NTBAS(0),
     &             NTBAS(0),NTBAS(0),-1,LUPRI)
             ELSE
              CALL OUTPAK(OP1INT,NTBAS(0),-1,LUPRI)
             ENDIF
          ENDIF
        ENDIF
C
C       Determine where integrals go and add them on
C
        IF(FIRST(IQP)) THEN
          CALL DCOPY(IDIMENS,OP1INT,1,PMAT(1,IQP),1)
          FIRST(IQP) = .FALSE.
        ELSE
          CALL DAXPY(IDIMENS,D1,OP1INT,1,PMAT(1,IQP),1)
        ENDIF
 11   CONTINUE
      IF(IPRINT.GE.8) THEN
        DO IZ = 1,NZ
          WRITE(LUPRI,'(A,I2,A)') '* Property matrix (IZ=',IZ,')'
          IF (FIRST(IZ)) THEN
             WRITE(LUPRI,'(/A)') '  --- zero matrix ---'
          ELSE
           IF (SQUARE) THEN
             CALL OUTPUT(PMAT(1,IZ),1,NTBAS(0),1,NTBAS(0),
     &                   NTBAS(0),NTBAS(0),-1,LUPRI)
           ELSE
             CALL OUTPAK(PMAT(1,IZ),NTBAS(0),-1,LUPRI)
           ENDIF
          END IF
        ENDDO
      ENDIF
      RETURN
!radovan: the rest of the routine is unreachable
C
C     Not able to read IDOINT information from PDOINT
C
      WRITE(LUPRI,'(A,A)') 'PRPMAO: Class label:',PRPLBL(KLBL)
      WRITE(LUPRI,'(A,A)') 'PRPMAO ERROR: Not able to read COMB =',
     &                     PDOINT(KLBL)
      CALL QUIT('PRPMAO: Not able to read COMB')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRPMSAO(INDXPR,IDIMENS,PTRI,PMAT,
     &              OP1INT,FIRST,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     PURPOSE: Get 2c/4c SA-AO matrix of one-electron operator
C     corresponding to property INDXPR of COMMON block XCBPRP
C
C     On input:  
C    =============
C                INDXPR - index of the operator in the list
C                IDIMENS - dimension of arrays (triangular/square)
C                PTRI,OP1INT,FIRST - arrays to be used inside
C    
C         PTRI, OP1INT have size of IDIMENS.
C
C    On output:
C   =============
C                PMAT - 2c/4c property matrix in SA-AO basis
C                PTRI(IDIMENS,NZ) - property matrix in BU AO basis
C
C NOTE: In the TWOCOM.and.BSS=.true. case one needs only PMAT array !
C       All other input parameters in this case can be dummy ...
C
C     Written by Miro ILIAS/Strasbourg, March 2006
C
C*****************************************************************************
      use x2c_fio
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0, PTOL = 1.0D-12)
C
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dummy.h"
C
      LOGICAL FIRST(NZ),IS_SQUARE, YESPCTRA, MAKEPCTRA
      DIMENSION PTRI(IDIMENS,NZ),PMAT(N2BBASX,NZ),
     &          OP1INT(IDIMENS),WORK(LWORK)
      CHARACTER PRPLABEL*8, prplab_x2c*12, NUMSTR*4
C

      CALL QENTER('PRPMSAO')
#include "memint.h"

CMI ... first MEMCHK does not work !!!
C      CALL MEMCHK('in PRPMSAO',WORK,1)

C
C     Extract symmetry information - ISYM, IREP, ITIM
C     =============================
      ISYM = IPRPSYM(INDXPR)
      IREP = ISYM - 1
      ITIM = IPRPTIM(INDXPR)

      IF (IDIMENS.EQ.N2BBASX.AND.ITIM.EQ.0) THEN
        IS_SQUARE = .TRUE.
      ELSE IF (IDIMENS.EQ.NNBBASX.AND.(ITIM.EQ.1.OR.ITIM.EQ.-1)) THEN
        IS_SQUARE = .FALSE.
      ELSE
        write(lupri,*) 'PRPMSAO: NNBBASX=',NNBBASX,' N2BBASX=',N2BBASX
        write(lupri,*) 'PRPMSAO: entering IDIMENS=',IDIMENS
        CALL QUIT('PRPMSAO: wrong value of IDIMENS !')
      ENDIF

      IPRINT2 = MAX(IPRPRP,IPRINT)
      IF (IPRINT2.GE.6) THEN
       CALL HEADER('*** Output from PRPMSAO ***',-1)
       write(lupri,'(2X,A,I3,A,A)')
     & '* dealing with the operator with INDXPR=',
     & INDXPR,' name=',PRPNAM(INDXPR)
       write(lupri,'(2X,A,I1,A,I1,A,I2)')
     & '* operator ISYM=',ISYM,' IREP=',
     & IREP,' ITIM=',ITIM
       IF (IS_SQUARE) THEN
        write(lupri,'(2X,A)') '* this is a SQUARE operator'
       ELSE
        write(lupri,'(2X,a)') '* this is a TRIANGULAR operator'
       ENDIF
      ENDIF

C ========================================================================================
C    ... in this case read picture change transformed property integrals from BSSMAT file
C      which are there under the name  'P2C_'//NUMSTR 
C
C     X2C module: picture change transformed property integrals reside
C                 on file X2CMAT
C ========================================================================================
      IF((TWOCOMP.AND.BSS).or.x2c) THEN

        IPRINT2 = MAX(IPRHAM,IPRINT)

        IF (N2BBASXQ.NE.(NTBAS(1)*NTBAS(1)*NZ)) THEN
          write(lupri,*) 'PRPMSAO error output:'
          write(lupri,*) 'N2BBASXQ=',N2BBASXQ
          write(lupri,*) 'NTBAS(1)=',NTBAS(1)
          write(lupri,*) 'NTBAS(1)*NTBAS(1)*NZ=',
     &                    NTBAS(1)*NTBAS(1)*NZ
          CALL QUIT('PRPMSAO: N2BBASXQ.ne.N2BBASXQ_L!!!')
       ENDIF

CMI    ...the label will contain the INDXPR operator index in the string form !
       CALL NUM2STR(INDXPR,NUMSTR)

       if(x2c)then
         write(prplab_x2c,'(a8,a4)') 'prpint2c',numstr
         open(lux2c,file='X2CMAT',status='old',form='unformatted',
     &        access='sequential',action="read",position='rewind')
         call x2c_read(prplab_x2c,pmat,n2bbasxq,lux2c)
         close(lux2c,status="keep")
       else
         PRPLABEL = 'P2C_'//NUMSTR
         CALL RFBSSMAT(PRPLABEL,.FALSE.,'PRPMSAO',
     &                 PMAT,N2BBASXQ,LUBSS,IPRINT2)
         IF (IPRINT2.GE.3) THEN  
           write(lupri,'(2x,A,A,A)')   
     &     'PRPMSAO: 2comp. property record ',PRPLABEL,
     &     ' was read from BSSMAT file.'      
         ENDIF
       end if


       CALL LSET(NZ,.FALSE.,FIRST)

      ELSE
CMI    -----  get the full four-component operator -------

        CALL LSET(NZ,.TRUE.,FIRST)

C     ... do resort PTRI to DIRAC sorting
        IPRINT2 = MAX(IPRPRP,IPRINT)
        IF (IS_SQUARE) THEN
          CALL PRPMAO(LU1INT,INDXPR,IDIMENS,PMAT,OP1INT,FIRST,
     &                IPRINT2)   
          CALL SORT2SAO(IDIMENS,DUMMY,PMAT,FIRST,
     &               IREP,ITIM,WORK(KFREE),LFREE,IPRINT2)
        ELSE
          CALL PRPMAO(LU1INT,INDXPR,IDIMENS,PTRI,OP1INT,FIRST,
     &                       IPRINT2)

          CALL SORT2SAO(IDIMENS,PTRI,PMAT,FIRST,
     &               IREP,ITIM,WORK(KFREE),LFREE,IPRINT2)
        ENDIF

CMI   ----  Do the picture change transformation -----                    
        IF (BSS) THEN
         IPRINT1 = MAX(IPRHAM,IPRPRP)
C     ... decide if the picture change transformation is needed
         CALL DECIDEPCTR(YESPCTRA,PRPNAM(INDXPR),IPRINT1)
C     ... if needed, do the picture change transformation of the PMAT
         IF (YESPCTRA) THEN
          MAKEPCTRA = DOPCT .AND. .NOT.NOPCT

C  ... allocate array for the picture change transformation matrix
          ISIZE = NFBAS(1,0)*NFBAS(1,1)*NZ
          IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NFBAS(2,1)*NZ
          CALL MEMGET('REAL',KPCTM,ISIZE,WORK,KFREE,LFREE)
C     ... array for the 2-component picture change transformed operator
          CALL MEMGET('REAL',KPCTO,N2BBASXQ,WORK,KFREE,LFREE)

          CALL DOPCTRA(PMAT,WORK(KPCTM),WORK(KPCTO),ISIZE,MAKEPCTRA,
     &                  PRPNAM(INDXPR),ISYM,IREP,ITIM,
     &                 'S',WORK(KFREE),LFREE,IPRINT1)

          CALL MEMREL('DOPCTRA',WORK,KWORK,KPCTO,KFREE,LFREE)
          CALL MEMREL('DOPCTRA',WORK,KWORK,KPCTM,KFREE,LFREE)
CMI       ... this is needed for the picture change transformed operator
          CALL LSET(NZ,.FALSE.,FIRST)
         ENDIF
        ENDIF

      ENDIF

      CALL QEXIT('PRPMSAO')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck decidepctr */
      SUBROUTINE DECIDEPCTR(DOTRANS,PROPNAME,IPRINT) 
C*****************************************************************************
C
C  Purpose: Decides if the property operator (PROPNAME from the operator list)
C           has to be picture change transformed (returns logical DOTRANS).
C
C  Called from: PRPMSAO
C               PROP2BSS
C
C  Written by Miro ILIAS/Strasbourg/March 2006
C  Last modfications: Miro Ilias/Tel Aviv/November 2007
C
C*****************************************************************************
#include "implicit.h"      
#include "priunit.h"                                        
                                                            
      LOGICAL DOTRANS     
      CHARACTER*16  PROPNAME
                                  
      CALL QENTER('DECIDEPCTR')                             
                                                            
      DOTRANS = .TRUE. 

C ... here is list of operators which DO not have to be picture change tranformed
C ... LL-forms of some operators (Beta matrix, Kinetic energy...) are not usable, however,
C ... some are usable (ZANGMOM) 

      IF ( PROPNAME.EQ.'Overlap matrix  ' .OR.  !                
     &     PROPNAME.EQ.'z-spin'           .OR.  ! overlap * 1/2
     &     PROPNAME.EQ.'Nuc. attraction ' .OR.
     &     PROPNAME.EQ.'Beta matrix     ' .OR.
     &     PROPNAME.EQ.'Kinetic energy  ' .OR.
     &     PROPNAME.EQ.'Orbital z-moment' .OR.  ! ... ZANGMOM
     &     PROPNAME.EQ.'ZANGMOM'          .OR.  ! ... ZANGMOM
     &     PROPNAME.EQ.'z-angular mom'    .OR.  ! ... ZANGMOM
     &     PROPNAME.EQ.'d|S>/dBX        ' .OR.
     &     PROPNAME.EQ.'d|S>/dBY        ' .OR.
     &     PROPNAME.EQ.'d|S>/dBZ        ' .OR.
     &     PROPNAME.EQ.'dS/dBX          ' .OR.
     &     PROPNAME.EQ.'dS/dBY          ' .OR.
     &     PROPNAME.EQ.'dS/dBZ          ' .OR.
     &     PROPNAME.EQ.'dS/dB2XX        ' .OR.
     &     PROPNAME.EQ.'dS/dB2XY        ' .OR.
     &     PROPNAME.EQ.'dS/dB2XZ        ' .OR.
     &     PROPNAME.EQ.'dS/dB2YY        ' .OR.
     &     PROPNAME.EQ.'dS/dB2YZ        ' .OR.
     &     PROPNAME.EQ.'dS/dB2ZZ        ' 
     &      ) DOTRANS = .FALSE.

      IF (IPRINT.GE.5) THEN
       WRITE(LUPRI,*)
       IF (DOTRANS) THEN
         WRITE(lupri,'(2X,A,A,A)')
     &   'DECIDEPCTR: operator of entering name ',PROPNAME,
     &   ' HAS to be picture change transformed.'
       ELSE
         WRITE(lupri,'(2X,A,A,A)')
     &   'DECIDEPCTR: operator of entering name ',PROPNAME,
     &   ' DOES NOT to be picture change transformed.'
       ENDIF
       CALL FLSHFO(LUPRI)
      ENDIF
                                                            
      CALL QEXIT('DECIDEPCTR')                              
      END                                                   

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sort2sao */
      SUBROUTINE SORT2SAO(IDIMENS,PTRI,PMAT,FIRST,
     &                    IREP,ITIM,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C    PURPOSE: Sort the general 4comp. property operator
C              to the symmetry adapted AO-basis
C
C    On input:  PTRI - matrix of 4c prop. operator in unsorted basis
C               IDIMENS - tells whether it is triangular or square matrix
C
C               IREP, ITIM
C
C               IPRINT
C               FIRST 
C
C  NOTE: In the case of SQUARE matrix only the PMAT is entering ! All others
C        parameters in this case are dummy (PTRI only...)
C
C    On output: PMAT
C
C    Call from: PRPMSAO
C
C    Written by Miro ILIAS/Strasbourg/March 2006
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbbas.h"
      DIMENSION PMAT(N2BBASX,NZ),PTRI(IDIMENS,NZ),WORK(LWORK) 
      LOGICAL   FIRST(NZ), IS_SQUARE
C
      CALL QENTER('SORT2SAO')

      IF (IDIMENS.EQ.N2BBASX.AND.ITIM.EQ.0) THEN
        IS_SQUARE = .TRUE.
      ELSE IF (IDIMENS.EQ.NNBBASX.AND.(ITIM.EQ.1.OR.ITIM.EQ.-1)) THEN
        IS_SQUARE = .FALSE.
      ELSE
        write(lupri,*) 'SORT2SAO: N2BBASX,NNBBASX,IDIMENS=',
     &  N2BBASX,NNBBASX,IDIMENS
        CALL QUIT('SORT2SAO: wrong value of IDIMENS !')
      ENDIF

      DO IZ = 1,NZ
         IF(.NOT.FIRST(IZ)) THEN
           IQ = IPQTOQ(IZ,IREP)
           IH = IHQMAT(IQ,ITIM)
           IF(IH.EQ.1) THEN
             CALL DSPTSI(NTBAS(0),PTRI(1,IZ),PMAT(1,IZ))
             IF (IS_SQUARE) THEN
              CALL QUIT('SORT2SAO: wrong value of IDIMENS !') 
             ENDIF
           ELSEIF(IH.EQ.2) THEN
             CALL DAPTGE(NTBAS(0),PTRI(1,IZ),PMAT(1,IZ))
             IF (IS_SQUARE) THEN
              CALL QUIT('SORT2SAO: wrong value of IDIMENS !') 
             ENDIF
           ELSE
CMI   ... in this case PTRI is dummy, only PMAT exists
            IF (.NOT.IS_SQUARE) THEN
              CALL QUIT('SORT2SAO: wrong value of IDIMENS !') 
            ENDIF
           ENDIF
             CALL BUTOBS(PMAT(1,IZ),1,WORK,LWORK)
         ELSE 
CMI ... zeroed out non-used quaternionic blocks of property matrixes...
           CALL DZERO(PMAT(1,IZ),N2BBASX)  
         ENDIF
      END DO

      CALL QEXIT('SORT2SAO')
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sort2sao */
      SUBROUTINE SORT2SAO_no_work(IDIMENS,PTRI,PMAT,FIRST,
     &                            IREP,ITIM,IPRINT)
C*****************************************************************************
C
C    PURPOSE: Sort the general 4comp. property operator
C             to the symmetry adapted AO-basis
C
C    On input:  PTRI - matrix of 4c prop. operator in unsorted basis
C               IDIMENS - tells whether it is triangular or square matrix
C
C               IREP, ITIM
C
C               IPRINT
C               FIRST 
C
C  NOTE: In the case of SQUARE matrix only the PMAT is entering ! All others
C        parameters in this case are dummy (PTRI only...)
C
C    On output: PMAT
C
C    Call from: PRPMSAO
C
C    Written by Miro ILIAS/Strasbourg/March 2006
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbbas.h"
      DIMENSION PMAT(N2BBASX,NZ),PTRI(IDIMENS,NZ)
      LOGICAL   FIRST(NZ), IS_SQUARE
C
      CALL QENTER('SORT2SAO_no_work')

      IF (IDIMENS.EQ.N2BBASX.AND.ITIM.EQ.0) THEN
        IS_SQUARE = .TRUE.
      ELSE IF (IDIMENS.EQ.NNBBASX.AND.(ITIM.EQ.1.OR.ITIM.EQ.-1)) THEN
        IS_SQUARE = .FALSE.
      ELSE
        write(lupri,*) 'SORT2SAO_no_work: N2BBASX,NNBBASX,IDIMENS=',
     &  N2BBASX,NNBBASX,IDIMENS
        CALL QUIT('SORT2SAO_no_work: wrong value of IDIMENS !')
      ENDIF

      DO IZ = 1,NZ
        IF(.NOT.FIRST(IZ)) THEN
          IQ = IPQTOQ(IZ,IREP)
          IH = IHQMAT(IQ,ITIM)
          IF(IH.EQ.1) THEN
            CALL DSPTSI(NTBAS(0),PTRI(1,IZ),PMAT(1,IZ))
            IF (IS_SQUARE) THEN
             CALL QUIT('SORT2SAO_no_work: wrong value of IDIMENS !') 
            ENDIF
          ELSEIF(IH.EQ.2) THEN
            CALL DAPTGE(NTBAS(0),PTRI(1,IZ),PMAT(1,IZ))
            IF (IS_SQUARE) THEN
             CALL QUIT('SORT2SAO_no_work: wrong value of IDIMENS !') 
            ENDIF
          ELSE
CMI   ... in this case PTRI is dummy, only PMAT exists
            IF (.NOT.IS_SQUARE) THEN
              CALL QUIT('SORT2SAO_no_work: wrong value of IDIMENS !') 
            ENDIF
          ENDIF
            CALL BUTOBS_no_work(PMAT(1,IZ),1)
        ELSE 
CMI ... zeroed out non-used quaternionic blocks of property matrixes...
          CALL DZERO(PMAT(1,IZ),N2BBASX)  
        ENDIF
      END DO

      CALL QEXIT('SORT2SAO_no_work')
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lowdin */
      SUBROUTINE LOWDIN(VMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     LOWDIN and LOWDI1 routines
C     ===========================
C
C     PURPOSE:
C        Construction of the transformation matrix for Lowdin's
C        canonical orthonormalization
C        NB! Not the symmetrical orthonormalization, which
C            often has been called the "Lowdin orthonormalization"
C
C     METHOD:
C        In this subroutine V+ is computed.  As the basis set is
C        real it reduces to the transpose:
C
C           V+ = TRN(V) = t(s**(-0.5))
C
C     Written by T.Saue
C     Last revision: Nov 23 1994 - tsaue
C     Last revision: Jan.2006 modified for pure nonrel.case...
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
C
#include "ccom.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "primit.h"
      LOGICAL LBIT
      DIMENSION VMAT(*),WORK(LWORK)
C
      CALL QENTER('LOWDIN')
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET('REAL',KSMAT,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,MXBBAS ,WORK,KFREE,LFREE)
C.....generate IRKB array if needed
      IF(.NOT.TWOCOMP.AND.LBIT(ISPHTR,2)) THEN
        CALL MEMGET('INTE',KIRKB,NPSHEL,WORK,KFREE,LFREE)
      ELSE
        KIRKB = KFREE
      ENDIF
C
C     Make Lowdin canonical matrix
C     ============================
C
      CALL LOWDI1(VMAT,WORK(KSMAT),WORK(KEIG),WORK(KIRKB),
     &            WORK(KFREE),LFREE,IPRINT)
C
C     Memory deallocation
C
      CALL MEMREL('LOWDIN.LOWDI1',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
      CALL QEXIT('LOWDIN')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lowdi1 */
      SUBROUTINE LOWDI1(VMAT,SMAT,EIG,IRKB,WORK,LWORK,IPRINT)
C*****************************************************************************
C     See SUBROUTINE LOWDIN for explanation
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
#include "ccom.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
      LOGICAL LBIT
      CHARACTER COMP(2)*1
      DIMENSION VMAT(*),SMAT(*),EIG(*),IRKB(*),WORK(LWORK)
C
      CALL QENTER('LOWDI1')
      IF(IPRINT.GE.5) CALL TITLER('Output from LOWDI1','*',103)
      COMP(1) = 'L'
      COMP(2) = 'S'
C
C     Get overlap matrix packed on fermion irreps
C     ===========================================
C
      CALL GTOVLT(SMAT,D1,IPRINT)
C
C     Generate IRKB array for use in RKBLOW
C     =====================================
C
      IF(LBIT(ISPHTR,2).AND..NOT.(TWOCOMP.OR.MDIRAC)) THEN
        CALL RKBREL(IRKB,IPRINT)
      ENDIF
C
C     Initialize Lowdin canonical matrix
C     ========================
C
      CALL DZERO(VMAT,N2BBASX)
      NBRP = 4/NZ
      ISOFF = 1
      IVOFF = 1
      CALL HEADER('Generating Lowdin canonical matrix:',-1)
      DO IFRP = 1,NFSYM
        JSOFF = ISOFF
        JVOFF = IVOFF
C
C       Large components
C       ----------------
C
        IC = 1
        NFORB(IFRP,IC) = 0
        IF(MDIRAC) THEN
          IP = IFRP
        ELSE
          IP = MOD(IFRP+IC,NFSYM) + 1
        ENDIF
        DO JSYM = 1,NBRP
          NBORB(JSYM,IFRP,IC) = 0
          IREP = JFSYM(JSYM,IP) - 1
          IF(NBBAS(IREP,IC).GT.0) THEN
            IF(LBIT(ISPHTR,IC)) THEN
              CALL SPHLOW(IREP,IC,SMAT(JSOFF),VMAT(JVOFF),
     &                    NFBAS(IFRP,0),NEFF,STOL(IC),EIG,
     &                    WORK,LWORK,IPRINT)
              NDEL  = NBBAS(IREP,IC) - NEFF
            ELSE
              CALL LOWGEN(SMAT(JSOFF),NBBAS(IREP,IC),
     &                    VMAT(JVOFF),NEFF,NFBAS(IFRP,0),
     &                    STOL(IC),IPRINT,EIG,WORK,LWORK)
              NDEL  = NBBAS(IREP,IC) - NEFF
              WRITE(LUPRI,'(3X,A1,3X,A3,3X,A,I10)') 
     &          COMP(IC),REP(IREP),'* Deleted: ',NDEL
            ENDIF
C
C           Zero columns of Lowdin canonical matrix that are to be deleted
C
            IF(NDEL.GT.0) THEN
              J = JVOFF + NFBAS(IFRP,0)*NEFF 
              DO I = 1,NDEL
                CALL DZERO(VMAT(J),NBBAS(IREP,IC))
                J = J + NFBAS(IFRP,0)
              ENDDO
            ENDIF
            JSOFF = JSOFF + (NFBAS(IFRP,0)+1)*NBBAS(IREP,IC)
            JVOFF = JVOFF + NFBAS(IFRP,0)*NEFF + NBBAS(IREP,IC)
            NBORB(JSYM,IFRP,IC) = NEFF
            NFORB(IFRP,IC) = NFORB(IFRP,IC) + NEFF
          ENDIF
        ENDDO
C
C       Small components
C       ----------------
C
        IC = 2
        NFORB(IFRP,IC) = 0
        IF(MDIRAC) THEN 
          IP = IFRP
        ELSE
          IP = MOD(IFRP+IC,NFSYM) + 1
        ENDIF
        DO JSYM = 1,NBRP
          NBORB(JSYM,IFRP,IC) = 0
          IREP = JFSYM(JSYM,IP) - 1
          IF(NBBAS(IREP,IC).GT.0) THEN
            IF(LBIT(ISPHTR,IC)) THEN
              IF(MDIRAC) THEN
                CALL SPHLOW(IREP,IC,SMAT(JSOFF),VMAT(JVOFF),
     &                      NFBAS(IFRP,0),NEFF,STOL(IC),EIG,
     &                      WORK,LWORK,IPRINT)
              ELSE
                CALL RKBLOW(IREP,IC,SMAT(JSOFF),VMAT(JVOFF),
     &                      NFBAS(IFRP,0),NEFF,STOL(IC),EIG,
     &                      IRKB,WORK,LWORK,IPRINT)
              ENDIF
              NDEL  = NBBAS(IREP,IC) - NEFF
            ELSE
              CALL LOWGEN(SMAT(JSOFF),NBBAS(IREP,IC),
     &                    VMAT(JVOFF),NEFF,NFBAS(IFRP,0),
     &                    STOL(IC),IPRINT,EIG,WORK,LWORK)
              NDEL  = NBBAS(IREP,IC) - NEFF
              WRITE(LUPRI,'(3X,A1,3X,A3,3X,A,I10)') 
     &           COMP(IC),REP(IREP),'* Deleted: ',NDEL
            ENDIF
C
C           Zero columns of Lowdin canonical matrix that are to be deleted
C
            IF(NDEL.GT.0) THEN
              J = JVOFF + NFBAS(IFRP,0)*NEFF 
              DO I = 1,NDEL
                CALL DZERO(VMAT(J),NBBAS(IREP,IC))
                J = J + NFBAS(IFRP,0)
              ENDDO
            ENDIF
            NBORB(JSYM,IFRP,IC) = NEFF
            NFORB(IFRP,IC) = NFORB(IFRP,IC) + NEFF
            JSOFF = JSOFF + (NFBAS(IFRP,0)+1)*NBBAS(IREP,IC)
            JVOFF = JVOFF + NFBAS(IFRP,0)*NEFF+NBBAS(IREP,IC)
          ENDIF
          NBORB(JSYM,IFRP,0) = NBORB(JSYM,IFRP,1) + NBORB(JSYM,IFRP,2)
        ENDDO
        NFORB(IFRP,0) = NFORB(IFRP,1) + NFORB(IFRP,2)
!       IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &    'LOWDI1: Lowdin canonical matrix (in ferm.symm.blocks) ',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** fermion corep ',IFRP,'/',NFSYM
           CALL PRQMAT(VMAT(IVOFF),NFBAS(IFRP,0),NFORB(IFRP,0),
     &           NFBAS(IFRP,0),NFORB(IFRP,0),1,IPQTOQ(1,0),LUPRI)
!       END IF
        ISOFF = ISOFF + NFBAS(IFRP,0)*NFBAS(IFRP,0)
        IVOFF = IVOFF + NFBAS(IFRP,0)*NFORB(IFRP,0)
      ENDDO

      CALL QEXIT('LOWDI1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgetint */
      SUBROUTINE RGETINT(WORD,OP1INT,DOINT,NCOMP,IORDER,NPQUAD,
     &                  MXFORM,TOFILE,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for getting one-electron property integrals
C     Written by T.Saue October 1994
C     Last revision : Nov 15 1994
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "shells.h"
#include "symmet.h"
#include "nuclei.h"
      DIMENSION OP1INT(*),WORK(LWORK)
      CHARACTER WORD*7,MXFORM*6,LABINT(3*MXCOOR)*8
      LOGICAL DOINT(2,2),PROPRI,TRIANG,TOFILE
C
      CALL QENTER('RGETINT')
#include "memint.h"
C
      MAXTYP = 3*MXCOOR
      IF ((WORD.EQ.'ELFGRDC').OR.(WORD.EQ.'ELFGRDS')) THEN
         NINTAD = 9*NUCIND*(MAXREP+1)
      ELSE
         NINTAD = MAXTYP
      END IF
      PROPRI   = .FALSE.
      TRIANG   = .TRUE.
C
C     Memory allocation
C
      CALL MEMGET('INTE',KINTRP,MAXTYP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINTAD,NINTAD,WORK,KFREE,LFREE)
C
C     Calculate integrals in abacus part
C
      CALL PR1INT_1(WORK(KFREE),LFREE,WORK(KINTRP),WORK(KINTAD),
     &            LABINT,WORD,IORDER,NPQUAD,TRIANG,
     &            PROPRI,IPRINT,OP1INT,NCOMP,TOFILE,MXFORM,DOINT)
C
C     Memory deallocation
C
      CALL MEMREL('RGETINT.PR1INT_1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('RGETINT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vfsort */
      SUBROUTINE VFSORT(NLEFF,NSEFF,VMAT,VBUF,WORK,LWORK)
C*****************************************************************************
C
C     Sort Lowdin canonical matrix on fermion symmetries
C
C     Written by T.Saue Sep 12 1995
C     Last revision: Sep 12 1995
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C Used from COMMON blocks:
C  DCBBAS: NTBAS(0)
C
#include "dcbbas.h"
      DIMENSION VMAT(*),VBUF(*),WORK(LWORK)
      CALL QENTER('VFSORT')
#include "memint.h"
C
C     Memory allocation
      CALL MEMGET('INTE',KSORT,NTBAS(0),WORK,KFREE,LFREE)
C
      CALL VFSOR1(NLEFF,NSEFF,VMAT,VBUF,WORK(KSORT))
C     Memory release
      CALL MEMREL('VFSORT.VFSOR1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('VFSORT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vfsor1 */
      SUBROUTINE VFSOR1(NLEFF,NSEFF,VMAT,VBUF,ISORT)
C*****************************************************************************
C
C     Sort Lowdin canonical matrix on fermion symmetries
C
C     Written by T.Saue Sep 12 1995
C     Last revision: Sep 12 1995
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcborb.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION VBUF(NBAST,*)
#else
      DIMENSION VBUF(NTBAS(0),*)
#endif
      DIMENSION VMAT(*),ISORT(*),NVEC(2)
C
C
C     Fermion symmetry 1
C     ==================
C
      CALL IZERO(NFORB,6)
      IOFF    = 1
      N2      = 0
      NVEC(1) = NLEFF
      NVEC(2) = NSEFF
      JJ = 0
      DO IC = 1,2
        DO J = 1,NVEC(IC)
          JJ = JJ + 1
          IND = IDAMAX(NTBAS(0),VBUF(1,JJ),1)
          IF(IND.LE.NFBAS(1,0)) THEN
            NFORB(1,IC) = NFORB(1,IC) + 1
            CALL DCOPY(NFBAS(1,0),VBUF(1,JJ),1,VMAT(IOFF),1)
            IOFF = IOFF + NFBAS(1,0)
          ELSE
            N2          = N2 + 1
            NFORB(2,IC) = NFORB(2,IC) + 1
            ISORT(N2)   = JJ
          ENDIF
        ENDDO
      ENDDO
C
C     Fermion  irrep 2
C     ================
C
      DO J = 1,N2
        CALL DCOPY(NFBAS(2,0),VBUF(NFBAS(1,0)+1,ISORT(J)),1,
     &             VMAT(IOFF),1)
        IOFF = IOFF + NFBAS(2,0)
      ENDDO
      NFORB(1,0) = NFORB(1,1)+NFORB(1,2)
      NFORB(2,0) = NFORB(2,1)+NFORB(2,2)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck freemt */
      SUBROUTINE FREEMT(FOCK,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for the construction of the one-electron free particle matrix
C     in SO-basis
C     Trond Saue, October 1995 - Oslo
C     Last revision: May 21 1996 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION FOCK(*),WORK(LWORK)
#include "dcbgen.h"
#include "cbihr1.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      CALL QENTER('FREEMT')
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET('REAL',KFTRI ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN,   NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST,NZ        ,WORK,KFREE,LFREE)
C
C     ****** Construct Fock matrix ******
C
      CALL FREEM1(FOCK,WORK(KFTRI),WORK(KONEIN),WORK(KFIRST),
     &              WORK(KFREE),LFREE,IPRINT)

C     Memory deallocation
      CALL MEMREL('FREEMT.FREEM1',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('FREEMT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck freem1 */
      SUBROUTINE FREEM1(FOCK,FTRI,OP1INT,FIRST,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Driver for the construction of the one-electron free particle matrix
C     in SO-basis
C     Trond Saue, October 1995 - Oslo
C
C-----------------------------------------------------------------------------
C
C     The free particle matrix has the following form in the C1 case:
C
C             La      Sa      Lb      Sb
C     La  |   0    | -icDz |   0   | -icD- |
C     Sa  | -icDz  |  -C2  | -icD- |   0   |
C     Lb  |   0    | -icD+ |   0   | -icDz |
C     Sb  | -icD+  |   0   | -icDz |  -C2  |
C
C
C     The integrals:
C        S  - overlap integrals                  symmetric
C        C2 - (2c^2)S                            symmetric
C        D  - dipole velocity integrals          antisymmetric
C            D+ = Dx + iDy
C            D- = Dx - iDy
C     The ordering of components is chosen to reveal the following
C     matrix structure which is due to time reversal symmetry:
C
C     |  A  |  B  |
C     | -B* |  A* |
C
C     In this routine the full matrices A and B are constructed:
C
C            A(real)   |   A(imag)   |    B(real)   |  B(imag)
C     ------------------------------------------------------------
C           La      Sa     La    Sa  |  Lb     Sb     Lb     Sb
C     La |   0  |   0  |   0  | -cDz |    0  | -cDy |   0  | -cDx |
C     Sa |   0  |  -C2 | -cDz |   0  |  -cDy |   0  | -cDx |   0  |
C
C     It should be noted that Ar is symmetric, whereas Ai,Br and Bi
C     are antisymmetric.
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbgen.h"
C
      LOGICAL FIRST(NZ)
      DIMENSION FOCK(N2BBASX,NZ),FTRI(NNBBASX,NZ),OP1INT(NNBBASX),
     &          WORK(LWORK)
      CALL QENTER('FREEM1')

C
C     Initialize
C
#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
      CALL LSET(NZ,.TRUE.,FIRST)
C
C     Beta matrix
C     ===========
C
      IOP = IPBETAMT
      CALL PRPMAO(LU1INT,IOP,NNBBASX,FTRI,
     &            OP1INT,FIRST,IPRHAM)
C
C     Kinetic energy
C     ==============
C
      IOP = IPKINERG
      CALL PRPMAO(LU1INT,IOP,NNBBASX,FTRI,
     &            OP1INT,FIRST,IPRHAM)
C
      CLOSE(LU1INT,STATUS = 'KEEP')
C
C     Make full matrix
C     ================
C
      DO IZ = 1,NZ
        IF(FIRST(IZ)) THEN
          CALL DZERO(FOCK(1,IZ),N2BBASX)
        ELSE
          IQ = IPQTOQ(IZ,0)
          IH = IHQMAT(IQ,1)
          IF(IH.EQ.1) THEN        
            CALL DSPTSI(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DAPTGE(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ENDIF
        ENDIF
      ENDDO
      IF(IPRHAM.GE.5) THEN
        CALL HEADER(
     &   'FREEM1:Unsorted total one-electron Fock matrix',-1)
        CALL PRQMAT(FOCK,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &              IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Reindex to sorted basis
C     =======================
C
      CALL BUTOBS(FOCK,NZ,WORK,LWORK)
C
C     Output section
C     ==============
C
      IF(IPRHAM.GE.4) THEN
        CALL HEADER(
     &  'FREEM1: One-electron (symm.blocked) free particle matrix',-1)
        DO 10 I = 1,NFSYM
          IF(NFBAS(I,0).EQ.0) GOTO 10
          WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(FOCK(I2BASX(I,I)+1,1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   10   CONTINUE
      ENDIF
C
      CALL QEXIT('FREEM1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gmotra */
      SUBROUTINE GMOTRA(GOFAST)
C*****************************************************************************
C
C PURPOSE:
C     This matrix provides initial transformation matrix from
C     SO- to MO-basis. With no restrictions on kinetic balance
C     the Lowdin canonical matrix is used. With restricted balance, the
C     Lowdin canonical matrix is used to diagonalize the free-particle
C     matrix. Unphysical solutions are then deleted from the set
C     of vectors and the remaining vectors used as transformation matrix.
C
C MI,HJAaJ/2002-3
C     For the two-component relativistic calculations
C     apply the Barysz-Sadlej-Snijders transformation
C     on the Dirac one-electron Hamiltonian.
C     Save its elements and the picture change 
C     transformation matrix as well.
C MI/2006
C     For the two-component nonrelativistic 
C     calculations run the "URKBAL" branch as only
C     Lowdin canonical matrix is needed. It's assumed that small
C     components were eliminated in HERINP.
C
C     Called from: PSIDHF (when BSS2DC or DC2BSS mode is activated)
C                  PAMSET (at the begining of the DIRAC run, and second time
C                          when pure two-comp. scheme+BSS is wished)
C
C     Written by T.Saue, October 1995
C     Last revision May 21 1996 - tsaue
C                    2006 - M.Ilias
C
C*****************************************************************************
      use memory_allocator
      use x2cmod_cfg
      use dirx2c
      use x2c_pct_ao,       only:
     &    pctrafo_all_op_driver
      use x2c_fio
      use x2c_pct_mo_coefficients, only: 
     &    x2c_molecular_mean_field_mos
!
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
C Used from COMMON blocks:
C   DCBBAS: N2BBASX
C   PAMIOU: LUTMAT
#include "dcbgen.h"
#include "dcbcls.h"
#include "dcbham.h"
#include "cbihr1.h"
#include "cbirea.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
! needed for restart from X2CMAT
#include "dcbxpr.h"
#include "maxaqn.h"
#include "ccom.h"
#include "dgroup.h"
      real(8), allocatable :: WORK(:)
      logical, intent(in)  :: GOFAST
!     local scratch
      logical              :: TOBE, FNDLAB, FNDLAB12, gmotra_case3
      logical              :: is_2c_rel_run
      integer              :: is_defining_h1mat
      integer              :: ndim_lowmat
      integer              :: lulmat = 99

      CALL QENTER('GMOTRA')

      IF(GOFAST.AND.(.NOT.(BSS.or.x2c))) GOTO 10

      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in GMOTRA')
      KFRSAV = KFREE
C
!     set file unit number for X2CMAT
      lux2c = 40

!     set local 2c-relativistic run flag
      is_2c_rel_run = bss.or.x2c

!     logical flag for keeping track of whether case 3 in gmotra has
!     been entered
      gmotra_case3 = .false.

C     By default we do not partition the space. This may be
C     modified if we have linear symmetry or approximate the
C     Hamiltonian : Levy-Leblond, Spinfree, ..
C
      SUB_BL = .FALSE.
C
C     Make one-center projection matrix
C     =================================
C
      IF (ONECAP.AND.(INTV1C.EQ.1.AND.ICTLV1C(2).EQ.1)) 
     &    CALL MK1CEX (WORK,KFREE,LFREE,IPRGEN)

!     we use N2BBASXQ instead of N2BBASX for BSS...to be checked
      CALL MEMGET('REAL',KVMAT,N2BBASXQ,WORK,KFREE,LFREE)
C
C     Get Lowdin canonical orthonormalization matrix
C     ==============================================
 
      inquire(exist=tobe,file='LOWDMAT')
!     delete existing LOWDIN matrix file when performing a geometry
!     optimization - it could be the one from a previous geometry.
!     stefan - june 2011
      if(tobe.and.optimi)then
        call opnfil(lulmat,'LOWDMAT','OLD','gmodel')
        close(lulmat,status='delete')
        tobe = .false.
      end if

      call opnfil(lulmat,'LOWDMAT','UNKNOWN','gmotra')
      rewind lulmat

      if(tobe)then
        call readi(lulmat,1,ndim_lowmat)
        call readt(lulmat,ndim_lowmat,work(kvmat))
      else
        CALL TIMER('START ',TIMSTR,TIMEND)
        CALL LOWDIN(WORK(KVMAT),WORK(KFREE),LFREE,IPRHAM)
        CALL TIMER('Lwdn_a',TIMSTR,TIMEND)

!       save LOWDIN matrix on file
        ndim_lowmat = 0
        do i = 1, nfsym
          ndim_lowmat = ndim_lowmat + nforb(i,0)*nfbas(i,0)
        end do
        call writi(lulmat,1,ndim_lowmat)
        call writt(lulmat,ndim_lowmat,work(kvmat))
      end if

      close(lulmat,status='keep')

!     release memory when we invoke the x2c module
      if(x2c)then 
       call memrel('x2cfre',work,kwork,kwork,kfree,lfree)
       call dealloc(work)
      end if

C     Unrestricted kinetic balance; use Lowdin canonical matrix
C     =========================================================
C
C     LV : The unrestricted procedure should be revised so that
C     the large components are also sorted on boson irreps !!!
C
      IF(.NOT.URKBAL.AND.NTBAS(2).EQ.0.AND.
     &   .NOT.(TWOCOMP.OR.is_2c_rel_run)) THEN
        URKBAL = .TRUE.
        WRITE(LUPRI,'(A)') 
     &  '*** WARNING *** : No small bases. Setting URKBAL = T !'
      ENDIF
C
C     ** Branching point **
C=======================================================================
C
!       ---------------------------------------------------
!       Case 1: Unrestricted kinetic balance 
!               or 2-component non-relativistic calculation
!       ---------------------------------------------------
      IF(URKBAL.OR.(TWOCOMP.AND..NOT.is_2c_rel_run))then
        NZT = 1
        KTMAT = KVMAT
        DO IFRP = 1,NFSYM
          NESH(IFRP)  = NFORB(IFRP,1)
          NPSH(IFRP)  = NFORB(IFRP,2)
          NORB(IFRP)  = NESH(IFRP) + NPSH(IFRP)
          NTMO(IFRP)  = NORB(IFRP)
        ENDDO
        IF (TWOCOMP) THEN
         SSMTRC = D0
         MC = 1
        ENDIF
        IF (IPRHAM.GE.5) THEN
          write(lupri,*) 'GMOTRA: branch URKBAL'
          write(lupri,*) 'NESH(IFRP),NPSH(IFRP)',
     &     (NESH(IFRP),NPSH(IFRP),IFRP=1,NFSYM)
        ENDIF
!       --------------------------------------------
!       Case 2: relativistic 2-component Hamiltonian
!       --------------------------------------------
      ELSE IF(is_2c_rel_run)THEN
!
        CALL TIMER('START ',TIMSTR,TIMEND)

        TOBE=.FALSE.
        IF(GOFAST) THEN
          if(bss)then
            inquire(file='BSSMAT',exist=tobe)
          else if(x2c)then
            inquire(file='X2CMAT',exist=tobe)
          end if
        ENDIF

!       case 2.1 - step a: restart 2c-run from file
!       -------------------------------------------
        IF(TOBE) THEN

C         BSSMAT or X2CMAT available:
C         1) open BSSMAT/X2CMAT and remove all property integrals
C         2) do picture change transformation of integrals on AOPROPER

          if(bss)then

            WRITE(LUPRI,'(A,A)') '* Info: DIRAC found BSSMAT ',
     &           '(BSS stuff) ... skipping generation of H2C'
            LUBSS = 40
            CALL OPNFIL(LUBSS,'BSSMAT','OLD','GMOTRA')
            IF (FNDLAB('P2C_0001',LUBSS)) THEN
              BACKSPACE(LUBSS)                   
              CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
            ENDIF
            CLOSE(LUBSS,STATUS='KEEP')
            CALL PROP2BSS(WORK(KVMAT),IPRHAM,WORK(KFREE),LFREE)

          else if(x2c)then

            write(lupri,'(a,a)') '* Info: DIRAC found X2CMAT ',
     &           '(X2C stuff) ... skipping generation of H2C'
            open(lux2c,file='X2CMAT',status='old',form='unformatted',
     &           access='sequential',action='readwrite',
     &           position='rewind')

!           find Hamiltonian on file property integral... and step back
            if(.not.x2cmod_fragment_x2c)then
              i2cofk =  0
!             check for hamiltonian integrals on file X2CMAT and set i2cofk
              if(fndlab12('h12cAOA   00',lux2c))then
                if(noamfi)then
                  write(lupri,'(a)') 
     &' *** restart error: try to read h1+amfi and noamfi == .true. ***'
                  call quit(
     &'*** error in gmotra: mismatch between X2CMAT data and input ***')
                end if
!               set i2cofk to |3|: since we restart, it does not matter
!               how the AMFI part was derived.
                i2cofk = 3
              else if(fndlab12('h12cAOn   00',lux2c))then
                i2cofk = 1
              else
                  write(lupri,'(a)') 
     &'*** restart error: no hamiltonian integrals stored on X2CMAT ***'
                  call quit(
     &'*** error in gmotra: no h1 data on file X2CMAT ***')
              end if
            else
!             nothing to do here for now: x2cmod_fragment_x2c is always false.
            end if
 
!           find first property integral... and step back
            if(fndlab12('prpint2c0001',lux2c))then
              backspace(lux2c)
              call x2c_write('thisisadummy',work,-1,lux2c)
            end if

            call pctrafo_all_op_driver(ntbas(0),
     &                                 ntbas(1),
     &                                 nfbas(1,0),
     &                                 nfbas(1,1),
     &                                 nfsym,
     &                                 nz,
     &                                 ipqtoq,
     &                                 i2basx,
     &                                 lux2c,
     &                                 lu1int,
     &                                 jbtof,
     &                                 nprps,
     &                                 lwork,
     &                                 dosphe,
     &                                 .false.,
     &                                 onesys, 
     &                                 -1,
     &                                 -1,
     &                                 'undef ',
     &                                 'undef ',
     &                                 iprham)

          close(lux2c,status='keep')

          end if

!       case 2.2 - step a: generate 2c-Hamiltonian
!       ------------------------------------------
        ELSE

          if(bss)then
!           bss module (if asked for it can also perform an X2C transformation)
            CALL MAKE_H2C(WORK(KVMAT),WORK(KFREE),LFREE)

          else if(x2c)then

!           put in an interface routine to be called inside x2c_driver
            if(mdirac)then
              NZT = 1
              DO IFRP = 1,NFSYM
                NESH(IFRP)  = NFORB(IFRP,1)
                NPSH(IFRP)  = NFORB(IFRP,1)
                NORB(IFRP)  = NESH(IFRP)    * 2
                NTMO(IFRP)  = NORB(IFRP)
              END DO
            end if

!           optimized x2c module
            call x2c_driver(x2cmod_fragment_x2c,.true.)
!           re-allocate 
            call alloc(WORK,LWORK,id='WORK in GMOTRA - part 2')
#include "memint.h"

          end if

        ENDIF ! restart X2C or re-do X2C switch

        CALL TIMER('mk_h2c',TIMSTR,TIMEND)

!       
!       case 2 - step b: generate Lowdin canonical orthonormalization matrix in 2c-mode
!       -------------------------------------------------------------------------------
        if(.not.(twocompbss.or.x2c).or.x2cmod_mmf)then

          if(x2cmod_mmf)then
            CALL MEMGET('REAL',KTMAT,N2BBASXQ,WORK,KFREE,LFREE)
            call dzero(work(KTMAT),N2BBASXQ)
            open(lux2c,file='X2CMAT',status='old',form='unformatted',
     &           access='sequential',action='readwrite',
     &           position='rewind')
            ndim_lowdinmat   = 0
            do i = 1, nfsym
              ndim_lowdinmat = ndim_lowdinmat + nfbas(i,0)*nesh(i)
            end do

            if(linear)then
              gmotra_case3 = .true.
            else
              call x2c_read('2cLWAOg   10',work(KTMAT),
     &                      ndim_lowdinmat,lux2c)
            end if
            
            close(lux2c,status='keep')
            sub_bl = .false.
          else
!           WORK(KTMAT) carries already the large-electronic part of the Lowdin canonical matrix
            KTMAT = KVMAT
          end if
        else
C
          if(.not.x2c)
     &    CALL MEMREL('BSS.2c',WORK,KWORK,KWORK,KFREE,LFREE)

          TWOCOMP = .TRUE.
          RDINPC = .FALSE.
          IPREAD = IPREAD - 5
          CALL READIN(.FALSE.)
          CALL SETDC1(IPREAD)
          IPREAD = IPREAD + 5
C.........reset NEWPRP to re-initialize AOPROPER
          NEWPRP = .TRUE.
          CALL LSET(NPRPCLS,.TRUE.,CLSCAL)
          CALL ONEGEN(WORK,LWORK)

          CALL MEMGET('REAL',KTMAT,N2BBASXQ,WORK,KFREE,LFREE)
C.........do canonical orthonormalization
          CALL TIMER('START ',TIMSTR,TIMEND)
          CALL LOWDIN(WORK(KTMAT),WORK(KFREE),LFREE,IPRHAM)
          CALL TIMER('Lwdn_b',TIMSTR,TIMEND)

!         save new 2c-LOWDIN matrix on file
          call opnfil(lulmat,'LOWDMAT','UNKNOWN','gmotra')
          rewind lulmat
          ndim_lowmat = 0
          do i = 1, nfsym
            ndim_lowmat = ndim_lowmat + nforb(i,0)*nfbas(i,0)
          end do
          call writi(lulmat,1,ndim_lowmat)
          call writt(lulmat,ndim_lowmat,work(kvmat))
          close(lulmat,status='keep')

!         reset orbital dimensions
          NZT = 1
          DO IFRP = 1,NFSYM
            NESH(IFRP)  = NFORB(IFRP,1)
            NPSH(IFRP)  = NFORB(IFRP,2)
            NORB(IFRP)  = NESH(IFRP) + NPSH(IFRP)
            NTMO(IFRP)  = NORB(IFRP)
          ENDDO
          SSMTRC = D0
          MC = 1

        end if

        IF (IPRHAM.GE.5) THEN
          write(lupri,*) 'GMOTRA: relativistic 2c-branch'
          write(lupri,*) 'NESH(IFRP),NPSH(IFRP)',
     &     (NESH(IFRP),NPSH(IFRP),IFRP=1,NFSYM)
        ENDIF

!
!     -------------------------------------------------------------------
!     Case 3: Transformation matrix with modified Dirac equation embedded
!     -------------------------------------------------------------------
      ELSE
        CALL MEMGET('REAL',KTMAT,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIG ,NTBAS(0),WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KTBUF,N2BBASXQ,WORK,KFREE,LFREE)

!       defining h1 matrix control variable (in the X2C module) -
!       set here for consistency to the maximum value
        is_defining_h1mat = 4

!        in: LOWDIN matrix                   --> KVMAT

!       out: AO2MO                           --> KTMAT
!       out: AO2MO SL resorted               --> KVMAT
!       out: MO2MO linear symmetry adapted   --> KTBUF (if linear == .true.; nothing otherwise)
        CALL MODHAM(WORK(KTMAT),WORK(KVMAT),WORK(KTBUF),WORK(KEIG),
     &              .not.(VEXTPJ.OR.FREEPJ.OR.LEVYLE),linear,
     &              is_defining_h1mat,WORK(KFREE),LFREE)

        CALL MEMREL('MODHAM',WORK,KEIG,KEIG,KFREE,LFREE)

        IF (IPRHAM.GE.5) THEN
          write(lupri,*) 'GMOTRA: branch MODHAM'
          write(lupri,*) 'NESH(IFRP),NPSH(IFRP)',
     &     (NESH(IFRP),NPSH(IFRP),IFRP=1,NFSYM)
        ENDIF

        gmotra_case3 = .true.

      ENDIF
C
C     ** End of branching **
C=======================================================================
C
      CALL MEMGET('REAL',KTBUF,N2BBASXQ,WORK,KFREE,LFREE)
      N2TMT   = 0
      N2TMOTQ = 0
      DO IFRP = 1,NFSYM
        NESHMF(IFRP) = NESH(IFRP)
        NPSHMF(IFRP) = NPSH(IFRP)
        NISHMF(IFRP) = NISH(IFRP)
        NOCCMF(IFRP) = NOCC(IFRP)
        I2TMT(IFRP)  = N2TMT
        N2TMT        = N2TMT + NFBAS(IFRP,0)*NORB(IFRP)*NZT
        I2TMOT(IFRP) = N2TMOTQ
        N2TMO(IFRP)  = NTMO(IFRP)*NTMO(IFRP)
        N2TMOTQ      = N2TMOTQ + N2TMO(IFRP)*NZ
      ENDDO
C
C     Copy information from NBORB to the subblock array if we have not
C     partioned the orbital space otherwise.
C
      IF (.NOT.SUB_BL) CALL INISUB
      IF (SPINFR.or.(x2c.and.x2cmod_mmf)) THEN
         SUB_BL = .TRUE.
      ENDIF

      if(linear.and..not.gmotra_case3)then
C.......SUB_BL will be set true if linear symmetry works
        CALL DCOPY(N2TMT,WORK(KTMAT),1,WORK(KTBUF),1)
!       in:  KTMAT + KTBUF
!       out: KTMAT
        CALL MEMGET('REAL',KEIG ,NTBAS(0),WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KTBUF2,N2BBASXQ,WORK,KFREE,LFREE)

!        in: AO2MO                         --> KTBUF

!       out: AO2MO linear symmetry adapted --> KTMAT
!       out: AO2MO                         --> KTBUF
!       out: MO2MO linear symmetry adapted --> KTBUF2 (discarded if not x2c-mmf)
        CALL GLINSM(WORK(KTMAT),WORK(KTBUF),WORK(KTBUF2),WORK(KEIG),
     &              WORK(KFREE),LFREE)

        CALL MEMREL('MODHAM',WORK,KEIG,KEIG,KFREE,LFREE)
      else if(linear)then
        if(x2c.and.x2cmod_mmf)then
          open(lux2c,file='X2CMAT',status='old',form='unformatted',
     &         access='sequential',action='readwrite',
     &         position='rewind')
          call x2c_read('2cLWAOl   10',work(KTMAT),ndim_lowdinmat,lux2c)
          close(lux2c,status='keep')
        else
          open(99,file='AOMOlin',status='old',form='unformatted',
     &         access='sequential',action='readwrite',
     &         position='rewind')
!         skip rec # 1
          read(99)
          call readt(99,N2TMT,work(KTMAT))
!         file is needed in the 4c-2c transformation after the 4c-SCF
!         in order to ensure numerical stability (avoiding double
!         diagonalization, etc.)
          if(x2cmod_mmf)then
            close(99,status='keep')
          else
            close(99,status='delete')
          end if
        end if
      end if
C
C     Construct right index transformed overlap matrix Saomo = Sao T
C     Since T(dagger) Sao T = I, Saomo is the inverse of T(dagger)
C
!      in: AO2MO (linear symmetry adapted if linsym = T)--> KTMAT
!     out: right index transformed overlap matrix Saomo --> KTBUF
      CALL MKSAOMO(WORK(KTBUF),WORK(KTMAT),.FALSE.,
     &              IPRHAM,WORK(KFREE),LFREE)
C
C     If requested, reduce transformation matrix by projection with
C     respect to a set of fragment orbitals
C
      IF(PROJEC) THEN
        CALL PSIPRJ(WORK(KTMAT),WORK(KTBUF),IPRONE,WORK,KFREE,LFREE)
        N2TMOTQ = 0
        DO IFRP = 1,NFSYM
          I2TMOT(IFRP) = N2TMOTQ
          N2TMO(IFRP) = NTMO(IFRP)*NTMO(IFRP)
          N2TMOTQ = N2TMOTQ + N2TMO(IFRP)*NZ
          NORB(IFRP) = NTMO(IFRP) + NFRO(IFRP)
        ENDDO
        CALL SETDC2(0)
      ENDIF
C
C     Write transformation matrices
C
      CALL OPNFIL(LUTMAT,'AOMOMAT','UNKNOWN','GMOTRA')
!     AO2MO transformation matrix (linear symmetry adapted if linsym == .true.)
      CALL WRITT(LUTMAT,N2TMT,WORK(KTMAT))
!     right index transformed overlap matrix Saomo
      CALL WRITT(LUTMAT,N2TMT,WORK(KTBUF))
!     orbital dimensions, etc --> "runtime" data
      CALL ORBDAT

      CLOSE(LUTMAT,STATUS='KEEP')
!
      if(x2c.and.x2cmod_mmf)then

!       perform the final step for the molecular mean-field procedure within the X2C-module
!       ------------------------------------------------------------------------------------

!       purpose: transform the 2c-MO-coefficients from the orthonormal reference basis to AO basis
!                and renew the coefficient file DFCOEF
!                a potential issue could be frozen orbitals - test me!
        call memget('REAL',kcmomat,ncmotq,work,kfree,lfree)
        call memget('REAL',kscratch,norbt,work,kfree,lfree)
        call dzero(work(ktbuf),   n2bbasxq)
        call dzero(work(kcmomat), ncmotq)
        call dzero(work(kscratch),norbt)
        call x2c_molecular_mean_field_mos(work(ktmat),work(ktbuf),
     &                                    work(kcmomat),work(kscratch),
     &                                    i2tmt,i2tmot,icmoq,ntmo,
     &                                    nfbas,nesh,nz,nzt,nfsym,
     &                                    nesht,ipqtoq,iprham)
      end if ! x2cmod_mmf .equiv. .true.

!     final output
!     ------------
      CALL TITLER('Orbital dimensions','*',125)

      IF ( NFSYM .EQ. 1 ) THEN
         WRITE(LUPRI,'(A,I5)')
     &        'No. of positive energy orbitals (NESH): ',NESH(1)
         WRITE(LUPRI,'(A,I5)')
     &        'No. of negative energy orbitals (NPSH): ',NPSH(1)
         WRITE(LUPRI,'(A,I5)')
     &        'Total no. of orbitals           (NORB): ',NORB(1)
      ELSE
         WRITE(LUPRI,'(35X,5(A5,1X,I1,1X))') 
     &        ('Irrep',I,I=1,NFSYM),'Sum'
         J = 0
         DO I = 1, NFSYM
            J = J + NESH(I)
         END DO
         WRITE(LUPRI,'(A,5(I5,3X))')
     &        'No. of electronic orbitals (NESH): ',
     &        (NESH(I),I=1,NFSYM),J
         J = 0
         DO I = 1, NFSYM
            J = J + NPSH(I)
         END DO
         WRITE(LUPRI,'(A,5(I5,3X))')
     &        'No. of positronic orbitals (NPSH): ',
     &        (NPSH(I),I=1,NFSYM),J
         J = 0
         DO I = 1, NFSYM
            J = J + NORB(I)
         END DO
         WRITE(LUPRI,'(A,5(I5,3X)/)')
     &        'Total no. of orbitals      (NORB): ',
     &        (NORB(I),I=1,NFSYM),J
      END IF
C
C     Memory deallocation
      CALL MEMREL('GMOTRA',WORK,KWORK,KWORK,KFREE,LFREE)
      if(allocated(work)) call dealloc(work)
 10   CONTINUE

      CALL QEXIT('GMOTRA')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck modham */
      SUBROUTINE MODHAM(TMAT,VMAT,TBUF,EIG,SLGROUP,LSYMAD,
     &                  is_defining_h1mat,WORK,LWORK)
C*****************************************************************************
C
C     Generate transformation matrix to MO-basis with the
C     quaternion modified Dirac equation embedded.
C     Ref: Visscher and Saue, JCP 113 (2000) 3996.
C
C     and for S-L resorting of columns (flag SLGROUP).
C
C     On input:
C   ============= 
C
C     VMAT - 4comp Lowdin transformation matrix for canonical orthonormalization
C            possibly with sperical and RKB transformations embedded.
C            (is destroyed!)
C
C     SLGROUP - if .true. (=Trond's 2comp.approach-one step BSS transformation)
C                resort columns of the transformation matrix into "S-L" order)
C
C     LSYMAD  - if .true., apply the LINSYM procedure for adaptation
C               of the transformation matrix for the linear symmetry - needed in connection
C               to 4c->2c transition (.DO4C2C keyword)
C
C     is_defining_h1mat - control parameter for the defining h1 matrix 
C                         in the X2C module wrt which the decoupling will be performed.
C        
C     On output:
C   ==============
C
C  TMAT - total 4c AO -> MO transformation matrix
C
C  VMAT - if SLGROUP=.true., then it contains the standard (ALWAYS non-linear symmetry) transformation
C         matrix with SL-resorted columns
C
C  TBUF - if linear symmetry (LSYMAD) is on, it contains the transformation matrix from non-linear MO
C         to the linear MO
C     
C   Written by T.Saue, October 1995
C   Modified by L. Visscher, Jan 20, 1998 to include Levy-Leblond 
C            M. Ilias, 2005, to include LSYMAD, SLGROUP
C            S. Knecht, 2010, to include is_defining_h1mat
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, DM2 = -2.0D0, D0 = 0.0D0)
C
C Used from COMMON blocks
C  CBIPAM:
C  DCBBAS:
C  DCBORB:
C  DGROUP:
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "pgroup.h"
#include "cbihr1.h"
C
      DIMENSION TMAT(*),TBUF(*),VMAT(*),EIG(*),WORK(LWORK)
      LOGICAL SLGROUP,LSYMAD
      integer is_defining_h1mat
      CALL QENTER('MODHAM')
#include "memint.h" 
C
CMI  ... check of keywords..
      IF (SLGROUP.AND.(VEXTPJ.OR.FREEPJ.OR.LEVYLE)) 
     &  CALL QUIT('MODHAM: SLGROUP does not correspond with'//
     &  ' VEXTPJ/FREEPJ/LEVYLE flags!')
!
      FAC = DM2*CVAL*CVAL
      NSYM = 4/NZ
      print *, 'slgroup', slgroup
      if(.not.MDIRAC)THEN
C
C       Generate free particle (or bare nucleus) matrix in Lowdin-basis
C       ===============================================================
C
        print *, 'blubb',vextpj
        IF(VEXTPJ) THEN
          CALL ONEFCK(TBUF,IPRONE,WORK(KFREE),LFREE)
        ELSE
          CALL FREEMT(TBUF,IPRONE,WORK(KFREE),LFREE)
        ENDIF
        IOFT = 1
        IOFV = 1
        DO 10 I = 1,NFSYM
          IF(NFORB(I,0).EQ.0) GOTO 10
CMI        ... transform to the Lowdin basis "Theta"
          CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                              NFORB(I,0),NFORB(I,0),
     &         TBUF(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &         TMAT(IOFT),NFORB(I,0),NFORB(I,0),NZ,IPQTOQ(1,0),
     &         VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &         VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &         WORK,LWORK,IPRGEN)
          IF (IPRHAM.GE.6) THEN
             IF(VEXTPJ) THEN
               CALL HEADER('MODHAM: '//
     &'Bare nucleus Dirac matrix in orthonormal canonical Lowdin basis'
     &,-1)
             ELSE
               CALL HEADER('MODHAM: '//
     &'Free particle Dirac matrix in orthonormal canonical Lowdin basis'
     &,-1)
             ENDIF
             WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
             CALL PRQMAT(TMAT(IOFT),NFORB(I,0),NFORB(I,0),NFORB(I,0),
     &                   NFORB(I,0),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
          IOFT = IOFT + NFORB(I,0)*NFORB(I,0)*NZ
          IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)
   10   CONTINUE
C
        NZT = NZ
      END IF

      IF(LEVYLE) THEN
C
C       Generate transformation matrix for the Levy-Leblond equation
C       ============================================================
C
        CALL LL_EQ(TMAT,TBUF,VMAT)

      ELSE IF(MDIRAC)THEN      

C       The modified Dirac equation at the integral level (InteRest):
C       =============================================================
        NZT = 1

        imat = 1
        DO IFRP = 1,NFSYM
          NESH(IFRP)  = NFORB(IFRP,1)
          NPSH(IFRP)  = NFORB(IFRP,1)
          NORB(IFRP)  = NESH(IFRP)    * 2
          NTMO(IFRP)  = NORB(IFRP)
          lMAT        = NFBAS(ifrp,0)*NORB(ifrp)*NZ
          call resort_mdirac_vmat(tmat(imat),vmat(imat),
     &                            nfbas(ifrp,0),norb(ifrp),
     &                            nforb(ifrp,1),nforb(ifrp,1),nz)
          iMAT        = imat + NFBAS(ifrp,0)*NORB(ifrp)*NZ
        END DO

      ELSE      
C
C       The modified Dirac equation at the matrix level (Hermit):
C       =========================================================

        CALL MOD_DIREQ(TMAT,TBUF,VMAT,EIG,WORK(KFREE),LFREE)
      ENDIF   
C
C     Print section
C
      IF(LEVYLE.OR.VEXTPJ.OR.FREEPJ) THEN
        CALL HEADER('Output from MODHAM',-1)
        WRITE(LUPRI,'(A)') '* All positronic solutions deleted !'
      ELSE IF(.not.mdirac)then
        CALL HEADER('Output from MODHAM',-1)
        WRITE(LUPRI,'(A)') '* Applied strict kinetic balance !'
      ENDIF
!     IF(IPRHAM.GE.5) THEN
        CALL HEADER(
     &    'MODHAM: Actual transformation matrix 4c AO->MO ',-1)
        IMAT = 1
        DO 80 I = 1,NFSYM
        IF(NORB(I).GT.0) THEN
          WRITE(LUPRI,'(A,I1,A,I1)')
     &    '* Fermion ircop no.',I,'/',NFSYM
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
        ENDIF
 80     CONTINUE
!     ENDIF

CMI ...  For the BSS.RKB (ie 1step approach) do resort columns of transformation matrix so that
CMI     they are in the order S1 S2 ... S(nbsym) L1 L2 ... L(nbsym) instead of former S1 L1 S2 L2 ...
      IF(SLGROUP)THEN
        ISIZE = NFBAS(1,0)*NORB(1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NORB(2)*NZ 
!       is_defining_h1mat is the control variable for the defining h1
!       matrix wrt which the decoupling is performed in the X2C module.
!       see the following lsymad procedure for more comments - SK Aug 2010
!       print *,'is_defining_h1mat is...',is_defining_h1mat
        open(99,file='AOMOSLR',status='replace',form='unformatted',
     &       access='sequential',action='readwrite',
     &       position='rewind')
        if(is_defining_h1mat.gt.3.or.
     &    (is_defining_h1mat.eq.2.and.LSYMAD))then
!       if(is_defining_h1mat.gt.3)then
          CALL MEMGET('REAL',KBUF,ISIZE,WORK,KFREE,LFREE)
!         in:  AO2MO transformation matrix             --> TMAT
!         out: AO2MO transformation matrix             --> TMAT
!              AO2MO transformation matrix SL resorted --> WORK(KBUF)
          CALL SLSORT(TMAT,WORK(KBUF),IPRHAM,WORK(KFREE),LFREE)
          call writt(99,ISIZE, WORK(KBUF))
!         print *,'SLSORT branch 1...'
        else
          call dzero(VMAT,n2bbasxq)
!         in:  AO2MO transformation matrix             --> TMAT
!         out: AO2MO transformation matrix             --> TMAT
!              AO2MO transformation matrix SL resorted --> VMAT
          CALL SLSORT(TMAT,VMAT,IPRHAM,WORK(KFREE),LFREE)
          call writt(99,ISIZE, vmat)
!         print *,'SLSORT branch 2...'
        end if
        close(99,status='keep')
      ENDIF

!     adapt the MO part of the AO2MO transformation matrix wrt linear symmetry by 
!     using the eigenvectors of <j_z> in MO basis as transformation matrix
!     ---------------------------------------------------------------------------
      open(99,file='AOMOlin',status='replace',form='unformatted',
     &     access='sequential',action='readwrite',
     &     position='rewind')
!     remember: in linear symmetry NZ ==> 1
      ISIZE  = NFBAS(1,0)*NORB(1)*NZ
      ISIZE2 = NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2)then
        ISIZE  = ISIZE  + NFBAS(2,0) * NORB(2)*NZ
        ISIZE2 = ISIZE2 + NORB(2)    * NORB(2)*NZ
      END IF

      IF(LSYMAD)THEN
        CALL MEMGET('REAL',KBUF1,ISIZE,WORK,KFREE,LFREE)
        CALL DCOPY(ISIZE,TMAT,1,VMAT,1)
        CALL DCOPY(ISIZE,TMAT,1,WORK(KBUF1),1)
        IF(IPRHAM.GE.5) THEN
          CALL HEADER(
     &    'MODHAM: BEFORE LINSYM final transformation'//
     &    ' matrix 4c AO->MO, VMAT(copyied from TMAT) ',-1)
          IMAT = 1
          DO I = 1,NFSYM
            IF(NORB(I).GT.0) THEN
              WRITE(LUPRI,'(/A,I1,A,I1)')
     &        '* Fermion ircop no.',I,'/',NFSYM
              CALL PRQMAT(VMAT(IMAT),NFBAS(I,0),NORB(I),
     &                    NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
              IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
            ENDIF
          ENDDO
        ENDIF

        CALL DZERO(TMAT,N2BBASXQ)
        CALL DZERO(TBUF,N2BBASXQ)
        IF (.NOT.SUB_BL) CALL INISUB
        IF (SPINFR) THEN
          SUB_BL = .TRUE.
        ENDIF

!       in:  VMAT

!       out: AO2MO linear symmetry adapted --> TMAT
!       out: AO2MO                         --> VMAT
!       out: MO2MO linear symmetry adapted --> TBUF (discarded if not x2c-mmf)
        CALL GLINSM(TMAT,VMAT,TBUF,EIG,WORK(KFREE),LFREE)

!       save matrices on file to be re-used in the x2c-module
        if(.not.spinfr)then
!         MO2MO lin
          call writt(99,ISIZE2,TBUF)
!         AO2MO lin
          call writt(99,ISIZE, TMAT)
        end if
!       AO2MO
        call writt(99,ISIZE, WORK(KBUF1))

!       by default: return the unmodified AO2MO transformation matrix
        CALL DCOPY(ISIZE,WORK(KBUF1),1,TMAT,1)

        IF(IPRHAM.GE.5) THEN
          CALL HEADER(
     &    'MODHAM: AFTER LINSYM final transformation'//
     &    ' matrix 4c AO->LINSYM MO, TMAT(output)',-1)
          IMAT = 1
          IBUF = 1
          DO I = 1,NFSYM
            IF(NORB(I).GT.0) THEN
              WRITE(LUPRI,'(A,I1,A,I1)')
     &        '* Fermion ircop no.',I,'/',NFSYM
              CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                    NFBAS(I,0),NORB(I),1,IPQTOQ(1,0),LUPRI)
              WRITE(LUPRI,'(/5X,A/)')
     &    'MODHAM: 4c ON  RKB "theta" -> LSYM "theta" transformation'//
     &        '  matrix (TBUF).'
              WRITE(LUPRI,'(A,I1,A,I1)')
     &        '* Fermion ircop no.',I,'/',NFSYM
              CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),
     &                    NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
              IMAT = IMAT + NFBAS(I,0)*NORB(I)*1
              IBUF = IBUF + NORB(I)*NORB(I)*1
            ENDIF
          ENDDO
        ENDIF
        CALL MEMREL('SLGROU1',WORK,KBUF1,KBUF1,KFREE,LFREE)
      else
!       AO2MO transformation matrix
        call writt(99,ISIZE, TMAT)
      endif
      close(99,status='keep')

      IF(SLGROUP.and.(is_defining_h1mat.gt.3.or.
     &    (is_defining_h1mat.eq.2.and.LSYMAD)))THEN
C  ...  copy the SL resorted transformation matrix to be exported into VMAT ...
        CALL DCOPY(ISIZE,WORK(KBUF),1,VMAT,1)
        IF(IPRHAM.GE.5) THEN
          CALL HEADER(
     &    'MODHAM: "SL" resorted transformation matrix'//
     &    ' 4c AO->MO, VMAT(ouput) ',-1)
          IBUF = 1
          DO I = 1, NFSYM
            IF(NORB(I).GT.0) THEN
              WRITE(LUPRI,'(A,I1,A,I1)')
     &        '* Fermion ircop no.',I,'/',NFSYM
              CALL PRQMAT(VMAT(IBUF),NFBAS(I,0),NORB(I),
     &                    NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
              IBUF = IBUF + (NFBAS(I,0)*NORB(I)*NZ)
            ENDIF
          ENDDO
        ENDIF
        CALL MEMREL('SLGROU2',WORK,KBUF,KBUF,KFREE,LFREE)
      ENDIF
 
      CALL FLSHFO(LUPRI)
      CALL QEXIT('MODHAM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prprea */
      SUBROUTINE PRPREA(LU1INT,PRPLBL,RTNLBL,PCOMB,PRPINT,NDIM,IPRINT)
C***********************************************************************
C
C     Read property integrals from file with the label PRPLBL
C     and class blocks defined by COMB 
C     Return integrals in PRPINT and return labels in RTNLBL
C
C     Called from: PRPMAO - general routine to get AO integrals from file
C                  GTOVLX - to get overlap integrals
C                  GTOVLT - to get symmetry packed overlap integrals 
C
C     Written by T.Saue, January 1995
C     Last revision: May 19 1996 - tsaue
C                    MI/March 2006 - more print-out
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      LOGICAL FNDLB3, FINDLABEL
      CHARACTER PRPLBL*8, RTNLBL(2)*8, PCOMB*4, COMB*4
      DIMENSION PRPINT(NDIM)

      CALL QENTER('PRPREA')
C
Chj   PDOINT is now '+', '-', or '0',
Chj   corresponding to a factor of +1, -1, or 0 on that block.
C     The four blocks (AO basis): LL, SL, LS, SS
      DO J = 1,4
        IF (PCOMB(J:J) .EQ. '0' .OR. PCOMB(J:J) .EQ. 'F') THEN
           COMB(J:J) = 'F'
        ELSE
           COMB(J:J) = 'T'
        END IF
      END DO
C
      REWIND LU1INT

      IF (IPRINT.GE.12) THEN
       WRITE(LUPRI,'(A,I4)')
     & 'PRPREA: >>> Going to read the property '//PRPLBL//
     & ' with COMB '//COMB//
     & ' from the file with LU1INT = ',LU1INT
      ENDIF
      IF(FNDLB3(PRPLBL,COMB,RTNLBL,LU1INT,IPRINT)) THEN
        CALL READT(LU1INT,NDIM,PRPINT)
        IF (IPRINT.GE.17) THEN
           WRITE(LUPRI,'(A,I4)')
     &    'PRPREA: >>> Property '//PRPLBL//' '//COMB//
     &    'succesfully read from the file with LU1INT = ',LU1INT
        ENDIF
      ELSE
        WRITE(LUPRI,'(A,I4)')
     &    'PRPREA: >>> Property '//PRPLBL//' '//COMB//
     &                     ' not found on LU1INT = ',LU1INT
        WRITE(LUPRI,*) 'label RTNLBL=',RTNLBL

        !MI: more printout for testing purposes
        !IF (IPRINT.GE.20) THEN
          REWIND LU1INT
          FINDLABEL=FNDLB3(PRPLBL,COMB,RTNLBL,LU1INT,32)
        !ENDIF
        CALL QUIT('PRPREA:Property '//PRPLBL//
     &   ' of COMB '//COMB//' not found on AOPROPER file !')
      ENDIF

      CALL FLSHFO(LUPRI)
      CALL QEXIT('PRPREA')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fndlb3 */
      LOGICAL FUNCTION FNDLB3(SRCLBL,COMB,RTNLBL,LU,IPRINT)
C***********************************************************************
C
C  18-Jan-1995 tsaue
C  1. Find MOLECULE label and block COMB on unformatted file LU
C       a)Cannot find label: return FALSE
C       b)Finds label      : continues

CMI: 12/12/2008: added printout for betted debugg
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      LOGICAL FNDLB2
      CHARACTER SRCLBL*8, RTNLBL(2)*8, COMB*4
      FNDLB3 = .FALSE.
   10 CONTINUE
      IF(FNDLB2(SRCLBL,RTNLBL,LU)) THEN
C     Found label SRCLBL on LU; check COMB
        IF(RTNLBL(2)(5:8).NE.COMB) THEN
CMI       ... control print out ...
          IF (IPRINT.GE.30) THEN
            WRITE(LUPRI,'(2X,A,A)') 'FNDLB3: searching for SRCLBL=',
     &      SRCLBL
            WRITE(LUPRI,'(2X,A,A4,A,A4)') 'FNDLB3: '//
     &      'RTNLBL(2)(5:8)=',RTNLBL(2)(5:8),'.ne. COMB=',COMB
          ENDIF
          GOTO 10
        ENDIF
        FNDLB3 = .TRUE.
      ENDIF
      RETURN
      END

      SUBROUTINE GTOVLX(SMAT,SSMETR)
C*****************************************************************************
C
C     Get overlap matrix from file AOPROPER
C     The SS block is scaled by the factor SSMETR
C     No symmetry packing
C
C     Written by T.Saue May 26 1996
C     Last revision: Jan 20 1998 - L. Visscher
C
C*****************************************************************************
      use memory_allocator
#ifdef MOD_INTEREST
      use module_interest_interface
#endif
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
#include "dcbham.h"
#include "dcbgen.h"
#include "cbihr1.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION SMAT(N2BBASX)
      CHARACTER*8 RTNLBL(2)
      real(8), allocatable :: sbuf(:)
      CALL QENTER('GTOVLX')
C
      IF(MDIRAC) THEN
        IF(NBSYM.GT.1) THEN
          CALL QUIT('GTOVLX: Do not know how to handle symmetry yet...')
        ENDIF
#ifdef MOD_INTEREST
        call process_overlap_integrals( smat )
#endif
      ELSE
C      
C       Allocate integral buffer
C       ========================
C
        call alloc(sbuf, NNBBASX)
C       
C       Read integral from file
C       =======================
C       
        OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &           FILE='AOPROPER')
        CALL PRPREA(LU1INT,'OVERLAP ',RTNLBL,'TFFT',sbuf,NNBBASX,
     &              IPRONE)
        CLOSE(LU1INT,STATUS = 'KEEP')
C       
C       Make full matrix
C       ================
C       
        CALL DSPTSI(NTBAS(0),SBUF,SMAT)
C       
C       Memory deallocation
C       ===================
C       
        call dealloc(sbuf)
C       
C       Reindex to sorted basis
C       =======================
C       
        CALL BUTOBS_no_work(SMAT,1)
C       
C       Scale SS-block
C       ==============
C       
        ssmetr_save = ssmetr
!       ssmetr      = 0.0d0
        IF(SSMETR.NE.D1) THEN
          DO I = 1,NFSYM
             JJ = I2BASX(I,I)+(NFBAS(I,1)*NTBAS(0))+NFBAS(I,1)+1
             DO J = 1, NFBAS(I,2)
                CALL DSCAL (NFBAS(I,2),SSMETR,SMAT(JJ),1)
                JJ = JJ + NTBAS(0)
             ENDDO
          ENDDO
        ENDIF
        ssmetr= ssmetr_save
      ENDIF
C
      CALL QEXIT('GTOVLX')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gtovlt */
      SUBROUTINE GTOVLT(SMAT,SSMETR,IPRINT)
C*****************************************************************************
C
C     Get overlap matrix from file AOPROPER
C     The SS block is scaled by the factor SSMETR
C     Packed on fermion irreps
C
C     Written by T.Saue May 26 1996
C     Last revision: Jan 20 1998 - L. Visscher
C
CMI/Sept2005 The SMAT must be allocated at least to the N2BAST size !
C
C*****************************************************************************
      use memory_allocator
#ifdef MOD_INTEREST
      use module_interest_interface
#endif
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "cbihr1.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
      DIMENSION SMAT(*)
      CHARACTER*8 RTNLBL(2)
      real(8), allocatable :: stri(:)
      real(8), allocatable :: sbuf(:)
      CALL QENTER('GTOVLT')
C      
C     Allocate integral buffer
C     ========================
C
      CALL DZERO(SMAT,N2BAST)
      NBRP = 4/NZ

      IF(MDIRAC) THEN
        IF(NBSYM.GT.1) THEN
          CALL QUIT('GTOVLT: Do not know how to handle symmetry yet...')
        ENDIF
#ifdef MOD_INTEREST
        call process_overlap_integrals( smat )
#endif
      ELSE
!       call alloc(sbuf, n2bbasx)
!       call alloc(stri, nnbbasx)
        allocate(sbuf(n2bbasx))
        allocate(stri(nnbbasx))
C
C       Read integral from file
C       =======================
C
        OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &           FILE='AOPROPER')
        CALL PRPREA(LU1INT,'OVERLAP ',RTNLBL,'TFFT',STRI,NNBBASX,
     &              IPRONE)
        CLOSE(LU1INT,STATUS = 'KEEP')
C
C       Make full matrix
C       ================
C
        CALL DSPTSI(NTBAS(0),STRI,SBUF)
C
        IOFF = 0
        DO I = 1,NFSYM
          DO IC = 1,2
            IP    = MOD(I+IC,NFSYM) + 1
            DO JSYM = 1,NBRP
              ISYM           = JFSYM(JSYM,IP)
              IREP           = ISYM - 1
!fixme irep can be < 0 in NONREL runs
              IF(NBBAS(IREP,IC).GT.0) THEN
                ISOFF = IOFF + (IBBAS(IREP,IC)-IBAS(I))*(NFBAS(I,0)+1)
                IUOFF = ICOS(ISYM,IC)*(NTBAS(0)+1)
                DO J = 1,NBBAS(IREP,IC)
                  CALL DCOPY(NBBAS(IREP,IC),SBUF(1+IUOFF),1,
     &                                      SMAT(1+ISOFF)    ,1)
                  IUOFF = IUOFF + NTBAS(0)
                  ISOFF = ISOFF + NFBAS(I,0)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
          IOFF = IOFF + N2BAS(I)
        ENDDO
        deallocate(stri)
        deallocate(sbuf)
!       call dealloc(stri)
!       call dealloc(sbuf)
      ENDIF
C
C     Scale SS-block
C     ==============
C
      IF(SSMETR.NE.D1) THEN
        IOFF = 1
        DO I = 1,NFSYM
          JJ = IOFF + NFBAS(I,1)*(NFBAS(I,0)+1)
          DO J = 1, NFBAS(I,2)
            CALL DSCAL (NFBAS(I,2),SSMETR,SMAT(JJ),1)
            JJ = JJ + NFBAS(I,0)
          ENDDO
          IOFF = IOFF + N2BAS(I)
        ENDDO
      ENDIF
C
C     Print section
C     =============
C
      IF(IPRINT.GE.5) THEN
        CALL HEADER('GTOVLT: Overlap (metric) matrix:',-1)
        IOFF = 1
        DO I = 1,NFSYM
        IF(NFBAS(I,0).GT.0) THEN
          WRITE(LUPRI,'(A,I1,A,I1)')
     &     '* Fermion ircop no.',I,'/',NFSYM   
          CALL PRQMAT(SMAT(IOFF),NFBAS(I,0),NFBAS(I,0),
     &                NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),LUPRI)
          IOFF = IOFF + N2BAS(I)
        ENDIF
        ENDDO
      ENDIF    
C
      CALL QEXIT('GTOVLT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reno */
      SUBROUTINE RENO(LDM,NDIM,SMAT)
#include "implicit.h"
#include "priunit.h"
      DIMENSION SMAT(LDM,NDIM)
      DO J = 1,NDIM
        DO I = 1,NDIM
          SMAT(I,J) = SMAT(I,J)/SQRT(SMAT(I,I)*SMAT(J,J))
        ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck sphlow */
      SUBROUTINE SPHLOW(IREP,IC,SMAT,VMAT,LDM,NEFF,STOL,
     &                  EIG,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate Lowdin's canoncial orthonormalization matrix with 
C     cartesian to spherical transformation embedded.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),
     &          WORK(*)
C
      CALL QENTER('SPHLOW')
#include "memint.h"
      NCAR = NBBAS(IREP,IC)
      NMAT = NCAR*NCAR
      CALL MEMGET('REAL',KTMAT,NMAT,WORK,KFREE,LFREE)
      CALL SPHLO1(IREP,IC,NCAR,SMAT,VMAT,LDM,NEFF,STOL,EIG,
     &            WORK(KTMAT),WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('SPHLOW',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('SPHLOW')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck sphlo1 */
      SUBROUTINE SPHLO1(IREP,IC,NCAR,SMAT,VMAT,LDM,NEFF,STOL,EIG,
     &                  TMAT,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate Lowdin's canoncial orthonormalization matrix with 
C     cartesian to spherical transformation embedded.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
      LOGICAL DONRM
      CHARACTER COMP(2)*1
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),TMAT(NCAR,NCAR),WORK(*)
      DATA COMP/'L','S'/
#include "ibtfun.h"
C
C**************************************************
C***  Generate spherical transformation matrix  ***
C**************************************************
C
      DONRM = .TRUE.
      CALL SPHCAR(TMAT,IREP,IC,NCAR,NSPH,
     &            DONRM,IPRINT,WORK,LWORK)
C
C*********************************************************
C***  Transform overlap matrix to spherical basis and
C***  generate Lowdin canonical orthonromalization matrix
C*********************************************************
C
      CALL QTRANS('AOMO','S',D0,NCAR,NCAR,NSPH,NSPH,
     &            SMAT,LDM ,NCAR,1,IPQTOQ(1,0),
     &            VMAT,LDM ,NSPH,1,IPQTOQ(1,0),
     &            TMAT,NCAR,NSPH,1,IPQTOQ(1,0),
     &            TMAT,NCAR,NSPH,1,IPQTOQ(1,0),
     &            WORK,LWORK,IPRINT)
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A,1X,A3,2X,A1)') 
     &    'SPHLOW: Overlap matrix in solid harmonics for',
     &    REP(IREP),COMP(IC)
        CALL OUTPUT(VMAT,1,NSPH,1,NSPH,LDM,NSPH,-1,LUPRI)
      ENDIF
      CALL LOWGEN(VMAT,NSPH,SMAT,NEFF,LDM,STOL,IPRINT,
     &            EIG,WORK,LWORK)
      CALL DGEMM('N','N',NCAR,NEFF,NSPH,D1,TMAT,NCAR,
     &           SMAT,LDM,D0,VMAT,LDM)
      NA = NCAR - NSPH
      NB = NSPH - NEFF 
      ND = NA + NB
      WRITE(LUPRI,'(3X,A1,3X,A3,3X,3(A,I10),A)') 
     &           COMP(IC),REP(IREP),'* Deleted: ',ND,
     &           '(Proj: ',NA,', Lindep: ',NB,')'
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck rkblow */
      SUBROUTINE RKBLOW(IREP,IC,SMAT,VMAT,LDM,NEFF,STOL,
     &                  EIG,IRKB,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate Lowdin's canoncial orthonormalization matrix with 
C     the restricted kinetic balance (RKB) transformation embedded.
C     See SUBROUTINE RKBLO1 for further explanations.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),IRKB(*),WORK(*)
C
      CALL QENTER('RKBLOW')
#include "memint.h"
      NCAR = NBBAS(IREP,IC)
      NMAT = NCAR*NCAR
      CALL MEMGET('REAL',KTMAT,NMAT,WORK,KFREE,LFREE)
      CALL RKBLO1(IREP,IC,NCAR,SMAT,VMAT,LDM,NEFF,STOL,EIG,
     &            WORK(KTMAT),IRKB,WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('RKBLOW',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('RKBLOW')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck rkblo1 */
      SUBROUTINE RKBLO1(IREP,IC,NCAR,SMAT,VMAT,LDM,NEFF,STOL,EIG,
     &                  TMAT,IRKB,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate Lowdin's canoncial orthonormalization matrix with 
C     the restricted kinetic balance (RKB) transformation embedded.
C
C     For each center and L value a transformation matrix is extracted
C     from spherical transformation matric CSP according to allowed symmetries.
C       IFUN(*,1) - points to Cartesian component in full list
C       IFUN(*,2) - points to spherical component in full list
C       IFUN(*,3) - points to Cartesian comp. in symmetry-reduced list
C     The symmetry-reduced spherical transformation matrix is
C     stored in CRED
C     For each L value
C       IRED(1,*) - gives number of symmetry-reduced Cartesians KXYZ
C       IRED(2,*) - gives number of symmetry-reduced sphericals KLM
C       IRED(3,*) - gives offsets to each L value
C     For each L value CRED has the dimension CRED(KXYZ,KLM)
C     The RKB part is adjoined to CRED and has the dimension
C        CRED(KRKB,KLMM2) where KRKB = KXYZ+KXYZM2
C     KXYZM2 and KLMM2 are the number of symmetry-reduced Cartesian
C     and spherical components, respectively, of the L-2 partner.
C
C      

C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
C
#include "dgroup.h"
#include "pgroup.h"
      LOGICAL DONRM
      CHARACTER COMP(2)*1
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),TMAT(NCAR,NCAR),
     &          IRKB(*),WORK(*)
      DATA COMP/'L','S'/
C
C**************************************************
C***  Generate spherical transformation matrix  ***
C**************************************************
C
      DONRM = .TRUE.
      CALL RKBCAR(TMAT,IREP,IC,NCAR,NSPH,
     &            DONRM,IRKB,IPRINT,WORK,LWORK)
C
C*********************************************************
C***  Transform overlap matrix to spherical basis and
C***  generate Lowdin canonical orthonormalization matrix
C*********************************************************
C
      CALL QTRANS('AOMO','S',D0,NCAR,NCAR,NSPH,NSPH,
     &            SMAT,LDM ,NCAR,1,IPQTOQ(1,0),
     &            VMAT,LDM ,NSPH,1,IPQTOQ(1,0),
     &            TMAT,NCAR,NSPH,1,IPQTOQ(1,0),
     &            TMAT,NCAR,NSPH,1,IPQTOQ(1,0),
     &            WORK,LWORK,IPRINT)
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A,1X,A3,2X,A1)') 
     &    'RKBLOW: Overlap matrix in modified solid harmonics for',
     &    REP(IREP),COMP(IC)
        CALL OUTPUT(VMAT,1,NSPH,1,NSPH,LDM,NSPH,-1,LUPRI)
      ENDIF
      CALL LOWGEN(VMAT,NSPH,SMAT,NEFF,LDM,STOL,IPRINT,
     &     EIG,WORK,LWORK)
      CALL DGEMM('N','N',NCAR,NEFF,NSPH,D1,TMAT,NCAR,
     &     SMAT,LDM,D0,VMAT,LDM)
      NA = NCAR - NSPH
      NB = NSPH - NEFF 
      ND = NA + NB
      WRITE(LUPRI,'(3X,A1,3X,A3,3X,3(A,I10),A)') 
     &     COMP(IC),REP(IREP),'* Deleted: ',ND,
     &     '(Proj: ',NA,', Lindep: ',NB,')'
C     
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lowgen */
      SUBROUTINE LOWGEN(SMAT,NS,VMAT,NV,LDM,SSTOL,IPRINT,
     &                  EIG,WORK,LWORK)
C*****************************************************************************
C
C     Generate Lowdin matrix
C     (Lowdin's canonical orthonormalization)
C
C     V_ij = U_ij/SQRT(s_j) where U_ij is the matrix that diagonalizes
C     the overlap matrix S and s_j is eigenvalue j.
C     Linear dependence is removed by eliminating columns of V
C     corresponding to eigenvalues of the overlap matrix below
C     the given threshold.
C
C     Written by T.Saue Sep 11 1995
C     Last revision Apr 1 1998 - tsaue
C                    Febr2006 - MI
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dummy.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),WORK(LWORK)
C
#include "memint.h"
      CALL QENTER('LOWGEN')

CMI .. mover before diagonalization which can destroy SMAT
      IF(IPRINT.GE.6) THEN
        CALL HEADER('LOWGEN: Lowdin input (overlap) matrix:',-1)
        WRITE(LUPRI,'(3X,A,2I5)') 'size :',NS,NS
        CALL OUTPUT(SMAT,1,NS,1,NS,LDM,NS,-1,LUPRI)
      ENDIF

#ifdef MOD_MATLAB_LOG
      call matexport_text('% Entering LOWGEN')
      call matexport_double2('S',SMAT,LDM,NS,NS)
#endif
C
C     Diagonalize overlap matrix
C
#ifdef MOD_UNRELEASED
      IF (.NOT.DOQJACO) THEN
#endif
         IF (LOWJACO) THEN
            CALL RSJACO(LDM,NS,NS,SMAT,EIG,1,-1,0,VMAT,WORK,LWORK)
         ELSE
            CALL QDIAG(1,NS,SMAT,LDM,NS,EIG,1,VMAT,LDM,NS,WORK,LWORK
     $           ,IERR)
c    QDIAG gives the reverse value ordering from RSJACO, so fix that here.
c    It would be better to rewrite the code below instead.. 
            CALL ORDER3(VMAT,EIG,LDM,NS,NS,-1)
            IF (IERR.NE.0) THEN
               CALL QUIT('Diagonalization of overlap matrix failed')
            ENDIF
         ENDIF
#ifdef MOD_UNRELEASED
      ELSE
CMI   ... (quaternion) Jacobi diagonalization routine
        NSS=NS*NS
        CALL MEMGET('REAL',KSMAT,NSS,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KVMAT,NSS,WORK,KFREE,LFREE)
        IPOS=0
        DO II=1,NS
        DO JJ=1,NS
          WORK(KSMAT+IPOS)=SMAT(JJ,II) ! needed NxN matrix for input
          IPOS=IPOS+1
        ENDDO
        ENDDO
        CALL QJACOBI(WORK(KSMAT),WORK(KVMAT),NS,1,0,IDUMMY,.TRUE.,2)
        IPOS=0
        DO II=1,NS
          EIG(II)=WORK(KSMAT+(NS*(II-1)+II-1)) ! extract eigenvalues
        DO JJ=1,NS
          VMAT(JJ,II)=WORK(KVMAT+IPOS) ! extract eigenvectors
          IPOS=IPOS+1
        ENDDO
        ENDDO
        CALL ORDER3(VMAT,EIG,LDM,NS,NS,-1) ! resort
        CALL MEMREL('LOWGEN',WORK,KSMAT,KSMAT,KFREE,LFREE)
      ENDIF
#endif
!    ... ensure V_ij = U_ij/SQRT(s_j)
      NV = 0
      DO I = 1,NS
        IF(EIG(I).LE.SSTOL) GOTO 10
        NV  = NV + 1
        FAC = D1/SQRT(EIG(I))
        CALL DSCAL(NS,FAC,VMAT(1,I),1)
      ENDDO
 10   CONTINUE
C
C     Print section
C
      IF(IPRINT.GE.2) THEN
          WRITE(LUPRI,'(7X,A,1P,D8.1)')
     &   'LOWGEN: Smallest ".LINDEP" test value of a kept orbital:',
     &   EIG(NV)
      ENDIF

C
      IF(IPRINT.GE.3) THEN
        CALL HEADER('LOWGEN: Eigenvalues in Lowdin diagonalization:',-1)
        CALL OUTPUT(EIG,1,NS,1,1,NS,1,-1,LUPRI)
      ENDIF
C
      IF(IPRINT.GE.6) THEN
        CALL HEADER(
     &  'LOWGEN: Canonical Lowdin matrix (i.e. eigenvectors):',-1)
        CALL OUTPUT(VMAT,1,NS,1,NV,LDM,NS,-1,LUPRI)
      ENDIF

#ifdef MOD_MATLAB_LOG
      call matexport_double2('V',VMAT,LDM,NS,NV)
      call matexport_text('% Leaving LOWGEN')
#endif

      CALL QEXIT('LOWGEN')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tscorr */
      SUBROUTINE TSCORR(FMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Generate SS nuclear attraction integrals to mimic
C     two-electron SS-integrals
C
C     Written by T.Saue July 1997
C     Last revision July 15 1997 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(D1 = 1.0D0)
C
      DIMENSION WORK(LWORK),CBUF(2),GBUF(2),FMAT(*)
      LOGICAL DOINT(2,2)
      CHARACTER OMITVNUC(2)*4
#include "dgroup.h"
#include "dcbbas.h"
#include "nuclei.h"
C
      CALL QENTER('TSCORR')
#include "memint.h"
C     Memory allocation
      CALL MEMGET('REAL',KONEMT,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN,NNBBASX,WORK,KFREE,LFREE)
C
C     Buffer info
C
      DO I = 1,NUCIND
        CBUF(I) = CHARGE(I)
        GBUF(I) = GNUEXP(I)
        GNUEXP(I) = 0.0D0
      ENDDO
      call quit('TSCORR not properly implemented in this version')
C     HJAaJ: this is only for a specific datomic!!! (TeSr ??)
C
      CHARGE(1) = -0.1812D0
      CHARGE(2) = -0.0832D0
C
C     Generate integrals
C
      DOINT(1,1) = .FALSE.
      DOINT(2,2) = .TRUE.
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.

      OMITVNUC(1) = 'TTTT'
      OMITVNUC(2) = 'FFFF'
      CALL NUCAT1(WORK(KONEIN),WORK(KFREE),LFREE,
     &            DOINT,OMITVNUC,IPRINT)
C
      CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
      CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
      IF(IPRINT.GE.6) THEN
        CALL HEADER('SS nuclear attraction integrals',-1)
          CALL PRQMAT(WORK(KONEMT),NTBAS(0),NTBAS(0),
     &                NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
      ENDIF
      CALL DAXPY(N2BBASX,D1,WORK(KONEMT),1,FMAT,1)
C
C     Memory deallocation
      CALL MEMREL('TSCORR',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Buffer retrieve
C
      DO I = 1,NUCIND
        CHARGE(I) = CBUF(I) 
        GNUEXP(I) = GBUF(I) 
      ENDDO
C
      CALL QEXIT('TSCORR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bncorr */
      SUBROUTINE BNCORR(FMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     BNCORR: Generate screened bare nucleus matrix
C     Jun 2003 Jesper Kielberg Pedersen (based on OLDBNCR)
C
C     Description : 
C    ===============
C     Extension of BNCORR that in a simple way takes the shell-structure
C     into account. For each center in a molecule we estimate the screening
C     neglected in the bare nucleus approximation by a sum of the 
C     contributions from each shells (n) in the atoms :
C           sum_(n) <X_A|(-Z(n)*G(n,a,r_C))/r_C|X_B> 
C           X=L,S ; G(a,r)=exp(-a*r^2)
C
C     On input/ouput:     FMAT - assumed one-electron (2c/4c bare nucleus) Fock matrix (input),
C    =================    screening bare nucleus added (ouput)
C
C     Called from: PREDHF/dirscf.F
C                  RH1DIAG/../krmc/krmcopt.F
C
C      In the case of two-component mode do the picture change transformation
C      of the four-compo. correction matrix. In pure two-component mode add only
C      LL block.
C
C    Last modifications: March 2006/MI,Strasbourg 
C
C*****************************************************************************
      use x2c_fio, only:
     &    x2c_read
      use x2c_utils, only:
     &    print_x2cmat
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0)
C
      LOGICAL DOINT(2,2), FINISH, ISONBSSMAT, WEHAVE 
      CHARACTER*4 OMITVNUC
      CHARACTER*8 BNCRLAB
      DIMENSION WORK(LWORK),FMAT(*)
      logical fndlab12
      character*12 flabel
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbiher.h"
#include "cbihr1.h"
#include "cbisol.h"
C
      CALL QENTER('BNCORR')
c
c     Start timer
c
      TIMSTR = SECOND() 
c
#include "memint.h"

      WEHAVE = .FALSE.

      if(bss)then
        WEHAVE = ISONBSSMAT('BNCR2C_2','BNCORR',LUBSS,IPRHAM).OR.
     &           ISONBSSMAT('BNCR2C_4','BNCORR',LUBSS,IPRHAM)

      else if(x2c)then

       open(lux2c,file='X2CMAT',status='old',form='unformatted',
     &      access='sequential',action='read',position='rewind')

       write(flabel,'(a11,i1)') 'h1bncc2cAOl',0
       if(fndlab12(flabel,lux2c)) then 
         wehave = .true.
       else
         close(lux2c,status='keep')
       end if
      end if

      IF (WEHAVE) THEN
CMI   ... read 2c screened bare nucleus from the BSSMAT file in 4c or 2c metric
        if(bss)then
          if(twocomp)then
            BNCRLAB='BNCR2C_2'
          else
            BNCRLAB='BNCR2C_4'
          end if
          CALL RFBSSMAT(BNCRLAB,.FALSE.,'BNCORR',
     &                  FMAT,N2BBASXQ,LUBSS,IPRHAM)
!         call print_x2cmat(FMAT,NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
!    &         'bss - h1barenuc',LUPRI)
        else if(x2c)then
          print *, 'blubb enter' 
          call print_x2cmat(FMAT,NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &         'x2c - enter fmat',LUPRI)
          call x2c_read(flabel,fmat,n2bbasx * nz,lux2c)
          print *, 'blubb leaving' 
          call print_x2cmat(FMAT,NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &         'x2c - h1barenuc',LUPRI)
          close(lux2c,status='keep')
        end if
      ELSE
CMI   ... continue the standard way
C     Memory allocation
      CALL MEMGET('REAL',KONEMT, N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN, NNBBASX,WORK,KFREE,LFREE)
      LFCBA = NUCDEP**3
      CALL MEMGET('REAL',KFCBA ,  LFCBA ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXGEXP,  NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCHRG ,  NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KESHEL,7*NUCIND,WORK,KFREE,LFREE)
!     allocate(eshell_vec(7*NUCIND))
      LEN_WRK   = LWORK - KFREE + 1
      CALL DZERO(WORK(KCHRG),NUCIND)
      CALL DZERO(WORK(KONEMT),N2BBASX)
      CALL DZERO(WORK(KESHEL),7*NUCIND)
C
C     Initialization (calculate LL and SS)
C
      DOINT(1,1) = .TRUE.
      IF (NOSMLV.OR.BSS.OR.x2c) THEN
C     IF (NOSMLV) THEN
         DOINT(2,2) = .FALSE.
         IF (IPRINT.GE.1) WRITE(LUPRI,'(/A)')
     &     '  BNCORR: SS contribution to the bare nuclei omitted !'
      ELSE
         DOINT(2,2) = .TRUE.
      END IF
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.
c
c   ... find max number of occupied shells
c
      Z = D0
      DO I = 1,NUCIND
         Z = MAX(Z,CHARGE(I))
      ENDDO

      CALL BNCOR1(GNU,ZEFF,0,Z,NSHELL,0,DUMMY,NUCIND,.TRUE.)

c   ... loop over shells
c       in the first shell (N=1) we only place one electron to avoid
c       self-interaction. This electron is delocalized (to some degree)
c       in the zeroth run (N=0) by using a diffuse Gaussian (see BNCOR1)
c
      IF (IPRINT .GE. 1) THEN
         WRITE(LUPRI,'(/4X,A/4X,A//14X,A/14X,A)')
     &' *** BNCORR for < X_A | (-Zeff*G(a,r_C))/r_C | X_B > ***',
     &'     X=L,S ; G(a,r_C)=exp(-a*r_C^2)',
     &      ' Nucleus    Charge=(-Zeff)   Exponent(a)',
     &      '-----------------------------------------'
      ENDIF
      DO N = 0,NSHELL
c     ... In Slater's rules we write the electron configuration as
c         [1s] [2s2p] [3s3p][3d] [4s4p][4d][4f]... so for each shell we have
c          1     1        2            3           loops
         IF (N.LE.1) THEN
            NSLOOP = 1
         ELSE
            NSLOOP = N - 1
         ENDIF
c        ... we will not go beyond f-orbitals
         NSLOOP = MIN(NSLOOP,3)
         DO L = 1,NSLOOP
c
         ZTEST = D0
         DO I = 1,NUCIND
c      ... Get charges and exponents for the electrons in this group
           CALL BNCOR1(GNU,ZEFF,N,CHARGE(I),L,I,WORK(KESHEL),
     &          NUCIND,.FALSE.)
           ZTEST = MAX(ZTEST,ZEFF)
           WORK(KXGEXP-1+I) =  GNU
           WORK(KCHRG -1+I) = -ZEFF
         ENDDO
c
         IF(IPRINT .GE. 1) THEN
            IF (N .LE. 1) THEN
               WRITE(LUPRI,'(14X,A,I2,A)')
     &         '-- Subshell n =',N,', l = 0'
            ELSE IF (L .EQ. 1) THEN
               WRITE(LUPRI,'(14X,A,I2,A)')
     &         '-- Subshell n =',N,', l = 0,1'
            ELSE
               WRITE(LUPRI,'(14X,A,I2,A,I2)')
     &         '-- Subshell n =',N,', l =',L
            END IF
            DO I = 1,NUCIND
               IF (WORK(KCHRG-1+I) .NE. D0)
     &         WRITE(LUPRI,'(15X,A5,F15.2,F15.2)')
     &         NAMN(I),WORK(KCHRG-1+I),WORK(KXGEXP-1+I)
            ENDDO
         ENDIF
c
         IF (ZTEST .NE. D0) THEN
C
C     Omit selected terms depending on OMITVNUC (See CATNUC for details)
C     In this version one-center terms are calculated with the same
C     charges/exponents as multi-center terms.
C
         OMITVNUC = 'FFFF'
c
c     ... Make list of nuclei.
c
         CALL CATNUC(NUCDEP,NCENTC,WORK(KCHRG),WORK(KFCBA),WORK(KGEXP),
     &            WORK(KXGEXP),WORK(KCOORC), WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),OMITVNUC,IPRONE)
c
c    ... Calculate integrals (both LL and SS)
c
         CALL CATDR1(WORK(KONEMT),DUMMY,NCENTC,WORK(KFCBA),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KFREE),LEN_WRK,IPRONE,
     &            .FALSE.,0,NNBBASX,DOINT)
c
         END IF ! if (ztest .ne. d0)
c  ... end loop within shells
         ENDDO ! DO L = 1,NSLOOP
c  ... end loop over shells
      ENDDO ! DO N = 0,NSHELL
      IF (IPRINT .GE. 1) WRITE(LUPRI,'(14X,A/)')
     &      '-----------------------------------------'
C
C     Symmetry unpack integrals
C
      CALL SYMUPK(WORK(KONEMT),WORK(KONEIN),1,NNBBASX)
      CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
      CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
C
      IF(IPRINT.GE.11) THEN
         CALL HEADER('BNCORR correction for DC FMAT',-1)
         CALL OUTPUT(WORK(KONEMT),1,NTBAS(0),1,NTBAS(0),
     &        NTBAS(0),NTBAS(0),-1,LUPRI)
      END IF

C
C     Add correction to (entering) Fock-matrix, FMAT
C
      CALL DAXPY(N2BBASX,D1,WORK(KONEMT),1,FMAT,1)
      IF(IPRINT.GE.11) THEN
         CALL HEADER('BNCORR: screening corrected bare nucleus FMAT',-1)
         CALL PRQMAT(FMAT,NTBAS(0),NTBAS(0),
     &               NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
      END IF
C
C     Memory deallocation
      CALL MEMREL('BNCORR',WORK,KWORK,KWORK,KFREE,LFREE)

      ENDIF    

C
c     Stop timer
c
      TIMEND = SECOND()
      TIME   = TIMEND - TIMSTR
      IF (IPRINT.GE.2)
     &  CALL TIMTXT('>>> Time used in BNCORR is',TIME,LUPRI)
      CALL QEXIT('BNCORR')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bncor1 */
      SUBROUTINE BNCOR1(GNU,ZEFF,N,Z,L,INUC,SHELE,NNUC,MAXSHELL)
C*****************************************************************************
C
C     BNCOR1 :  Find effective charges and Gaussian exponent 
C               (from Slater rules) for a given shell (N) for
C               a given atom (Z)
C
C     Jun 2003 Jesper Kielberg Pedersen (used for BNCORR)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0,
     &          D6 = 6.0D0, D8 = 8.0D0, D10=10.0D0, D14=14.0D0,
     &          D18=18.0D0)
      LOGICAL MAXSHELL
      DIMENSION STO1G(7), XN_EFF(7), SHELE(NNUC,7)
      SAVE STO1G
C     STO1G(n) is the exponent of the 1s-Gaussian that
C     best fits the ns Slater-type atomic orbitals of exponent 1.
C     (Robert F. Stewart, JCP 52(1),431 (Jan 1970))
C     This article provides exponent up to the n=5 shell. The exponents
C     for the last two shells have been generated from a fit :
C     (y = 0.28196 x^-1.5532)
      DATA STO1G /0.2709D0, 0.1012D0, 0.05297D0,0.03265D0,0.02217D0,
     &            0.01744D0,0.01373D0/
!     DATA XN_EFF/1.0D0, 2.0D0, 3.0D0, 3.7D0, 4.0D0, 4.2D0, 4.4D0/
      DATA XN_EFF/1.0D0, 2.0D0, 3.0D0, 3.7D0, 4.4D0, 5.1D0, 5.8D0/
C     ... Slater's n_eff according to wikipedia for n=1:6,
C         I extrapolated to n=7 /hjaaj Sep 2011.
C     ... Empirical modification Feb. 2012 hjaaj: XeF6 did not
C         converge with the Sep 2011 values; on the other hand
C         we had an example in Sep 2011 which benefitted from the
C         n_eff values.
C
C     Table of s,p,d,f orbital occupations of the elements
C     (from which the electron-configuration can be deduced)
      DIMENSION NOCC(118*4)
      SAVE NOCC
      DATA NOCC 
C     Periode 1
     &          / 1, 0, 0, 0, 2, 0, 0, 0,
C     Periode 2
     &            3, 0, 0, 0, 4, 0, 0, 0,
     &            4, 1, 0, 0, 4, 2, 0, 0, 4, 3, 0, 0,
     &            4, 4, 0, 0, 4, 5, 0, 0, 4, 6, 0, 0,
C     Periode 3
     &            5, 6, 0, 0, 6, 6, 0, 0,
     &            6, 7, 0, 0, 6, 8, 0, 0, 6, 9, 0, 0,
     &            6,10, 0, 0, 6,11, 0, 0, 6,12, 0, 0,
C     Periode 4
     &            7,12, 0, 0, 8,12, 0, 0,
     &            8,12, 1, 0, 8,12, 2, 0, 8,12, 3, 0,
     &            8,12, 4, 0, 8,12, 5, 0, 8,12, 6, 0,
     &            8,12, 7, 0, 8,12, 8, 0, 8,12, 9, 0,
     &            8,12,10, 0,
     &            8,13,10, 0, 8,14,10, 0, 8,15,10, 0,
     &            8,16,10, 0, 8,17,10, 0, 8,18,10, 0,
C     Periode 5 (5s block)
     &            9,18,10, 0, 10,18,10, 0,
C     (4d block)
     &           10,18,11, 0, 10,18,12, 0, 10,18,13, 0,
     &           10,18,14, 0, 10,18,15, 0, 10,18,16, 0,
     &           10,18,17, 0, 10,18,18, 0, 10,18,19, 0,
     &           10,18,20, 0,
C     (5p block)
     &           10,19,20, 0, 10,20,20, 0, 10,21,20, 0,
     &           10,22,20, 0, 10,23,20, 0, 10,24,20, 0,
C     Periode 6 (6s block)
     &           11,24,20, 0, 12,24,20, 0, 12,24,21, 0,
C    --> Lanthanides (4f block)
     &           12,24,21, 1, 12,24,20, 3, 12,24,20, 4,
     &           12,24,20, 5, 12,24,20, 6, 12,24,20, 7,
     &           12,24,20, 8, 12,24,20, 9, 12,24,20,10,
     &           12,24,20,11, 12,24,20,12, 12,24,20,13,
     &           12,24,20,14, 12,24,21,14,
C    <-- Lanthanides; Hf -> Hg (5d block)
     &           12,24,22,14, 12,24,23,14, 12,24,24,14,
     &           12,24,25,14, 12,24,26,14, 12,24,27,14,
     &           12,24,28,14, 12,24,29,14, 12,24,30,14,
C     Tl -> Rn (6p block)
     &           12,25,30,14, 12,26,30,14, 12,27,30,14,
     &           12,28,30,14, 12,29,30,14, 12,30,30,14,
C     Periode 7 (7s block)
     &           13,30,30,14,
     &           14,30,30,14,
     &           14,30,31,14,
C    --> Actinides
     &           14,30,32,14, 14,30,31,16, 14,30,31,17,
     &           14,30,31,18, 14,30,30,20, 14,30,30,21,
     &           14,30,31,21, 14,30,30,23, 14,30,30,24,
     &           14,30,30,25, 14,30,30,26, 14,30,30,27,
     &           14,30,30,28, 14,30,31,28,
C    <-- Actinides; (6d block and 7p block)
     &           14,30,32,28, 14,30,33,28, 14,30,34,28,
     &           14,30,35,28, 14,30,36,28, 14,30,37,28,
     &           13,30,39,28, 13,30,40,28, 14,30,40,28,
     &           14,31,40,28, 14,32,40,28, 14,33,40,28,
     &           14,34,40,28, 14,35,40,28, 14,36,40,28/
      CALL QENTER('BNCOR1')
c
      IF (NINT(Z) .GT. 118) THEN
         WRITE(LUPRI,'(1X,2(A/))') 'FATAL : The bare nucleus '//
     &   'correction is only implemented for elemtents 1-118',
     &   '         You can disable BNCORR with the .NOBNCR keyword '//
     &   'in *SCF'
         CALL QUIT('BNCORR : Unsupported element!')
      ENDIF

      ZEFF = D0
      GNU  = D0

      IF (MAXSHELL) THEN
c    ... just find maximum number of shells we need to run over 
c        and return number in L-variable (used in BNCORR)
         NST = NOCC((NINT(Z)-1)*4 + 1)
         L = NST/2
         IF(MOD(NST,2).NE.0) L = L + 1
         GOTO 10
      ENDIF
C
c     ... If this is a ghost-atom (f.ex for Counter-Poise 
c         correction) we return (charge is zero in BNCORR)
      IF (Z .EQ. D0) GOTO 10
      IF (N.EQ.0) THEN
c     ...Always make 1 electron semi-delocalized
         ZEFF = D1
         GNU  = 0.4D0
         SHELE(INUC,1) = SHELE(INUC,1) + D1
         GOTO 10
      ENDIF
c     ... If Hydrogen we are done
      IF (Z.EQ.D1) GOTO 10
c     ... special case for the (2-1) 1s electron
      IF (N.EQ.1) THEN
         ZEFF = D1
         SEXP = Z - 0.35
         IF (Z.EQ.2) SEXP = SEXP + 0.05
         SHELE(INUC,1) = SHELE(INUC,1) + D1
         GOTO 20
      ENDIF
c     ... Generel case. We are past H and He and have thus accounted 
c         for two electrons.
c     ... Find numbers of s,p,d,f-electrons for element (Z)
         NST=0
         NPT=0
         NDT=0
         NFT=0
      IF (Z.GE.1.0D0) THEN ! miro: prevent NOCC of out-of-bounds 
         NST = NOCC((NINT(Z)-1)*4 + 1)
         NPT = NOCC((NINT(Z)-1)*4 + 2)
         NDT = NOCC((NINT(Z)-1)*4 + 3)
         NFT = NOCC((NINT(Z)-1)*4 + 4)
      ENDIF
C
c     ... Now figure out the occupation used for the Slater rules (N >= 2) :
c      [1s][2s2p][3s3p][3d][4s4p][4d][4f]...
c
c     =================================
c     How many s-electrons in this shell
c     =================================
c
      NS = NST - (N-1)*2
c
      IF (NS .LE. 0) THEN
         NS = 0
      ELSEIF (NS.GE.2) THEN
         NS = 2
      ENDIF
c
c     =================================
c     How many p-electrons in this shell
c     =================================
c
      NP = NPT - (N-2)*6
c
      IF (NP .LE. 0) THEN
         NP = 0
      ELSEIF (NP.GE.6) THEN
         NP = 6
      ENDIF
c
c     =================================
c     How many d-electrons in this shell
c     =================================
c
      ND = NDT - (N-3)*10
c
      IF (ND .LE. 0) THEN
         ND = 0
      ELSEIF (ND.GE.10) THEN
         ND = 10
      ENDIF
c
c     =================================
c     How many f-electron in this shell
c     =================================
c
      NF = NFT - (N-4)*14
c
      IF (NF .LE. 0) THEN
         NF = 0
      ELSEIF (NF.GE.14) THEN
         NF = 14
      ENDIF
      IF (L.EQ.1) ZEFF = NS + NP
      IF (L.EQ.2) ZEFF = ND
      IF (L.EQ.3) ZEFF = NF
      SHELE(INUC,N) = SHELE(INUC,N) + ZEFF
      IF (ZEFF.EQ.D0) GO TO 10
C
C     ... Finally use Slater's rules to estimate screening 
C         (and orbital exponent)
C
      IF (L.EQ.1) THEN
c        [nsnp]
         EMN = D0
         DO I = 1,N-2
            EMN = EMN + SHELE(INUC,I)
         ENDDO
         SEXP = (Z - ((ZEFF-1.0)*0.35D0 + SHELE(INUC,N-1)*0.85D0 + EMN))
         SEXP = SEXP / XN_EFF(N)
      ELSE
c        [nd],[nf]
         EMN = D0
         DO I = 1,N-1
            EMN = EMN + SHELE(INUC,I)
         ENDDO
         SEXP = (Z - ((ZEFF-1.0D0)*0.35D0 + EMN)) / XN_EFF(N)
      ENDIF
c
   20 CONTINUE
c     ... Convert STO to GTO using 1 Gaussian (STO1G)
      GNU = D2*SEXP*SEXP*STO1G(N)
!     write (lupri,*) 'BNCOR1: N,L,STO1G(N),SEXP,GNU',
!    &                         N,L,STO1G(N),SEXP,GNU
C     ... Make sure we dont get more delocalized than the first 
C         1s electron.
C     GNU = MAX(GNU,0.50D0)
   10 CONTINUE
      CALL QEXIT ('BNCOR1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck eonecap */
      SUBROUTINE EONECAP(FMAT,DMAT,WORK,LWORK,IQ)
C****************************************************************************
C     PURPOSE:
C        Calculate correction to the nuclear attraction energy 
C        from the one-center models
C
C      E = D(IJ)(*)F(IJ)
C
C     Sep 2001 : Jesper Kielberg Pedersen
C     Based on ERGCAL
C     Revised Juli 2002 : jkp
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D2=2.0D0)
C
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
      CHARACTER FMT*6,MXFORM*6
      DIMENSION FMAT(*), DMAT(*),WORK(LWORK)
C
      CALL QENTER('EONECAP')
C
      ECORR  = D0
C
C        Factor D2 comes from Time-Reversal Symmetry.
C
      IF(INTV1C.EQ.3.OR.(INTV1C.EQ.2.AND.ICTLV1C(2).EQ.0)) THEN
         IF (NOPEN .GE. 1) THEN
            CALL QUIT('EONECAP not implemented for open shell '//
     &                'for this ONECAP model yet.')
         END IF
      END IF
      ECORR = D2*DDOT(N2BBASX,DMAT,1,FMAT,1)
C
      CALL HEADER('Output from EONECAP',-1)
      FMT = MXFORM(ABS(ECORR),20)
      IF (IQ .EQ. 1) THEN
         WRITE(LUPRI,'(A,'//FMT//')')
     &  'ONECAP correction to LL-Nuc att.  :  ',ECORR
      ELSE IF (IQ .EQ. 2) THEN
         WRITE(LUPRI,'(A,'//FMT//')')
     &  'ONECAP correction to SS-Nuc att.  :  ',ECORR
      ELSE IF (IQ .EQ. 3) THEN
         WRITE(LUPRI,'(A,'//FMT//')')
     &  'Estimation of error on the energy in this model :  ',ECORR
      END IF
C
      CALL QEXIT('EONECAP')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiprj */
      SUBROUTINE PSIPRJ(TMAT,TINV,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Eliminate a set of fragment orbitals from the variational space;
C     this is done by projecting them out of the MO transformation matrix
C     Written by T.Saue March 26 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "dcbgen.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dummy.h"
C
      LOGICAL TOBE
      DIMENSION TMAT(*),TINV(*),WORK(*)
      DIMENSION NSTR(2,0:2,MFRAG),KVEC(2,MFRAG),NVECS(2)
      DIMENSION IOFF(2)
C
      KFRSAV = KFREE
      CALL TITLER('Reduction of variational space by projection',
     &     '*',125)
C
C     Set dimensions for coefficients
C     ===============================
C
      NORBT = 0
      NCMOT  = 0
      DO IFRP = 1,NFSYM
        NORBT       = NORBT  + NORB(IFRP)
        ICMO(IFRP)  = NCMOT
        ICMOQ(IFRP) = NCMOT*NZ
        NCMO(IFRP)  = NFBAS(IFRP,0)*NORB(IFRP)
        NCMOQ(IFRP) = NCMO(IFRP)*NZ
        NCMOT       = NCMOT + NCMO(IFRP)
      ENDDO
      NCMOTQ  = NCMOT*NZ
C
C     Find total number of fragment orbitals
C     =======================================
C
C     NFRAG - number of fragments
C     NSTR(IFRP,IC,IFRAG) 
C     - number of fragment orbitals of fermion ircop IFRP for fragment IFRAG
C        IC = 0 total
C        IC = 1 positive energy
C        IC = 2 negative energy
C     NVECS(IFRP) - total number of fragment orbitals in fermion ircop IFRP
C     NVECT       - total number of fragment orbtials
C        
      NVECT = 0
C.....loop over fermion ircops
      DO IFRP = 1,NFSYM
        NVECS(IFRP) = 0
C.......loop over fragments
        DO IFRAG = 1,NFRAG
          NSTR(IFRP,1,IFRAG) = NESH(IFRP)
          NSTR(IFRP,2,IFRAG) = NPSH(IFRP)
          CALL ORBNUM(VCPROJ(IFRP,IFRAG),IFRP,KVEC(IFRP,IFRAG),
     &                NSTR(1,0,IFRAG),WORK,KFREE,LFREE)
          NVECS(IFRP) = NVECS(IFRP) + NSTR(IFRP,0,IFRAG)
        ENDDO
        NVECT = NVECT + NVECS(IFRP)
      ENDDO
      IF(NVECT.EQ.0) GOTO 30
C
C     Calculate dimensions of fragment coefficients array
C     ====================================================
C
C     KQ - start adress for coefficient array
C     KE - start adress for eigenvalue array
C     KI - start adress for boson irreps info
C
      NCDIM = 0
      NEDIM = 0
      DO IFRP = 1,NFSYM
        NCDIM = NCDIM + NFBAS(IFRP,0)*NVECS(IFRP)*NZ
        NEDIM = NEDIM + NVECS(IFRP)
      ENDDO       
      CALL MEMGET('REAL',KQ,NCDIM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KE,NEDIM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KI,NEDIM,WORK,KFREE,LFREE)
      CALL IZERO(WORK(KI),NEDIM)
C
C     Allocate buffer space for coefficients, eigenvalues
C     and boson irrep information
C     ===================================================
C
      CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT   ,WORK,KFREE,LFREE)
      IF(SUB_BL) THEN
        IOPT = 15
        CALL MEMGET('REAL',KIBE,NORBT,WORK,KFREE,LFREE)
      ELSE
        IOPT = 7
        KIBE = KFREE
      ENDIF
C
C     Select fragment orbitals
C     =========================
C
      DO IFRP = 1,NFSYM
        IOFF(IFRP) = 0
      ENDDO
      CALL SELFRAG(WORK(KQ),NVECS,WORK(KE),WORK(KI),
     &             NFRAG,NSTR,PRJFIL,PROOWN,NPRJNUC,IOPT,KVEC,
     &             WORK(KCMO),WORK(KEIG),WORK(KIBE),IOFF,
     &             WORK,KFREE,LFREE)
C
C     Throw away the full set of coefficients and eigenvalues
C     =======================================================
C
      CALL MEMREL('PSIPRJ.cfs',WORK,KFRSAV,KCMO,KFREE,LFREE)
C
C     Reduce variational space by projection using the fragment orbitals
C     ===================================================================
C     If more than one set of fragment orbitals, they must be orthogonalized !
C
C     Do projection
C
      CALL PSIPR1(NVECS,WORK(KQ),WORK(KE),WORK(KI),TMAT,TINV,
     &            NSTR,IPRINT,WORK,KFREE,LFREE)
 30   CONTINUE
      CALL MEMREL('PSIPRJ.final',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PRJORT */
      SUBROUTINE PRJORT(VECS,IBVEC,VREF,IBREF,IREF,NREF,
     &                  THRLDP,THRMAX,NVDIM,NVECS,NVREF,NDEL,IP,IPRINT)
C*****************************************************************************
C
C     On input:
C        VECS  - unit matrix (transformation matrix in MO-basis)
C        IBVEC - integer array with boson irrep info for transformation matrix
C        VREF  - fragment vectors in MO-basis
C        IBREF - integer array with boson irrep info for fragment vectors
C        NREF  - number of fragment vectors to project out
C        IREF  - starting adress of fragment vectors to project out
C        
C     Written by T. Saue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,DM1=-1.0D0)
      PARAMETER (T1MIN = 1.0D-8, THRRND = 1.0D-4)
C
#include "dgroup.h"
      DIMENSION VECS(NVDIM,NVECS,NZ),VREF(NVDIM,NVREF,NZ),
     &          IBVEC(NVECS),IBREF(NVREF)
      DIMENSION P(4)
C
CTROND
      WRITE(6,*) 'TEST, IPRINT = ', IPRINT
      IF(IPRINT.GE.7) THEN
         CALL HEADER('PRJORT: MO transformation matrix on input',-1)
         CALL PRQMAT(VECS,NVDIM,NVECS,NVDIM,NVECS,1,IPQTOQ(1,0),LUPRI)
      ENDIF
      THRMAX = D0
      NDEL   = 0
      N2VECS = NVDIM*NVECS
      N2VREF = NVDIM*NVREF
      IROUND = 0
      ITURN  = 0
 10   CONTINUE
      ITURN  = ITURN + 1
      IF(IPRINT.GE.2) WRITE(6,'(A,I5)') '* Projection turn :',ITURN
C
C     First project out fragment vectors
C     ===================================
C
      JTOP = IREF+NREF-1
      DO JREF = IREF,JTOP
        IF(IPRINT.GE.3) THEN
          CALL QDOT(VNRM,IPQTOQ(1,0),1,NVDIM,
     &           'H',IPQTOQ(1,0),VREF(1,JREF,1),N2VREF,NZ,1,
     &           'N',IPQTOQ(1,0),VREF(1,JREF,1),N2VREF,NZ,1)
          VNRM=SQRT(VNRM)
          WRITE(LUPRI,'(A,I3,A,1P,D12.5,4X,A,I5)')
     &      'Fragment vector no.',JREF,' Norm: ',VNRM,
     &      'Bsym: ',IBREF(JREF)
        ENDIF
        DO 100 IVEC = 1,NVECS
          IF(IBVEC(IVEC).NE.IBREF(JREF)) GOTO 100
C
C         Form -<VREF|IVEC>
C
          CALL QDOT(P,IPQTOQ(1,0),NZ,NVDIM,
     &         'H',IPQTOQ(1,0),VREF(1,JREF,1),N2VREF,NZ,1,
     &         'N',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1)
          DO IZ = 1,NZ
            P(IZ) = -P(IZ)
          ENDDO
C
C         Form VECS = VECS -VREF<VREF|VECS>
C
          CALL  QXAPY(NVDIM,'N','N',P,IPQTOQ(1,0),NZ,
     &                VREF(1,JREF,1),N2VREF,1,IPQTOQ(1,0),NZ,
     &                VECS(1,IVEC,1),N2VECS,1,IPQTOQ(1,0),NZ)
 100    CONTINUE
      ENDDO
C
C     Next commence internal orthogonalization
C     =========================================
C
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,'(A)')
     &      '* MO-vectors after projecting out reference vectors:'
      ENDIF
      DO 200 IVEC = 1,NVECS
CTROND: This is a bit dangerous; should be fixed !
        IF(IBVEC(IVEC).LT.-100) GOTO 200
C
C       Check size of norm
C       ==================
C
        CALL QDOT(VNRM,IPQTOQ(1,0),1,NVDIM,
     &           'H',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1,
     &           'N',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1)
        VNRM=SQRT(VNRM)

        THRMAX=MAX(THRMAX,VNRM)
        IF(VNRM.LT.THRLDP) THEN
          IF(IPRINT.GE.1) THEN
            WRITE(LUPRI,'(3X,A,I3,A,1P,D12.5,3X,A,I5)')
     &        'Vector no.',IVEC,' removed. Linear dependence: ',VNRM,
     &        ' Bsym: ',IBVEC(IVEC)
          ENDIF
C         Store information about irrep and neg/pos solution
C         of deleted vector
          IBVEC(IVEC) = -(IP*1000+IBVEC(IVEC))
          NDEL = NDEL + 1
          GOTO 200
        ELSE
          IF(IPRINT.GE.2) WRITE(LUPRI,'(3X,A,I3,A,1P,D12.5,3X,A,I5)')
     &      'Vector no.',IVEC,' Norm: ',VNRM,
     &      'Bsym:',IBVEC(IVEC)
          IF(VNRM.LT.THRRND.AND.ITURN.EQ.1) IROUND=IROUND+1
          IF(VNRM.LT.T1MIN) THEN
            VNRM = D1/VNRM
            DO IZ = 1,NZ
              CALL DSCAL(NVDIM,VNRM,VECS(1,IVEC,IZ),1)
            ENDDO
            CALL QDOT(VNRM,IPQTOQ(1,0),1,NVDIM,
     &           'H',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1,
     &           'N',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1)
          END IF
          VNRM = D1 / VNRM
          DO IZ = 1,NZ
            CALL DSCAL(NVDIM,VNRM,VECS(1,IVEC,IZ),1)
          ENDDO
        ENDIF
        DO 300 JVEC = (IVEC+1),NVECS
          IF(IBVEC(JVEC).NE.IBVEC(IVEC)) GOTO 300
C
C         Form -<IVEC|JVEC>
C
          CALL QDOT(P,IPQTOQ(1,0),NZ,NVDIM,
     &         'H',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1,
     &         'N',IPQTOQ(1,0),VECS(1,JVEC,1),N2VECS,NZ,1)
          DO IZ = 1,NZ
            P(IZ) = -P(IZ)
          ENDDO
C
C         Form JVEC = JVEC -IVEC<IVEC|JVEC>
C
          CALL  QXAPY(NVDIM,'N','N',P,IPQTOQ(1,0),NZ,
     &                VECS(1,IVEC,1),N2VECS,1,IPQTOQ(1,0),NZ,
     &                VECS(1,JVEC,1),N2VECS,1,IPQTOQ(1,0),NZ)
 300    CONTINUE
 200  CONTINUE
      IF (ITURN.LT.4) THEN
        GO TO 10
      ENDIF
      IF(IPRINT.GE.7) THEN
         CALL HEADER('PRJORT: MO transformation matrix on output',-1)
         CALL PRQMAT(VECS,NVDIM,NVECS,NVDIM,NVECS,1,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck prjcnt */
      SUBROUTINE PRJCNT(VMO,VAO,IFRP,NMO,NAO,NEFF,IVEC,TMAT)
C***********************************************************************
C
C     After a projection step, backtransform to SO-basis
C     VMO is Mo transformation matrix after transformation.
C     The columns for which IVEC=1 are backtransformed to
C     VAO, the transformation matrix in AO-basis, using
C     the transformation matrix TMAT
C
C     Written by T.Saue Jan 31 2001
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
      DIMENSION VMO(NMO,NMO,NZ),VAO(NAO,NEFF,NZ),IVEC(*),
     &          TMAT(NAO,NMO,NZT)
      DIMENSION IBOS(0:7)
      NBRP = 4/NZ
      DO ISYM = 1, NBRP
        IREP = JFSYM(ISYM,IFRP) - 1
        IBOS(IREP) = ISYM
      ENDDO
      II = 0
      DO I = 1,NMO
      IF(IVEC(I).GE.-100) THEN
        II = II + 1
        CALL QGEMM(NAO,1,NMO,D1,
     &          'N','N',IPQTOQ(1,0),TMAT,NAO,NMO,NZT,
     &          'N','N',IPQTOQ(1,0),VMO(1,I,1), NMO,NMO,NZ,
     &          D0,IPQTOQ(1,0),VAO(1,II,1),NAO,NEFF,NZ)
      ELSE
        IF(SUB_BL) THEN
          ID = -MOD(IVEC(I),1000)
          IP   = -(IVEC(I)+ID)/1000
          DO IB = 1, N_SUB_BL(IFRP)
            IF(ID.eq.ID_SUB_BL(IB,IFRP)) THEN
              NORB_SUB(IB,IFRP,IP) = NORB_SUB(IB,IFRP,IP) - 1
              GOTO 10
            ENDIF
          ENDDO
 10       CONTINUE
        ENDIF
      ENDIF
      ENDDO
      DO IB = 1, N_SUB_BL(IFRP)
        NORB_SUB(IB,IFRP,0) = NORB_SUB(IB,IFRP,1)+NORB_SUB(IB,IFRP,2)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck prjcnt2 */
      SUBROUTINE PRJCNT2(IVEC,IFRP,NBAS,TMAT,NMAT,TBUF,NBUF)
C***********************************************************************
C
C     After a projection step, update dimension arrays in 
C     COMMON blocks DCBBAS and DCBORB
C     Written by T.Saue Jan 31 2001
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION IVEC(*),TMAT(NBAS,NMAT,NZ),TBUF(NBAS,NBUF,NZ)
      NSYM = 4/NZ
      JVEC = 0
      JBUF = 0
      DO ISYM = 1,NSYM
        DO IC = 2,1,-1
          NVEC = NBORB(ISYM,IFRP,IC)
          IF(NVEC.GT.0) THEN
            NBORB(ISYM,IFRP,IC)=ISUM(NVEC,IVEC(JVEC+1),1)
            DO II = 1,NVEC
              JVEC = JVEC + 1
              IF(IVEC(JVEC).EQ.1) THEN
                JBUF = JBUF + 1
                DO IZ = 1,NZ
                  CALL DCOPY(NBAS,TMAT(1,JVEC,IZ),1,TBUF(1,JBUF,IZ),1)
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        NBORB(ISYM,IFRP,0) = NBORB(ISYM,IFRP,1) + NBORB(ISYM,IFRP,2)
      ENDDO
      NESH(IFRP) = ISUM(NSYM,NBORB(1,IFRP,1),1)
      NPSH(IFRP) = ISUM(NSYM,NBORB(1,IFRP,2),1)
      NORB(IFRP) = NESH(IFRP) + NPSH(IFRP)
      RETURN
      END
C
C*****************************************************************************
C  /* Deck mk1cex */
      SUBROUTINE MK1CEX (WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C     PURPOSE:
C        Construct one-center expansion coefficients
C        
C        Written by L.Visscher
C        Last revision: Aug 16 2001 - lv
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
#include "ccom.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "nuclei.h"
      DIMENSION WORK(*)
C
      CALL QENTER('MK1CEX')
      KSAVE = KFREE
C
      LDM = NTBAS(0)
      LDM2 = LDM * LDM
      CALL MEMGET('REAL',KSMAT1,LDM2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT2,LDM2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDAO,LDM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDAOS,LDM,WORK,KFREE,LFREE)
      NIND = 2 * NUCDEP + 2
      CALL MEMGET('INTE',KINDNAO,NIND,WORK,KFREE,LFREE)
C
C     Get overlap matrix
C     ==================
C
      CALL GTOVLX(WORK(KSMAT1),D1)
C
C     Backtransform overlap matrix to AO-basis
C     ========================================
C
      CALL MTBSBU(WORK(KSMAT1),WORK(KSMAT2))
      CALL MTSOAO(WORK(KSMAT2),WORK(KSMAT1),NTBAS(0),0,IPRINT)
C
C     Get index arrays to sort on expansion centers (nuclei)
C     ======================================================
C
      CALL INDAOC (WORK(KINDAO),WORK(KINDAOS),WORK(KINDNAO))
C
C     Sort the overlap matrix to this order
C     =====================================
C
      CALL REORMT (WORK(KSMAT1),WORK(KSMAT2),LDM,LDM,
     &             WORK(KINDAO),WORK(KINDAO))
C
C     Construct the projection operator and put it in SMAT1
C     =====================================================
C
      CALL MK1CEX2 (LDM,WORK(KSMAT2),WORK(KSMAT1),WORK(KINDNAO),
     &              WORK,KFREE,LFREE,IPRINT)
C
C     Sort the projection operator back to the standard AO order
C     ========================================================
C
      CALL REORMT (WORK(KSMAT1),WORK(KSMAT2),LDM,LDM,
     &             WORK(KINDAOS),WORK(KINDAOS))
C
C     Write the projection operator to file for later use
C     ===================================================
C
      CALL OPNFIL(LUPMAT,'DFPMAT','UNKNOWN','MK1CEX')
      CALL WRITT(LUPMAT,N2BBASX,WORK(KSMAT2))
      CLOSE(LUPMAT,STATUS='KEEP')
C
      CALL MEMREL('MK1CEX',WORK,KSAVE,KSAVE,KFREE,LFREE)
C
      CALL QEXIT('MK1CEX')
      RETURN
      END
C*****************************************************************************
C  /* Deck mk1cex2 */
      SUBROUTINE MK1CEX2 (LDM,STOT,PTOT,INDNAO,WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C     PURPOSE:
C        Construct one-center expansion coefficients
C        
C        Written by L.Visscher
C        Last revision: Aug 16 2001 - lv
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
#include "ccom.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "nuclei.h"
      DIMENSION STOT(LDM,LDM),PTOT(LDM,LDM)
      DIMENSION WORK(*)
      DIMENSION INDNAO(2,*)
C
      CALL QENTER('MK1CEX2')
      KSAVE = KFREE
      CALL DZERO (PTOT,LDM*LDM)
C
C     Loop over large/small component
C     Only the small part is considered because the
C     one-center approximation is not used in the 
C     large component.
C 
      DO IC = 2, 2
         NB0 = INDNAO(IC,1)+1
C
C        Loop over expansion centers
C
         DO IEXPC = 1, NUCDEP
           NB1 = INDNAO(IC,IEXPC)+1
           NB2 = INDNAO(IC,IEXPC+1)
           NLOC = NB2-NB1+1
           NV = NLOC*LDM
           NW = NLOC*NLOC
           CALL MEMGET('REAL',KEIG,NLOC,WORK,KFREE,LFREE)
           CALL MEMGET('REAL',KV,NV,WORK,KFREE,LFREE)
           CALL MEMGET('REAL',KW,NW,WORK,KFREE,LFREE)
C
C          Make the local projection operator
C
           CALL MK1CEX3 (STOT(NB1,NB0),STOT(NB1,NB1),PTOT(NB1,NB0),LDM,
     &                   NTBAS(IC),NLOC,
     &                   WORK(KV),WORK(KW),WORK(KEIG),STOL(IC),
     &                   IPRINT,WORK,KFREE,LFREE)
           CALL MEMREL('MK1CEX2',WORK,KEIG,KEIG,KFREE,LFREE)
C
C          This operator contains also the projection on the functions of the same
C          center which should be (close to) the unit matrix. Make this part zero.
C
C          TODO : check that this part is indeed the unit matrix.
C
           DO I = NB1, NB2
              CALL DZERO(PTOT(NB1,I),NLOC)
           ENDDO
         ENDDO
      ENDDO
C
      IF(IPRINT.GE.5) THEN
        CALL HEADER('Projection operator in AOC basis',-1)
        CALL OUTPUT(PTOT,1,LDM,1,LDM,LDM,LDM,-1,LUPRI)
      ENDIF
C
      CALL MEMREL('MK1CEX2',WORK,KSAVE,KSAVE,KFREE,LFREE)
      CALL QEXIT('MK1CEX2')
      RETURN
      END
C*****************************************************************************
C  /* Deck mk1cex3 */
      SUBROUTINE MK1CEX3 (STOT,SLOC,PTOT,LDM,NTOT,NLOC,V,W,EIG,STOL,
     &                    IPRINT,WORK,KFREE,LFREE)
C*****************************************************************************
C     PURPOSE:
C        Make projection operator for this center. 
C        We reduce the basis if the overlap matrix on this
C        center is nearly singular. This may happen as we work in
C        the untransformed non-orthogonal AO basis. STOL should be
C        consistent with the values used in LOWDIN !!
C        
C        On input :
C           STOT : The complete overlap matrix
C           SLOC : The overlap matrix for the current center
C           LDM : Leading dimension (should be NBAST)
C           NTOT : The total number of functions to be projected 
C                  (now all large or small functions in the molecule
C                   but generalization to domains is easy)
C           NLOC : The number of functions on the current center
C           STOL : Linear dependency treshold
C
C        On output :
C           PTOT : The updated projection operator matrix (added the
C                  rows for the current center)
C
C        NB : STOT and PTOT should be positioned such that the current
C             center is the first row !
C 
C        Scratch arrays : V, W, EIG
C 
C        Written by L.Visscher
C        Last revision: Aug 16 2001 - lv
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
      DIMENSION WORK(*)
      DIMENSION STOT(LDM,*),PTOT(LDM,*)
      DIMENSION SLOC(LDM,NLOC),V(LDM,NLOC),W(NLOC,NLOC),EIG(NLOC)
C
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Local overlap matrix',-1)
        CALL OUTPUT(SLOC,1,NLOC,1,NLOC,LDM,NLOC,-1,LUPRI)
      ENDIF
C
C     Form (Sloc)**-1/2 * U
C
      CALL LOWGEN(SLOC,NLOC,V,NV,LDM,STOL,IPRINT,EIG,WORK(KFREE),LFREE)
C
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Local Lowdin canonical matrix',-1)
        CALL OUTPUT(V,1,NLOC,1,NV,LDM,NLOC,-1,LUPRI)
      ENDIF
C
C     Form (Sloc)**-1/2 * U * U+ * (Sloc)**-1/2 = (Sloc)**-1
C
      CALL DGEMM('N','T',NLOC,NLOC,NV,D1,V,LDM,V,LDM,D0,W,NLOC)
C
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Local inverse overlap matrix',-1)
        CALL OUTPUT(W,1,NLOC,1,NLOC,NLOC,NLOC,-1,LUPRI)
      ENDIF
C
C     Make (Sloc)**-1  * Stot : the projection operator  
C
      CALL DGEMM('N','N',NLOC,NTOT,NLOC,D1,W,NLOC,STOT,LDM,D0,PTOT,LDM)
C
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Local projection operator',-1)
        CALL OUTPUT(PTOT,1,NLOC,1,NTOT,LDM,NLOC,-1,LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PR1CEX1 */
      SUBROUTINE PR1CEX1 (FTRI,FIRST,FOCK,WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C     Correct for neglect of multicenter integrals in symmetry packed SO matrix
C
C     Written by L. Visscher, august 2001
C     Last revision: Aug 16, 2001
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbham.h"
C
      LOGICAL FIRST(NZ)
      DIMENSION FOCK(N2BBASX,NZ),FTRI(NNBBASX,NZ),
     &          WORK(*)
C
C     Make full matrix
C     ================
C
      DO IZ = 1,NZ
        IF(FIRST(IZ)) THEN
          CALL DZERO(FOCK(1,IZ),N2BBASX)
        ELSE
          IQ = IPQTOQ(IZ,0)
          IH = IHQMAT(IQ,1)
          IF(IH.EQ.1) THEN        
            CALL DSPTSI(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DAPTGE(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ENDIF
        ENDIF
      ENDDO
C
C     Allocate work space to store the projection operator and
C     a work array.
C
      CALL MEMGET('REAL',KPROJ,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF1,N2BBASX*4,WORK,KFREE,LFREE)
C
C     Get the projection operator
C
      CALL OPNFIL(LUPMAT,'DFPMAT','OLD','PR1CEX1')
      CALL READT(LUPMAT,N2BBASX,WORK(KPROJ))
      CLOSE(LUPMAT,STATUS='KEEP')
C
C     Add contributions from missing multicenter integrals
C     ====================================================
C
      IREP = 0
      CALL PR1CEX2 (NTBAS(0),IREP,FOCK,WORK(KPROJ),WORK(KBUF1),IPRINT)
C
C     Release the memory
C
      CALL MEMREL('PR1CEX1',WORK,KPROJ,KPROJ,KFREE,LFREE)
C
C     Pack matrix
C     ===========
C
      DO IZ = 1,NZ
        IF(.NOT.FIRST(IZ)) THEN
          IQ = IPQTOQ(IZ,0)
          IH = IHQMAT(IQ,1)
          IF(IH.EQ.1) THEN        
            CALL DSITSP(NTBAS(0),FOCK(1,IZ),FTRI(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DGETAP(NTBAS(0),FOCK(1,IZ),FTRI(1,IZ))
          ENDIF
        ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PR1CEX2 */
      SUBROUTINE PR1CEX2 (NB,IREP,FOCK,PROJ,BUF1,IPRINT)
C***********************************************************************
C     Correct for neglect of multicenter integrals in SO matrix
C
C     A basis function |q> can be projected on the set of functions
C     |r> on center A using the projection operators P(r,q) 
C     made in MK1CEX1. 
C
C     We use this to correct the input matrix F(p,q) that contains
C     only the one-center contributions (F1) because the multicenter
C     contributions (F2) were neglected.
C 
C     F(p,q) = F1(p,q) + F2(p,q)
C
C     FR(p,q) = F1(p,r) * P(r,q)
C     FL(p,q) = P+(p,r) * F1(r,q) 
C
C     F(p,q) is approximated by F1(p,q) + 1/2 (FL(p,q) + FR(p,q)) 
C     We use the average of the left and right projected matrices to
C     preserve the permutation symmetry.
C
C     Written by L. Visscher, august 2001
C     Last revision: Aug 16, 2001
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
C
      PARAMETER (D0=0.D0, D1=1.D0, HALF=0.5D0)
      DIMENSION FOCK(NB,NB,NZ),PROJ(NB,NB),BUF1(NB,NB,4)
C
C     Transform Fock matrix to AO-basis
C     ==================================
C
      IF (IPRINT.GE.5) THEN
         CALL HEADER('PR1CEX2: Fock matrix in SO basis',-1)
         CALL PRQMAT(FOCK,NB,NB,NB,NB,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      IPAR = JBTOF(IREP,1)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL MTSOAO(FOCK(1,1,IPQ),BUF1(1,1,IZ),NB,IREPD,IPRINT)
      ENDDO
C
      DO IZ = 1, 4
         IF (IPRINT.GE.8) THEN
            CALL HEADER('PR1CEX2: Fock matrix in AO basis',-1)
            CALL OUTPUT(BUF1(1,1,IZ),1,NB,1,NB,NB,NB,-1,LUPRI)
         ENDIF
C
C        Construct the right projected fock matrix and store it in FOCK
C
         CALL DGEMM('N','N',NB,NB,NB,D1,BUF1(1,1,IZ),NB,PROJ,NB,
     &               D0,FOCK,NB)
C
C        Construct and add the left projected fock matrix
C
         CALL DGEMM('T','N',NB,NB,NB,D1,PROJ,NB,BUF1(1,1,IZ),NB,
     &               D1,FOCK,NB)
C
C        Add the correction to the Fock matrix
C 
         IF (IPRINT.GE.10) THEN
            CALL HEADER('PR1CEX2 : Projection operator',-1)
            CALL OUTPUT(PROJ,1,NB,1,NB,NB,NB,-1,LUPRI)
            CALL HEADER('PR1CEX2: Correction to Fock matrix (*2)',-1)
            CALL OUTPUT(FOCK,1,NB,1,NB,NB,NB,-1,LUPRI)
         ENDIF
C
         NB2 = NB * NB
         CALL DAXPY (NB2,HALF,FOCK,1,BUF1(1,1,IZ),1)
C
         IF (IPRINT.GE.8) THEN
            CALL HEADER('PR1CEX2: Corrected Fock matrix',-1)
            CALL OUTPUT(BUF1(1,1,IZ),1,NB,1,NB,NB,NB,-1,LUPRI)
         ENDIF
      ENDDO
C     
C     Transform Fock matrix to SO-basis
C     =================================
C
      CALL DZERO(FOCK,NB2*NZ)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL MTAOSO(BUF1(1,1,IZ),FOCK(1,1,IPQ),NB,IREPD,IPRINT)
      ENDDO
C
      IF (IPRINT.GE.5) THEN
         CALL HEADER('PR1CEX2: Corrected Fock matrix in SO basis',-1)
         CALL PRQMAT(FOCK,NB,NB,NB,NB,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PR1CEX3 */
      SUBROUTINE PR1CEX3 (NB,IREP,DENS,PROJ,BUF1,IPRINT)
C***********************************************************************
C     Correct for neglect of multicenter integrals in SO matrix
C     Note that this routines resembles PR1CEX2 that is used for
C     Fock matrix (contributions). This one is good for densities
C     where the transpose of P needs to be used and the transformation
C     to AO basis requires a different routine.
C
C     A basis function |q> can be projected on the set of functions
C     |r> on center A using the projection operators P(r,q) 
C     made in MK1CEX1. 
C
C     We use this to correct the input matrix D(q,p) to map multicenter
C     contributions to the calculated one-center integrals
C 
C     D(q,p) = D1(q,p) + D2(q,p)
C
C     DR(q,p) = D1(q,r) * P+(r,p)
C     DL(q,p) = P(q,r) * D1(r,p) 
C
C     D(q,p) is approximated by F1(q,p) + 1/2 (DL(q,p) + DR(q,p)) 
C     We use the average of the left and right projected matrices to
C     preserve the permutation symmetry.
C
C     Written by L. Visscher, august 2001
C     Last revision: Aug 17, 2001
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
C
      PARAMETER (D0=0.D0, D1=1.D0, HALF=0.5D0)
      DIMENSION DENS(NB,NB,NZ),PROJ(NB,NB),BUF1(NB,NB,4)
C
C     Transform density matrix to AO-basis
C     =====================================
C
      IF (IPRINT.GE.5) THEN
         CALL HEADER('PR1CEX3: Density matrix in SO basis',-1)
         CALL PRQMAT(DENS,NB,NB,NB,NB,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      IPAR = JBTOF(IREP,1)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL DTSOAO(DENS(1,1,IPQ),BUF1(1,1,IZ),NB,IREPD,IPRINT)
      ENDDO
C
      DO IZ = 1, 4
         IF (IPRINT.GE.8) THEN
            CALL HEADER('PR1CEX3: Density matrix in AO basis',-1)
            CALL OUTPUT(BUF1(1,1,IZ),1,NB,1,NB,NB,NB,-1,LUPRI)
         ENDIF
C
C        Construct the right projected fock matrix and store it in DENS
C
         CALL DGEMM('N','T',NB,NB,NB,D1,BUF1(1,1,IZ),NB,PROJ,NB,
     &               D0,DENS,NB)
C
C        Construct and add the left projected fock matrix
C
         CALL DGEMM('N','N',NB,NB,NB,D1,PROJ,NB,BUF1(1,1,IZ),NB,
     &               D1,DENS,NB)
C
C        Add the correction to the density matrix
C 
         IF (IPRINT.GE.10) THEN
            CALL HEADER('PR1CEX3 : Projection operator',-1)
            CALL OUTPUT(PROJ,1,NB,1,NB,NB,NB,-1,LUPRI)
            CALL HEADER('PR1CEX3: Correction to density matrix (*2)',-1)
            CALL OUTPUT(DENS,1,NB,1,NB,NB,NB,-1,LUPRI)
         ENDIF
C
         NB2 = NB * NB
         CALL DAXPY (NB2,HALF,DENS,1,BUF1(1,1,IZ),1)
C
         IF (IPRINT.GE.8) THEN
            CALL HEADER('PR1CEX3: Corrected density matrix',-1)
            CALL OUTPUT(BUF1(1,1,IZ),1,NB,1,NB,NB,NB,-1,LUPRI)
         ENDIF
      ENDDO
C     
C     Transform density matrix to SO-basis
C     ====================================
C
      CALL DZERO(DENS,NB2*NZ)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL DTAOSO(BUF1(1,1,IZ),DENS(1,1,IPQ),NB,IREPD,IPRINT)
      ENDDO
C
      IF (IPRINT.GE.5) THEN
         CALL HEADER('PR1CEX3: Corrected density matrix in SO basis',-1)
         CALL PRQMAT(DENS,NB,NB,NB,NB,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck glinsm */
      SUBROUTINE GLINSM(TMAT,VMAT,TBUF,EIG,WORK,LWORK)
C*****************************************************************************
C
C     This routine makes an additional transformation from the
C     matrix supplied by GMOTRA to a basis that is diagonal in
C     one component of angular momentum (j in general, l for spinfree
C     and Levy-Leblond). This can be used for systems that posess
C     linear symmetry.
C
C     Called from GMOTRA
C
C     Written by L. Visscher october/november 2003
C     Modified by T. Saue Jan 12 2007
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "cbihr1.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION WORK(LWORK),TMAT(*),VMAT(*), TBUF(*), EIG(*)
      CALL QENTER('GLINSM')

#include "memint.h"
      IF (URKBAL) THEN
         WRITE (LUPRI,1000) "unrestricted kinetic balance"
         LINEAR = .FALSE.
         CALL QEXIT('GLINSM')
         RETURN
      ENDIF
C
      IF (SPINFR) THEN
         WRITE (LUPRI,1000) "spinfree"
         LINEAR = .FALSE.
         CALL QEXIT('GLINSM')
         RETURN
      ENDIF
C
      IF (ZORA) THEN
         WRITE (LUPRI,1000) "ZORA"
         LINEAR = .FALSE.
         CALL QEXIT('GLINSM')
         RETURN
      ENDIF
C
      CALL LINSYM(TMAT,VMAT,TBUF,EIG,WORK(KFREE),LFREE)
      SUB_BL = .TRUE.

      IF (IPRHAM.GE.2) THEN
       WRITE(LUPRI,'(/2X,A/)')
     &'GLINSM: TMAT and Saomo were adapted for the linear symmetry.'
      ENDIF
      CALL FLSHFO(LUPRI)

      CALL QEXIT('GLINSM')
      RETURN
C
 1000 FORMAT (//" WARNING : linear symmetry is not implemented for ",
     &A," yet" /"           continue calculation in lower symmetry")
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck linsym */
      SUBROUTINE LINSYM(TMAT,VMAT,TBUF,EIG,WORK,LWORK)
C*****************************************************************************
C
C     Construct matrix Jz(x,y,z) and transform to eigenvectors of this operator
C
C     Routine based on MODHAM and parts of the ZORA diagonalization routine
C
C     Called from GLINSM, LOWD2C
C
C     On input:  VMAT 4c AO-MO RKB transformation matrix
C
C     On output: TMAT - new 4c AO-MO RKB transf.matrix adapted for the linear symmetry
C                TBUF - MO to MO transformation matrix adapted for linear symmetry...
C
C     Written by L.Visscher, october 2003
C     Last changes: M.Ilias, Strasbourg, 2005 
C                  M.Ilias, Tel Aviv, April 2009 - extended with QJACOBI
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0, D1=1.D0, D2=2.D0)
C
C Used from COMMON blocks
C  CBIPAM:
C  DCBBAS:
C  DCBORB:
C  DGROUP:
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "cbihr1.h"
#include "dummy.h"
C
      DIMENSION TMAT(*),TBUF(*),VMAT(*),EIG(*),WORK(LWORK)
      DIMENSION IQBF(4),IQW(4)

      CALL QENTER('LINSYM')

#include "memint.h"
C
!     iprham_save = iprham
!     iprham      = 15
      CALL HEADER('Output from LINSYM',-1)
      NBRP = 4 / NZ
C
C     Construct mj matrix in AO basis
C
      CALL AMJMAT(TBUF,IPRHAM,WORK,LWORK)
C
C     Transform to MO-basis
C     =====================
C
C Note : It's tempting to use I2ORBT to find the start of a block, but this is not
C possible as this array is not yet set (is done in dirset). Took me half an hour to
C figure this one out (LV).
C
      IOFT = 1
      IOFV = 1
      DO 10 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 10
        IOFTB = I2BASX(I,I) + 1
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &            NORB(I),NORB(I),
     &            TBUF(IOFTB),NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &            TMAT(IOFT),NORB(I),NORB(I),1,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(I,0),NORB(I),1,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(I,0),NORB(I),1,IPQTOQ(1,0),
     &            WORK,LWORK,IPRHAM)
!    &            IPRHAM)

        IF (IPRHAM.GE.15) THEN

         CALL HEADER(
     &    'LINSYM: Entering AO->MO transformation matrix',-1) 
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IOFV),NFBAS(I,0),NORB(I),
     &   NFBAS(I,0),NORB(I),1,IPQTOQ(1,0),LUPRI)

         CALL HEADER('LINSYM: MJ matrix in AO basis',-1) 
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IOFTB),NFBAS(I,0),NFBAS(I,0),
     &   NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)

         CALL HEADER('LINSYM: MJ matrix in MO basis',-1) 
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IOFT),NORB(I),NORB(I),
     &   NORB(I),NORB(I),1,IPQTOQ(1,0),LUPRI)

         CALL FLSHFO(LUPRI)

        ENDIF

        IOFT = IOFT + NORB(I)*NORB(I)
        IOFV = IOFV + NORB(I)*NFBAS(I,0)
   10 CONTINUE

C
C     Initialize TBUF to zero (we will fill in parts later and its easier to initialize
C     the rest here). Also initialize the array that counts the number of functions for 
C     each Mj value.
C
      CALL DZERO (TBUF,N2TMT)

      CALL IZERO (NORB_SUB,MAX_SUB_BL*3*2)
      WRITE (LUPRI,1000)

!     stefan - june 2011: save LL block of <j_z>-eigenvectors (temporary) on file
!                         this block is needed to transform the 2c-LOWDIN matrix 
!                         in the x2c-mmf approach (thus ensuring numerical stability by not re-doing this step).
      if(do4c2c)then
        open(13,file='JzMOLLb',status='replace',form='unformatted',
     &       access='sequential',action='readwrite',
     &       position='rewind')
      end if
C
      IOFT = 1
      DO 20 I = 1,NFSYM
C
        IF(NORB(I).EQ.0) GOTO 20
C         
C       Set dimensions
C
        NL = NESH(I) 
        NS = NPSH(I)
        NT = NORB(I)
        NLL = NL * NL
        NSS = NS * NS
        NTT = NT * NT
        IEIGP = IORB(I) + 1
        IEIGE = IEIGP + NS
C
C       Allocate memory for the diagonalization. We want to
C       have the LL and SS blocks in separate blocks.
C
        CALL MEMGET('REAL',KFLL,NLL,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFSS,NSS,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KVLL,NLL,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KVSS,NSS,WORK,KFREE,LFREE)
C
C       Extract LL block of Jz matrix
C
        CALL EXTRSB (1,I,NBORB,4,NBRP,NBORB,4,NBRP,0,0,1,1,1,NZ,
     &               WORK(KFLL),NL,NL,1,
     &               TMAT(IOFT),NT,NT,NZ)
C
C       Find eigenvalues of Jz
C
        IF (NL.GT.0) THEN
          IF (IPRHAM.GE.15) THEN
           CALL HEADER('LINSYM: Jz_LL to be diagonalized',-1) 
           WRITE(LUPRI,'(2X,A,L1)') 'DOJACO=',DOJACO
           WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(WORK(KFLL),NL,NL,NL,NL,1,
     &          IPQTOQ(1,0),LUPRI)
          ENDIF

#ifdef MOD_UNRELEASED
          IF (.NOT.DOQJACO) THEN
#endif
            IF(DOJACO) THEN
              CALL RSJACO(NL,NL,NL,WORK(KFLL),EIG(IEIGE),1,1,0,
     &                  WORK(KVLL),WORK(KFREE),LFREE)
            ELSE
              CALL QDIAG(1,NL,WORK(KFLL),NL,NL,
     &                EIG(IEIGE),1,
     &                WORK(KVLL),NL,NL,
     &                WORK(KFREE),LFREE,IERR)
              IF (IERR.NE.0) THEN
                WRITE(LUPRI,'(2X,A)')
     &         'LINSYM: Errorneous QDIAG diagonalization of Jz_LL!'
               CALL QUIT(
     &        'LINSYM: Errorneous QDIAG diagonalization of Jz_LL !')
              ENDIF
            ENDIF
#ifdef MOD_UNRELEASED
          ELSE
            CALL QJACOBI(WORK(KFLL),WORK(KVLL),NL,1,0,IDUMMY,.TRUE.,2)
            DO II=1,NL
              EIG(IEIGE+II-1)=WORK(KFLL+NL*(II-1)+II-1) ! store eigenvalues
            ENDDO
          ENDIF
#endif
          IF (IPRHAM.GE.15) THEN
           CALL HEADER('LINSYM: Jz_LL eigenvalues / eigenvectors',-1) 
           WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,'/',NFSYM
           write(lupri,*) (EIG(IEIGE+II),II=0,NL-1)
           CALL PRQMAT(WORK(KVLL),NL,NL,NL,NL,1,
     &          IPQTOQ(1,0),LUPRI)
           CALL FLSHFO(LUPRI)
           !call quit('in linsym')
          ENDIF
        ENDIF
C
C       Extract SS block of Jz matrix
C
        CALL EXTRSB (1,I,NBORB,4,NBRP,NBORB,4,NBRP,0,0,2,2,1,NZ,
     &               WORK(KFSS),NS,NS,1,
     &               TMAT(IOFT),NT,NT,NZ)
C
        IF (NS.GT.0) THEN
          IF (IPRHAM.GE.15) THEN
           CALL HEADER('LINSYM: Jz_SS to be diagonalized',-1) 
           WRITE(LUPRI,'(2X,A,L1)') 'DOJACO=',DOJACO
           WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(WORK(KFSS),NS,NS,NS,NS,1,
     &          IPQTOQ(1,0),LUPRI)
          ENDIF
#ifdef MOD_UNRELEASED
          IF (.NOT.DOQJACO) THEN
#endif
            IF(DOJACO) THEN
              CALL RSJACO(NS,NS,NS,WORK(KFSS),EIG(IEIGP),1,1,0,
     &                  WORK(KVSS),WORK(KFREE),LFREE)
            ELSE
              CALL QDIAG(1,NS,WORK(KFSS),NS,NS,
     &                EIG(IEIGP),1,
     &                WORK(KVSS),NS,NS,
     &                WORK(KFREE),LFREE,IERR)
              IF (IERR.NE.0) THEN
                WRITE(LUPRI,'(2X,A)')
     &         'LINSYM: Errorneous QDIAG diagonalization of Jz_SS!'
               CALL QUIT(
     &        'LINSYM: Errorneous QDIAG diagonalization of Jz_SS !')
              ENDIF
            ENDIF
#ifdef MOD_UNRELEASED
          ELSE
            CALL QJACOBI(WORK(KFSS),WORK(KVSS),NS,1,0,IDUMMY,.TRUE.,2)
            DO II=1,NS
              EIG(IEIGP+II-1)=WORK(KFSS+NS*(II-1)+II-1) ! store eigenvalues
            ENDDO
          ENDIF
#endif
          IF (IPRHAM.GE.15) THEN
           CALL HEADER('LINSYM: Jz_SS eigenvalues/eigenvectors',-1) 
           WRITE(LUPRI,'(3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,'/',NFSYM
           write(lupri,*) (EIG(IEIGP+II),II=0,NS-1)
           CALL PRQMAT(WORK(KVSS),NS,NS,NS,NS,1,
     &          IPQTOQ(1,0),LUPRI)
           CALL FLSHFO(LUPRI)
          ENDIF
        ENDIF
C
C       Count the Mj-values.
C       We will either find a negative or a positive value for each Mj.
C
        MJMAX = 0
        DO IX = IEIGE, IEIGE+NL-1
           MJ = NINT(D2*EIG(IX))
           MJABS = (ABS(MJ)+1)/2
           IF (MJABS.GT.MAX_SUB_BL) THEN
             write(lupri,*) 'MAX_SUB_BL:',MAX_SUB_BL,
     &       ' actual MJ:',MJABS,' max value of MJ:',
     &       (ABS( NINT(D2*EIG(IEIGE+NL-1)))+1)/2
             CALL QUIT ("increase max_sub_bl")
           ENDIF
!        miro: for spin-free ensure the MJABS index not to be out-of-bounds
           IF (MJABS.GT.0) THEN
             ID_SUB_BL(MJABS,I)  = MJ
             NORB_SUB(MJABS,I,1) = NORB_SUB(MJABS,I,1) + 1
             NORB_SUB(MJABS,I,0) = NORB_SUB(MJABS,I,0) + 1
             MJMAX = MAX(MJABS,MJMAX)
           ENDIF
        ENDDO
C
        DO IX = IEIGP, IEIGP+NS-1
           MJ = NINT(D2*EIG(IX))
           MJABS = (ABS(MJ)+1)/2
           IF (MJABS.GT.MAX_SUB_BL) THEN 
             write(lupri,*) 'MAX_SUB_BL:',MAX_SUB_BL,
     &        ' max value of MJ:',
     &       (ABS( NINT(D2*EIG(IEIGE+NL-1)))+1)/2
             CALL QUIT ("increase max_sub_bl")
           ENDIF
!        miro: for spin-free avoid MJABS index out-of-bounds
           IF (MJABS.GT.0) THEN
             ID_SUB_BL(MJABS,I)  = MJ
             NORB_SUB(MJABS,I,2) = NORB_SUB(MJABS,I,2) + 1
             NORB_SUB(MJABS,I,0) = NORB_SUB(MJABS,I,0) + 1
             MJMAX = MAX(MJABS,MJMAX)
           ENDIF
        ENDDO
C
        N_SUB_BL(I) = MJMAX
C
C       Write information to output
C
        DO MJABS = 1, MJMAX
           WRITE (LUPRI,1010) 3-2*I,ID_SUB_BL(MJABS,I),
     &                       (NORB_SUB(MJABS,I,IC),IC=0,2)
        ENDDO
C
C       We need to be careful now. In making NORB_SUB we chose to
C       use the absolute value of MJ because this gives the most logical
C       block structure. The eigenvectors and solutions are, however,
C       still ordered by the signed value. Reorder them first before
C       continuing.
C
        DO IX = IEIGP, IEIGP+NS-1
           EIG(IX) = ABS(EIG(IX))
        ENDDO
C
C       The SELBOS routine does just the reordering that we want,
C       use work(kfss) as scratch, and allocate a few other small
C       scratch arrays. The IBEIG part will not be used but we will
C       initialize it to prevent problems with uninitialized data while
C       sorting it.
C
        CALL MEMGET('INTE',KINDX,NS,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIDUM,NS,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBEI,NS,WORK,KFREE,LFREE)
        MJ = 0
        CALL ICOPY (NS,MJ,0,WORK(KIBEI),1)
        CALL SELBOS (I,NS,WORK(KINDX),EIG(IEIGP),
     &               WORK(KIBEI),WORK(KVSS),
     &               WORK(KFSS),WORK(KIDUM),IPRHAM)
        CALL MEMREL('LINSYM',WORK,KINDX,KINDX,KFREE,LFREE)
C
C       Same thing for the positive energy solutions
C
        DO IX = IEIGE, IEIGE+NL-1
           EIG(IX) = ABS(EIG(IX))
        ENDDO
        CALL MEMGET('INTE',KINDX,NL,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIDUM,NL,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBEI,NL,WORK,KFREE,LFREE)
        MJ = 0
        CALL ICOPY (NL,MJ,0,WORK(KIBEI),1)
        CALL SELBOS (I,NL,WORK(KINDX),EIG(IEIGE),
     &               WORK(KIBEI),WORK(KVLL),
     &               WORK(KFLL),WORK(KIDUM),IPRHAM)
        CALL MEMREL('LINSYM',WORK,KINDX,KINDX,KFREE,LFREE)
C
C       Put the MO vectors in TBUF
C
C       We now have the solutions ordered on on the right index 
C       by absolute mj-value and at the left index on boson irreps.
C       Keep this ordering when inserting into the full matrix so
C       that the new transformation matrix (generated below) will
C       give an ordering on MJ values.
C

        CALL DZERO(TBUF(IOFT),NT*NT)

        CALL EXTRSB (2,I,NBORB,4,NBRP,
     &               NORB_SUB,MAX_SUB_BL,N_SUB_BL(I),0,0,1,1,1,1,
     &               WORK(KVLL),NL,NL,1,
     &               TBUF(IOFT),NT,NT,1)

        if(do4c2c)then ! x2c-mmf approach
          call writt(13,NL**2,WORK(KVLL))
        end if

        IF (IPRHAM.GE.15) THEN
         CALL HEADER('LINSYM: MO vectors in TBUF,LL',-1) 
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IOFT),NT,NT,NT,NT,1,
     &   IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)

         CALL HEADER('LINSYM: MO vectors in buffer space,LL',-1) 
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(WORK(KVLL),NL,NL,NL,NL,1,
     &   IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
        ENDIF

        CALL EXTRSB (2,I,NBORB,4,NBRP,
     &               NORB_SUB,MAX_SUB_BL,N_SUB_BL(I),0,0,2,2,1,1,
     &               WORK(KVSS),NS,NS,1,
     &               TBUF(IOFT),NT,NT,1)

        IF (IPRHAM.GE.15) THEN
         CALL HEADER('LINSYM: MO vectors in TBUF,SS',-1) 
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IOFT),NT,NT,NT,NT,1,
     &   IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
        ENDIF
C
C       Done; release the extra memory used for diagonalization
C
        CALL MEMREL('LINSYM',WORK,KFLL,KFLL,KFREE,LFREE)
C
C       Update offset
C
        IOFT = IOFT + NORB(I)*NORB(I)
C
   20 CONTINUE

      if(do4c2c)then ! x2c-mmf approach
        close(13,status='keep')
      end if
C
C     Generate transformation matrix in non-orthogonal basis
C     ======================================================
C
      IOFT = 1
      IOFU = 1
      IOFV = 1
      DO 30 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 30
        IF (IPRHAM.GE.5) THEN
         CALL HEADER(
     &'LINSYM: Transformation matrix in orthon. basis (MO->MO),TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IOFU),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),1,IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
        ENDIF
C
C       Backtransform all solutions
C       ===========================
C
        CALL BCKTRA(TMAT(IOFT),NFBAS(I,0),NORB(I),
     &              TBUF(IOFU),NORB(I),NORB(I),
     &              NORB(I),1,
     &              NORB(I),1,NFBAS(I,0),
     &              VMAT(IOFV),NFBAS(I,0),NORB(I),1,
     &              IPRHAM)

        IF (IPRHAM.GE.5) THEN
         CALL HEADER(
     &'LINSYM: Final (symm.blocked) transformation'//
     &' matrix (AO->MO),TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IOFT),NFBAS(I,0),NORB(I),
     &         NFBAS(I,0),NORB(I),1,IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
        ENDIF

        IOFT = IOFT + NFBAS(I,0)*NORB(I)
        IOFU = IOFU + NORB(I)*NORB(I)
        IOFV = IOFV + NORB(I)*NFBAS(I,0)
   30 CONTINUE
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('LINSYM')
!     iprham = iprham_save
      RETURN
 1000 FORMAT 
     & (/' Parity   MJ  Functions(total) Functions(LC) Functions(SC)')
 1010 FORMAT (T5,I2,T9,I3,'/2',T20,I5,T33,I5,T46,I5)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck amjmat */
      SUBROUTINE AMJMAT(FOCK,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Generate one-electron matrix of Jz
C     IPRINT = IPRHAM
C
C     Called from LINSYM 
C
C     Written by Luuk Visscher, oct 1999
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "cbihr1.h"
      DIMENSION FOCK(*),WORK(LWORK)
C
      CALL QENTER('AMJMAT')
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET('REAL',KFTRI ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN,   NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST,NZ        ,WORK,KFREE,LFREE)
C
C     *************************************************
C     ****** Construct one-electron matrix of Jz ******
C     *************************************************
C
      CALL AMJMA1(FOCK,WORK(KFTRI),WORK(KONEIN),WORK(KFIRST),
     &            WORK,KFREE,LFREE,IPRINT)
C
C     Memory deallocation
      CALL MEMREL('AMJMAT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('AMJMAT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck amjma1 */
      SUBROUTINE AMJMA1(FOCK,FTRI,OP1INT,FIRST,WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C
C   Generate one-electron matrix of Jz
C
C  Written by Luuk Visscher, oct 1999
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbham.h"
C
      PARAMETER (HALF=0.5D0, D1=1.0D0)
      LOGICAL FIRST
      DIMENSION FOCK(*),FTRI(NNBBASX),OP1INT(NNBBASX),
     &          WORK(*)
      CALL QENTER('AMJMA1')
C
C     Sanity check : linear symmetry requires use of real groups (C2v or D2h)
C
      IF (NZ.NE.1) CALL QUIT ("Linear symmetry requires NZ=1")
C
#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
C
C     Get integrals of component of the z-angular momentum that we want to
C     quantize on.
C
      FIRST = .TRUE.
      IOP =  IPANGMOM
      CALL PRPMAO(LU1INT,IOP,NNBBASX,FTRI,
     &         OP1INT,FIRST,IPRHAM)
      CALL DSPTSI(NTBAS(0),FTRI,FOCK)
C     
C     Get integrals of the spin-z (sigma_z) operator
C     ===============================================
C     
      FIRST = .TRUE.
      IOP =  IPSPNMOM
      CALL PRPMAO(LU1INT,IOP,NNBBASX,FTRI,
     &        OP1INT,FIRST,IPRHAM)
      CALL MEMGET('REAL',KBUF,N2BBASX,WORK,KFREE,LFREE)
      CALL DSPTSI(NTBAS(0),FTRI,WORK(KBUF))
C     
C     Make the Jz operator (in case of spinfree we make Lz)
C     =====================================================
C     
      IF (.NOT.SPINFR)
     &   CALL DAXPY (N2BBASX,D1,WORK(KBUF),1,FOCK,1)
C
C     Release memory & close AOPROPER
C
      CALL MEMREL('AMJMA1',WORK,KBUF,KBUF,KFREE,LFREE)
      CLOSE(LU1INT,STATUS = 'KEEP')
C
C     Reindex to sorted basis
C     =======================
C
      CALL BUTOBS(FOCK,1,WORK(KFREE),LFREE)
C
C     Output section
C     ==============
C
      IF(IPRHAM.GE.4) THEN
        CALL HEADER(
     &  'AMJMA1: One-electron MJ matrix in symm.blocked SA-AO basis',-1)
        DO 10 I = 1,NFSYM
          IF(NFBAS(I,0).EQ.0) GOTO 10
          WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(FOCK(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
   10   CONTINUE
      ENDIF
C     
      CALL QEXIT('AMJMA1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Wrifrz */
      SUBROUTINE WRIFRZ(IUNIT,IFRP,CSEL,ESEL,ISEL,IVEC,BUF,IBUF,
     &                  NVEC,NDEL)
C***********************************************************************
C
C     Dump frozen orbitals to file
C     Written by T. Saue May 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
      DIMENSION CSEL(*),ESEL(*),ISEL(*),
     &          BUF(*),IVEC(*),IBUF(*),NDEL(2)
C
C     Get index array for frozen orbitals;
C
      NTOT=1
      NESH(IFRP) = NESH(IFRP) + NDEL(1)
      NPSH(IFRP) = NPSH(IFRP) + NDEL(2)
      CALL NUMLS3(VCFROZ(IFRP),IVEC,NFBAS(IFRP,0),
     &            -NPSH(IFRP),NESH(IFRP),NTOT,NZERO)
C
C     Consistency check
C
      IF(NFRO(IFRP).NE.(NTOT-NZERO)) THEN
        CALL QUIT('WRIFRZ:Inconsistency !')
      ENDIF
C
C     Selected orbitals that are not frozen have index zero;
C     they must be taken out and pointers made to the remaining ones
C     Write compressed index array to file
C
C     We recompute NASHMFT as this value may change if there are frozen open shells
C
      NE = 0
      NP = 0
      NASHMFT = 0
      II = 0
      DO I = 1,NTOT
      IP = IVEC(I)
      IF(IP.NE.0) THEN
        II = II + 1
        IBUF(II)=I        
C       Positronic position
        IF(IP.LT.0.AND.IP.GE.(-NPSH(IFRP))) THEN
          IVEC(II)=1+IP
          NP = NP + 1
C       Electronic position
        ELSEIF (IP.GT.0.AND.IP.LE.NESH(IFRP)) THEN
          IVEC(II)=IP
          NE = NE + 1
C         Check if orbital belongs to closed or open shells.
C         If so, subtract from occupation number
          IF(IP.GT.0.AND.IP.LE.NISH(IFRP)) THEN
            NISHMF(IFRP) = NISHMF(IFRP) - 1
          ELSE
            I1 = NISH(IFRP)
            DO IOPEN = 1,NOPEN
              I2 = I1 + NACSH(IFRP,IOPEN)
              IF(IP.GT.I1.OR.IP.LE.I2) THEN
                NACSHMF(IFRP,IOPEN) = NACSHMF(IFRP,IOPEN) - 1
                GOTO 20
              ENDIF
              I1 = I2
            ENDDO
 20         CONTINUE
          ENDIF
        ENDIF
      ENDIF
      ENDDO
      NOCCMF(IFRP) = NISHMF(IFRP)
      DO IOPEN = 1,NOPEN
        NOCCMF(IFRP) = NOCCMF(IFRP) + NACSHMF(IFRP,IOPEN)
        NASHMFT = NASHMFT + NACSHMF(IFRP,IOPEN)
      ENDDO
C
C     Final adjustments
C
      NPSHMF(IFRP) = NPSH(IFRP) - NP
      NESHMF(IFRP) = NESH(IFRP) - NE
      DO I = 1,NFRO(IFRP)
        IVEC(I) = IVEC(I) + NPSH(IFRP) - NDEL(2)
      ENDDO
C
C     Write index array to file
C
      CALL WRITI(IUNIT,NFRO(IFRP),IVEC)
C
C     Write coefficients and eigenvalues to file
C
      IDIM = NFBAS(IFRP,0)*NVEC
      JDIM = NFBAS(IFRP,0)*NZ
      DO II = 1,NFRO(IFRP)
        I = IBUF(II)
        IOFF = NFBAS(IFRP,0)*(I-1)+1
        JOFF = 1
        DO IZ = 1,NZ
          CALL DCOPY(NFBAS(IFRP,0),CSEL(IOFF),1,BUF(JOFF),1)
          IOFF = IOFF + IDIM
          JOFF = JOFF + NFBAS(IFRP,0)
        ENDDO
        CALL WRITT(IUNIT,JDIM,BUF)
        WRITE(IUNIT) ESEL(I),ISEL(I)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Inisub */
      SUBROUTINE INISUB
C***********************************************************************
C
C
C     Copy information from NBORB to the subblock array if we have not
C     partioned the orbital space otherwise.
C
C     Written by Luuk and copied by Trond Oct 19 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dgroup.h"
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
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck DELPOS */
      SUBROUTINE DELPOS(TMAT,TMOD,FMAT,EIG,IBEIG,VMAT,WORK,KFREE,LFREE)
C***********************************************************************
C
C
C     Deleting negative-energy solutions of:
C     ======================================
C       FREEPJ - free particle matrix
C       VEXTPJ - bare nucleus matrix          
C
C     On input:
C       TMAT(NFORB,NFORB,IFRP) 
C            - free particle/bare nucleus matrix
C              in orthonormal canonical Lowdin basis
C       TMOD(NFORB,NORB,IFRP) 
C            - transformation matrix with modified Dirac equation
C              embedded, in orthonormal basis
C       FMAT - empty
C     On output:
C       TMAT - transformation matrix with negative-energy solutions deleted
C     Written by T.Saue, October 1995
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "dgroup.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
      DIMENSION TMAT(*),TMOD(*),FMAT(*),EIG(*),IBEIG(*),VMAT(*),WORK(*)
C
      NSYM = 4/NZ
C
C     Transform free particle/bare nucleus matrix to modified Dirac form; 
C
      IOFT = 1
      IOFV = 1
      DO 10 I = 1,NFSYM
        IF(NTMO(I).EQ.0) GOTO 10
        CALL QTRANS('AOMO','S',D0,NFORB(I,0),NFORB(I,0),
     &       NTMO(I),NTMO(I),
     &       TMAT(IOFT),NFORB(I,0),NFORB(I,0),NZ,IPQTOQ(1,0),
     &       FMAT(1+I2TMOT(I)),NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),
     &       TMOD(1+I2TMT(I)),NFORB(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &       TMOD(1+I2TMT(I)),NFORB(I,0),NTMO(I),NZT,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)
        IF (IPRHAM.GE.6) THEN
          IF(VEXTPJ) THEN
            CALL HEADER('DELPOS: '//
     &      'Bare nucleus Dirac matrix in mod. Dirac eq. MO-basis',-1)
          ELSE
            CALL HEADER('DELPOS: '//
     &      'Free particle Dirac matrix in mod. Dirac eq. MO-basis',-1)
          ENDIF
          CALL PRQMAT(FMAT(1+I2TMOT(I)),NTMO(I),NTMO(I),
     &              NTMO(I),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
C       get full transformation as well !
        CALL BCKTR1(TMAT(1+I2ORBT(I)),NFBAS(I,0),NTMO(I),
     &              TMOD(1+I2TMT(I)),NFORB(I,0),NTMO(I),
     &              NFORB(I,0),NZ,NTMO(I),1,NFBAS(I,0),
     &              VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,
     &              IPRHAM)
        IF (IPRHAM.GE.6) THEN
          CALL HEADER('DELPOS: Transformation matrix',-1)
          CALL PRQMAT(TMAT(1+I2ORBT(I)),NFBAS(I,0),NTMO(I),
     &                NFBAS(I,0),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
        IOFT = IOFT + NFORB(I,0)*NFORB(I,0)*NZ
        IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)
 10   CONTINUE
C
C     Diagonalize: eigenvectors are now stored in TMOD
C
      CALL DFDIAG (FMAT,EIG,IBEIG,TMOD,.FALSE.,WORK,KFREE,LFREE)
C
C     Backtransform only electronic solutions
C
      IBUF = 1
      DO 20 I = 1,NFSYM
        IF(NESH(I).EQ.0) GOTO 20
          IF(IPRHAM.GE.5) THEN
            CALL HEADER('DELPOS: Eigenvalues',-1)
            WRITE(LUPRI,'(I5,E16.8)') (J,EIG(IORB(I)+J),J=1,NTMO(I))
          ENDIF
          CALL BCKTR1(FMAT(IBUF),NFBAS(I,0),NESH(I),
     &                TMOD(1+I2TMOT(I)),NTMO(I),NTMO(I),NTMO(I),NZ,
     &                NESH(I),(NPSH(I)+1),NFBAS(I,0),
     &                TMAT(1+I2ORBT(I)),NFBAS(I,0),NTMO(I),NZT,IPRHAM)
          NPSH(I) = 0
          NORB(I) = NESH(I)
          NBUF = NFBAS(I,0)*NESH(I)*NZ
          CALL DCOPY(NBUF,FMAT(IBUF),1,TMAT(IBUF),1)
          IBUF = IBUF + NBUF
          NTMO(I) = NORB(I)
C         Finish by adjusting the number of positron solutions to zero
          DO ISYM = 1, NSYM
            NBORB(ISYM,I,2) = 0
            NBORB(ISYM,I,0) = NBORB(ISYM,I,1) + NBORB(ISYM,I,2)
          ENDDO
 20     CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck LL_EQ */
      SUBROUTINE LL_EQ(TMAT,TBUF,VMAT)
C***********************************************************************
C
C     The Levy-Leblond equation:
C     ==========================
C
C       Set up basis for Levy-Leblond equation as:
C
C              [                (Psi(L), 0 ]
C        Psi = [                           ]
C              [(2mc)^-1 sigma.p Psi(L), 0 ]
C
C     Reference:
C       L. Visscher and T. Saue, J. Chem. Phys. 113(2000) 3996
C       "Approximate relativistic electronic structure methods based on 
C       the quaternion modified Dirac equation"
C
C     On input:
C       TMAT - free particle matrix in Lowdin orthonormal basis     
C       VMAT - Lowdin canonical transformation matrix
C
C     On output:
C       TMAT - transformation matrix with Levy-Leblond equation embedded
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, D2 = 2.0D0, D0 = 0.0D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "cbihr1.h"
C
      DIMENSION TMAT(*),TBUF(*),VMAT(*)
C      
      FAC  = D2*CVAL*CVAL
      FAC  = D1/FAC
      NSYM = 4/NZ
C
      IBUF  = 1
      IMAT  = 1
      DO 10 IFRP = 1,NFSYM
        NESH(IFRP)  = NFORB(IFRP,1)
        NPSH(IFRP)  = 0
        NORB(IFRP)  = NESH(IFRP) + NPSH(IFRP)
        NTMO(IFRP)  = NORB(IFRP)
        IF(NORB(IFRP).EQ.0) GOTO 10
        CALL DZERO(TBUF(IBUF),NFORB(IFRP,0)*NORB(IFRP)*NZT)
C.......set LL block to unit matrix
        JJ = IBUF
        DO J = 1, NFORB(IFRP,1)
          TBUF(JJ) = D1
          JJ = JJ + NFORB(IFRP,0) + 1
        ENDDO
C.......set SL block to (2mc)^-1 sigma.p
        DO IZ = 1, NZ
          JSL = IBUF + NFORB(IFRP,1)
          ISL = IMAT + NFORB(IFRP,1)
          DO J = 1, NFORB(IFRP,1)
            CALL DCOPY (NFORB(IFRP,2),TMAT(ISL),1,TBUF(JSL),1)
            CALL DSCAL (NFORB(IFRP,2),FAC,TBUF(JSL),1)
            ISL = ISL + NFORB(IFRP,0)
            JSL = JSL + NFORB(IFRP,0)
          ENDDO
          IF (IPRHAM.GE.6) THEN
           CALL HEADER('Transformation matrix in orthon. basis',-1)
           CALL PRQMAT(TBUF(IBUF),NFORB(IFRP,0),NORB(IFRP),
     &                 NFORB(IFRP,0),NORB(IFRP),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
          IMAT = IMAT + NFORB(IFRP,0) * NFORB(IFRP,0)
          IBUF = IBUF + NFORB(IFRP,0) * NORB(IFRP)
        ENDDO
C.......finish by adjusting the number of positron solutions to zero
        DO ISYM = 1, NSYM
          NBORB(ISYM,IFRP,2) = 0
          NBORB(ISYM,IFRP,0) = NBORB(ISYM,IFRP,1) + NBORB(ISYM,IFRP,2)
        ENDDO
   10 CONTINUE
C
C    Generate transformation matrix
C    ==============================
C
      IBUF = 1
      IMAT = 1
      IOFV = 1
      DO 20 IFRP = 1,NFSYM
        IF(NORB(IFRP).EQ.0) GOTO 20
        CALL BCKTR1(TMAT(IMAT),NFBAS(IFRP,0),NFORB(IFRP,1),
     &              TBUF(IBUF),NFORB(IFRP,0),NFORB(IFRP,1),
     &              NFORB(IFRP,0),NZ,
     &              NFORB(IFRP,1),1,NFBAS(IFRP,0),
     &              VMAT(IOFV),NFBAS(IFRP,0),NFORB(IFRP,0),1,
     &              IPRGEN)
        IMAT = IMAT + NFBAS(IFRP,0)*NORB(IFRP)*NZ
        IBUF = IBUF + NFORB(IFRP,0)*NORB(IFRP)*NZ
        IOFV = IOFV + NFORB(IFRP,0)*NFBAS(IFRP,0)
   20   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MOD_DIREQ */
      SUBROUTINE MOD_DIREQ(TMAT,TBUF,VMAT,EIG,WORK,LWORK)
C***********************************************************************
C
C       The modified Dirac eqation:
C       ===========================
C
C       Basis as in Levy-Leblond but without fixed ratio
C
C
C               [                (Psi(L),                      0  ]
C         Psi = [                                                 ]
C               [                      0, (2mc)^-1 sigma.p Psi(L) ]
C
C     Reference:
C       L. Visscher and T. Saue, J. Chem. Phys. 113(2000) 3996
C       "Approximate relativistic electronic structure methods based on 
C       the quaternion modified Dirac equation"
C
C     On input:
C       TMAT - free particle matrix in Lowdin orthonormal basis     
C       VMAT - Lowdin canonical transformation matrix
C
C     On output:
C       TMAT - transformation matrix with modified Dirac equation embedded
C
C     Note that there is a difference in boson irrep ordering 
C     between TMAT and the Lowdin canonical transformation matrix VMAT:
!
!
C     After transformation with VMAT the boson irrep ordering
C     in a given fermion ircop is:
C       L1 - L2 - L3 .. S1 - S2 - S3 - ..
C     whereas with TMAT the ordering is:
C       S1 - L1 - S2 - L2 - S3 - L3 ..
C
C     where L and S refer to the large and small components,
C     respectively.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, DM2 = -2.0D0, D0 = 0.0D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "cbihr1.h"
      DIMENSION TMAT(*),TBUF(*),VMAT(*),EIG(*),WORK(*)
#include "memint.h"

      NSYM = 4/NZ
      IMAT = 1
      IBUF = 1
      DO 10 I = 1,NFSYM
        NESH(I)  = NFORB(I,1)
        NPSH(I)  = NFORB(I,1)
        NORB(I)  = NESH(I) + NPSH(I)
        NTMO(I)  = NORB(I)
        IF(NORB(I).EQ.0) GOTO 10
C
C       Prepare small component part
C       ----------------------------
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
C         Make (sigma.p)(sigma.p) = p2
C
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
            CALL HEADER('Kinetic energy matrix in orthon. basis',-1)
            CALL PRQMAT(WORK(KSMT),NBORB(ISYM,I,1),NBORB(ISYM,I,1),
     &                  NBORB(ISYM,I,1),NBORB(ISYM,I,1),
     &                  1,IPQTOQ(1,0),LUPRI)
          ENDIF
          CALL LOWGEN(WORK(KSMT),NBORB(ISYM,I,1),WORK(IVMT),NEFFI,
     &                NBORB(ISYM,I,1),STOL(2),IPRONE,EIG,
     &                WORK(KFREE),LFREE)
          IF (NEFFI.LT.NBORB(ISYM,I,1)) WRITE(LUPRI,'(A)') 
     &            '*** WARNING *** Linear dependence in RKB basis'
          NEFF = NEFF + NEFFI
          NBORB(ISYM,I,2) = NEFFI
          IF (IPRHAM.GE.7) THEN
            CALL HEADER('Ext.transf. matrix in orthon. basis',-1)
            CALL PRQMAT(WORK(IVMT),NBORB(ISYM,I,1),NBORB(ISYM,I,2),
     &                  NBORB(ISYM,I,1),NBORB(ISYM,I,2),
     &                  1,IPQTOQ(1,0),LUPRI)
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
        NTMO(I)  = NORB(I)
C
        CALL DZERO(TBUF(IBUF),NFORB(I,0)*NORB(I)*NZ)
        IBO  = 0
        IBOL = 0
        IVMT = KVMT
        DO ISYM = 1,NSYM
        IF(NBORB(ISYM,I,1).GT.0) THEN
          IREP = JFSYM(ISYM,I) - 1
C
C         Small component part
C
          JBUF = IBUF + NFORB(I,0)*IBO + NFORB(I,1)
          JSL  = IMAT + IBOL * NFORB(I,0) + NFORB(I,1)
          CALL QGEMM(NFORB(I,2),NBORB(ISYM,I,2),NBORB(ISYM,I,1),D1,
     &         'N','N',IPQTOQ(1,0),TMAT(JSL),NFORB(I,0),NFORB(I,0),NZ,
     &         'N','N',IPQTOQ(1,0),WORK(IVMT),
     &                 NBORB(ISYM,I,1),NBORB(ISYM,I,2),1,
     &         D0,IPQTOQ(1,0),TBUF(JBUF),NFORB(I,0),NORB(I),NZ)
          IBO = IBO + NBORB(ISYM,I,2)
C
C         Large component part
C
          JBUF = IBUF + NFORB(I,0)*IBO + IBOL
          DO J = 1, NBORB(ISYM,I,1)
            TBUF(JBUF) = D1
            JBUF = JBUF + NFORB(I,0) + 1
          ENDDO
          IBO = IBO + NBORB(ISYM,I,1)
C
C         Update offsets
C
          IBOL = IBOL + NBORB(ISYM,I,1)
          IVMT = IVMT + NBORB(ISYM,I,1)*NBORB(ISYM,I,2)
        ENDIF
        ENDDO
        IF (IPRHAM.GE.6) THEN
         CALL HEADER('Transformation matrix in orthon. basis',-1)
         CALL PRQMAT(TBUF(IBUF),NFORB(I,0),NORB(I),
     &                 NFORB(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
        IMAT = IMAT + NFORB(I,0)*NFORB(I,0)*NZ
        IBUF = IBUF + NFORB(I,0)*NORB(I)*NZ
   10 CONTINUE
      CALL MEMREL('MOD_DIREQ',WORK,KSMT,KSMT,KFREE,LFREE)
C
      IF(VEXTPJ.OR.FREEPJ) THEN
C     These variables are later set by SETDC2, but needed here...
        N2TMT  = 0
        NORBT  = 0
        N2ORBT = 0
        N2TMOTQ = 0
        DO IFRP = 1,NFSYM
          I2TMT(IFRP)  = N2TMT
          N2TMT        = N2TMT + NFORB(IFRP,0)*NORB(IFRP)*NZT
          IORB(IFRP)   = NORBT
          NORBT        = NORBT  + NORB(IFRP)
          I2ORBT(IFRP) = N2ORBT*NZ
          N2ORBT       = N2ORBT + NFBAS(IFRP,0)*NORB(IFRP)
          I2TMOT(IFRP) = N2TMOTQ
          N2TMO(IFRP)  = NTMO(IFRP)*NTMO(IFRP)
          N2TMOTQ      = N2TMOTQ + NTMO(IFRP)*NTMO(IFRP)*NZ
        ENDDO
        CALL MEMGET('REAL',KBUF,N2ORBT*NZ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KBEIG,NTBAS(0),WORK,KFREE,LFREE)
        CALL DELPOS(TMAT,TBUF,WORK(KBUF),EIG,WORK(KBEIG),
     &              VMAT,WORK,KFREE,LFREE)
        CALL MEMREL('MOD_DIREQ.delpos',WORK,KBUF,KBUF,KFREE,LFREE)
      ELSE
C
C     Generate transformation matrix
C     ==============================
C
        IMAT = 1
        IBUF = 1
        IOFV = 1
        DO 20 I = 1,NFSYM
          IF(NFORB(I,0).EQ.0) GOTO 20
          CALL BCKTR1(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                TBUF(IBUF),NFORB(I,0),NORB(I),
     &                NFORB(I,0),NZ,
     &                NORB(I),1,NFBAS(I,0),
     &                VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,
     &                IPRGEN)
          IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
          IBUF = IBUF + NFORB(I,0)*NORB(I)*NZ
          IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)
   20   CONTINUE
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck psipr1 */
      SUBROUTINE PSIPR1(NVECS,CSEL,ESEL,ISEL,TMAT,TINV,NSTR,
     &                  IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
C
      LOGICAL TOBE
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbdhf.h"
      DIMENSION CSEL(*),ESEL(*),ISEL(*),TMAT(*),TINV(*),WORK(*)
      DIMENSION NVECS(2),NQ(2),NSTR(2,0:2,MFRAG)
C
      KFRSAV = KFREE
      NFROT = 0
      DO IFRP = 1,NFSYM
        NFROT = NFROT + NFRO(IFRP)
      ENDDO
      IF(NFROT.GT.0) CALL OPNFIL(LUCOEF,'DFFROZ','UNKNOWN','PSIPR1')
      TOBE = .FALSE.
      N2TMT = 0
      ICOFF = 1
      IEOFF = 1
      DO 10 IFRP = 1,NFSYM
        WRITE(LUPRI,'(A,A3,A2)') '* Fermion ircop ',FREP(IFRP),' :'
C       Transform fragment coefficients to orthonormal basis
        NDEL = 0
        NQ(1) = 0
        NQ(2) = 0
        IF(NVECS(IFRP).EQ.0) THEN
           WRITE(LUPRI,'(2X,A)') '- No orbitals deleted.'
          GOTO 11
        ENDIF
        NREF = NORB(IFRP)*NVECS(IFRP)*NZ
        NPRJ = NORB(IFRP)*NORB(IFRP)*NZ
        CALL MEMGET('REAL',KBUF,NPRJ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBOS,NORB(IFRP),WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KREF,NREF,WORK,KFREE,LFREE)
        CALL QGEMM(NORB(IFRP),NVECS(IFRP),NFBAS(IFRP,0),D1,
     &            'H','N',IPQTOQ(1,0),TINV(I2TMT(IFRP)+1),
     &                    NFBAS(IFRP,0),NORB(IFRP),NZT,
     &            'N','N',IPQTOQ(1,0),CSEL(ICOFF),
     &                    NFBAS(IFRP,0),NVECS(IFRP),NZ,
     &            D0,IPQTOQ(1,0),WORK(KREF),NORB(IFRP),NVECS(IFRP),NZ)
        IF(IPRINT.GE.5) THEN
          CALL HEADER('PSIPR1: Fragment coefs in orthonormal basis',-1)
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',IFRP
          CALL PRQMAT(WORK(KREF),NORB(IFRP),NVECS(IFRP),
     &                NORB(IFRP),NVECS(IFRP),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
C       Initialize unit matrix (transformation matrix in MO-basis)
        CALL DUNIT2(WORK(KBUF),NORB(IFRP),NORB(IFRP),NORB(IFRP),NZ)
        IF(SUB_BL) THEN
C         Spin-free or linear case: Fill array with subblock info
          CALL INIBOS(WORK(KIBOS),IFRP,.FALSE.,0)
        ELSE
C         Spin-orbit coupling: Fill array with zeros..
          CALL ICOPY (NORB(IFRP),0,0,WORK(KIBOS),1) 
        ENDIF
C       Do projection
        IREF  = 1
        DO IFRAG = 1,NFRAG
          DO IP = 2,1,-1
            NREF = NSTR(IFRP,IP,IFRAG)        
            IF(NREF.GT.0) THEN
              CALL PRJORT(WORK(KBUF),WORK(KIBOS),
     &              WORK(KREF),ISEL(IEOFF),IREF,NREF,
     &              PRJTHR,THRMAX,NORB(IFRP),NORB(IFRP),
     &              NVECS(IFRP),NDEL,IP,IPRHAM)
              IREF = IREF + NREF
              NQ(IP) = NQ(IP) + NDEL
            ENDIF
          ENDDO
        ENDDO
        CALL MEMREL('PSIPRJ.prjort',WORK,KFRSAV,KREF,KFREE,LFREE)
        NDEL = NQ(1) + NQ(2)
        WRITE(LUPRI,'(2X,A,I6)') 
     &     '- Total number of orbitals deleted: ',NDEL,
     &     '       - positive energy solutions: ',NQ(1),
     &     '       - negative energy solutions: ',NQ(2)
 11     CONTINUE
        NTMO(IFRP)  = NORB(IFRP)-NDEL
        NESH(IFRP)  = NESH(IFRP)-NQ(1)
        NPSH(IFRP)  = NPSH(IFRP)-NQ(2)
        IEFF        = I2TMT(IFRP)+1
        I2TMT(IFRP) = N2TMT
        N2TMT = N2TMT + NFBAS(IFRP,0)*NTMO(IFRP)*NZT
        IF(NFRO(IFRP).GT.0) THEN
          IF(NDEL.LT.NFRO(IFRP)) THEN
            WRITE(LUPRI,'(2X,A,I6,A)') 
     &      'is smaller than number of frozen orbitals',
     &      NFRO(IFRP),'!'
            WRITE(LUPRI,'(2X,A)')
     &      'You may rerun with PRINT.GE.2 under **HAMILTONIAN',
     &      'to see the norms of the remaining vectors.',
     &      'Based on this, you may set the projection threshold'
            WRITE(LUPRI,'(2X,A,E12.4)')
     &      'PRJTHR lower. Current value is ',PRJTHR
            CALL QUIT
     &      ('PSIPRJ:Fewer deleted orbitals than frozen ones!')
          ENDIF
C         Backtransform to SO-basis
          NDIM = NFBAS(IFRP,0)*NORB(IFRP)*NZT
          CALL DCOPY(NDIM,TMAT(IEFF),1,TINV(IEFF),1)
          CALL PRJCNT(WORK(KBUF),TMAT(I2TMT(IFRP)+1),
     &                IFRP,NORB(IFRP),NFBAS(IFRP,0),
     &                NTMO(IFRP),WORK(KIBOS),TINV(IEFF))
          CALL MEMREL('PSIPRJ.prj',WORK,KFRSAV,KBUF,KFREE,LFREE)
C
C         Save frozen orbitals to file
C         
C
          NBUF = NFBAS(IFRP,0)*NZ
          CALL MEMGET('INTE',KSEL,NVECS(IFRP),WORK,KFREE,LFREE)
          CALL MEMGET('INTE',KIND,NFRO(IFRP),WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KBUF,NBUF,WORK,KFREE,LFREE)
          CALL WRIFRZ(LUCOEF,IFRP,CSEL(ICOFF),ESEL(IEOFF),ISEL(IEOFF),
     &                  WORK(KSEL),WORK(KBUF),WORK(KIND),NVECS(IFRP),NQ)
          CALL MEMREL('PSIPRJ.froz',WORK,KFRSAV,KBUF,KFREE,LFREE)       
          TOBE = .TRUE.       
        ELSEIF((I2TMT(IFRP)+1).NE.IEFF) THEN
C         Move transformation matrix if columns was deleted in previous
C         irrep ...
          NDIM = NFBAS(IFRP,0)*NTMO(IFRP)*NZT
          CALL DCOPY(NDIM,TMAT(IEFF),1,TMAT(I2TMT(IFRP)+1),1)
        ENDIF
        ICOFF = ICOFF + NFBAS(IFRP,0)*NVECS(IFRP)*NZ
        IEOFF = IEOFF + NVECS(IFRP)
 10   CONTINUE
      IF(NFROT.GT.0) CLOSE(LUCOEF,STATUS='KEEP')
      IF(IPRHAM.GE.5) THEN
        CALL HEADER('PSIPRJ: Final transformation matrix ',-1)
        DO 20 I = 1,NFSYM
        IF(NTMO(I).GT.0) THEN
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',I
          CALL PRQMAT(TMAT(I2TMT(I)+1),NFBAS(I,0),NTMO(I),
     &                NFBAS(I,0),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
 20     CONTINUE
      ENDIF
      IF(TOBE) THEN
C
C       Construct right index transformed overlap matrix Saomo = Sao T
C       Since T(dagger) Sao T = I, Saomo is the inverse of T(dagger)
C
        CALL MKSAOMO(TINV,TMAT,.FALSE.,IPRINT,WORK(KFREE),LFREE)
      ENDIF
C
      RETURN
      END
! -- end of dirone.F --
