      SUBROUTINE WRT_CC_VEC2(CC,LU,CCTYPE)
*
* Print vector of CC amplitudes, type is governed by CCTYPE
*     
#include "implicit.inc"
*. General input 
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "crun.inc"
#include "ctcc.inc"
#include "wrkspc.inc"

      CHARACTER*6 CCTYPE
      DIMENSION CC(*) 
*     
      IF(CCTYPE(1:2).EQ.'CC') THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' Single excitation part '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        CALL WRT_SX(CC,1)
*
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' Double excitation part '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        CALL WRT_DX1(CC(1+NSXE),1)
      ELSE IF (CCTYPE(1:6).EQ.'GEN_CC') THEN
        CALL WRTBLKN_CC(CC,NSPOBEX_TP,WORK(KLLSOBEX))
      END IF
*
      RETURN
      END
*
      SUBROUTINE WRT_CC_VEC_REL(CC,LU)
*
* Print vector of CC amplitudes 
* Timo Fleig, Winter 2009/2010
*
#include "implicit.inc"
*. Generel input 
#include "mxpdim.inc"
#include "crun.inc"
      DIMENSION CC(*)
*
      WRITE(6,*)
      WRITE(6,*) ' ====================== '
      WRITE(6,*) ' Single excitation part '
      WRITE(6,*) ' CAUTION : This is not yet adapted to the rel. case !'
      WRITE(6,*) ' ====================== '
      WRITE(6,*)
      CALL WRT_SX_REL(CC,1)
*
      WRITE(6,*)
      WRITE(6,*) ' ====================== '
      WRITE(6,*) ' Double excitation part '
      WRITE(6,*) ' CAUTION : This is not yet adapted to the rel. case !'
      WRITE(6,*) ' ====================== '
      WRITE(6,*)
      CALL WRT_DX_REL(CC(1+NSXE),1)
*
      RETURN
      END
*
      SUBROUTINE WRT_SX_REL(C,ISX_SYM)
*
* Write list of single excitations given in  compact form
*
* Jeppe Olsen, Summer of 98
* Timo Fleig, Winter of 2010
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "multd2h.inc"
*. Input 
      DIMENSION C(*)
* 
      WRITE(6,*)
      WRITE(6,*) ' List of single excitations in compressed form '
      WRITE(6,*) ' ============================================= '
      WRITE(6,*)
C
      IJ_OFF = 1
      print*,'NSMOB ', NSMOB



      DO ISM =1, NSMOB
        JSM = MULTD2H(ISM,ISX_SYM)
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
C                   I_SX_CCACT(IGAS,JGAS)
            IJACT = I_SX_CCACT(IGAS,JGAS)
            IF(IJACT.EQ.1) THEN
*
              NI = NOBPTS(IGAS,ISM)
              NJ = NOBPTS(JGAS,JSM)
*
*
              IF(NI*NJ.GT.0) THEN
                WRITE(6,'(A,A, 4I3)')
     &          ' Block of single excitations with sym and type ',
     &          ' (ism,itp,jsm,jtp)', ISM,IGAS,JSM,JGAS
                CALL WRTMAT(C(IJ_OFF),NI,NJ,NI,NJ)
                IJ_OFF = IJ_OFF + NI*NJ
              END IF
*             ^ End if nonvanishing block
            END IF
*           ^ End if active block
          END DO
        END DO
      END DO
*
      RETURN
      END
*
      SUBROUTINE WRT_DX_REL(C,IDX_SYM)
*
*
* Print list of double excitations in compressed form
* - without singlet-singlet, triplet-triplet separation
*
* Jeppe Olsen, summer of 98
* Timo Fleig, Winter of 2010
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "lucinp.inc"
*. Input 
      DIMENSION C(*)
*
      IJBL_GT_KLBL = 0
*
C     WRITE(6,*)
C     WRITE(6,*) ' ================================================'
C     WRITE(6,*) ' List of double excitations without S/T splitting'
C     WRITE(6,*) ' ================================================'
C     WRITE(6,*) 
      WRITE(6,*)
     & '  (Blocks for E(IJ) E(KL) written as matrices C(IJ,KL) )'
      IJKL_OFF = 1
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
      DO KSM = 1, ISM
       IJSM = MULTD2H(ISM,JSM)
       IJKSM = MULTD2H(IJSM,KSM)
       LSM   = MULTD2H(IJKSM,IDX_SYM)
       IF(ISM.GT.KSM.OR.(ISM.EQ.KSM.AND.JSM.GT.LSM)) THEN
         IJSM_GT_KLSM = 1
       ELSE IF( ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
         IJSM_GT_KLSM = 0
       ELSE
         IJSM_GT_KLSM = -1
       END IF
       IF( IJSM_GT_KLSM.GE.0) THEN
         DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
         DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
          IJKL_ACT = I_DX_CCACT(IGAS,KGAS,JGAS,LGAS)
*. Check of block fulfills (IJ.GE.KL)
          IF( IJSM_GT_KLSM .EQ. 1 ) THEN
            IJBL_GT_KLBL = 1
          ELSE IF ( IJSM_GT_KLSM .EQ. 0 ) THEN
            IF(IGAS.GT.KGAS.OR.(IGAS.EQ.KGAS.AND.JGAS.GT.LGAS)) THEN
              IJBL_GT_KLBL = 1
            ELSE IF(IGAS.EQ.KGAS.AND.JGAS.EQ.LGAS) THEN
              IJBL_GT_KLBL = 0
            ELSE
              IJBL_GT_KLBL = -1
            END IF
          END IF
          IF(IJKL_ACT.EQ.1 .AND. IJBL_GT_KLBL.GE.0 ) THEN
*
            NI = NOBPTS(IGAS,ISM)
            I_OFF = IOBPTS(IGAS,ISM)
*
            NJ = NOBPTS(JGAS,JSM)
            J_OFF = IOBPTS(JGAS,JSM)
*
            NK = NOBPTS(KGAS,KSM)
            K_OFF = IOBPTS(KGAS,KSM)
*
            NL = NOBPTS(LGAS,LSM)
            L_OFF = IOBPTS(LGAS,LSM)
*
            IF(NI*NJ*NK*NL.GT.0) THEN
              WRITE(6,'(A,8I3)')
     &        ' Orbital indeces I,J,K,L (type and sym)',
     &          IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM
              IF(IJBL_GT_KLBL.EQ.1) THEN
                CALL WRTMAT(C(IJKL_OFF),NI*NJ,NK*NL,NI*NJ,NL*NL)
                IJKL_OFF = IJKL_OFF + NI*NJ*NK*NL
              ELSE IF (IJBL_GT_KLBL .EQ. 0) THEN
                CALL PRSM2(C(IJKL_OFF),NI*NJ)
                IJKL_OFF = IJKL_OFF + NI*NJ*(NI*NJ+1)/2
              END IF
            END IF
*           ^ End if nonvanishing block
          END IF
*         ^ End if allowed block
         END DO
         END DO
         END DO
         END DO
*        ^ End of loop over gasspaces
       END IF
*      ^ End if IJ_SM .GT. KL_SM
      END DO
      END DO
      END DO
*     ^ End of loop over orbital symmetries
*
      RETURN
      END
*
      SUBROUTINE WRTBLKN_CC(VEC,NBLOCK,LBLOCK)
*     
* Write the NBLOCK blocks of VEC
*  Offset 1 for skipping Fermi vacuum state
*
#include "implicit.inc"
      DIMENSION VEC(*)
      DIMENSION LBLOCK(NBLOCK)
*     
      IOFF = 1
      DO IBLOCK = 1, NBLOCK
        IF( LBLOCK(IBLOCK).GT.0) THEN
          WRITE(6,*) ' Block : ', IBLOCK
          WRITE(6,*) ' ==================='
          WRITE(6,*)
          WRITE(6,*) ' Length : ', LBLOCK(IBLOCK)
          WRITE(6,*)
          CALL WRITVE(VEC(IOFF),LBLOCK(IBLOCK))
          IOFF = IOFF + LBLOCK(IBLOCK)
          WRITE(6,*)
        END IF
      END DO
*
      RETURN
      END

