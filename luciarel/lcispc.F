      SUBROUTINE LCISPC(IPRNT)
*
* Number of dets and combinations
* per symmetry for each type of internal space
*
* Jeppe Olsen , Winter 1994/1995 ( woops !)
*               July 1997, multiple Ms spaces allowed
*
* GAS VERSION
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* ===================
*.Input common blocks
* ===================
*
*./LUCINP/
#include "mxpdim.inc"
#include "lucinp.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "strbas.inc"
*./CSM/ : NSMST is used
#include "csm.inc"
#include "stinf.inc"
*./CSTATE/ : IDC is used
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
*
* ====================
*. Output common block : XISPSM is calculated
* ====================
*
#include "cicisp.inc"
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
*
      CALL QENTER('LCISP')
*. Number of spaces
      NICISP = NCMBSPC
*
*.Local memory
*
*.  Largest NOCTPA*NOCTPB
      MXNOCAB = 0
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
        NOCAB = NOCTYP(IATP)*NOCTYP(IBTP)
        MXNOCAB = MAX(MXNOCAB,NOCAB)
      END DO
C     WRITE(6,*) ' MXNOCAB = ',MXNOCAB
      CALL MEMMAN(KLIOIO,MXNOCAB,   'ADDL  ',2,'KLIOIO')
*
      CALL MEMMAN(KLBLTP,NSMST,'ADDL  ',2,'KLBLTP')
*. Obtain array giving symmetry of sigma v reflection times string
*. symmetry.
      IF(IDC.EQ.3 .OR. IDC .EQ. 4 )
     &CALL MEMMAN(KLCVST,NSMST,'ADDL  ',2,'KLCVST')
*. Array defining symmetry combinations of internal strings
*. Number of internal dets for each symmetry
      CALL SMOST(NSMST,NSMCI,MXPCSM,ISMOST)
*
      MXSB = 0
      MXSOOB = 0
*. Loop over MS2 spaces
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
*
        NOCTPA =  NOCTYP(IATP)
        NOCTPB =  NOCTYP(IBTP)
*
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
*
C?      WRITE(6,*) ' IMS2 IATP IBTP NOCTPA NOCTPB '
C?      WRITE(6,*)   IMS2,IATP,IBTP,NOCTPA,NOCTPB
*
        DO ICI = 1, NICISP
*. allowed combination of types
          CALL IAIBCM_GAS(LCMBSPC(ICI),ICMBSPC(1,ICI),
     &         IGSOCCX,NOCTPA,
     &         NOCTPB,ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &         NELFGP,MXPNGAS,NGAS,WORK(KLIOIO),IPRDIA)
CFEB19
          KLCVST = 1
          DO  ISYM = 1, NSMCI
            CALL ZBLTP(ISMOST(1,ISYM),NSMST,IDC,WORK(KLBLTP))
            CALL NGASDT(IGSOCCX(1,1,ICI),IGSOCCX(1,2,ICI),NGAS,ISYM,
     &           NSMST,NOCTPA,NOCTPB,WORK(KNSTSO(IATP)),
     &           WORK(KNSTSO2(IBTP)),
     &           ISPGPFTP(1,IBSPGPFTP(IATP)),
     &           ISPGPFTP(1,IBSPGPFTP(IBTP)),MXPNGAS,NELFGP,
     &           NCOMB,XNCOMB,MXS,MXSOO,WORK(KLBLTP),NTTSBL,
     &           LCOL,WORK(KLIOIO))

            IF(IMS2.EQ.1) THEN
              XISPSM(ISYM,ICI) = XNCOMB
              NBLKIC(ISYM,ICI) = NTTSBL
              LCOLIC(ISYM,ICI) = LCOL
            ELSE
              XISPSM(ISYM,ICI) = XISPSM(ISYM,ICI) + XNCOMB
              NBLKIC(ISYM,ICI) = NBLKIC(ISYM,ICI) + NTTSBL
              LCOLIC(ISYM,ICI) = LCOLIC(ISYM,ICI) + LCOL
            END IF
            MXSOOB = MAX(MXSOOB,MXSOO)
            MXSB = MAX(MXSB,MXS)
          END  DO
*         ^ End of loop over symmetries
        END DO
*       ^ End of loop over CI spaces
      END DO
*     ^ End of loop over MS2 spaces
*
      WRITE(6,*)
     &' Number of internal combinations per symmetry '
      WRITE(6,*)
     & ' =========================================== '
*
      DO 200 ICI = 1, NCMBSPC
          WRITE(6,*) ' CI space ', ICI
          CALL WRTMAT(XISPSM(1,ICI),1,NSMCI,1,NSMCI)
  200 CONTINUE
      WRITE(6,*) ' Largest symmetry block           ',MXSB
      WRITE(6,*) ' Largest Symmetry-type-type block ',MXSOOB
*
      WRITE(6,*)
     &' Number of TTS subblocks per CI expansion '
      WRITE(6,*)
     & ' ======================================== '
*
      DO  ICI = 1,  NCMBSPC
          WRITE(6,*) ' Internal CI space ', ICI
          CALL IWRTMA(NBLKIC(1,ICI),1,NSMCI,1,NSMCI)
      END DO
*. Largest number of BLOCKS in a CI expansion
      MXNTTS = 0
      DO ICI = 1,NCMBSPC
       DO ISM =1, NSMCI
        MXNTTS = MAX(MXNTTS,NBLKIC(ISM,ICI))
       END DO
      END DO
*
      WRITE(6,*) ' Largest number of blocks in CI expansion',
     &MXNTTS
*
      IF(NTEST.GE.5) THEN
      WRITE(6,*)
     &' Number of columns per CI expansion '
      WRITE(6,*)
     & ' =================================== '
*
      DO  ICI = 1,  NCMBSPC
          WRITE(6,*) ' Internal CI space ', ICI
          CALL IWRTMA(LCOLIC(1,ICI),1,NSMCI,1,NSMCI)
      END DO
      END IF
*
      CALL QEXIT('LCISP')
*
      RETURN
      END
