      SUBROUTINE IAIBCM_REL(ICISPC,IATP_TF,IBTP_TF,IAIB)
*
* obtain allowed combination of alpha- and beta- supergroups
* for CI space ICISPC
*
* Master for IAIBCM_GAS
*
*      Jeppe Olsen, august 1995
*
*  Corrected for rel. GASCI
*       Timo Fleig, 1999
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
#include "stinf.inc"
#include "strinp.inc"
*. Output
      INTEGER IAIB(*)
*
      NAEL = NELEC(IATP_TF)
      NBEL = NELEC(IBTP_TF)
*
      NOCTPA = NOCTYP(IATP_TF)
      NOCTPB = NOCTYP(IBTP_TF)
*
      IOCTPA = IBSPGPFTP(IATP_TF)
      IOCTPB = IBSPGPFTP(IBTP_TF)
*
C?    write(6,*) ' IAIB ::::::'
C?    write(6,*) ' LCMBSPC, ICISPC, ICMBSPC '
C?    WRITE(6,*) ICISPC,  LCMBSPC(ICISPC)
C?    WRITE(6,*) (ICMBSPC(II,ICISPC),II=1, LCMBSPC(ICISPC))

      CALL IAIBCM_GAS(LCMBSPC(ICISPC),ICMBSPC(1,ICISPC),
     &            IGSOCCX,NOCTPA,
     &            NOCTPB,ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &            NELFGP,MXPNGAS,NGAS,IAIB,IPRDIA)
*
      RETURN
      END
