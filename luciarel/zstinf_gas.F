      SUBROUTINE ZSTINF_GAS(IPRNT)
*
* Set up common block /STINF/ from information in /STINP/
*
*=========
* Input
*=========
* Information in /CGAS/ and /GASSTR/
*
*======================
* Output ( in /STINF/ )
*======================
* ISTAC (MXPSTT,2) : string type obtained by creating (ISTAC(ITYP,2))
*                    or annihilating (ISTAC(ITYP,1)) an electron
*                    from a string of type  ITYP . A zero indicates
*                    that this mapping is not included
*                    Only strings belonging to the same
*                    Orbital group are mapped
*                    mapped
*. Input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
C     COMMON/CGAS/IDOGAS,NGAS,NGSSH(MXPIRR,MXPNGAS),
C    &            NGSOB(MXPOBS,MXPNGAS),
C    &            NGSOBT(MXPNGAS),IGSOCC(2,MXPNGAS),IGSINA,IGSDEL
C     COMMON/GASSTR/MNGSOC(MXPNGAS),MXGSOC(MXPNGAS),NGPSTR(MXPNGAS),
C    &              IBGPSTR(MXPNGAS),NELFGP(MXPSTT),IGSFGP(MXPSTT),
C    &              NSTFGP(MXPSTT),MNELFGP(MXPNGAS),MXELFGP(MXPNGAS)
*. Output
#include "stinf.inc"
*./STINF/
C     COMMON/STINF/ISTAC(MXPSTT,2),NOCTYP(MXPSTT),NSTFTP(MXPSTT),
C    &             INUMAP(MXPSTT),INDMAP(MXPSTT)
*. Only the first element, i.e. ISTAC  is defined

*
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
* ******************************************************************
* Mappings between strings with the same type ISTTP index , +/- 1 el
* ******************************************************************
      CALL ISETVC(ISTAC,0,2*MXPSTT)
      DO  IGAS = 1, NGAS
*. groups for a given gas spaces goes with increasing number of orbitals,
*  so the first space does not have any creation mapping
*  and the last space does not have any annihilation mapping
*
        MGRP = NGPSTR(IGAS)
        DO IGRP = 1, MGRP
          IIGRP = IGRP + IBGPSTR(IGAS) -1
          IF(IGRP.NE.1) THEN
*. Annihilation map is present : IIGRP => IIGRP - 1
            ISTAC(IIGRP,1) = IIGRP -1
          END IF
          IF(IGRP.NE.MGRP) THEN
*. Creation map is present : IIGRP => IIGRP + 1
             ISTAC(IIGRP,2) = IIGRP + 1
          END IF
        END DO
      END DO
*
      IF(NTEST .GE. 10 ) THEN
        WRITE(6,*) ' Type - type mapping array ISTAC '
        WRITE(6,*) ' =============================== '
        CALL IWRTMA(ISTAC,NGRP  ,2,MXPSTT,2)
      END IF
*
      RETURN
      END
