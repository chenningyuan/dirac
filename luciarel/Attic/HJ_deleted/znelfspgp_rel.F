      SUBROUTINE ZNELFSPGP_REL(NTESTG)
*
* Generate for each supergroup the number of electrons in each active
* orbital space and store in NELFSPGP
*
* Jeppe Olsen, July 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "cgas.inc"
*. Input and Output ( NELFSPGP(MXPNGAS,MXPSTT) )
#include "gasstr.inc"
*
      NTESTL = 0
      NTEST = MAX(NTESTG,NTESTL)
*
      DO ITP = 1, NSTTP
        NSPGP = NSPGPFTP(ITP)
        IBSPGP = IBSPGPFTP(ITP)
        DO ISPGP = IBSPGP,IBSPGP + NSPGP - 1
          DO IGAS = 1, NGAS
            NELFSPGP(IGAS,ISPGP) = NELFGP(ISPGPFTP(IGAS,ISPGP))
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Distribution of electrons in Active spaces '
        DO ITP = 1, NSTTP
          WRITE(6,*) ' String type ', ITP
          WRITE(6,*) ' Row : active space, Column: supergroup '
          NSPGP = NSPGPFTP(ITP)
          IBSPGP = IBSPGPFTP(ITP)
          CALL IWRTMA(NELFSPGP(1,IBSPGP),NGAS,NSPGP,MXPNGAS,NSPGP)
        END DO
      END IF
*
      RETURN
      END
