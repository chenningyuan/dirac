      SUBROUTINE LEN_GENOP_STR_MAP(NGENOP,IGENOP,NSPGRP,ISPGRP,NOBPT,
     &                             NGAS,MAXLEN)
*
* A set of general operators, IGENOP, and a set of supergroups, ISPGRP, 
* are defined. Obtain max.length of mappings from resolution strings 
* to istrings. 
*
* The resolution strings KSTR are assumed to be inserted between the 
* creation and annihilation operators 
*
*   <ISPGRP!ICREA!KSTR> <KSTR!IANNI!ISPGRP'>
* So the mappings from KSTR to ISTR are always creation mappings 
*
* Jeppe Olsen, July 2000 ( At summerschool, HNIE)
*
#include "implicit.inc"
#include "mxpdim.inc"
*
*. Specific input
      INTEGER IGENOP(2*NGAS,*), ISPGRP(MXPNGAS,*), NOBPT(NGAS)
*. Local scratch
      INTEGER KSPGRP(MXPNGAS), IGENOP_EXP(MXPLCCOP)
*
      MAXLEN = 0
      IONE = 1
      IMONE = -1
      DO ICA = 1, 2
      DO JSPGRP = 1, NSPGRP
*. ICREA/IANNI(dag)*!KSTR> = !ISTR>
       DO JGENOP = 1, NGENOP
*. Occupation of KSPGRP
         IF(ICA.EQ.1) THEN
           CALL IVCSUM(KSPGRP,ISPGRP(1,JSPGRP),IGENOP(1,JGENOP),
     &                 IONE,IMONE,NGAS)
         ELSE
           CALL IVCSUM(KSPGRP,ISPGRP(1,JSPGRP),IGENOP(NGAS+1,JGENOP),
     &                 IONE,IMONE,NGAS)
         END IF
*. Is KSPGRP a correct supergroup( all occ larger than zero) 
         I_AM_OKAY = 1
         DO IGAS = 1, NGAS
           IF(KSPGRP(IGAS).LT.0) I_AM_OKAY = 0
         END DO
         IF(I_AM_OKAY.EQ.1) THEN
*. IGENOP in expanded form
C  REF_OP(IOPGAS,IOP,NOP,NGAS,IWAY)
          LEN = 0
          IF(ICA.EQ.1) THEN
            CALL REF_OP(IGENOP(1,JGENOP),IGENOP_EXP,NOP,
     &                  NGAS,1,MXPLCCOP)
          ELSE
            CALL REF_OP(IGENOP(NGAS+1,JGENOP),IGENOP_EXP,NOP,
     &                  NGAS,1,MXPLCCOP)
          END IF
          DO IOP = 1, NOP
            IOPTP = IGENOP_EXP(IOP)
*. Number of Kstrings
C                  NST_FOR_OCC(NELEC_PER_GAS,NORB_PER_GAS,NGAS)
            NSTR = NST_FOR_OCC(KSPGRP,NOBPT,NGAS)
            NORB = NOBPT(IOPTP)
            LEN = LEN + NSTR*NORB 
*. Update KSPGRP 
            KSPGRP(IOPTP) = KSPGRP(IOPTP) + 1
          END DO
          MAXLEN = MAX(MAXLEN,LEN)
         END IF
       END DO
      END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Max length of KSTR => ISTR map = ', MAXLEN
      END IF
*
      RETURN
      END
