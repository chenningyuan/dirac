      SUBROUTINE MAP_EXSTR(IOP,IAC,NOP,IREFOC,IX,SX,NK,IBX,SCLFAC,IAB)

* Information for mapping |Kstr> = IOP |Irefoc>
*
* Mapping for each creation/annihilation operator is given
*
*     Jeppe Olsen, Summer of 99
*     Timo Fleig, summer of 2001, mod. for rel. case
*                 IAB added for alpha (unbarred) and beta (barred)
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "glbbas.inc"
#include "cgas.inc"
#include "csm.inc"
#include "symm.inc"
#include "ctcc.inc"
*. Input
      INTEGER IOP(NGAS),IAC(NGAS),IREFOC(NGAS)
*. Output
      INTEGER IX(*),IBX(MXNDGIRR,MXNDGIRR,NOP)
      INTEGER NK(MXPOBS,NOP)
      DIMENSION SX(*)
*. Local Scratch
      INTEGER KOC(MXPNGAS),KGRP(MXPNGAS)
*. IBX(ISM,JSM,JOP) : Start of map in JX for orbitals of symmetry ISM,
*.                    input strings of symmetry JSM for operator JOP
*. Notice at each level N the mapping is from !K(N)> to !K(N-1)> 
*  where !K(0)> is input type and !K(NOP)> is output type
*
C     MAP_EXSTR(IOP,IAC,NOP,IREFOC,IX,SX,NK,IB,SCLFAC)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' MAP_EXSTR speaking' 
         WRITE(6,*) ' =================='
        WRITE(6,*) ' Reference occupation '
        CALL IWRTMA(IREFOC,1,NGAS,1,NGAS) 
        WRITE(6,*) ' NOP = ', NOP
        WRITE(6,*) ' Operator string, Orbspaces and C/A '
        CALL IWRTMA(IOP,1,NOP,1,NOP)
        CALL IWRTMA(IAC,1,NOP,1,NOP)
      END IF
*
      CALL ICOPVE(IREFOC,KOC,NGAS)
      CALL OCC_TO_GRP(KOC,KGRP,1)
      IOFF = 1
      IF(NTEST.GE.1000) WRITE(6,*) ' NOP = ', NOP
      DO JOP = NOP,1,-1
       IF(NTEST.GE.1000) WRITE(6,*) ' JOP = ', JOP
*
       JGAS = IOP(JOP)
       JAC  = IAC(JOP)
*
       IF(JAC.EQ.1) THEN 
         JAC_ADJ = 2 
       ELSE 
         JAC_ADJ = 1 
         WRITE(6,*)' Map <I|a|K> not debugged for DBG symmetry'
         STOP' Map <I|a|K> not debugged for DBG symmetry'
       END IF
*. 
COLD
CTF Depending on the evaluation "direction" of the coupling coefficient,
CTF the symmetry handling changes: If |K> is on the right, the inverse
CTF element is correct (INVELM(JSMOB)), if |J> is on the right, the
CTF element is NOT inversed!
CTF
       DO JSMOB = 1, NSMOB
        DO JSMSPGP = 1, NSMST
*. Symmmetry of |Kgrp> = a(JSMOB) |JSMSPGP>
*. = symmetry obtained from |JSMSPGP> = a+(JSMOB) |Kgrp>
         KSMSPGP = IDBGMULT(INVELM(JSMOB),JSMSPGP)
*
         IBX(JSMOB,KSMSPGP,JOP) = IOFF
         CALL ADAST_GAS(JSMOB,JGAS,NGAS,KGRP,JSMSPGP,
     &                  IX(IOFF),SX(IOFF),NKSTR,IEND,
     &                  IFRST,KFRST,KACT,SCLFAC,JAC_ADJ,IAB)
C             ADAST_GAS(IOBSM,IOBTP,NIGRP,IGRP,ISPGPSM,
C    &                  I1,XI1S,NKSTR,IEND,
C    &                  IFRST,KFRST,KACT,SCLFAC,IAC,IAB)
         NK(KSMSPGP,JOP) = NKSTR
         if (IAB.eq.1) then
           IOFF = IOFF + NOBPTS(JGAS,JSMOB)*NKSTR
         else if (IAB.eq.2) then
           IOFF = IOFF + NOBPTS2(JGAS,JSMOB)*NKSTR
         end if
         IF(IOFF-1.GT. MAXLEN_I1) THEN
           WRITE(6,*) ' MAP... MAXLEN_I1 too small '
           WRITE(6,*) ' IOFF and MAXLEN_I1 = ', IOFF, MAXLEN_I1
           STOP       ' MAP... MAXLEN_I1 too small '
         END IF
*
        END DO
       END DO
*. Updated Kgroup  
*. Info on !K(JOP)> 
       IF(JAC.EQ.1) THEN 
         JDELTA = -1 
       ELSE 
         JDELTA = +1
       END IF
       KOC(JGAS) = KOC(JGAS) + JDELTA
       CALL OCC_TO_GRP(KOC,KGRP,1)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' KOC and KGRP '
         CALL IWRTMA(KOC,1,NGAS,1,NGAS)
         CALL IWRTMA(KGRP,1,NGAS,1,NGAS)
       END IF
*
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from MAP_EXSTR '
        WRITE(6,*) ' Mappings from a+/a !K(J)> to !K(J-1)> '
*
        DO JOP = 1, NOP
         IF(IAC(JOP).EQ.1) THEN
           WRITE(6,*) 
     &     ' Operator ', JOP, ' is annihilation operator of type',
     &     IOP(NOP) 
         ELSE
           WRITE(6,*) 
     &     ' Operator ', JOP, ' is creation     operator of type',
     &     IOP(NOP) 
         END IF
         DO JSM = 1, NSMOB
           DO KNSM = 1, NSMST
             WRITE(6,*) ' Sym of operator and !K(N)> ',
     &       JSM,KNSM 
             LK = NK(KNSM,JOP)
             WRITE(6,*) ' Number of strings !K(N)> ',LK
             if (IAB.eq.1) then
               LJ = NOBPTS(IOP(JOP),JSM) 
             else if (IAB.eq.2) then
               LJ = NOBPTS2(IOP(JOP),JSM) 
             end if
             DO JJ = 1, LJ
              WRITE(6,*) ' Info for orbital ', JJ
              WRITE(6,*) ' Excited strings and sign '
              IOFF2 = IBX(JSM,KNSM,JOP) + (JJ-1)*LK
              WRITE(6,*) ' JOP JSM KNSM LK ',JOP,JSM,KNSM,LK
              WRITE(6,*) ' IOFF2 ', IOFF2
              CALL IWRTMA(IX(IOFF2),1,LK,1,LK)        
              CALL WRTMAT(SX(IOFF2),1,LK,1,LK)        
             END DO
           END DO
         END DO
        END DO
      END IF
*
      RETURN
      END 
