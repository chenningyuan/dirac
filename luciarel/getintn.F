      SUBROUTINE GETINTN(XINT,ITP,ISMGET,JTP,JSMGET,KTP,KSMGET,
     &                   LTP,LSMGET,
     &                   IXCHNG,IKSMGET,JLSMGET,ICOUL,IRI)
*
* Interface to DIRAC relativistic integrals
* (in spinfree and full spin-dependent formalism) implemented
*
*---------------------------------------------
* Timo Fleig, October 2000
*---------------------------------------------
*
* Outer routine for accessing integral block
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "cgas.inc"
#include "glbbas.inc"
#include "symm.inc"
*
      CALL QENTER('GETIN')
      NTEST = 000
*
      if (INTIMP.eq.1.or.INTIMP.eq.5.or.INTIMP.eq.6) then
*. Interface to MOTRA or SIRIUS integrals
*  or to DIRAC (spinfree) integrals
        if (INCORE.eq.0) then
          if (INTIMP.eq.6) then
            write(6,*) 'Out of core integral handling'
            write(6,*) 'not implemented yet.'
            call abend('Quitting.')
          end if
          call getinm(XINT,ITP,ISMGET,JTP,JSMGET,KTP,KSMGET,LTP,LSMGET,
     &                IXCHNG,IKSMGET,JLSMGET,SYDI,NIRR_PN)
        else if (INCORE.ne.0) then
* Integrals are held in core
          if (INTIMP.eq.6) then
* Read two-electron integrals from DIRAC file MDCINT
            call getind(XINT,ITP,ISMGET,JTP,JSMGET,KTP,KSMGET,
     &                  LTP,LSMGET,
     &                  WORK(KINVREO),WORK(KKR),WORK(KNUMINT),
     &                  WORK(KINDK),WORK(KINDL),WORK(KIOFF12),
     &                  WORK(KDITWOE+(IRI-1)*IMOFFSET))
          else
*. Interface to MOTRA or SIRIUS integrals
            call getincn(XINT,ITP,ISMGET,JTP,JSMGET,KTP,KSMGET,
     &                   LTP,LSMGET,
     &                   IXCHNG,IKSMGET,JLSMGET,WORK(KINT2),
     &                   WORK(KPINT2),NSMSH,WORK(KINH1),ICOUL)
          end if
        end if
*
      else if (INTIMP.EQ.2) THEN
C. Interface to LUCAS integrals
        CALL GETINL(XINT,ITP,ISMGET,JTP,JSMGET,KTP,KSMGET,LTP,LSMGET,
     &                  IXCHNG,IKSMGET,JLSMGET)
*
      else if (INTIMP.eq.8) then
*
* Molecule Sweden interface
*
        call GETINCN(XINT,ITP,ISMGET,JTP,JSMGET,KTP,KSMGET,LTP,LSMGET,
     &               IXCHNG,IKSMGET,JLSMGET,WORK(KINT2),
     &               WORK(KPINT2),NSMSH,WORK(KINH1),ICOUL)
      END IF
*
      IF(NTEST.NE.0) THEN
        NI = NGSSH(ISMGET,ITP)
        NK = NGSSH(KSMGET,KTP)
        IF(IKSMGET.EQ.0) THEN
          NIK = NI * NK
        ELSE
          NIK = NI*(NI+1)/2
        END IF
        NJ = NGSSH(JSMGET,JTP)
        NL = NGSSH(LSMGET,LTP)
        IF(JLSMGET.EQ.0) THEN
          NJL = NJ * NL
        ELSE
          NJL = NJ*(NJ+1)/2
        END IF
        WRITE(6,*) ' 2 electron integral block for TS blocks '
        WRITE(6,*) ' Ixchng :', IXCHNG
        WRITE(6,*) ' After GETINC '
        WRITE(6,'(1H ,4(A,I2,A,I2,A))')
     &  '(',ITP,',',ISMGET,')','(',JTP,',',JSMGET,')',
     &  '(',KTP,',',KSMGET,')','(',LTP,',',LSMGET,')'
        CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
      END IF
*
      CALL QEXIT('GETIN')
*
      RETURN
      END
