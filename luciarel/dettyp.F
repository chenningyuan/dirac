      SUBROUTINE DETTYP(IPRNT)
*
* Allowed set of determinants : Distribution of
* electrons in alpha(up) -and beta(down) spinors
*
* Defined in terms of
* NMS2VAL : Number of MS2 values (MK2 values )
* MS2VAL  : The allowed MS2 values
*
* Jeppe Olsen, July 97
* Timo Fleig, August '97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
#include "mxpdim.inc"
#include "crun.inc"
#include "lucinp.inc"
*. Output
#include "cstate.inc"
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
*
      if (NTEST.ge.10) then
        write(6,*)
        write(6,*) ' -------------------------------- '
        write(6,*) '   DETTYP speaking  '
        write(6,*) ' -------------------------------- '
        write(6,*)
        write(6,'(A17,I4,I4)') 'MK2REF,MK2DEL are',MK2REF,MK2DEL
      end if
*
      IF(IDBGCI.EQ.0.AND.IRELACI.EQ.0) THEN
*. Single MS2 space
        NMS2VAL = 1
        MS2VAL(1) = MK2REF
      ELSE
*. Potentially several MS2/MK2 values
        MK2MAX = MK2REF + MK2DEL
*. MK2MAX gives highest MK2 value for possibly coupled strings (?)
        if (MK2MAX.gt.NACTEL) then
          MK2MAX = NACTEL
        else
          call evenodd(IEVODMX,MK2MAX)
          call evenodd(IEVODN,NACTEL)
          if (IEVODN.eq.1.and.IEVODMX.eq.2) MK2MAX = MK2MAX - 1
        end if
        MK2MIN = MK2REF - MK2DEL
        IF(MK2MIN.LT.-NACTEL) MK2MIN = -NACTEL
* Pairing ? Ensure symmetry between +MK and - MK spaces
        IPAIR = 1
        call evenodd(IEVOD,NACTEL)
        IF(IPAIR.EQ.1) THEN
          IF(MK2MIN.GT.0) THEN
            MK2MINP = -MK2MIN
          ELSE
            if (IEVOD.eq.1) then
               MK2MIN = 1
               MK2MINP = -1
            else if (IEVOD.eq.2) then
               MK2MIN = 0
               MK2MINP = -2
            else
               write(6,*) 'Problem in evenodd'
               Call Abend2( 'Quitting in dettyp.' )
            end if
          END IF
          MK2MAXP = -MK2MAX
        ELSE
          if (IEVOD.eq.1) then
             MK2MINP = -1
             MK2MAXP = 1
          else if (IEVOD.eq.2) then
             MK2MINP = -2
             MK2MAXP = 0
          else
             write(6,*) 'Problem in evenodd'
             Call Abend2( 'Quitting in dettyp.' )
          end if
        END IF
*. And construct the MS2 spaces
        NMS2VAL = 0
        DO MS2 = MK2MIN,MK2MAX,2
          NMS2VAL = NMS2VAL + 1
          MS2VAL(NMS2VAL) = MS2
        END DO
        DO MS2 = MK2MAXP,MK2MINP,2                           !sic
          NMS2VAL = NMS2VAL+1
          MS2VAL(NMS2VAL) = MS2
        END DO
      END IF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,'(A,I2)') ' Number of MS2/MK2 values', NMS2VAL
        WRITE(6,*)
        WRITE(6,'(A)') ' Included MS2/MK2 values '
        WRITE(6,'(A)') ' ======================= '
        DO IMS2 = 1, NMS2VAL
          WRITE(6,'(1H ,I3)') MS2VAL(IMS2)
        END DO
      END IF
*
      RETURN
      END
