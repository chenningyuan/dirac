      SUBROUTINE TS_SYM_PNT2(IGRP,NIGRP,MAXVAL,MINVAL,ISYM,
     &                       IPNT,LPNT,IAB)
*
* Construct pointers to start of symmetry distributions
* for supergroup of strings with given symmetry
*
* The start of symmetry block ISYM1 ISYM2 ISYM3 .... ISYMN
* is given as
*     1
*     + (ISM1-MINVAL(1))
*     + (ISM2-MINVAL(2))*(MAXVAL(1)-MINVAL(1)+1)
*     + (ISM3-MINVAL(3))*(MAXVAL(1)-MINVAL(1)+1)*(MAXVAL(2)-MINVAL(2)+1)
*     +
*     +
*     +
*     + (ISM L-1-MINVAL(L-1))*Prod(i=1,L-2)(MAXVAL(i)-MINVAL(i)+1)
*
* Where L is the last group of strings with nonvanishing occupation
*
* Jeppe Olsen, September 1997
* Same kind of bug fix inserted as in getstr_totsm_spgp. (July 2001)
* IAB inserted for relativistic calculations. (August 2001)
*   Timo Fleig
*
* Version 2 : Uses IGRP and NIGRP to define supergroup
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "csm.inc"
*. Specific Input
      INTEGER IGRP(NIGRP)
*. Local scratch
      INTEGER NELFGS(MXPNGAS),ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
*. Output
      INTEGER MINVAL(*),MAXVAL(*),IPNT(*)
*
      NTEST = 000
*. Info on groups of strings in supergroup
      NGASL = 1
      DO IGAS = 1, NIGRP
       ITPFGS(IGAS) = IGRP(IGAS)
       NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
       IF(NELFGS(IGAS).GT.0) NGASL = IGAS
*. Number of strings per symmetry in each gasspace
        if (IAB.eq.1) then
          CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                 NNSTSGP(1,IGAS))
        else if (IAB.eq.2) then
          CALL ICOPVE2(WORK(KNSTSGP2(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &                 NNSTSGP(1,IGAS))
        end if
      END DO
*
      IZERO = 0
      DO IGAS = 1, NIGRP
        MINVAL(IGAS) = IZERO
        MAXVAL(IGAS) = IZERO
        DO ISMST = 1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
CTF
* If no strings are found in a given GAS, we set the symmetry
* value in MINVAL to 1, just to prevent crashes (only a test,
* waiting for more elegant solution ...)
      IONE = 1
      do IGAS = 1,NGAS-1,1
        if (MINVAL(IGAS).eq.0) then
          MINVAL(IGAS) = IONE
        end if
      end do
CTF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)  ' MINVAL and MAXVAL '
        CALL IWRTMA(MINVAL,1,NIGRP,1,NIGRP)
        CALL IWRTMA(MAXVAL,1,NIGRP,1,NIGRP)
        WRITE(6,*) ' NIGRP = ', NIGRP
      END IF

*. Total number of strings that will be generated
      NBLKS = 1
      DO IGAS = 1, NGASL-1
       NBLKS = NBLKS*(MAXVAL(IGAS)-MINVAL(IGAS)+1)
      END DO
      IF(NBLKS.GT.LPNT) THEN
        WRITE(6,*) ' Problem in TS_SYM_PNT'
        WRITE(6,*) ' Dimension of IPNT too small'
        WRITE(6,*) ' Actual and required length',NBLKS,LPNT
        WRITE(6,*)
        WRITE(6,*) ' I will Stop and wait for instructions'
        Call Abend2( ' TS_SYM_PNT too small' )
      END IF
*. Loop over symmetry blocks in standard order
      IFIRST = 1
      ISTRBS = 1
      NSTRINT = 0
 2000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 2001
        END IF
        IFIRST = 0
*. Symmetry of NGASL -1 spaces given, symmetry of full space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*.  sym of SPACE NGASL
        CALL SYMCOM(2,1,ISTSMM1,ISMGSN,ISYM)
        ISMFGS(NGASL) = ISMGSN
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' next symmetry of NGASL spaces '
          CALL IWRTMA(ISMFGS,1,NGASL,1,NGASL)
        END IF
*. Number of strings with this symmetry combination
        NSTRII = 1
        DO IGAS = 1, NGASL
          NSTRII = NSTRII*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for this symmetry distribution in IOFFI
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGASL-1
          IOFF = IOFF + (ISMFGS(IGAS)-MINVAL(IGAS))*MULT
          MULT = MULT * (MAXVAL(IGAS)-MINVAL(IGAS)+1)
        END DO
*
        IPNT(IOFF) = NSTRINT + 1
        NSTRINT = NSTRINT + NSTRII
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IOFF, IPNT(IOFF) NSTRII ',
     &                 IOFF, IPNT(IOFF),NSTRII
        END IF
*
      IF(NGASL-1.GT.0) GOTO 2000
 2001 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from TS_SYM_PNT'
        WRITE(6,*) ' Required total symmetry',ISYM
        WRITE(6,*) ' Number of symmetry blocks ', NBLKS
        WRITE(6,*)
        WRITE(6,*) ' Offset array  for symmetry blocks'
        CALL IWRTMA(IPNT,1,NBLKS,1,NBLKS)
      END IF
*
      RETURN
      END
