      SUBROUTINE EXPCIVS_REL(LUIN,VECIN,NCOMBIN,NTTSIN,ITTSIN,IABIN,
     &                       NOCTPA,NOCTPB,IBLTPIN,
     &                       LUUT,VECUT,NCOMBUT,NTTSUT,ITTSUT,IABUT,
     &                       IBLTPUT,
     &                       ICISTR,IDC,NSMST,LBLK,ITTSS_ORD)
*
* Obtain those part of vector in cispace UT ,
* that can be obtained from terms in cispace IN
*
* Input vector on LUIN, Output vector in LUUT
* Output vector is supposed on start of vector
*
* LUIN is assumed to be single vector file,
* so rewinding will place vector on start of vector
*
* Both files are assumed on start of vector
*
* Jeppe Olsen, February 1994
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
      INTEGER IABIN(NOCTPA,NOCTPB),IABUT(NOCTPA,NOCTPB)
      INTEGER IBLTPIN(NSMST),IBLTPUT(NSMST)
      INTEGER NTTSIN(NOCTPA,NOCTPB,NSMST),ITTSIN(NOCTPA,NOCTPB,NSMST)
      INTEGER NTTSUT(NOCTPA,NOCTPB,NSMST),ITTSUT(NOCTPA,NOCTPB,NSMST)
      DIMENSION VECIN(*)
C     DIMENSION IATPUI(*),IBTPUI(*)
*. Output
      DIMENSION VECUT(*)
*
      WRITE(6,*) ' EXPCIVS in action'
      WRITE(6,*) '        Number of input  parameters',NCOMBIN
      WRITE(6,*) '        Number of output parameters',NCOMBUT
*
      IF(ICISTR.EQ.1) THEN
        CALL FRMDSC(VECIN,NCOMBIN,-1,LUIN,IMZERO,IAMPACK)
C       WRITE(6,*) ' Input vector '
C       CALL WRTMAT(VECIN,1,NCOMBIN,1,NCOMBIN)
      END IF
*
*. Loop over TTS blocks of output vector
*
      IATPIN = 1
      IBTPIN = 0
      IASMIN = 1
*
      IATPUT = 1
      IBTPUT = 0
      IASMUT = 1
*
 1000 CONTINUE
*. Next output block
        CALL NXTBLK(IATPUT,IBTPUT,IASMUT,NOCTPA,NOCTPB,NSMST,
     &              IBLTPUT,IDC,NONEWUT,IABUT,ITTSS_ORD)
*. Corresponding input TTS block
        JATPIN = IATPUT
        JBTPIN = IBTPUT
        IF(IABIN(JATPIN,JBTPIN).EQ.0) THEN
          IZERO = 1
        ELSE
          IZERO = 0
        END IF
*
        NELMNT = NTTSUT(IATPUT,IBTPUT,IASMUT)
        IF(ICISTR.EQ.1) THEN
          IF(IZERO.EQ.0)
     &    IOFFIN = ITTSIN(JATPIN,JBTPIN,IASMUT)
          IOFFUT = ITTSUT(IATPUT,IBTPUT,IASMUT)
          IF(IZERO.EQ.0) THEN
            CALL COPVEC(VECIN(IOFFIN),VECUT(IOFFUT),NELMNT)
          ELSE
            ZERO = 0.0D0
            CALL SETVEC(VECUT(IOFFUT),ZERO,NELMNT)
          END IF
        END IF
*
        IF(ICISTR.NE.1.AND.NONEWUT.EQ.0) THEN
          IF(IZERO.EQ.1) THEN
C           ZERO = 0.0D0
C           CALL SETVEC(VECUT,ZERO,NELMNT)
            CALL ITODS(NELMNT,1,-1,LUUT)
            CALL ZERORC_REL(-1,LUUT)
C           CALL TODSC(VECUT,NELMNT,-1,LUUT)
          ELSE
*. Obtain input block
  999      CONTINUE
           CALL NXTBLK(IATPIN,IBTPIN,IASMIN,NOCTPA,NOCTPB,NSMST,
     &                 IBLTPIN,IDC,NONEWIN,IABIN,ITTSS_ORD)
           IF(NONEWIN.NE.0) THEN
             CALL REWINO(LUIN)
             IATPIN = 1
             IBTPIN = 0
             IASMIN = 1
             GOTO 999
           ELSE IF (NONEWIN.EQ. 0 ) THEN
             CALL IFRMDS(LENGTH,1,-1,LUIN)
             CALL FRMDSC(VECIN,LENGTH,-1,LUIN,IMZERO,IAMPACK)
             IF(IATPIN.EQ.JATPIN.AND.IBTPIN.EQ.JBTPIN.AND.
     &          IASMIN.EQ.IASMUT) THEN
*. Correct block, save it
                  CALL ITODS(LENGTH,1,-1,LUUT)
                  IF(IMZERO.EQ.1) THEN
                    CALL ZERORC_REL(-1,LUUT)
                  ELSE
                    IF(IAMPACK.EQ.0) THEN
                      CALL TODSC(VECIN,LENGTH,-1,LUUT)
                    ELSE
                      CALL TODSCP(VECIN,LENGTH,-1,LUUT)
                    END IF
                  END IF
             ELSE
               GOTO 999
             END IF
           END IF
          END IF
        END IF
      IF(NONEWUT.EQ.0) GOTO 1000
*. End of loop over output blocks
      IF(ICISTR.EQ.1) THEN
        CALL TODSC(VECUT,NCOMBUT,-1,LUUT)
      ELSE
        CALL ITODS(-1,1,-1,LUUT)
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' EXPCIVS Speaking '
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' ============ '
        WRITE(6,*) ' Input Vector '
        WRITE(6,*) ' ============ '
        WRITE(6,*)
        IF(ICISTR.EQ.1) THEN
          CALL WRTMAT(VECIN,1,NCOMBIN,1,NCOMBIN)
        ELSE
          CALL WRTVCD(VECIN,LUIN,1,LBLK)
        END IF
        WRITE(6,*)
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' Output Vector '
        WRITE(6,*) ' =============== '
        WRITE(6,*)
        IF(ICISTR.EQ.1) THEN
          CALL WRTMAT(VECUT,1,NCOMBUT,1,NCOMBUT)
        ELSE
          CALL WRTVCD(VECUT,LUUT,1,LBLK)
        END IF
      END IF
*
      RETURN
      END
