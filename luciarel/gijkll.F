      FUNCTION GIJKLL(IREOTS,IPNTR,ISL,XINT,ISMFTO,IBSO,NACOB,NSMOB,
     &         NOCOBS,I,J,K,L)
*
* Obtain (IJ!KL), Lucas order
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IREOTS(*),IPNTR(NSMOB,NSMOB,NSMOB)
      DIMENSION ISL(NSMOB,NSMOB,NSMOB)
      DIMENSION IBSO(*),NOCOBS(*) ,ISMFTO(*)
      DIMENSION XINT(*)
*
C?    write(6,*) ' Hi from GIJKLL'
      II = IREOTS(I)
C?    write(6,*) ' II ',II
      ISM = ISMFTO(I)
C?    write(6,*) ' ISM ',ISM
      NI = NOCOBS(ISM)
C?    write(6,*) ' NI ',NI
      IREL = II - IBSO(ISM) + 1
C?    write(6,*) ' IREL ',IREL
      JJ = IREOTS(J)
      JSM = ISMFTO(J)
      JREL = JJ - IBSO(JSM) + 1
      NJ = NOCOBS(JSM)
      IJ = (IREL-1)*NJ + JREL
      JI = (JREL-1)*NI + IREL
      NJI = NI * NJ
      IJSM = (ISM-1)*NSMOB + JSM
*
      KK = IREOTS(K)
      KSM = ISMFTO(K)
      KREL = KK - IBSO(KSM) + 1
      NK = NOCOBS(KSM)
      LL = IREOTS(L)
      LSM = ISMFTO(L)
      LREL = LL - IBSO(LSM) + 1
      NL = NOCOBS(LSM)
      LK = (LREL-1)*NK + KREL
      KL = (KREL-1)*NL + LREL
      NLK = NK * NL
      KLSM = (KSM-1)*NSMOB + LSM
C?    WRITE(6,*) ' IJSM KLSM ', IJSM,KLSM
C?    WRITE(6,*) ' ISM JSM KSM LSM ',ISM,JSM,KSM,LSM

      IF(  (IJSM.GE.KLSM.AND.LSM.NE.ISL(ISM,JSM,KSM))
     &.OR. (IJSM.LT.KLSM.AND.JSM.NE.ISL(KSM,LSM,ISM)) )   THEN
        GIJKLL = 0.0D0
        write(6,*) ' Symmetry zero returned '
      ELSE
*
        IF(IJSM.GT.KLSM) THEN
C         IJKLO = (IJ-1)*NKL + KL + IPNTR(ISM,JSM,KSM)-1
          IJKLO = (LK-1)*NJI + JI + IPNTR(ISM,JSM,KSM)-1
        ELSE IF(IJSM.LT.KLSM) THEN
C         IJKLO = (KL-1)*NIJ + IJ + IPNTR(KSM,LSM,ISM)-1
          IJKLO = (JI-1)*NLK + LK + IPNTR(KSM,LSM,ISM)-1
        ELSE IF( IJSM.EQ.KLSM) THEN
C         IF(IJ.GE.KL) THEN
          IF(JI.GE.LK) THEN
C           IJKLO = IJ*(IJ-1)/2+KL + IPNTR(ISM,JSM,KSM)-1
            IJKLO = JI*(JI-1)/2+LK + IPNTR(ISM,JSM,KSM)-1
          ELSE
C           IJKLO = KL*(KL-1)/2+IJ + IPNTR(ISM,JSM,KSM)-1
            IJKLO = LK*(LK-1)/2+JI + IPNTR(ISM,JSM,KSM)-1
          END IF
        END IF
        GIJKLL = XINT(IJKLO)
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,5I4,3X,E18.12)')
     &  ' GIJKLL I J K L ,IJKLO,(IJ!KL) ', I,J,K,L,IJKLO,GIJKLL
      END IF
*
      RETURN
      END
