      SUBROUTINE IAIBCM_GAS(LCMBSPC,ICMBSPC,
     &                      MNMXOC,NOCTPA,NOCTPB,IOCA,IOCB,NELFTP,
     &                      MXPNGAS,NGAS,IOCOC,IPRNT)
*
* Allowed combinations of alpha and beta types, GAS version
*
*
* =====
*.Input
* =====
*
* LCMBSPC : Number of GAS spaces included in this expnasion
* ICMBSPC : Gas spaces included in this expansion
*
* MXMNOC(IGAS,1,IGASSPC) : Min accumulated occ for AS 1-IGAS for space IGASSPC
* MXMNOC(IGAS,2,IGASSPC) : Max accumulated occ for AS 1-IGAS for space IGASSPC
*
* NOCTPA : Number of alpha types
* NOCTPB : Number of beta types
*
* IOCA(IGAS,ISTR) occupation of AS IGAS for alpha string type ISTR
* IOCB(IGAS,ISTR) occupation of AS IGAS for beta  string type ISTR
*
* MXPNGAS : Largest allowed number of gas spaces
* NGAS    : Actual number of gas spaces

*
* ======
*.Output
* ======
*
* IOCOC(IATP,IBTP)  = 1 =>      allowed combination
* IOCOC(IATP,IBTP)  = 0 => not allowed combination
*
*.Input
      INTEGER ICMBSPC(LCMBSPC)
      INTEGER MNMXOC(MXPNGAS,2,*)
C     INTEGER MNOCC(NGAS),MXOCC(NGAS)
      INTEGER IOCA(MXPNGAS,NOCTPA),IOCB(MXPNGAS,NOCTPB)
      INTEGER NELFTP(*)
*.Output
      INTEGER IOCOC(NOCTPA,NOCTPB)
*
      NTEST = 00
C     NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' IAICBM_GAS entered '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        WRITE(6,*) ' Number of GAS spaces included ', LCMBSPC
        WRITE(6,*) ' GAS spaces included ',(ICMBSPC(II),II=1,LCMBSPC)
        WRITE(6,*)
        IF(NTEST.GE.20) THEN
          WRITE(6,*) ' IOCA and IOCB '
          CALL IWRTMA(IOCA,NGAS,NOCTPA,MXPNGAS,NOCTPA)
          CALL IWRTMA(IOCB,NGAS,NOCTPB,MXPNGAS,NOCTPB)
        END IF
      END IF
*
      CALL ISETVC(IOCOC,0,NOCTPA*NOCTPB)
      DO 100 IATP = 1, NOCTPA
         DO 90 IBTP = 1, NOCTPB
*. is this combination allowed in any of the GAS spaces included
           INCLUDE = 0
           DO JJCMBSPC = 1, LCMBSPC
             JCMBSPC = ICMBSPC(JJCMBSPC)
             IEL = 0
             IAMOKAY = 1
             DO IGAS = 1, NGAS
               IEL = IEL
     &             + NELFTP(IOCA(IGAS,IATP))+NELFTP(IOCB(IGAS,IBTP))
               IF(IEL.LT.MNMXOC(IGAS,1,JCMBSPC).OR.
     &            IEL.GT.MNMXOC(IGAS,2,JCMBSPC))
     &         IAMOKAY = 0
             END DO
             IF(IAMOKAY.EQ.1) INCLUDE = 1
           END DO
           IF(INCLUDE .EQ.1) THEN
*. Congratulations , you are allowed
              IOCOC(IATP,IBTP) = 1
          END IF
   90   CONTINUE
  100 CONTINUE
*
      IF ( NTEST .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' Matrix giving allowed combinations of types '
        WRITE(6,*)
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      END IF
*
      RETURN
      END
