      SUBROUTINE H1TERMS(SIGMA,C,NR,NK,NI,NJ,MAXK,HIJ,
     &                  KI,XKI,KJ,XKJ,SCLFAC,NCS,NCC,ITRANS)
*
* Update Sigma with one-electron contributions
*
* Sigma(Ir,Ic) = Sigma(Ir,Ic) + Sclfac*H_ij <Ic|a+ ic a jc!Jc> C(Ir,Jc)
*
* If ITRANS.EQ.1, then Sigma and C are transposed
*
* corresponding to current set of intermediate strings
*
* Jeppe Olsen, Fall of 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION C(*),HIJ(NI,NJ)
*. single operator connections
      DIMENSION KI(MAXK,*),XKI(MAXK,*)
      DIMENSION KJ(MAXK,*),XKJ(MAXK,*)
*. Input and output
      DIMENSION SIGMA(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' H1TERMS speaking '
        WRITE(6,*) ' HIJ matrix '
        CALL WRTMAT(HIJ,NI,NJ,NI,NJ)
        IF(ITRANS.EQ.1) THEN
          WRITE(6,*) ' Initial sigma block (transposed )'
          CALL WRTMAT(SIGMA,NCS,NR,NCS,NR)
        END IF
      END IF
*. Loop over intermediate strings
C?    WRITE(6,*) ' NK = ', NK
      DO KC = 1, NK
*
*. Number of nonvanishing connections from KC
*
       JJ = 0
       DO JORB = 1, NJ
         IF(KJ(KC,JORB).NE.0) JJ = JJ + 1
       END DO
       II = 0
       DO IORB = 1, NI
         IF(KI(KC,IORB).NE.0) II = II + 1
       END DO
C?     WRITE(6,*) ' II, JJ = ', II,JJ
*
       IF(II.NE.0.AND.JJ.NE.0) THEN
         DO IORB = 1, NI
           IC = KI(KC,IORB)
           IF(IC.NE.0) THEN
             SGNI = XKI(KC,IORB)
             DO JORB = 1, NJ
               JC = KJ(KC,JORB)
C?             write(6,*) ' IORB JORB', IORB,JORB
               IF(JC.NE.0) THEN
                 SGNJ = XKJ(KC,JORB)
C?               WRITE(6,*) ' SGNI SGNJ SCLFAC',SGNI,SGNJ,SCLFAC
                 FACTOR = SGNI*SGNJ*HIJ(IORB,JORB)*SCLFAC
C?               WRITE(6,*) ' IC JC FACTOR ', IC,JC,FACTOR
                 ONE = 1.0D0
                 IF(ITRANS.EQ.0) THEN
                   CALL VECSUM(SIGMA(1+(IC-1)*NIA),SIGMA(1+(IC-1)*NIA),
     &                  C(1+(JC-1)*NIA),ONE,FACTOR,NR)
                 ELSE
*. Not pretty, should be done somewhere else or on transposed matrices
                   DO IR = 1, NR
                     SIGMA(IC+(IR-1)*NCS) = SIGMA(IC+(IR-1)*NCS)
     &              +FACTOR * C(JC+(IR-1)*NCC)
                   END DO
                 END IF
*                ^ End of transpose switch
               END IF
             END DO
           END IF
         END DO
       END IF
      END DO
*
      IF(NTEST.GE.100) THEN
        IF(ITRANS.EQ.1) THEN
          WRITE(6,*) ' Updated sigma block (transposed )'
          CALL WRTMAT(SIGMA,NCS,NR,NCS,NR)
        END IF
      END IF

      RETURN
      END
