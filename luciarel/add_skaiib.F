      SUBROUTINE ADD_SKAIIB(SB,NI,NIA,SKAIIB,NKA,NIB,
     &                          I,ISCA,SSCA)
*
* Update Transposed sigma block with contributions for given orbital index j
* from the matrix S(Ka,i,Ib)
*
* S(Ib,Isca(Ka)) =  S(Ib,Isca(Ka)) + Ssca(Ka)*S(Ka,I,Ib)
*
*
* For efficient processing of alpha-beta loop
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
       DIMENSION SKAIIB(*),SSCA(*),ISCA(*)
*. Input and Output
       DIMENSION SB(NIB,NIA)
*
C     LBLK = 100
      LBLK = 40
      NBLK = NIB/LBLK
      IF(LBLK*NBLK.LT.NIB) NBLK = NBLK + 1
      DO ICBL = 1, NBLK
        IF(ICBL.EQ.1) THEN
          ICOFF = 1
        ELSE
          ICOFF = ICOFF + LBLK
        END IF
        ICEND = MIN(ICOFF+LBLK-1,NIB)
        ICONST = NKA*NI
        IADR0 =  (I-1)*NKA+(ICOFF-1-1)*NKA*NI
        IF(ICEND.GT.ICOFF) THEN
*. Use form with Inner loop over IB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
C             IADR = KA + (I-1)*NKA+(ICOFF-1-1)*NKA*NI
              IADR = IADR0 + KA
              DO IB = ICOFF,ICEND
*. Adress of S(Ka,i,Ib)
                IADR = IADR + ICONST
                SB(Ib,IROW) = SB(Ib,IROW)+S*SKAIIB(IADR)
              END DO
            END IF
          END DO
        ELSE
*. Form with no loop over IB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
              IADR = IADR0 + KA + ICONST
C             DO IB = ICOFF,ICEND
*. Adress of S(Ka,i,Ib)
C               IADR = IADR + ICONST
                SB(ICOFF,IROW) = SB(ICOFF,IROW)+S*SKAIIB(IADR)
C             END DO
            END IF
          END DO
        END IF
*       ^ End of test of ICOFF=ICEND
      END DO
*
      RETURN
      END
