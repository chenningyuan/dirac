      SUBROUTINE SORLOW(WRK,STVAL,ISTART,KZVAR,KEXST2,JEXSTV,IPRT)
C
C PURPOSE: FIND THE KEXSTV LOWEST VALUES IN WRK(KZVAR)
C          CHECK IF THE LAST ELEMENTS ARE DEGENERATE
C          JEXSTV IS THE NUMBER OF SORTED ELEMENTS WHERE
C          NO DEGENERACIES OCCUR AMONG THE HIGHEST ELEMENTS
* INPUT
*======
* WRK : ARRAY TO BE SORTED
* KEXST2 : NUMBER OF ELEMENTS TO BE OBTAINED  + 1
* KZVAR : LENGTH OG WRK
* OUTPUT
*=======
* STVAL  : STVAL(I) IS VALUES OF SORTED ELEMENT I
* ISTART : SCATTER POINTER ISTART(I) IS ADRESS IN FULL LIST OF
*          SORTED ELEMENT I
* KZVAR : LENGTH OG WRK
* JEXSTV : NUMBER OF ELEMENTS RETURNED , CAN BE LESS THAN
*          KEXSTV-1 IF DEGENERNCIES OCCURS  AMONG THE LAST ELEMENTS
*
C
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION WRK(*),STVAL(*),ISTART(*)
C
      PARAMETER ( TOLEQL=1.0D-6 )
      PARAMETER ( D0=0.0D0 , D1=1.0D0 , DM1 = -1.0D0 )
      LOGICAL FULL
C
      IF(KEXST2.GE.KZVAR) THEN
         FULL = .TRUE.
         KEXSTV = KZVAR
      ELSE
         FULL = .FALSE.
         KEXSTV = KEXST2+1
      END IF
*
      DO 100 K=1,KEXSTV
         ISTART(K) = K
         STVAL(K)  = WRK(K)
 100  CONTINUE
      KK=KEXSTV
C
      DO 210 I=1,KEXSTV
         DO 220 J=I+1,KEXSTV
            IF ((STVAL(J)-STVAL(I)) .LT. D0) THEN
               X =STVAL(I)
               II=ISTART(I)
               STVAL(I) =STVAL(J)
               ISTART(I)=ISTART(J)
               STVAL(J) =X
               ISTART(J)=II
            ENDIF
 220     CONTINUE
 210  CONTINUE
      GO TO 115
C     REPEAT UNTIL ...
 105  CONTINUE
         DO 110 I = KEXSTV,2,-1
            J = I - 1
            IF ((STVAL(J)-STVAL(I)) .GT. D0) THEN
               X =STVAL(I)
               II=ISTART(I)
               STVAL(I) =STVAL(J)
               ISTART(I)=ISTART(J)
               STVAL(J) =X
               ISTART(J)=II
            ELSE
               GO TO 115
            ENDIF
 110     CONTINUE
 115     CONTINUE
         STMAX=STVAL(KEXSTV)
C
 125     CONTINUE
            KK=KK+1
            IF (KK.LE.KZVAR) THEN
               IF ((STMAX-WRK(KK)).GT.D0) THEN
                  ISTART(KEXSTV) = KK
                  STVAL(KEXSTV)  = WRK(KK)
                  GO TO 105
C     ^--------------------
               ENDIF
               GO TO 125
C        ^--------------
            END IF
C
C     Check for degeneracy among diagonal elements
C
      JEXSTV = KEXSTV
      IF(.NOT.FULL) THEN

  160   JEXSTV = JEXSTV - 1
        IF (
     *  (ABS(STVAL(JEXSTV+1)-STVAL(JEXSTV))).LE.TOLEQL) GO TO 160
      END IF
      IF((IPRT.GT.8).AND.(KEXST2.NE.JEXSTV)) WRITE(6,1600)KEXST2,JEXSTV

1600  FORMAT(/' NUMBER OF START VECTORS IS DIMINISHED TO',I5,' FROM',I5)
      IF (IPRT.GE.100) THEN
         WRITE(6,*) '(I,(ISTART(I),STVAL(I)),I=1,JEXSTV)'
      DO 170 I = 1,JEXSTV
            WRITE(6,*) I,ISTART(I),STVAL(I)
  170    CONTINUE
         WRITE(6,*) 'THE FIRST',JEXSTV,' ELEMENTS ARE SELECTED.'
      END IF
C
C     END OF SORVAL
C
      RETURN
      END
