      SUBROUTINE DMTVDS(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,
     &                  ISCAT,XSCAT,NSCAT,LBLK,XINOUT)
*
* Multiply/divide elements of two vectors residing on disc
* Elements corresponding to absolute adresses in ISCAT
* are set to the elements of XSCAT
*
* For elements not in ISCAT the operation is thus :
*
* For INV.NE. 0 :
*
*    V3(I) = (V1(I)+FAC)-1 * V2(I)
*    LU3      LU1            LU2
*
* For INV.NE. 0 :
*
*    V3(I) = (V1(I)+FAC) * V2(I)
*    LU3         LU1        LU2
*
* LBLK defines structure of files
*
* Overlap between input and output vector is also calculated
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  VEC1(*),VEC2(*)
      DIMENSION  XSCAT(NSCAT),ISCAT(NSCAT)
*
      REAL*8 INPROD
*
      XINOUT = 0.0D0
*
      IF ( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU1
          REWIND LU2
          REWIND LU3
        ELSE
          CALL REWINE( LU1,LBLK)
          CALL REWINE( LU2,LBLK)
          CALL REWINE( LU3,LBLK)
         END IF
      END IF
*
* Loop over blocks
*
        IOFF = 1
 1000 CONTINUE
        IF (LBLK .GT. 0 ) THEN
          LBL1 = LBLK
          LBL2 = LBLK
        ELSE IF( LBLK .EQ. 0 ) THEN
          READ(LU1) LBL1
          READ(LU2) LBL2
          WRITE(LU3) LBL1
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL1,1,-1,LU1)
          CALL IFRMDS(LBL2,1,-1,LU2)
          CALL ITODS (LBL1,1,-1,LU3)
        END IF
        IF(LBL1 .NE. LBL2 ) THEN
          WRITE(6,'(A,2I3)') ' DIFFERENT BLOCKSIZES IN DMTVSD : '
     &                     , LBL1,LBL2
          Call Abend2( ' DIFFERENT BLOCKSIZES IN DMTVSD ' )
        END IF
        IF(LBL1 .GE. 0 ) THEN
          IF(      LBLK .GE.0 ) THEN
            KBLK = LBL1
          ELSE
            KBLK = -1
          END IF
          LENGTH = LBL1
          CALL FRMDSC(VEC1,LENGTH,KBLK,LU1,IMZERO,IAMPACK)
          CALL FRMDSC(VEC2,LENGTH,KBLK,LU2,IMZERO,IAMPACK)
          IF( LBL1 .GT. 0 )THEN
            IF(INV .NE. 0 ) THEN
              CALL DIAVC2(VEC1,VEC2,VEC1,FAC,LENGTH)
            ELSE
              CALL VVTOV(VEC1,VEC2,VEC1,LBL1)
              CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FAC,LENGTH)
            END IF
          END IF
*
          IF(NSCAT.NE.0) THEN
            IFIRST = IOFF
            ILAST = IOFF + LENGTH - 1
            DO 100 I = 1, NSCAT
              IF(IFIRST .LE. ISCAT(I) .AND. ISCAT(I) .LE. ILAST )
     &        VEC1(ISCAT(I)-IOFF+1) = XSCAT(I)
  100       CONTINUE
          END IF
*
          XINOUT = XINOUT + INPROD(VEC1,VEC2,LENGTH)
          CALL TODSC(VEC1,LENGTH,KBLK,LU3)
          IOFF = IOFF + LENGTH
        END IF

*
      IF( LBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
*
      RETURN
      END
