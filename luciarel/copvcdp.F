      SUBROUTINE COPVCDP(LUIN,LUOUT,SEGMNT,IREW,LBLK)
C
C COPY VECTOR ON FILE LUIN TO FILE LUOUT
*
* Packed version
C
C
C LBLK DEFINES STRUCTURE OF FILE
C Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        CALL REWINE( LUIN ,LBLK)
        CALL REWINE( LUOUT ,LBLK)
      END IF

C
C LOOP OVER BLOCKS
C
C?      write(6,*) ' COPVCD LBLK : ', LBLK
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
          WRITE(LUOUT) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          CALL ITODS (LBL,1,-1,LUOUT)
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
