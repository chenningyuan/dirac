      SUBROUTINE MVCSMD2(LUIN,FAC,FACLUOUT,LUOUT,LUSCR,
     &           VEC1,VEC2,NVEC,IREW,LBLK)
*
* LUOUT = Factor * LUOUT + sum(I=1,nvec)fac(I) LUIN(I)
*
* LUOUT and LUSCR are rewinded and input data
* on these files are lost.
*
* MVCSMD2 is identical MVCSMD, except that
* the input vector on LUOUT can be included
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(1)
*
      IF(IREW .NE. 0 ) CALL REWINE(LUIN,LBLK)
*
      LLUOUT = LUOUT
      LLUSCR = LUSCR
      DO IVEC = 1, NVEC
        CALL REWINE(LLUSCR,LBLK)
        CALL REWINE(LLUOUT,LBLK)
        IF( IVEC .EQ. 1 ) THEN
          CALL VECSMD(VEC1,VEC2,FAC(1),FACLUOUT,LUIN,
     &                LUOUT,LUSCR,0,LBLK)
        ELSE
          ONE = 1.0D0
          CALL VECSMD(VEC1,VEC2,FAC(IVEC),ONE,LUIN,LLUOUT,LLUSCR,
     &                0,LBLK)
        END IF
*
        LBUF = LLUOUT
        LLUOUT = LLUSCR
        LLUSCR = LBUF
      END DO
*
      IF(LLUOUT.NE.LUOUT) THEN
C            COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
        CALL COPVCD(LLUOUT,LUOUT,VEC1,1,LBLK)
      END IF
*
      RETURN
      END
