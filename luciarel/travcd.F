      SUBROUTINE TRAVCD(VEC1,VEC2,X,NVECIN,NVECOUT,LUIN,LUOUT,
     &                  ICOPY,LBLK,LUSCR1,LUSCR2)
*
* NVECIN vectors reside on LU1, Transform these vectors,
*  using LUSCR1 and  LUSCR2  as
* scratch files,  with matrix X to produce output file LUOUT.
*
* Since LUIN is accessed several times it is always
* REWINDED. LUOUT is written to from current start.
*
* I ICOPY .ne. 0 the transformed vectors are written back to LUIN
*
* Jeppe Olsen, April 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION X(NVECIN,NVECOUT)
*. Scratch
      DIMENSION VEC1(*),VEC2(*)
*
C             MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
      DO IVECOUT = 1, NVECOUT
         CALL MVCSMD(LUIN,X(1,IVECOUT), LUSCR1,LUSCR2,VEC1,VEC2,
     &               NVECIN,1,LBLK)
         CALL REWINO(LUSCR1)
         CALL COPVCD(LUSCR1,LUOUT,VEC1,0,LBLK)
      END DO
*
      IF(ICOPY.EQ.1) THEN
        CALL REWINO(LUIN)
        CALL REWINO(LUOUT)
        DO IVECOUT = 1, NVECOUT
          CALL COPVCD(LUOUT,LUIN,VEC1,0,LBLK)
        END DO
      END IF
*
      RETURN
      END
