      SUBROUTINE MERGE_VEC(LU,VEC)
*
* A vector is defined as a multiple record file on LU
* Assemble vector
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
      DIMENSION VEC(*)
*
      NTEST = 000
*
*. Loop over records
      IBASE = 1
      IBLK = 0
 1000 CONTINUE
      CALL IFRMDS(LBLK,1,-1,LU)
      IBLK = IBLK + 1
       IF(LBLK .GE. 0 ) THEN
         CALL FRMDSC(VEC(IBASE),LBLK,-1,LU,IMZERO,IAMPACK)
         IBASE = IBASE + LBLK
       END IF
      if (NTEST.ge.100) then
        WRITE(6,*) ' IBLK,LBLK,IBASE',IBLK,LBLK,IBASE
      end if
      IF(LBLK.GE.0) GOTO 1000
*
      LENGTH = IBASE - 1
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Vector assembled from file LU = ', LU
        WRITE(6,*) ' Length = LENGTH '
        CALL WRTMAT(VEC,1,LENGTH,1,LENGTH)
      END IF
*
      RETURN
      END
