      SUBROUTINE SHTOOB2(NSHPIR,NIRREP,MXPOBS,NSMOB,NOSPIR,IOSPIR,
     &           NOBPS,NOB)
*
* Number of shells per irrep => Number of orbitals per symmetry
*
* =====
* Input
* =====
*
*  NSHPIR : Number of shells per irrep
*  NIRREP : Number of irreps
*  MXPOBS : Largest allowed number of orbitals symmetries
*  NSMOB  : Number of orbital symmetries
*  NOSPIR : Number of orbital symmetries per irrep
*  IOSPIR : Orbital symmetries per irrep
*
* ======
* Output
* ======
*  NOBPS  : Number of orbitals per symmetry
*  NOB    : Number of orbitals
*
* Jeppe Olsen, Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION NSHPIR(*),NOSPIR(*),IOSPIR(*)
*. Output
      DIMENSION NOBPS(*)
      CALL ISETVC(NOBPS,0,NSMOB)
      NOB = 0
      JSM = 0
      DO 100 IRREP = 1, NIRREP
      C DO 90 ISM = 1, NOSPIR(IRREP)
C         IISM = IOSPIR(ISM,IRREP)
          IISM  = IOSPIR(IRREP)
          if (IISM.ne.0) then
             JSM = JSM + 1
             NOBPS(IISM) = NOBPS(IISM) + NSHPIR(IRREP)
          end if
          NOB = NOB + NSHPIR(IRREP)
   90 C CONTINUE
  100 CONTINUE
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
         WRITE(6,*) ' SHTOOB Speaking '
         WRITE(6,*) ' =============== '
         WRITE(6,*) ' Number of orbitals obtained ', NOB
         WRITE(6,*) ' Number of orbitals per symmetry '
         CALL IWRTMA(NOBPS,1,NSMOB,1,NSMOB)
      END IF
*
      RETURN
      END
