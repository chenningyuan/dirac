      SUBROUTINE PNT2DM2(I12SM,NSMOB,IPSM,JPSM,IJSM,IPNTR)
*
* Pointer to two dimensional array
*
* Jeppe Olsen, Jan. 98
*
* Modified from PNT2DM2 : Reference to symmetry products referred to
*                         SYMCOM
*                         Generation of symmetry of last index eliminated
*
* =====
* Input
* =====
* I12SM  : ne.0 => restrict to lower half
*          eq.0 => completer matrix
* NSMOB : Number of orbital symmetries
* IPSM : Number of orbitals per symmetry for index 1
* JPSM : Number of orbitals per symmetry for index 2
* IJSM  : Symmetry of two index array
*
* =======
* Output
* =======
* IPNTR : Pointer to block with first index of given symmetry
*         = 0 indicates forbidden block
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IPSM(*),JPSM(*)
*.Output
      DIMENSION IPNTR(*)
*
      IZERO = 0
      CALL ISETVC(IPNTR,IZERO,NSMOB)
COLD  CALL ISETVC(ISM2 ,IZERO,NSMOB)
*
      IOFF = 1
      DO 100 ISM = 1,NSMOB
C            SYMCOM(ITASK,IOBJ,I1,I2,I12)
        CALL SYMCOM(2,1,ISM,JSMA,IJSM)
        IF(JSMA.EQ.0) GOTO 100
        JSM = I_GET_ADJSYM(JSMA)
        IF(I12SM.EQ.0.OR.ISM.GE.JSM) THEN
*. Allowed block
          IPNTR(ISM) = IOFF
          IF(I12SM.GT.0.AND.ISM.EQ.JSM) THEN
            IOFF = IOFF + IPSM(ISM)*(IPSM(ISM)+1)/2
          ELSE
            IOFF = IOFF + IPSM(ISM)*JPSM(JSM)
          END IF
        END IF
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' dimension of two-dimensional array ',IOFF-1
      END IF
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Pointer '
        CALL IWRTMA(IPNTR,1,NSMOB,1,NSMOB)
      END IF
*
      RETURN
      END
