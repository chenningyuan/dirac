      SUBROUTINE SIGVST_REL(ISGVST,NSMST)
*
* Obtain ISGVST(ISM) : Symmetry of sigma v on string of symmetry ism
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER ISGVST(*)
*
      DO 100 ISM = 1, NSMST
C            MLSM_REL(IML,IPARI,ISM,TYPE,IWAY)
        CALL MLSM_REL(IML,IPARI,ISM,'ST',2)
        MIML = - IML
        CALL MLSM_REL(MIML,IPARI,MISM,'ST',1)
        ISGVST(ISM) = MISM
  100 CONTINUE
*
      NTEST = 1
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' ISGVST array '
        WRITE(6,*) ' ============ '
        CALL IWRTMA(ISGVST,1,NSMST,1,NSMST)
      END IF
*
      RETURN
      END
