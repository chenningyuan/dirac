      SUBROUTINE FIFAMS_REL(FIFA,RHO1,IOBSM,NSMOB,LOBSM,NACOB,
     &                      NORBT,ISTOB)
*
* Update inactive fock matrix with active contributions
*
*     FIFA(I,J) = FIFA(I,J) + sum(k,l) ((ij!ab)-0.5*(ib!ja))
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FIFA(*),RHO1(NACOB,NACOB)
      INTEGER IOBSM(*),LOBSM(*),ISTOB(*)
*
*.  Assume spatial symmetric fock matrix
      IJSM = 1
      IJ = 0
      DO ISM = 1, NSMOB
        CALL SYMCOM_REL(2,6,ISM,JSM,IJSM)
C            SYMCOM_REL(ITASK,IOBJ,I1,I2,I12)
        IF(JSM.NE.0) THEN
          DO I = IOBSM(ISM),IOBSM(ISM) + LOBSM(ISM)-1
            DO J = IOBSM(JSM),I                          
C?            write(6,*) ' I J ', I,J
              IP = ISTOB(I)
              JP = ISTOB(J)
C?            write(6,*) ' IP JP ', IP,JP
               IJ= IJ + 1
               DO IA = 1, NACOB
                 DO IB = 1, NACOB
                   FIFA(IJ) = FIFA(IJ) 
     &           + RHO1(IA,IB)
     &           *(GTIJKL_REL(IP,JP,IA,IB)-0.5*GTIJKL_REL(IP,IB,IA,JP))
                 END DO
               END DO
            END DO
          END DO
        END IF
      END DO
*
      NTEST = 1
      IF(NTEST.NE.0) THEN
*
       WRITE(6,*) ' FI + FA in Symmetry blocked form '
       WRITE(6,*) ' ================================='
       WRITE(6,*) 
       ISYM = 1
       CALL APRBLM2(FIFA,LOBSM,LOBSM,NSMOB,ISYM)
      END IF
* 
      RETURN
      END
