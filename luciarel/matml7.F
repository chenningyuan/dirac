      SUBROUTINE MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
     &                  NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
C
C MULTIPLY A AND B TO GIVE C
C
C     C =  FACTORC*C + FACTORAB* A * B             FOR ITRNSP = 0
C
C     C =  FACTORC*C + FACTORAB* A(TRANSPOSED) * B FOR ITRNSP = 1
C
C     C =  FACTORC*C + FACTORAB* A * B(TRANSPOSED) FOR ITRNSP = 2
C
C... JEPPE OLSEN,
C
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION A(NAROW,NACOL),B(NBROW,NBCOL)
      DIMENSION C(NCROW,NCCOL)
*
c     CALL QENTER('MATML')
*
      NTEST = 00
      IF ( NTEST .NE. 0 ) THEN
      WRITE(6,*)
      WRITE(6,*) ' A, B and C MATRIX FROM MATML7 '
      WRITE(6,*)
      CALL WRTMAT(A,NAROW,NACOL,NAROW,NACOL)
      CALL WRTMAT(B,NBROW,NBCOL,NBROW,NBCOL)
      CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
      WRITE(6,*)      ' NCROW NCCOL NAROW NACOL NBROW NBCOL '
      WRITE(6,'(6I6)')  NCROW,NCCOL,NAROW,NACOL,NBROW,NBCOL
      WRITE(6,*) ' ITRNSP : ', ITRNSP
      END IF
*
      IF(NAROW*NACOL*NBROW*NBCOL*NCROW*NCCOL .EQ. 0 ) THEN
        IZERO = 1
      ELSE
        IZERO = 0
      END IF

      IESSL = 0
      ICONVEX = 0
      IF ((ICONVEX.EQ.1 .OR.IESSL.EQ.1).AND. IZERO.EQ.0 ) THEN
*. DGEMM from CONVEX/ESSL  lib
        LDA = MAX(1,NAROW)
        LDB = MAX(1,NBROW)
*
        LDC = MAX(1,NCROW)
        IF(ITRNSP.EQ.0) THEN
        CALL DGEMM('N','N',NAROW,NBCOL,NACOL,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        ELSE IF (ITRNSP.EQ.1) THEN
        CALL DGEMM('T','N',NACOL,NBCOL,NAROW,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        ELSE IF(ITRNSP.EQ.2) THEN
        CALL DGEMM('N','T',NAROW,NBROW,NACOL,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
        END IF
      ELSE
* Use Jeppes version ( it should be working )
        IF( ITRNSP .EQ. 0 ) THEN
* ======
* C=A*B
* ======
*
C         CALL SCALVE(C,FACTORC,NCROW*NCCOL)
          DO J =1, NCCOL
*. Initialize with FACTORC*C(I,J) + FACTORAB*A(I,1)*B(1,J)
            B1J = FACTORAB*B(1,J)
            DO I = 1, NCROW
              C(I,J) = FACTORC*C(I,J) + B1J*A(I,1)
            END DO
*. and the major part
            DO K =2, NBROW
              BKJ = FACTORAB*B(K,J)
              DO I = 1, NCROW
                C(I,J) = C(I,J) + BKJ*A(I,K)
              END DO
            END DO
          END DO
*
        END IF
        IF ( ITRNSP .EQ. 1 ) THEN
*
* =========
* C=A(T)*B
* =========
*
          DO J = 1, NCCOL
            DO I = 1, NCROW
              T = 0.0D0
              DO K = 1, NBROW
                T = T  + A(K,I)*B(K,J)
              END DO
              C(I,J) = FACTORC*C(I,J) + FACTORAB*T
            END DO
          END DO
        END IF
C
        IF ( ITRNSP .EQ. 2 ) THEN
C ===========
C. C = A*B(T)
C ===========
C         CALL SCALVE(C,FACTORC,NCROW*NCCOL)
          DO J = 1,NCCOL
*. Initialization
            BJ1 = FACTORAB*B(J,1)
            DO I = 1, NCROW
              C(I,J) = FACTORC*C(I,J) + BJ1*A(I,1)
            END DO
*. And the rest
            DO K = 2,NBCOL
              BJK = FACTORAB*B(J,K)
              DO I = 1, NCROW
                C(I,J) = C(I,J) + BJK*A(I,K)
              END DO
            END DO
          END DO
        END IF
      END IF
C
      IF ( NTEST .NE. 0 ) THEN
      WRITE(6,*)
      WRITE(6,*) ' C MATRIX FROM MATML7 '
      WRITE(6,*)
      CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
      END IF
C
c     CALL QEXIT('MATML')
      RETURN
      END
