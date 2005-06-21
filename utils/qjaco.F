       Program QJACOBITEST

C       PARAMETER(N=2,NZ=4,MAXBLOCKS=3)
C       PARAMETER(N=4,NZ=4,MAXBLOCKS=3)
        PARAMETER(N=20,NZ=2,NBLOCKS=0)
C       PARAMETER(N=30,NZ=4,MAXBLOCKS=3)
C       PARAMETER(N=50,NZ=4,MAXBLOCKS=3)
#include"implicit.h"
#include"priunit.h"
         DIMENSION AMTX(N,N,NZ)
         DIMENSION EV(N,N,NZ)

C        DIMENSION IBLOCKS(NBLOCKS,1:4)
         LOGICAL DOEVE


C          ...  start - row
C        IBLOCKS(1,1)=1
C          ...  start - col.
C        IBLOCKS(1,2)=2
C          ...  end  - row
C        IBLOCKS(1,3)=1
C          ... c end - col.
C        IBLOCKS(1,4)=4

C          ...  start - row
C        IBLOCKS(2,1)=2
C          ...  start - col.
C        IBLOCKS(2,2)=3
C          ...  end  - row
C        IBLOCKS(2,3)=2
C          ... c end - col.
C        IBLOCKS(2,4)=4

C       ... prepare a hermitian (quaternion) matrix
        CALL PREPMATRIX(AMTX,N,NZ)

        IPRINT=2


        DOEVE=.TRUE.
C       DOEVE=.FALSE.
        MXB=0
        IBL=0 
C       CALL QJACOBI(AMTX,EVXDUMMY,N,NZ,NBLOCKS,MXB,
C       CALL QJACOBI(AMTX,EVXDUMMY,N,NZ,NBLOCKS,MAXBLOCKS,
        CALL QJACOBI(AMTX,EV,N,NZ,NBLOCKS,
     &                IBL,DOEVE,IPRINT) 
C    &                IBL,DOEVE,IPRINT) 

       END

       SUBROUTINE PREPMATRIX(A,N,NZ)
#include"implicit.h"
#include"priunit.h"
       DIMENSION A(N,N,NZ)
       LOGICAL DORAND

C      DORAND=.FALSE.
       DORAND=.TRUE.

       IF (DORAND) THEN
       DO I=1,N
       DO J=I,N
         IF (I.EQ.J) THEN
C          A(I,I,1)=RAND(0)
C          A(I,I,1)=DFLOAT(3*I+4*J)
           A(I,I,1)=DFLOAT(3*I+4*J)*RAND(0)
C          A(I,I,1)=DFLOAT(3*I+4*J)/DFLOAT(J)
C          A(I,I,1)=DFLOAT(3*I+4*J)+1.8
C          A(I,I,1)=DFLOAT(3*I+4*J)+1.8+RAND(0)
           DO IZ=2,NZ
            A(I,I,IZ)=0.0D0
           ENDDO
         ELSE
           A(I,J,1)=DFLOAT(I+J)*RAND(0)
C          A(I,J,1)=DFLOAT(I+J+J)
           DO IZ=2,NZ
           A(I,J,IZ)=RAND(0)*DFLOAT(J+10+I) 
C          A(I,J,IZ)=DFLOAT(J+10+I)
           ENDDO
           A(J,I,1)=A(I,J,1)
           DO IZ=2,NZ
           A(J,I,IZ)=-A(I,J,IZ)
           ENDDO
         ENDIF
       ENDDO 
       ENDDO 

      ELSE
C...  read matrix from the file...
       LU=12

       OPEN(LU,FILE="matica",STATUS="OLD")

       READ(LU,*) N1,NZ1
       IF (N.NE.N1.OR.NZ.NE.NZ1) THEN
           write(lupri,*) 'N1,NZ1:',N1,NZ1
          CALL QUIT('PREPMATRIX....')
       ENDIF

       IF (NZ.EQ.4) THEN
       DO I=1,N
       DO J=I,N
         
         read(lu,*) I1,J1,AR,AI,AJ,AK
         IF (I.NE.I1.OR.J.NE.J1) CALL QUIT('wrond order !')

          A(I,J,1)=AR
          A(I,J,2)=AI
          A(I,J,3)=AJ
          A(I,J,4)=AK

          A(J,I,1)=AR
          A(J,I,2)=-AI
          A(J,I,3)=-AJ
          A(J,I,4)=-AK

       ENDDO 
       ENDDO 

       ENDIF

       CLOSE(LU,status="keep")

      ENDIF

       RETURN
       END
