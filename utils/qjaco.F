!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C  Program for testing (and improving) the quaternion Jacobi (block) diagonalization procedure.
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!> Program for testing (and improving) the quaternion Jacobi (block)
!! diagonalization procedure
      Program QJACOBITEST
      PARAMETER(MAXMTX=900000,MAXBL=400)
#include"implicit.h"
#include"priunit.h"
      !> process also eigenvectors
      LOGICAL DOEVE 
      !> check the presence of the file
      LOGICAL ISFILE 
      DIMENSION A(MAXMTX),EV(MAXMTX),IBL(MAXBL)

     
C    ... get N,NZ from the file
      CALL GETNNZ(N,NZ,NBL,LUF,ISFILE,IPRINT)
      IF (MAXMTX.LT.N*N*NZ) CALL QUIT('too big N! increase the MAXMTX!')

C       ... prepare a hermitian (quaternion) matrix
      CALL PREPMATRIX(A,N,NZ,NBL,IBL,DOEVE,IPRINT,LUF,ISFILE)
 
C     ... call of the main routine ...
      CALL QJACOBI(A,EV,N,NZ,NBL,IBL,DOEVE,IPRINT)

      STOP "QJACOBI test finished"
      END

C *******************************************************************8
C
C  Get basic date: N,NZ,NBL,IPRINT
C
C *******************************************************************8
!> get the basic variables:  N,NZ,NBL,IPRINT
      SUBROUTINE GETNNZ(N,NZ,NBL,LUF,ISFILE,IPRINT)
#include"implicit.h"
      PARAMETER (D0=0.0D0)
#include"priunit.h"
      LOGICAL DOEVE,ISFILE,FROMFILE

      luf = 10
      INQUIRE(FILE="matrix_input",EXIST=ISFILE)
      IF (ISFILE) THEN
         write(lupri,*) '...going to read the "matrix_input" file' 
         open(luf,file="matrix_input",status="old")
         read(luf,*)
         read(luf,*) N,NZ
         read(luf,*)
         read(luf,*) IPRINT
         read(luf,*)
         read(luf,*) NBL
      ELSE
         write(lupri,*) 'No matrix input file !'
         write(lupri,*) 'Enter the N,NZ :'
         read(*,*) N,NZ
         NBL = 0
         write(lupri,*) 'Enter the NBL :'
         read(*,*) NBL
         write(lupri,*) 'Enter IPRINT :'
         read(*,*) IPRINT
      ENDIF

      RETURN
      END    
C ============================================================
C            ... prepare testing matrix
C ============================================================
!> ... prepare the testing matrix
      SUBROUTINE PREPMATRIX(A,N,NZ,NBL,IBL,DOEVE,IPRINT,LUF,ISFILE)
#include"implicit.h"
      PARAMETER (D0=0.0D0)
#include"priunit.h"
      !> matrix to be diagonalized
      DIMENSION A(N,N,NZ) 
      !> blocks positions in the matrix to be diagonalized
      INTEGER IBL(NBL,1:4)
      !> function for getting randon numbers
      REAL*8  RANDOMNUMBER  !< random number ...
      LOGICAL DOEVE, !< process eigenvactors
     &        ISFILE, !< check presence of file
     &        FROMFILE !< somethink...

      FROMFILE=.FALSE.
      DOEVE = .TRUE.
      IF (ISFILE) THEN
         write(lupri,*) 'PRINT level=',IPRINT
         write(lupri,*) 'Number of blocks=',NBL
         IF (NBL.GT.0) THEN
C         ... read the indexes of blocks
          DO I = 1, NBL
           read(luf,*) IBL1,IBL2,IBL3,IBL4
              IBL(I,1) = IBL1
              IBL(I,2) = IBL2 
              IBL(I,3) = IBL3 
              IBL(I,4) = IBL4
           IF (IPRINT.GE.2) THEN
            write(lupri,*) 'Read block coordinates:',
     &      IBL(I,1),IBL(I,2),IBL(I,3),IBL(I,4)
           ENDIF
          ENDDO
         ENDIF
         read(luf,*)
         read(luf,*) FROMFILE
         read(luf,*)
      ELSE
         IPRINT = 1
         IF (NBL.GT.0) THEN
          DO I=1,NBL
            write(*,*) 'Enter ',I,'. block coordinates:'
            read(*,*)  IBL(I,1),IBL(I,2),IBL(I,3),IBL(I,4)
          ENDDO
         ELSE
           write(*,*) 'No blocks....full diagonalization'
         ENDIF
      ENDIF

      DO I = 1, N
      DO J = I, N
          IF (FROMFILE) THEN
            AI = D0
            AJ = D0
            AK = D0
            IF (NZ.EQ.1) THEN
             read(luf,*) I1,J1,AR
            ELSE IF (NZ.EQ.2) THEN
             read(luf,*) I1,J1,AR,AI
            ELSE IF (NZ.EQ.4) THEN
             read(luf,*) I1,J1,AR,AI,AJ,AK
             If (IPRINT.GE.2) THEN
              write(lupri,*) 'read i1,j1,ar,ai,aj,ak:',I1,J1,AR,AI,AJ,AK
             EndIf
            ENDIF
            IF (I1.NE.I.OR.J1.NE.J) THEN
             write(lupri,*) '   read I,J:',I1,J1
             write(lupri,*) 'desired I,J:',I,J
             CALL QUIT('...wrong indexes !')
            ENDIF
            IF (I1.EQ.J1) THEN
             IF (AI.NE.D0.OR.AJ.NE.D0.OR.AK.NE.D0) THEN
              write(lupri,*) 'I1,J1 --> AI,AJ,AK:',I1,J1,AI,AJ,AK
              call QUIT('diagonal i,j,k elements must be zero !!!')
             ENDIF
            ENDIF
          ELSE
C          ... fill with random numbers
           IF (I.NE.J) THEN
            AR = RANDOMNUMBER(I,J)
            AI = RANDOMNUMBER(I,J)
            AJ = RANDOMNUMBER(I,J)
            AK = RANDOMNUMBER(I,J)
           ELSE
            AR = RANDOMNUMBER(I,J)
            AI = D0
            AJ = D0
            AK = D0
           ENDIF
          ENDIF

C       ....now fill the matrix 
          IF (NZ.EQ.1) THEN
            A(I,J,1) = AR
            A(J,I,1) = AR
          ELSE IF (NZ.EQ.2) THEN
            A(I,J,1) = AR
            A(J,I,1) = AR
            A(I,J,2) = AI
            A(J,I,2) = -AI
          ELSE IF (NZ.EQ.4) THEN
            A(I,J,1) =  AR
            A(J,I,1) =  AR
            A(I,J,2) =  AI
            A(J,I,2) = -AI
            A(I,J,3) =  AJ
            A(J,I,3) = -AJ
            A(I,J,4) =  AK
            A(J,I,4) = -AK
            if (IPRINT.GE.10) then
             write(lupri,*) 'A filled:',
     &       I,J,':',A(I,J,1),A(I,J,2),A(I,J,3),A(I,J,4)
            endif
          ENDIF
          IF (IPRINT.GE.2) THEN
            write(lupri,*)
     &     'filled ready matrix element A(',I,';',J,')=',
     &     (A(I,J,IZ),IZ=1,NZ)
          ENDIF
      ENDDO 
      ENDDO 

      IF (ISFILE) CLOSE(LUF,STATUS='KEEP')

      RETURN
      END

      REAL*8 FUNCTION RANDOMNUMBER(I,J)
!> Function used to generate numbers to filling the matrix with them
#include"implicit.h"
#include"priunit.h"
C     RANDOMNUMBER = DFLOAT(3*I+4*J)*RAND(0) 
      !> get the random number
      RANDOMNUMBER = DFLOAT(3*I+4*J)*RAND() !< restriction...
      RETURN
      END

