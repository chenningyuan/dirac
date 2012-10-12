!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

!> Program for testing (and improving) the quaternion matrix Jacobi (block) diagonalization procedure
!> written by Miro ILIAS, Strasbourg,2005, improved in Tel Aviv, 2012
module qjacobi_mod
      !> process also eigenvectors
      LOGICAL :: DOEVE 
      !> check the presence of the file with input data
      LOGICAL::  ISFILE 
      LOGICAL :: FROMFILE
      real(kind=8), allocatable :: A(:,:,:),EV(:,:,:)
      integer, allocatable :: IBL(:,:)
      integer :: N,NZ,NBL,IPRINT
      integer, parameter :: LUF=10
      real*8, PARAMETER :: D0=0.0D0
contains

subroutine GET_N_NZ
!> get in the basic variables:  N,NZ,NBL,IPRINT
implicit none
#include"priunit.h"
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
         write(lupri,*) 'No "matrix_input" file found, reading in from terminal.'
         write(lupri,*) 'Enter the N>0,NZ=1,2 or 4 :'
         read(*,*) N,NZ ! TODO: verify entering data
         NBL = 0
         write(lupri,*) 'Enter the NBL (zero for no blocks) :'
         read(*,*) NBL
         write(lupri,*) 'Enter IPRINT, print level :'
         read(*,*) IPRINT
      ENDIF
end subroutine

subroutine PREPMATRIX
!> ... prepare the testing matrix - allocate
implicit none
#include"priunit.h"
 integer :: I,J,IBL1,IBL2,IBL3,IBL4,I1,J1,IZ,istat
 real(kind=8) :: AR, AI,AJ,AK
      FROMFILE=.FALSE.
      DOEVE = .TRUE.
      IF (ISFILE) THEN
         write(lupri,*) 'PRINT level=',IPRINT
         write(lupri,*) 'Number of blocks=',NBL
         IF (NBL.GT.0) THEN
!         ... read the indexes of blocks
          DO I = 1, NBL
           read(luf,*) IBL1,IBL2,IBL3,IBL4
              IBL(I,1) = IBL1
              IBL(I,2) = IBL2 
              IBL(I,3) = IBL3 
              IBL(I,4) = IBL4
           IF (IPRINT.GE.2) THEN
            write(lupri,*) 'Read block coordinates:', & 
            IBL(I,1),IBL(I,2),IBL(I,3),IBL(I,4)
           ENDIF
          ENDDO
         ENDIF
         read(luf,*)
         read(luf,*) FROMFILE
         read(luf,*)
      ELSE
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
!          ... fill with random numbers
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

!       ....now fill the whole matrix 
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
             write(lupri,*) 'A filled:', & 
             I,J,':',A(I,J,1),A(I,J,2),A(I,J,3),A(I,J,4)
            endif
          ENDIF
          IF (IPRINT.GE.3) THEN
            write(lupri,'(a,i3,a,i3,a,4d12.5)')  & 
           'random number in matrix element A(',I,';',J,')=',  & 
           (A(I,J,IZ),IZ=1,NZ)
          ENDIF
      ENDDO 
      ENDDO 
      IF (ISFILE) CLOSE(LUF,STATUS='KEEP')
end subroutine

REAL*8 FUNCTION RANDOMNUMBER(I,J)
!> Function used to generate pseudorandom numbers, also according to 
!> entering indexes I,J
! See also http://gcc.gnu.org/onlinedocs/gfortran/RAND.html
implicit none
#include"priunit.h"
integer, intent(in) :: I,J
integer (kind=4) :: rand_param=1
real(kind=8), external :: RAND
  !print *,'rand=',RAND(rand_param)
  !RANDOMNUMBER = DFLOAT(3*I+4*J)*RAND(rand_param) 
 ! miro: deactivated RAND, as this is not universally working
  RANDOMNUMBER = DFLOAT(3*I+4*J)
END FUNCTION
end module qjacobi_mod

Program QJACOBITEST
use qjacobi_mod
implicit none
#include"priunit.h"
integer :: istat
  ! ... get N,NZ from the file
  CALL GET_N_NZ
  !> allocate matrix to be diagonalized
  allocate(A(N,N,NZ),stat=istat)
  allocate(EV(N,N,NZ),stat=istat) !> space for eigenvectors (if desired)
  !> allocate space for blocks positions in the matrix to be diagonalized
  allocate(IBL(NBL,1:4),stat=istat)
  ! ... fill hermitian (real, quaternion) matrix
  CALL PREPMATRIX
  ! ... call of the external diagonalization routine ...
#ifdef MOD_UNRELEASED
  CALL QJACOBI(A,EV,N,NZ,NBL,IBL,DOEVE,IPRINT)
#else
  CALL QUIT('QJACOI is not in this version !')
#endif
END PROGRAM
