!************************************************************************ 
! Standalone program for comprehensive testing of all diagonalization
! routines used by DIRAC. 

! These are Householder, RSJACOBI, QJACOBI and
! DSYEVR.
!
!
! The program either reads the matrix file, or generates random matrix.
!
!
! Written by Miro Ilias, July 2014 Banska Bystrica, Slovakia
!   with the help of P.S. Bagus, T.Saue
!
!
!************************************************************************ 

module time

  character*10 start_date, end_date ! format of date: "dd_mmm_yy "
  character*8  start_time, end_time ! format of time: "hh:mm:ss"

end module time

module  diagmod
  integer :: N=0,NZ=0
  integer :: lu=12
  real*8, allocatable :: Ar(:,:), ASr(:,:),Ai(:,:),ASi(:,:),  & 
                         Aj(:,:),ASj(:,:), Ak(:,:), ASk(:,:)
  character*7 :: matrixfile="hmatrix"

contains

  subroutine read_matrix
  ! allocates space and reads the matrix for the diagonalization testing from the file
    integer :: i,j, ii,jj
   
    open (lu,file=matrixfile, access="sequential", status="old")  
    read(lu,*) N,NZ
    print *,'read N=',N,' NZ=',NZ

    ! allocates space for matrices
    call allocate_space
    ! read matrix elements from its open file
    do i=1,N
    do j=1,N
     if (NZ==1) then
       read(lu,*) ii,jj,Ar(ii,jj)
       ASr(ii,jj)=Ar(ii,jj)
     else if (NZ==2) then
       read(lu,*) ii,jj,Ar(ii,jj),Ai(ii,jj)
       ASr(ii,jj)=Ar(ii,jj); ASi(ii,jj)=Ai(ii,jj)
     else if (NZ==4) then
       read(lu,*) ii,jj,Ar(ii,jj),Ai(ii,jj),Aj(ii,jj),Ak(ii,jj)
       ASr(ii,jj)=Ar(ii,jj); ASi(ii,jj)=Ai(ii,jj)
       ASj(ii,jj)=Aj(ii,jj); ASk(ii,jj)=Ak(ii,jj)
     endif
    enddo
    enddo
    close(lu,status="keep")
    print *,"matrix for diagonalization was read"
 
  end subroutine read_matrix

  subroutine allocate_space
    if (N<=1.or.NZ<=0) then 
       print *,'allocate_space: N=',N,' NZ=',NZ
       call quit('wrong allocation parameters N or NZ !')
    endif
    if (NZ==1) then
      allocate( Ar(N,N), stat=ierr )
      call alloc_stat_check(" Ar(N,N)",ierr)
      allocate( ASr(N,N), stat=ierr )
      call alloc_stat_check(" ASr(N,N)",ierr)
    else if (NZ==2) then
      allocate( Ar(N,N), Ai(N,N), stat=ierr )
      call alloc_stat_check(" Ar(N,N), Ai(N,N)",ierr)
      allocate( ASr(N,N), ASi(N,N), stat=ierr )
      call alloc_stat_check(" ASr(N,N), ASi(N,N)",ierr)
    else if (NZ==4) then
      allocate( Ar(N,N), Ai(N,N), Aj(N,N), Ak(N,N), stat=ierr )
      call alloc_stat_check(" Ar(N,N), Ai(N,N), Aj(N,N), Ak(N,N)",ierr)
      allocate( ASr(N,N), ASi(N,N), ASj(N,N), ASk(N,N), stat=ierr )
      call alloc_stat_check(" ASr(N,N), ASi(N,N), ASj(N,N), ASk(N,N)",ierr)
    else 
       print *,'wrong value of NZ=',NZ
       call quit('NZ must be either 1 or 2 or 4')
    endif
  end subroutine allocate_space

  subroutine generate_random_matrix
  ! allocates and generates random Hemitian matrix instead of reading it from a file
  end subroutine generate_random_matrix

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine test_lapack_dsyevr
 !
 ! tests lapack's routine DSYERV for real only symmetric matrices
 !
   use machine_parameters
   implicit none
#include "priunit.h"
   integer :: ierr, MATZ, LLWORK, KLWORK, LILWORK
   integer :: IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, LRA, LRV
   character*1 :: JOBZ, RANGE, UPLO
   real*8 :: ABSTOL, VL, VU
   integer, allocatable :: ISUPPZ( : ), IWORK( : )
   real*8, allocatable ::  W( : ), WORK( : ), Z( :,: ), EIG(:),VEC(:,:)
   real*8 :: DDUM, WORKDUM
   integer :: IWORKDUM,IDUM
   real*8, parameter :: D0=0.0D0

   ABSTOL=mp_sfmin; JOBZ='V'; RANGE = 'A'; UPLO = 'U'
   DDUM = D0; IDUM = 0; WORKDUM = D0; IWORKDUM = 0

   MATZ = 1 
   IF(MATZ.EQ.1) THEN
        JOBZ='V'
   ELSEIF(MATZ.EQ.0) THEN
        JOBZ='N'
   ELSE
     WRITE(LUPRI,*) 'QDIAG: Illegal value of MATZ: ',MATZ
     CALL QUIT('QDIAG: Illegal value of MATZ !')
   ENDIF

   allocate(EIG(N), stat=ierr)
   call alloc_stat_check("test_lapack_dsyerv EIG(N)",ierr)
    
   ! 1.call to determine size 
   IERR =  -1 ! on input: calculate temporary storage
   M = N ; LRA = N ;LRV = N
   CALL DSYEVR(JOBZ,'A','U',N,Ar,LRA,DDUM,DDUM,IDUM,IDUM,ABSTOL,  & 
               M,EIG,VEC,LRV,IDUM,EIG,-1,IWORKDUM,-1,IERR)
   LLWORK = NINT(EIG(1))
   LILWORK = IWORKDUM

   allocate(WORK(KLWORK), stat=ierr)
   call alloc_stat_check("test_lapack_dsyerv WORK(KLWORK)",ierr)

   allocate(IWORK(LILWORK), stat=ierr)
   call alloc_stat_check("test_lapack_dsyerv IWORK(LILWORK)",ierr)

   allocate(W(N),Z(N,N), stat=ierr)
   call alloc_stat_check("test_lapack_dsyerv W,Z",ierr)

   allocate(ISUPPZ(2*N),IWORK(N), stat=ierr)
   call alloc_stat_check("test_lapack_dsyerv ISUPPZ, IWORK",ierr)

   CALL DSYEVR( JOBZ, RANGE, UPLO, N, Ar, LDA, VL, VU, IL, IU,   & 
                ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,      & 
                IWORK, LIWORK, INFO )
    if (INFO /= 0) then
      print *,'The lapack_dsyerv routine ended with error ! INFO=',INFO
      !call quit('RS routine ended with error !')
    endif

  end subroutine test_lapack_dsyevr

  subroutine test_dirac_jacobi
  ! tests DIRAC's RSJACO for real only symmetric matrices
  end subroutine test_dirac_jacobi

  subroutine test_dirac_qjacobi
  ! tests DIRAC's QJACO for real, complex and quaternion Hermitian matrices
  end subroutine test_dirac_qjacobi


  subroutine test_dirac_householder
  ! tests DIRAC's Houeholder for real, complex and quaternion Hermitian matrices
  end subroutine test_dirac_householder

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  subroutine test_dirac_rs
! tests DIRAC's EISPACK routine RS real Hermitian matrix
    real*8, allocatable :: eigval(:),eigvect(:,:),temp1(:),temp2(:)
    integer :: ierr,i

! ... allocate arrays for the subroutine
    allocate(eigval(N),eigvect(N,N),temp1(N),temp2(N),stat=ierr)
    call alloc_stat_check("eigval(N),eigvect(N,N),temp1(N),temp2(N) ",ierr)
    print *,"arrays for RS subroutine allocated"

    matz=1 ! have both eigenvalues and eigenvectors
    call RS(N,N,Ar,eigval,matz,eigvect,temp1,temp2,ierr)
    if (ierr /= 0) then
      print *,'The RS routine ended with error !'
      !call quit('RS routine ended with error !')
    endif
! ... print out eigenvalues
    print *,"Eispack's RS eigenvalues:"
    do i=1,N
      print *,i,eigval(i)
    enddo
    deallocate(eigval,eigvect,temp1,temp2,stat=ierr)
  end subroutine test_dirac_rs

  subroutine test_dirac_ch
! tests DIRAC's EISPACK routine CH for complex Hermitian matrix
  end subroutine test_dirac_ch

end module diagmod


Program Test_DIRAC_Diagonalization_Routines
  use diagmod
 
  print *, "Welcome to the testing program for DIRAC's diagonalization routines !"

  call read_matrix

  call test_dirac_rs
  call test_lapack_dsyevr
  

End Program Test_DIRAC_Diagonalization_Routines




