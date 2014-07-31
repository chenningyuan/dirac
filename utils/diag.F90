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

  implicit none
  character*10 start_date, end_date ! format of date: "dd_mmm_yy "
  character*8  start_time, end_time ! format of time: "hh:mm:ss"

end module time

module  diagmod

  integer :: N=0,NZ=0, print_level
  integer :: lu=12
  real*8, allocatable :: Ar(:,:), ASr(:,:),Ai(:,:),ASi(:,:),  & 
                         Aj(:,:),ASj(:,:), Ak(:,:), ASk(:,:)

  character*21 :: inpfile="DIAGONALIZE_TESTS.INP"

  character*50 :: title_text
  character*50 :: matrix_file_name

  logical :: do_dsyevr = .false.
  logical :: do_rs = .false.
  logical :: do_rsjaco = .false.
  logical :: do_qjaco = .false.
  logical :: do_eigv_check = .false.
  
contains

  subroutine read_input_file
! reads the text input file "DIAGONALIZE_TESTS.INP"
  implicit none
#include "priunit.h"
    logical :: input_found

    print *,"inpfile=",inpfile
    call read_menu_input(inpfile, lucmd, 'ALL', input_found)
    print *,'input_found=',input_found
 
    !call quit('read_input_file')
    ! stop "end in read_input_file"
  end  subroutine read_input_file


  subroutine read_matrix
  ! allocates space and reads the matrix for the diagonalization testing from the file
    integer :: i,j, ii,jj
   
    open (lu,file=trim(matrix_file_name), access="sequential", status="old")  
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
   integer :: i
   integer :: IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, LRA, LRV
   character*1 :: JOBZ, RANGE, UPLO
   real*8 :: ABSTOL, VL, VU
   integer, allocatable :: ISUPPZ(:), IWORK(:)
   real*8, allocatable ::  WORK(:), EIG(:),VEC(:,:)
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

   allocate(EIG(N),VEC(N,N), stat=ierr)
   call alloc_stat_check("test_lapack_dsyerv EIG(N),VEC(N,N)",ierr)
    
   ! 1.call to determine size 
   IERR =  0 
   ! FYI: LLWORK=-1; LIWORK=-1 on input means calculate temporary storage
   M = N ; LRA = N ;LRV = N
   CALL DSYEVR(JOBZ,'A','U',N,Ar,LRA,DDUM,DDUM,IDUM,IDUM,ABSTOL,  & 
               M,EIG,VEC,LRV,IDUM,EIG,-1,IWORKDUM,-1,IERR)

   LLWORK = NINT(EIG(1))
   LILWORK = IWORKDUM 

   !write(LUPRI,*) "LLWORK=",LLWORK
   allocate(WORK(LLWORK), stat=ierr)
   call alloc_stat_check("test_lapack_dsyerv WORK(KLWORK)",ierr)

   !write(LUPRI,*) "LLILWORK=",LILWORK
   allocate(IWORK(LILWORK), stat=ierr)
   call alloc_stat_check("test_lapack_dsyerv IWORK(LILWORK)",ierr)

   !print *,"N=",N," 2*N=",(2*N)
   allocate(ISUPPZ(2*N),stat=ierr)
   call alloc_stat_check("test_lapack_dsyerv ISUPPZ",ierr)

   IERR = 0
   CALL DSYEVR(JOBZ,'A','U',N,Ar,LRA,DDUM,DDUM,IDUM,IDUM,ABSTOL,  &
               M,EIG,VEC,LRV,ISUPPZ,WORK,LLWORK,IWORK,LILWORK,IERR)     
   if (IERR /= 0) then
     print *,'The lapack_dsyerv routine ended with error ! ierr=',ierr
   endif

  ! ... print out eigenvalues
  if (print_level >= 2) then
     print *,"LAPACK DSYEVR eigenvalues:"
     do i=1,N
       print *,i,EIG(i)
     enddo
   endif

   if (do_eigv_check) call eigv_check(EIG,VEC,'**** LAPACK DSYEVR ****')

   deallocate(WORK,IWORK,EIG,VEC,ISUPPZ,stat=ierr)

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

  subroutine test_PSB_householder
  ! tests Pauls Bagus own Houeholder routine for real symmetric matrices
  end subroutine test_PSB_householder


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  subroutine test_dirac_rs
! tests DIRAC's EISPACK routine RS real Hermitian matrix
#include "priunit.h"
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
    if (print_level >= 2) then
      print *,"Eispack's RS eigenvalues:"
      do i=1,N
        print *,i,eigval(i)
      enddo
    endif

    if (do_eigv_check) call eigv_check(eigval,eigvect,'**** EISPACK RS ****')

    deallocate(eigval,eigvect,temp1,temp2,stat=ierr)
  end subroutine test_dirac_rs

  subroutine test_dirac_ch
! tests DIRAC's EISPACK routine CH for complex Hermitian matrix
  end subroutine test_dirac_ch

  subroutine eigv_check(eigval,eigvect,routine_name)
  !-----------------------------------------------------------------------------------
  !
  ! Do the diagonalization checks using original Hermitian matrix, eigenvalues and
  ! eigenvectors
  ! 
  !           eigvect^+ * AS * eigvect - eigval  = 0
  !
  !           eigvect^+ * eigvect = I
  !
  !           eigvect * eigvect^+ = I
  !
  !-----------------------------------------------------------------------------------
   implicit none
#include "priunit.h"
#include "dgroup.h"
    real*8, intent(in) :: eigval(N),eigvect(N,N)
    character(*), intent(in) :: routine_name
    real*8, allocatable :: AXr(:,:),AYr(:,:),C(:,:)
    integer :: i,j, ierr, iprint
    real*8 :: C_norm_diag, C_norm_offdiag
    real*8 :: AXr_norm_diag, AXr_norm_offdiag
    real*8 :: AYr_norm_diag, AYr_norm_offdiag
    integer(kind=8) :: maxsize = -1
    
    allocate( AXr(N,N),AYr(N,N),C(N,N),stat=ierr )
    call alloc_stat_check("alloc AXr(N,N),AYr(N,N),C(N,N)",ierr)
 
   ! write(LUPRI,*) 'Eigenvectors matrix:'
   ! call prqmat(eigvect,N,N,N,N,1,IPQTOQ(1,0),LUPRI)

   ! Do  eigvect^{+)*ASr = AXr
    call dgemm('T','N',N,N,N,1.0d0,eigvect,N,ASr,N,0.0d0,AXr,N)

   ! Do  AXr * eigvect = C
    call dgemm('N','N',N,N,N,1.0d0,AXr,N,eigvect,N,0.0d0,C,N)

   ! Another check:  eigvect^{+)*eigvect = I = AXr
    call dgemm('T','N',N,N,N,1.0d0,eigvect,N,eigvect,N,0.0d0,AXr,N)
    if (print_level>=1) then
      write(LUPRI,*) "Printout of eigvect^{+)*eigvect =? I "
      call prqmat(AXr,N,N,N,N,1,IPQTOQ(1,0),LUPRI)
    endif

   ! Another check:  eigvect*eigvect^{+} = I = AYr
    call dgemm('N','T',N,N,N,1.0d0,eigvect,N,eigvect,N,0.0d0,AYr,N)
    if (print_level>=1) then
      write(LUPRI,*) "Printout of eigvect*eigvect^{+} =? I "
      call prqmat(AYr,N,N,N,N,1,IPQTOQ(1,0),LUPRI)
    endif

    ! substract eigenvalues ...
    C_norm_diag = 0.0d0; AXr_norm_diag = 0.0d0; AYr_norm_diag = 0.0d0
    do i=1,N
      C(i,i) = C(i,i) -eigval(i); 
      AXr(i,i)=AXr(i,i)-1.0d0; AYr(i,i)=AYr(i,i)-1.0d0
      C_norm_diag = C_norm_diag + DABS(C(i,i))
      AXr_norm_diag = AXr_norm_diag  + DABS(AXr(i,i))
      AYr_norm_diag = AYr_norm_diag  + DABS(AYr(i,i))
    enddo
    C_norm_diag   = C_norm_diag/DFLOAT(N)
    AXr_norm_diag = AXr_norm_diag/DFLOAT(N)
    AYr_norm_diag = AYr_norm_diag/DFLOAT(N)

    C_norm_offdiag = 0.0d0; AXr_norm_offdiag = 0.0d0;  AYr_norm_offdiag = 0.0d0
    do i = 1, N
    do j = 1, N
     if (i /= j) then
       C_norm_offdiag = C_norm_offdiag + DABS(C(i,j))
       AXr_norm_offdiag = AXr_norm_offdiag + DABS(AXr(i,j))
       AYr_norm_offdiag = AYr_norm_offdiag + DABS(AYr(i,j))
     endif
    enddo
    enddo
    C_norm_offdiag = C_norm_offdiag/((DFLOAT(N)*DFLOAT(N))-DFLOAT(N))
    AXr_norm_offdiag = AXr_norm_offdiag/((DFLOAT(N)*DFLOAT(N))-DFLOAT(N))
    AYr_norm_offdiag = AYr_norm_offdiag/((DFLOAT(N)*DFLOAT(N))-DFLOAT(N))

    write(LUPRI,"(2X,A)") routine_name
    write(LUPRI,"(2X,A)") "Checking average differences U^{+}*A*U - eps =  0 = C"
    write(LUPRI,"(A,D10.4,2X,A,D10.4)") "C_norm_diag:",C_norm_diag,"C_norm_offdiag:",C_norm_offdiag
    write(LUPRI,"(2X,A)") "Checking average differences U^{+}*U -I =  0 = AXr"
    write(LUPRI,"(A,D10.4,2X,A,D10.4)") "AXr_norm_diag:",AXr_norm_diag,"AXr_norm_offdiag:",AXr_norm_offdiag
    write(LUPRI,"(2X,A)") "Checking average differences U*U^{+} -I =  0 = AYr"
    write(LUPRI,"(A,D10.4,2X,A,D10.4)") "AYr_norm_diag:",AYr_norm_diag,"AYr_norm_offdiag:",AYr_norm_offdiag

    deallocate( AXr, AYr, C, stat=ierr )
    call alloc_stat_check("dealloc AXr(N,N), AYr(N,N) C(N,N)",ierr)

  end subroutine eigv_check

end module diagmod


Program Test_DIRAC_Diagonalization_Routines
  use diagmod
#include "priunit.h"
 
   write(LUPRI,'(2x,a)') "Welcome to the testing program for DIRAC's diagonalization routines !"

  call DBLGRP

  call read_input_file

  call read_matrix

  if (do_rs)      call test_dirac_rs
  if (do_dsyevr)  call test_lapack_dsyevr
  

End Program Test_DIRAC_Diagonalization_Routines

 ! reading subroutine
     subroutine read_input_diagonalization_tests(word, kw_section)
      use input_reader
      use diagmod
      implicit none
!     ---------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     ---------------------------------------------------------------------------
      call reset_available_kw_list()

   !  if (kw_matches(word, '.RADINT')) then
   !     call kw_read(word, num_grid_cfg_radint)
   !  end if

     if (kw_matches(word, '.TITLE ')) then
       ! read(get_file_unit(),'(a)') title_text
         call kw_read(word, title_text)
        print *,'title_text=',title_text
     end if

     if (kw_matches(word, '.PRINT ')) then
         call kw_read(word, print_level)
        print *,'print_level=',print_level
     end if

     if (kw_matches(word, '.MTXFIL')) then
         call kw_read(word, matrix_file_name)
        print *,'matrix_file_name=',matrix_file_name
     end if

     if (kw_matches(word, '.DSYEVR')) then
        do_dsyevr = .true.
        print *,'doing dsyevr !'
     end if

     if (kw_matches(word, '.RS    ')) then
        do_rs = .true.
        print *,'doing RS !'
     end if

     if (kw_matches(word, '.RSJACO')) then
        do_rsjaco = .true.
        print *,'doing RSJACO !'
     end if

     if (kw_matches(word, '.QJACO ')) then
        do_qjaco = .true.
        print *,'doing QJACO !'
     end if

     if (kw_matches(word, '.EIGVTS')) then
        do_eigv_check = .true.
        print *,'doing EIGVTST !'
     end if

     call check_whether_kw_found(word, kw_section)

   end subroutine

