!----------------------------------------------
!  simple test of Fortran GPU CUDA libraries
!----------------------------------------------

module gpu_stuff

implicit none

integer, parameter :: size_of_real=4
integer, parameter :: size_of_double=8
integer :: AllocateStatus, stat, ios

external :: cublas_init, cublas_set_matrix, cublas_get_matrix,cublas_set_vector,cublas_get_vector
external :: cublas_shutdown, cublas_alloc
external :: cublas_sgemm
integer :: cublas_alloc, cublas_set_matrix, cublas_get_matrix, cublas_set_vector, cublas_get_vector, cublas_sgemm
external :: snrm2
real(kind=size_of_real) :: snrm2

CHARACTER*100  BUFFER
character*10   method

integer :: n

integer :: i,j

integer :: c1,c2,cr,cm

real(kind=size_of_real),   allocatable:: A_r(:,:),B_r(:,:),C_r(:,:)
real(kind=size_of_double), allocatable:: A_d(:,:),B_d(:,:),C_d(:,:)

!#if ARCH_64
! 8-byte pointer on 64-bit architectures
integer*8:: devPtrA, devPtrB, devPtrC
!#else
!integer:: devPtrA, devPtrB, devPtrC
!#endif

contains

subroutine read_arguments
CALL GETARG(1,BUFFER); READ(BUFFER,*,IOSTAT=IOS) method
IF (IOS.NE.0) THEN
    stop "error in reading online argument !"
ELSE
    print *,'online argument read ... method=',method
ENDIF
CALL GETARG(2,BUFFER); READ(BUFFER,*,IOSTAT=IOS) n
IF (IOS.NE.0) THEN
    stop "error in reading online argument !"
ELSE
    print *,'online argument read in ...n=',n
ENDIF
end subroutine read_arguments

subroutine sgemm_test

allocate (A_r(n,n),B_r(n,n),C_r(n,n),STAT=AllocateStatus)
if (AllocateStatus.ne.0) then
   stop "error in main matrix allocations !"
else
   print *,"A,B,C matrixes allocated in CPU memory"
endif

do i=1,n
do j=1,n
 if (i.eq.j) then
   A_r(i,j)=1.0e0
 else
   A_r(i,j)=0.0e0
 endif
 B_r(i,j)=float(i+j) ! same jednotky...
 C_r(i,j)=0.0e0
enddo
enddo


call system_clock( c1, cr, cm )
call cublas_init
stat= cublas_Alloc(n*n,size_of_real, devPtrA)
if (stat .NE. 0) then
write(*,*) "device memory allocation of devPtrA failed"
call cublas_shutdown
stop
endif
stat=cublas_Alloc(n*n,size_of_real, devPtrB)
if (stat .NE. 0) then
write(*,*) "device memory allocation of devPtrB failed"
call cublas_shutdown
stop
endif
stat=cublas_Alloc(n*n,size_of_real, devPtrC)
if (stat .NE. 0) then
write(*,*) "device memory allocation of devPtrC failed"
call cublas_shutdown
stop
endif

stat=cublas_Set_Matrix(n,n,size_of_real,A_r,n,devPtrA,n)
if (stat .NE. 0) then
call cublas_free (devPtrA)
write(*,*) "data upload mtx Ar(CPU) -> devPtrA(GPU) failed"
call cublas_shutdown
stop
else
write(*,*) "data upload mtx Ar(CPU) -> devPtrA(GPU) : OK"
endif
stat=cublas_Set_Matrix(n,n,size_of_real,B_r,n,devPtrB,n)
if (stat .NE. 0) then
call cublas_free (devPtrB)
write(*,*) "data upload mtx Br(CPU) -> devPtrA(GPU) failed"
call cublas_shutdown
stop
else
write(*,*) "data upload mtx Br(CPU) -> devPtrB(GPU) : OK"
endif
stat=cublas_Set_Matrix(n,n,size_of_real,C_r,n,devPtrC,n)
if (stat .NE. 0) then
call cublas_free (devPtrC)
write(*,*) "data upload mtx Cr(CPU) -> devPtrC(GPU) failed"
call cublas_shutdown
stop
else
write(*,*) "data upload mtx Cr(CPU) -> devPtrC(GPU) : OK"
endif

stat=cublas_SGEMM('n','n', n,n,n,1.0e0,devPtrA,n,devPtrB,n,1.0e0,devPtrC,n)
if (stat .NE. 0) then
 stop "cublas_SGEMM failed!"
endif

stat=cublas_Get_Matrix(n,n,size_of_real,devPtrC,n,C_r,n)
if (stat .NE. 0) then
call cublas_free (devPtrC)
write(*,*) "data upload failed"
call cublas_shutdown
stop
endif

call system_clock( count=c2 )
print *, 'GPU cblas_sgemm  runtime:', 1.e3*real(c2-c1) / real(cr), 'ms'
call mtx_check(C_r)

call cublas_Free(devPtrA)
call cublas_Free(devPtrB)
call cublas_Free(devPtrC)
call cublas_shutdown

! now CPU sgemm...
call system_clock( c1, cr, cm )

call sgemm('n','n',n,n,n,1.0e0,A_r,n,B_r,n,-2.0e0,C_r,n)

call system_clock( count=c2 )
print *, ' CPU sgemm  runtime:', 1.e3*real(c2-c1) / real(cr), 'ms'
call mtx_check(C_r)

call cublas_shutdown
end subroutine sgemm_test

subroutine mtx_check(Mtx)
real*4, intent(in) :: Mtx(n,n)
real*4 ::  diag,offdiag
diag=0.0e0; offdiag=0.0e0
do i=1,n
 diag = diag + Mtx(i,i)
 do j=1,i
  if (j.ne.i) then
    offdiag = offdiag + (Mtx(i,j)+Mtx(j,i))
  endif
  !print *,i,j,c(i,j)
 enddo
enddo
!offdiag_man=(float(n)*float(n-1))
!diag_man=float(n)
!write(*,*)  'matica: ',diag,'- man:',diag_man,'offdiag:',offdiag,'- man:',offdiag_man
write(*,*)  'sum of diag ',diag,'sum of offdiag:',offdiag,'1,n:',Mtx(1,n),'n,1:',Mtx(n,1)
!print *,'matx norm:',snrm2(n*n,Mtx(1,1))

end subroutine mtx_check

end module gpu_stuff

Program GPU_math_tests
 
  use gpu_stuff

  print *,'hello'
  call read_arguments

  select case (method)
    case("sgemm")
      call sgemm_test
    case default
      write(*,*) 'no method for GPU math testing selected'
  end select

End Program GPU_math_tests



