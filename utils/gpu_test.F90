!----------------------------------------------
!
!  simple test of Fortran GPU CUDA libraries
! in the DIRAC enviroment
!
!----------------------------------------------

module gpu_stuff

implicit none

integer, parameter :: size_of_real  =4
integer, parameter :: size_of_double=8
integer :: AllocateStatus, stat, ios

external ::  cublas_init, cublas_free, cublas_set_matrix
external ::  cublas_get_matrix,cublas_set_vector,cublas_get_vector
external ::  cublas_shutdown, cublas_alloc
external ::  cublas_sgemm, cublas_sgemv, cublas_dgemm
integer ::   cublas_alloc, cublas_set_matrix, cublas_get_matrix, cublas_free
integer ::  cublas_set_vector, cublas_get_vector, cublas_sgemm, cublas_sgemv,cublas_dgemm
external :: snrm2
real(kind=size_of_real) :: snrm2

CHARACTER*100  BUFFER
character*10   method

integer :: verb=0 ! default verbosity level

integer :: n
integer :: i,j
integer(kind=8) :: c1,c2
integer(kind=8) :: cr,cm ! accuracy in nanoseconds

real(kind=size_of_real),   allocatable:: A_r(:,:),B_r(:,:),C_r(:,:),x_r(:),y_r(:)
real(kind=size_of_double), allocatable:: A_d(:,:),B_d(:,:),C_d(:,:),x_d(:),y_d(:)

#if ARCH_64
! 8-byte pointers on 64-bit architectures
integer*8 :: devPtrA,devPtrB,devPtrC,devPtrx,devPtry
#else
integer   :: devPtrA,devPtrB,devPtrC,devPtrx,devPtry
#endif

contains

subroutine usage
write(*,*) 'usage: gpu_test.x  <sgemm/sgemv/dgemm>  <n> <verbosity level>'
end subroutine usage

subroutine  print_time(message)
! http://gcc.gnu.org/onlinedocs/gfortran/SYSTEM_005fCLOCK.html
! print spent time in routine ...
character(*), intent(in) :: message
!!print *, 'GPU cblas_sgemm  runtime:', 1.e3*real(c2-c1) / real(cr), 'ms'
write(*,'(2x,a,a,a,f13.8,a)') 'Routine <',message,'> spent runtime:', float(c2-c1) / float(cr), ' seconds.'
end subroutine  print_time

subroutine read_arguments
CALL GETARG(1,BUFFER); READ(BUFFER,*,IOSTAT=IOS) method
IF (IOS.NE.0) THEN
    call usage
    stop "error in reading online argument !"
ELSE
    continue
    !print *,'online argument read ... method=',method
ENDIF
CALL GETARG(2,BUFFER); READ(BUFFER,*,IOSTAT=IOS) n
IF (IOS.NE.0) THEN
  call usage
  stop "error in reading online argument !"
ELSE
  continue   
  !print *,'online argument read in ...n=',n
ENDIF
! read in verbosity level - optional argument
CALL GETARG(3,BUFFER); READ(BUFFER,*,IOSTAT=IOS) verb
IF (IOS.NE.0) THEN
  write(*,*) 'using default verbosity level=',verb
ELSE
  write(*,*) 'online argument read in - verbosity level =',verb
ENDIF
end subroutine read_arguments

subroutine sgemm_test
allocate (A_r(n,n),B_r(n,n),C_r(n,n),STAT=AllocateStatus)
if (AllocateStatus.ne.0) then
  stop "error in matrix allocations !"
else
  if (verb.ge.1) then
    write(*,*) "A_r,B_r,C_r matrixes allocated in CPU memory"
  endif
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
!call system_clock( c1, cr, cm )
call cublas_init
stat= cublas_Alloc(n*n,size_of_real, devPtrA)
if (stat .NE. 0) then
  write(*,*) "device memory allocation of devPtrA failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "device memory allocation of devPtrA OK"
  endif
endif
stat=cublas_Alloc(n*n,size_of_real, devPtrB)
if (stat .NE. 0) then
  write(*,*) "device memory allocation of devPtrB failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "device memory allocation of devPtrB OK"
  endif
endif
stat=cublas_Alloc(n*n,size_of_real, devPtrC)
if (stat .NE. 0) then
  write(*,*) "device memory allocation of devPtrC failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "device memory allocation of devPtrC OK"
  endif
endif
! ... start the time measurement
call system_clock( c1, cr, cm )
stat=cublas_Set_Matrix(n,n,size_of_real,A_r,n,devPtrA,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrA)
  write(*,*) "data upload mtx A_r(CPU) -> devPtrA(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx A_r(CPU) -> devPtrA(GPU) : OK"
  endif
endif
stat=cublas_Set_Matrix(n,n,size_of_real,B_r,n,devPtrB,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrB)
  write(*,*) "data upload mtx B_r(CPU) -> devPtrB(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx B_r(CPU) -> devPtrB(GPU) : OK"
  endif
endif
stat=cublas_Set_Matrix(n,n,size_of_real,C_r,n,devPtrC,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrC)
  write(*,*) "data upload mtx C_r(CPU) -> devPtrC(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx C_r(CPU) -> devPtrC(GPU) : OK"
  endif
endif
stat=cublas_SGEMM('n','n', n,n,n,1.0e0,devPtrA,n,devPtrB,n,1.0e0,devPtrC,n)
if (stat .NE. 0) then
  stop "cublas_SGEMM failed!"
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx C_r(CPU) -> devPtrC(GPU) : OK"
  endif
endif
stat=cublas_Get_Matrix(n,n,size_of_real,devPtrC,n,C_r,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrC)
write(*,*) "data upload failed"
call cublas_shutdown
stop
endif
call system_clock( count=c2 )
call print_time("GPU sgemm")
call mtx_check(C_r)
stat=cublas_free(devPtrA)
stat=cublas_free(devPtrB)
stat=cublas_free(devPtrC)
call cublas_shutdown
! now CPU sgemm...
call system_clock( c1, cr, cm )
call sgemm('n','n',n,n,n,1.0e0,A_r,n,B_r,n,-2.0e0,C_r,n)
call system_clock( count=c2 )
call print_time("CPU sgemm")
call mtx_check(C_r)
call cublas_shutdown
end subroutine sgemm_test

subroutine dgemm_test
allocate (A_d(n,n),B_d(n,n),C_d(n,n),STAT=AllocateStatus)
if (AllocateStatus.ne.0) then
   stop "error in main matrix allocations !"
else
   if (verb.ge.1) then
     write(*,*) "A_d,B_d,C_d matrixes allocated in CPU memory"
   endif
endif
do i=1,n
do j=1,n
 if (i.eq.j) then ! A_d is unit matrix
   A_d(i,j)=1.0e0
 else
   A_d(i,j)=0.0e0
 endif
 B_d(i,j)=dfloat(i+j) ! B_d is symmetric matrix
 C_d(i,j)=0.0e0
enddo
enddo
!call system_clock( c1, cr, cm )
call cublas_init
stat= cublas_Alloc(n*n,size_of_double, devPtrA)
if (stat .NE. 0) then
write(*,*) "device memory allocation of devPtrA failed"
call cublas_shutdown
stop
else
  if (verb.ge.1) then
    write(*,*) "device memory allocation of devPtrA done"
  endif
endif
stat=cublas_Alloc(n*n,size_of_double,devPtrB)
if (stat .NE. 0) then
write(*,*) "device memory allocation of devPtrB failed"
call cublas_shutdown
stop
else
  if (verb.ge.1) then
    write(*,*) "device memory allocation of devPtrB done"
  endif
endif
stat=cublas_Alloc(n*n,size_of_double, devPtrC)
if (stat .NE. 0) then
write(*,*) "device memory allocation of devPtrC failed"
call cublas_shutdown
stop
else
  if (verb.ge.1) then
    write(*,*) "device memory allocation of devPtrC done"
  endif
endif
! start the time measurement
call system_clock( c1, cr, cm )
stat=cublas_Set_Matrix(n,n,size_of_double,A_d,n,devPtrA,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrA)
  write(*,*) "data upload mtx A_d(CPU) -> devPtrA(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx A_d(CPU) -> devPtrA(GPU) : OK"
  endif
endif
stat=cublas_Set_Matrix(n,n,size_of_double,B_d,n,devPtrB,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrB)
  write(*,*) "data upload mtx B_d(CPU) -> devPtrA(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx B_d(CPU) -> devPtrB(GPU) : OK"
  endif
endif
stat=cublas_Set_Matrix(n,n,size_of_double,C_d,n,devPtrC,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrC)
  write(*,*) "data upload mtx C_d(CPU) -> devPtrC(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx C_d(CPU) -> devPtrC(GPU) : OK"
  endif
endif
stat=cublas_DGEMM('n','n', n,n,n,1.0d0,devPtrA,n,devPtrB,n,1.0d0,devPtrC,n)
if (stat .NE. 0) then
  stop "cublas_DGEMM failed!"
else
  if (verb.ge.1) then
    write(*,*) "cublas_DGEMM done OK"
  endif
endif
stat=cublas_Get_Matrix(n,n,size_of_double,devPtrC,n,C_d,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrC)
  write(*,*) "data download GPU devPtrC -> CPU mtx C_d failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data download GPU devPtrC -> CPU mtx C_d OK"
  endif
endif
call system_clock( count=c2 )
call print_time("GPU dgemm")
call mtx_check_d(C_d)
stat=cublas_free(devPtrA)
stat=cublas_free(devPtrB)
stat=cublas_free(devPtrC)
call cublas_shutdown
! now CPU sgemm...
call system_clock( c1, cr, cm )
call dgemm('n','n',n,n,n,1.0d0,A_d,n,B_d,n,-2.0e0,C_d,n)
call system_clock( count=c2 )
call print_time("CPU dgemm")
call mtx_check_d(C_d)
call cublas_shutdown
end subroutine dgemm_test

subroutine sgemv_test
allocate (A_r(n,n),x_r(n),y_r(n),STAT=AllocateStatus)
if (AllocateStatus.ne.0) then
  stop "error in allocations !"
else
  if (verb.ge.1) then
    write(*,*) "A_r,x_r,y_r objects allocated in CPU memory"
  endif
endif
! fill allocated objects
do i=1,n
do j=1,n
 if (i.eq.j) then
   A_r(i,j)=1.0e0
 else
   A_r(i,j)=0.0e0
 endif
 x_r(i)=float(1)
 y_r(i)=0.0e0
enddo
enddo
!call system_clock( c1, cr, cm )
call cublas_init
stat=cublas_Alloc(n*n,size_of_real, devPtrA)
if (stat .NE. 0) then
  write(*,*) "GPU device memory allocation of devPtrA failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "GPU device memory allocation of devPtrA OK"
  endif
endif
stat= cublas_Alloc(n,size_of_real,devPtrx)
if (stat .NE. 0) then
  write(*,*) "GPU device memory allocation of devPtrx failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "GPU device memory allocation of devPtrx OK"
  endif
endif
stat=cublas_Alloc(n,size_of_real,devPtry)
if (stat .NE. 0) then
  write(*,*) "GPU device memory allocation of devPtry failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "GPU device memory allocation of devPtry OK"
  endif
endif
! start the time measurement
call system_clock( c1, cr, cm )
stat=cublas_Set_Matrix(n,n,size_of_real,A_r,n,devPtrA,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrA)
  write(*,*) "data upload mtx A_r(CPU) -> devPtrA(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx A_r(CPU) -> devPtrA(GPU) : OK"
  endif
endif
stat=cublas_set_vector(n,size_of_real,x_r,1,devPtrx,1)
if (stat .NE. 0) then
  stat=cublas_free (devPtrx)
  write(*,*) "data upload vector x_r(CPU) -> devPtrx(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload vector x_r(CPU) -> devPtrx(GPU) : OK"
  endif
endif
stat=cublas_set_vector(n,size_of_real,y_r,1,devPtry,1)
if (stat .NE. 0) then
  stat=cublas_free (devPtry)
  write(*,*) "data upload vector y_r(CPU) -> devPtry(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload vector y_r(CPU) -> devPtry(GPU) : OK"
  endif
endif
!y = beta*y+alpha*A*x
stat=cublas_sgemv('n',n,n,1.0e0,devPtrA,n,devPtrx,1,0.0e0,devPtry,1)
if (stat .NE. 0) then
  stop "cublas_sgemv failed!"
else
  if (verb.ge.1) then
    write(*,*) 'cublas_sgemv done'
  endif
endif
! get the vector from GPU into y_r...
stat=cublas_get_vector(n,size_of_real,devPtry,1,y_r,1)
call system_clock( count=c2 )
call print_time("GPU sgemv")
write(*,*) 'CPU sgemv vector control output: sum,norm=',sum(y_r),snrm2(n,y_r,1)
stat=cublas_free(devPtrA)
stat=cublas_free(devPtrx)
stat=cublas_free(devPtry)
call system_clock( c1, cr, cm )
!now perform CPU lib routine:  y = beta*y+alpha*A*x
call sgemv('n',n,n,-2.0e0,A_r,n,x_r,1,1.0e0,y_r,1)
call system_clock( count=c2 )
call print_time("CPU sgemv")
write(*,*) 'CPU sgemv vector control output: sum,norm=',sum(y_r),snrm2(n,y_r,1)
end subroutine sgemv_test

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
  enddo
enddo
write(*,*)  'sum of diag ',diag,'sum of offdiag:',offdiag,'1,n:',Mtx(1,n),'n,1:',Mtx(n,1)
end subroutine mtx_check

subroutine mtx_check_d(Mtx)
real*8, intent(in) :: Mtx(n,n)
real*8 ::  diag,offdiag
diag=0.0d0; offdiag=0.0d0
do i=1,n
  diag = diag + Mtx(i,i)
  do j=1,i
    if (j.ne.i) then
      offdiag = offdiag + (Mtx(i,j)+Mtx(j,i))
    endif
  enddo
enddo
write(*,*)  'sum of diag ',diag,'sum of offdiag:',offdiag,'1,n:',Mtx(1,n),'n,1:',Mtx(n,1)
end subroutine mtx_check_d

end module gpu_stuff

Program GPU_math_tests
  use gpu_stuff
  call read_arguments
  select case (method)
    case("sgemm")
      write(*,'(2x,a)') 'sgemm GPU vs CPU test'
      call sgemm_test
    case("dgemm")
      write(*,'(2x,a)') 'dgemm GPU vs CPU test'
      call dgemm_test
    case("sgemv")
      write(*,'(2x,a)') 'sgemv GPU vs CPU test'
      call sgemv_test
    case default
      write(*,*) 'no method for GPU math testing selected'
      call usage
      stop
  end select
End Program GPU_math_tests

