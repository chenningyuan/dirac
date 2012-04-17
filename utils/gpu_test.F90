!---------------------------------------------------------------------------------------------
!
! Simple test of Fortran GPU CUDA & CULA libraries
! in the DIRAC environment. 
!
! This program measures performances of few  GPU CUBLAS/CULA
! mathematical routines against corresponding CPU BLAS/LAPACK counterparts.
!
! Note:  CULA Dense Free Edition implements popular single precision real and complex functions from LAPACK,
! full commercial version of CULA meets all the LAPACK/BLAS needs of the DIRAC suite,
! see http://www.culatools.com/downloads/dense/
!
! Written by Miro Ilias, Banska Bystrica, March/April 2012
!
!-----------------------------------------------------------------------------------------------

subroutine fortran_compiler_info
! see  http://fortranwiki.org/fortran/show/Predefined+preprocessor+macros ...
write(*,'(A)') '--- Fortran compiler info --- '
#if defined VAR_GFORTRAN
! http://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
write(*,'(A,3(I1,A))') 'GNU, gfortran, version  ',__GNUC__,'.',__GNUC_MINOR__,'.',__GNUC_PATCHLEVEL__
write(*,*) 'version:',__VERSION__
#endif
#if defined VAR_IFORT
!write(*,'(A)') 'Intel Fortran info :'
#if defined __ECC
write(*,*) __ECC
#endif
#if defined __INTEL_COMPILER
write(*,*) 'Intel Fortran, version:',__INTEL_COMPILER
#endif
write(*,*) 'intel compiler build date:', __INTEL_COMPILER_BUILD_DATE
!#include "mkl.fi"
!DEC$IF DEFINED INTEL_MKL_VERSION
!#pragma message "INTEL_MKL_VERSION=",INTEL_MKL_VERSION
!write(*,*) 'Intel MKL  version:',INTEL_MKL_VERSION

!!DEC$IF INTEL_MKL_VERSION .EQ. 100304
!!* 	Code to be conditionally compiled 
!!DEC$ENDIF
!DEC$ENDIF
#endif
end subroutine fortran_compiler_info

subroutine CulaAbouts
implicit none
#if defined USE_CULA
! get info about CULA ...
! see  nm -D /usr/local/cula/lib64/libcula_fortran.so | less
! CULA Programmer’s Guide,  www.culatools.com,  Release R14 (CUDA 4.1)
integer, external ::  CULA_GETVERSION, CULA_GETCUBLASRUNTIMEVERSION, CULA_GETCUDADRIVERVERSION
integer, external :: CULA_GETCUDARUNTIMEVERSION, CULA_GETDEVICEINFO, CULA_GETCUDAMINIMUMVERSION
write(*,'(A)') '--- CULA (and CUDA) GPU info ---'
  write(*,'(1X,A,I5)') 'CULA version=',CULA_GETVERSION()
  write(*,'(1X,A,I5)') 'CULA_CUBLAS runtime version :',CULA_GETCUBLASRUNTIMEVERSION()
  write(*,'(1X,A,I5)') 'CULA_CUDA   runtime version :',CULA_GETCUDARUNTIMEVERSION()
  write(*,'(1X,A,I5)') 'CULA_CUDA   driver version  :', CULA_GETCUDADRIVERVERSION()
 if (CULA_GETCUDARUNTIMEVERSION().LT.CULA_GETCUDAMINIMUMVERSION()) then
   write(*,*) ' CULA_GETCUDAMINIMUMVERSION=',CULA_GETCUDAMINIMUMVERSION()
   write(*,*) "warning: CUDA runtime version < cuda minumum required version !"
 endif
 if (CULA_GETCUBLASRUNTIMEVERSION().LT.CULA_GETCUDAMINIMUMVERSION()) then
   write(*,*) "warning: CUBLAS runtime version < cuda minumum required version !"
   write(*,*) '  cublasRuntimeVersion',cula_GetCublasRuntimeVersion()
 endif
 if (CULA_GETCUDADRIVERVERSION() .LT.CULA_GETCUDAMINIMUMVERSION()) then
   write(*,*) "CUDA driver version < cuda minumum required version !"
 endif
#endif
end subroutine CulaAbouts

subroutine CudaAbouts
! No fortran functions about CUDA version, obtainable only from CULA
! see less /usr/local/cuda/src/fortran.c
implicit none
external :: CudaAbouts_C
call CudaAbouts_C()
end subroutine CudaAbouts

module fortran_timing
! http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/fortran/lin/compiler_f/lref_for/source_files/rfsystmc.htm
! If the type is INTEGER(2), count_rate is 1000. If the type is INTEGER(4), count_rate is 10000. If the type is INTEGER(8), count_rate is 1000000.
! but gfortran 4.8.0 ([trunk revision 185723]) has COUNT_RATE=1000000000 for INTEGER(8)...
implicit none
integer(kind=8) :: c1,c2,cr,cm ! accuracy in nanoseconds
!integer(kind=4) :: c1,c2,cr,cm ! accuracy in miliseconds
! default control values - these differ from compiler to compiler !
integer, parameter :: count_rate_kind4=1000, count_rate_kind8=1000000
contains

subroutine print_time_info
call system_clock( count_rate=cr, count_max=cm )
write(*,'(1x,a,i8,a,i20)') "--- Fortran system_clock timing routine, COUNT_RATE=",cr," COUNT_MAX=",cm
!if (cr.ne.count_rate_kind8) then
!  write(*,'(2x,a,i7)') 'wrong system value of system COUNT_RATE parameter, must be (kind=8) ',count_rate_kind8
  !write(*,'(2x,a,i7)') 'wrong system value of system COUNT_RATE parameter, must be (kind=4) ',count_rate_kind4
!endif
end subroutine print_time_info

subroutine  print_time(message)
! print spent time in routine between c1 and c2 
! http://gcc.gnu.org/onlinedocs/gfortran/SYSTEM_005fCLOCK.html
! print spent time in routine ...
character(*), intent(in) :: message
write(*,'(2x,a,a,a,f13.8,a)') 'Routine <',message,'> spent runtime:', (dfloat(c2-c1) / dfloat(cr)), ' seconds.'
end subroutine  print_time

end module fortran_timing

module gpu_stuff

  use fortran_timing

implicit none

integer, parameter :: size_of_real  =4
integer, parameter :: size_of_double=8
integer :: AllocateStatus, stat, ios

! CUDA specifix library obejcts 
#if USE_CUBLAS
external ::  cublas_init, cublas_free, cublas_set_matrix
external ::  cublas_get_matrix,cublas_set_vector,cublas_get_vector
external ::  cublas_shutdown, cublas_alloc
external ::  cublas_sgemm, cublas_sgemv, cublas_dgemm
integer ::   cublas_alloc, cublas_set_matrix, cublas_get_matrix, cublas_free
integer ::   cublas_set_vector, cublas_get_vector, cublas_sgemm, cublas_sgemv,cublas_dgemm
#endif

! CULA specific library objects 
#if USE_CULA
!external  :: cula_initialize, cula_shutdown,cula_sgemv
integer, external  :: cula_initialize, cula_shutdown,cula_sgemv, cula_dgemv
integer, external  :: CULA_SGEQRF, cula_sgesv, cula_sgemm, cula_dgemm
#endif

real(kind=size_of_real),   external :: snrm2 ! for checking  the norm of a vector ..
real(kind=size_of_double), external :: dnrm2 ! for checking  the norm of a vector ..

CHARACTER*100  BUFFER

integer :: verb=0 ! default verbosity level

integer :: n
integer :: i,j

! common allocatble objects - both for CUDA (BLAS) and CULA (BLAS+LAPACK) routines...
real(kind=size_of_real),   allocatable:: A_r(:,:),B_r(:,:),C_r(:,:),x_r(:),y_r(:)
real(kind=size_of_double), allocatable:: A_d(:,:),B_d(:,:),C_d(:,:),x_d(:),y_d(:)
REAL(kind=size_of_real),   ALLOCATABLE :: TAU_r(:), TAU1_r(:)
REAL(kind=size_of_real),   ALLOCATABLE :: WORK_r(:)
integer, allocatable :: ipiv(:), ipiv1(:)

#if ARCH_64
! 8-byte pointers on 64-bit architectures
integer*8 :: devPtrA,devPtrB,devPtrC,devPtrx,devPtry
#else
integer   :: devPtrA,devPtrB,devPtrC,devPtrx,devPtry
#endif

integer, parameter  :: maxnu=20
integer :: no_rout=0
character*15 :: routine(maxnu), method

contains

subroutine usage
 integer :: no_rout1=0
routine(1)="sgemm_cuda"
routine(2)="sgemv_cuda"
routine(3)="dgemm_cuda" 
routine(4)="sgemm_cuda_c"
routine(5)="sgemv_cuda_c"
routine(6)="dgemm_cuda_c"
no_rout=6
#if USE_CULA
 ! cula lapack
 no_rout1=no_rout
routine(no_rout1+1)="sgemm_cula"
routine(no_rout1+2)="sgemv_cula"
routine(no_rout1+3)="sgesv_cula"
routine(no_rout1+4)="sgeqrf_cula"
routine(no_rout1+5)="sgemm_cula_c"
routine(no_rout1+6)="sgemv_cula_c"
routine(no_rout1+7)="dgemm_cula"
routine(no_rout1+8)="dgemv_cula"
routine(no_rout1+9)="sgesv_cula_c"
routine(no_rout1+10)="dgemm_cula_c"
routine(no_rout1+11)="sgeqrf_cula_c"
  no_rout=no_rout1+11
#endif
write(*,*) 'usage: gpu_test.x  <routine>  <n> <verbosity level>'
! using own C-function trim 
write(*,*) 'routine: ',(trim(routine(i)),',',i=1,no_rout)
end subroutine usage

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
CALL GETARG(3,BUFFER); READ(BUFFER,*,IOSTAT=ios) verb
IF (ios.NE.0) THEN
  write(*,'(1x,a,i1)') 'using default verbosity level:',verb
ELSE
  write(*,'(1x,a,i2)') 'online argument read in - verbosity level:',verb
ENDIF
end subroutine read_arguments

subroutine sgemm_cuda_test
write(*,'(2x,a)') 'cuda_sgemm GPU vs CPU test'
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
  B_r(i,j)=float(i+j) 
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
call print_time("GPU CUDA SGEMM")
call mtx_check(C_r)
stat=cublas_free(devPtrA)
stat=cublas_free(devPtrB)
stat=cublas_free(devPtrC)
call cublas_shutdown
! now CPU sgemm...
call system_clock( c1, cr, cm )
call sgemm('n','n',n,n,n,1.0e0,A_r,n,B_r,n,-2.0e0,C_r,n)
call system_clock( count=c2 )
call print_time("CPU BLAS SGEMM")
call mtx_check(C_r)
call cublas_shutdown
end subroutine sgemm_cuda_test

subroutine sgemm_cula_test
#if defined USE_CULA
write(*,'(2x,a)') 'cula_sgemm GPU vs CPU test'
allocate (A_r(n,n),B_r(n,n),C_r(n,n),STAT=AllocateStatus)
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
   A_r(i,j)=1.0e0
 else
   A_r(i,j)=0.0e0
 endif
 B_r(i,j)=dfloat(i+j) ! B_d is symmetric matrix
 C_r(i,j)=0.0e0
enddo
enddo
! start the time measurement
call system_clock( c1, cr, cm )
stat = CULA_INITIALIZE(); call check_cula_status(stat)
stat = CULA_SGEMM('n','n',n,n,n,1.0e0,A_r,n,B_r,n,1.0e0,C_r,n); call check_cula_status(stat)
call system_clock( count=c2 )
call print_time("GPU CULA_SGEMM")
stat = CULA_SHUTDOWN(); call check_cula_status(stat)
! ---
call mtx_check(C_r)
! now CPU sgemm...
call system_clock( c1, cr, cm )
call sgemm('n','n',n,n,n,1.0e0,A_r,n,B_r,n,-2.0e0,C_r,n)
call system_clock( count=c2 )
call print_time("CPU BLAS SGEMM")
call mtx_check(C_r)
#else
  print *,'CULA library not active !'
#endif
end subroutine sgemm_cula_test

subroutine dgemm_cula_test
#if defined USE_CULA
write(*,'(2x,a)') 'cula_dgemm GPU vs CPU test'
allocate (A_d(n,n),B_d(n,n),C_d(n,n),STAT=AllocateStatus)
if (AllocateStatus.ne.0) then
   stop "error in dgemm_cula_test allocations !"
else
   if (verb.ge.1) then
     write(*,*) "A_d,B_d,C_d matrixes allocated in CPU memory for dgemm_cula_test"
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
! start the time measurement
call system_clock( c1, cr, cm )
stat = CULA_INITIALIZE(); call check_cula_status(stat)
stat = CULA_DGEMM('n','n',n,n,n,1.0d0,A_d,n,B_d,n,1.0e0,C_d,n); call check_cula_status(stat)
call system_clock( count=c2 )
call print_time("GPU CULA_DGEMM")
stat = CULA_SHUTDOWN(); call check_cula_status(stat)
! ---
call mtx_check_d(C_d)
! now CPU sgemm...
call system_clock( c1, cr, cm )
call dgemm('n','n',n,n,n,1.0d0,A_d,n,B_d,n,-2.0e0,C_d,n)
call system_clock( count=c2 )
call print_time("CPU BLAS DGEMM")
call mtx_check_d(C_d)
#else
  print *,'CULA library not active !'
#endif
end subroutine dgemm_cula_test


subroutine dgemm_cuda_test
write(*,'(2x,a)') 'cuda_dgemm GPU vs CPU test'
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
call print_time("GPU CUDA DGEMM")
call mtx_check_d(C_d)
stat=cublas_free(devPtrA)
stat=cublas_free(devPtrB)
stat=cublas_free(devPtrC)
call cublas_shutdown
! now CPU sgemm...
call system_clock( c1, cr, cm )
call dgemm('n','n',n,n,n,1.0d0,A_d,n,B_d,n,-2.0e0,C_d,n)
call system_clock( count=c2 )
call print_time("CPU BLAS DGEMM")
call mtx_check_d(C_d)
call cublas_shutdown
end subroutine dgemm_cuda_test

subroutine sgemv_cuda_test
write(*,'(2x,a)') 'cuda_sgemv GPU vs CPU test'
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
  enddo
  x_r(i)=float(1)
  y_r(i)=0.0e0
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
! get the vector from GPU (devPtry) into y_r...
stat=cublas_get_vector(n,size_of_real,devPtry,1,y_r,1)
call system_clock( count=c2 )
call print_time("GPU sgemv")
write(*,*) 'GPU cublas_sgemv vector control output: sum,norm=',sum(y_r),snrm2(n,y_r,1)
stat=cublas_free(devPtrA)
stat=cublas_free(devPtrx)
stat=cublas_free(devPtry)
call system_clock( c1, cr, cm )
!now perform CPU lib routine:  y = beta*y+alpha*A*x
call sgemv('n',n,n,-2.0e0,A_r,n,x_r,1,1.0e0,y_r,1)
call system_clock( count=c2 )
call print_time("CPU sgemv")
write(*,*) 'CPU sgemv vector control output: sum,norm=',sum(y_r),snrm2(n,y_r,1)
end subroutine sgemv_cuda_test

subroutine sgemv_cula_test
#if defined USE_CULA
write(*,'(2x,a)') 'cula_sgemv GPU vs CPU test'
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
enddo
 x_r(i)=float(1)
 y_r(i)=0.0e0
enddo
call system_clock( c1, cr, cm )
stat = CULA_INITIALIZE(); call check_cula_status(stat)
stat=cula_sgemv('n',n,n,1.0e0,A_r,n,x_r,1,0.0e0,y_r,1)
call system_clock( count=c2 )
call print_time("GPU CULA SGEMV")
stat = CULA_SHUTDOWN(); call check_cula_status(stat)
write(*,*) 'GPU cula_sgemv vector control output: sum,norm=',sum(y_r),snrm2(n,y_r,1)
! ---
call system_clock( c1, cr, cm )
!now perform CPU lib routine:  y = beta*y+alpha*A*x
call sgemv('n',n,n,-2.0e0,A_r,n,x_r,1,1.0e0,y_r,1)
call system_clock( count=c2 )
call print_time("CPU BLAS SGEMV")
write(*,*) 'CPU sgemv vector control output: sum,norm=',sum(y_r),snrm2(n,y_r,1)
#else
 print *,'sorry, CULA library is not active'
#endif
end subroutine sgemv_cula_test

subroutine dgemv_cula_test
#if defined USE_CULA
write(*,'(2x,a)') 'cula_dgemv GPU vs CPU test'
allocate (A_d(n,n),x_d(n),y_d(n),STAT=AllocateStatus)
if (AllocateStatus.ne.0) then
  stop "error in matrix,vectors allocations in dgemv_cula_test !"
else
  if (verb.ge.1) then
    write(*,*) "A_d,x_d,y_d objects allocated in CPU memory"
  endif
endif
! fill allocated objects
do i=1,n
do j=1,n
 if (i.eq.j) then
   A_d(i,j)=1.0e0
 else
   A_d(i,j)=0.0e0
 endif
enddo
 x_d(i)=dfloat(1)
 y_d(i)=0.0d0
enddo
call system_clock( c1, cr, cm )
stat = CULA_INITIALIZE(); call check_cula_status(stat)
stat=cula_dgemv('n',n,n,1.0d0,A_d,n,x_d,1,0.0d0,y_d,1); call check_cula_status(stat)
call system_clock( count=c2 )
call print_time("GPU CULA SGEMV")
stat = CULA_SHUTDOWN(); call check_cula_status(stat)
write(*,*) 'GPU cula_dgemv vector control output: sum,norm=',sum(y_d),dnrm2(n,y_d,1)
! ---
call system_clock( c1, cr, cm )
!now perform CPU lib routine:  y = beta*y+alpha*A*x
call dgemv('n',n,n,-2.0d0,A_d,n,x_d,1,1.0d0,y_d,1)
call system_clock( count=c2 )
call print_time("CPU BLAS DGEMV")
write(*,*) 'CPU dgemv vector control output: sum,norm=',sum(y_d),dnrm2(n,y_d,1)
#else
 print *,'sorry, CULA library is not active'
#endif
end subroutine dgemv_cula_test




subroutine sgeqrf_test
#if defined USE_CULA
! test on matrix QR decomposition ....
  real(kind=size_of_real) :: diffs = 0.0e0
#if defined VAR_IFORT
  real(kind=size_of_real), external :: rand
#endif
  write(*,'(2x,a)') 'cula_sgeqrf GPU vs lapack CPU test'
  allocate (A_r(n,n),B_r(n,n),STAT=AllocateStatus)
  if (AllocateStatus.ne.0) then
    stop "error in allocations of A_r,B_r !"
  else
    if (verb.ge.1) then
      write(*,*) "A_r, B_r matrixes allocated in CPU memory"
    endif
  endif
  allocate(WORK_r(n),STAT=AllocateStatus)
  if (AllocateStatus.ne.0) then
    stop "error in WORK_r array allocation"
  else
    if (verb.ge.1) then
      print *,"allocated array, WORK_r(LWORK)"
    endif
  endif
  allocate(TAU_r(n),TAU1_r(n),STAT=AllocateStatus)
  if (AllocateStatus.ne.0) then
    stop "error in TAU_r, TAU1_r array allocation"
  else
    if (verb.ge.1) then
      print *,"allocated auxiliary array, TAU_r(N),N=",n
    endif
  endif
! fill allocated objects A_r, B_r with (random?) data for later checking
  do i=1,n
  do j=1,n
     !miro: CPU-lapack SGEQRF of MKL library is giving wrong results !!! (found on 31.3.2012)
     ! A_r(i,j) = rand(0)
     A_r(i,j) = float(i)+float(j) 
     B_r(i,j)=A_r(i,j)
  enddo
  enddo
  ! ... start the time measurement
  call system_clock( c1, cr, cm )
  stat = CULA_INITIALIZE(); call check_cula_status(stat)
  stat = CULA_SGEQRF(N, N, A_r, N, TAU_r); call check_cula_status(stat)

  call system_clock( count=c2 )
  call print_time("GPU CULA_SGEQRF")

  stat = CULA_SHUTDOWN(); call check_cula_status(stat)
! time for lapack SGEQRF
!see  http://www.netlib.org/lapack/single/sgeqrf.f
! SUBROUTINE SGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
  CALL SGEQRF(N, N, B_r, N, TAU1_r, WORK_r, N , stat)
  IF (stat.ne.0) then
    stop "error in CPU LAPACK SGEQRF..."
  else
    if (verb.ge.1) then
      print *,"LAPACK SGEQRF finished successfully ...."
    endif
  endif
call system_clock( count=c2 )
call print_time("CPU LAPACK SGEQRF")
! check numbers in A_r, b_r - must be the same...
  if (verb.ge.3) then
    print *,'          QR results:     GPU-SGEQRF    CPU-SGEQRF    diff'
  endif
  do i=1,n
  do j=1,n
    diffs = diffs + abs(A_r(i,j)-B_r(i,j))  
    if (verb.ge.3) then
      print *,i,j,A_r(i,j), B_r(i,j), A_r(i,j)-B_r(i,j)
    endif
  enddo
    diffs =  diffs + abs(TAU_r(i)-TAU1_r(i))
    if (verb.ge.3) then
      print *,i,TAU_r(i),TAU1_r(i),TAU_r(i)-TAU1_r(i)
    endif
  enddo
  write(*,*) 'CULA_SGEQRF vs lapack SGEQRF averaged diffs are ',diffs/float(n*n)
#else
 print *,'sorry, CULA library is not active'
#endif
end subroutine sgeqrf_test

subroutine sgesv_test
#if defined USE_CULA
#if defined VAR_IFORT
  real(kind=size_of_real), external :: rand
#endif
  real(kind=size_of_real) :: diffs = 0.0e0
  write(*,'(2x,a)') 'cula_sgesv GPU vs lapack CPU test'
  allocate (A_r(n,n),B_r(n,n),STAT=AllocateStatus)
  if (AllocateStatus.ne.0) then
    stop "error in allocations of A_r,B_r !"
  else
    if (verb.ge.1) then
      write(*,*) "A_r, B_r matrixes allocated in CPU memory"
    endif
  endif
  allocate(TAU_r(n),TAU1_r(n),STAT=AllocateStatus)
  if (AllocateStatus.ne.0) then
    stop "error in TAU_r, TAU1_r array allocation"
  else
    if (verb.ge.1) then
      print *,"allocated arrays for solving linear equations"
    endif
  endif
  allocate(ipiv(n),ipiv1(n),STAT=AllocateStatus)
  if (AllocateStatus.ne.0) then
    stop "error in ipiv, ipiv1 integer arrays allocation"
  else
    if (verb.ge.1) then
      print *,"allocated integer arrays ipiv, ipiv1 for solving linear equations"
    endif
  endif
! fill allocated objects A_r, B_r, TAU1_r, TAU_r with random data for later checking
  do i=1,n
  do j=1,n
     A_r(i,j) = rand(0);  B_r(i,j)=A_r(i,j)
  enddo
     TAU_r(i) = rand(0);  TAU1_r(i) = TAU_r(i) 
  enddo
  ! ... start the time measurement
  call system_clock( c1, cr, cm )
  stat = CULA_INITIALIZE(); call check_cula_status(stat)
  stat = cula_sgesv(n, 1 , A_r, n , ipiv, TAU_r, n, stat ); call check_cula_status(stat)

  call system_clock( count=c2 )
  call print_time("GPU CULA_SGESV")

  stat = CULA_SHUTDOWN(); call check_cula_status(stat)
  ! now lapack sgesv...
  call sgesv(n, 1 , B_r, n , ipiv1, TAU1_r, n, stat )
  IF (stat.ne.0) then
    stop "error in CPU LAPACK SGESV ..."
  else
    if (verb.ge.1) then
      print *,"LAPACK SGESV finished successfully ..."
    endif
  endif
  call system_clock( count=c2 )
  call print_time("CPU LAPACK SGESV")
  ! check  results - roots - must be identical 
  if (verb.ge.3) then
    print *," roots and their differences  "
  endif
  do i=1,n
    diffs = diffs + abs(TAU_r(i)-TAU1_r(i))  
    if (verb.ge.3) then
      print *,i,TAU_r(i),TAU1_r(i),TAU_r(i)-TAU1_r(i)
    endif
  enddo
  write(*,*) 'control of results - CULA_SGESV vs lapack SGESV averaged diffs are ',diffs/float(n)
#else
 print *,'sorry, CULA library is not active'
#endif
end subroutine sgesv_test

#if defined USE_CULA
SUBROUTINE check_cula_status(STATUS)
! TOOD: make it more general wrt cula&lapack routines exit codes
         INTEGER, intent(in) :: STATUS
         INTEGER INFO
         external :: CULA_GETERRORINFO
         INTEGER  :: CULA_GETERRORINFO
         INFO = CULA_GETERRORINFO()
         IF (STATUS .NE. 0) THEN
            write(*,'(2x,a,$)') "CULA system status: "
            IF (STATUS .EQ. 7) THEN
!*              culaArgumentError
               WRITE(*,*) 'Invalid value for parameter ', INFO
            ELSE IF (STATUS .EQ. 8) THEN
!*              culaDataError
               WRITE(*,*) 'Data error (', INFO ,')'
            ELSE IF (STATUS .EQ. 9) THEN
!*              culaBlasError
               WRITE(*,*) 'Blas error (', INFO ,')'
            ELSE IF (STATUS .EQ. 10) THEN
!*              culaRuntimeError
               WRITE(*,*) 'Runtime error (', INFO ,')'
            ELSE
!*              others
               call CULA_GETSTATUSSTRING(STATUS)
            ENDIF
            STOP 1
         END IF
END subroutine check_cula_status
#endif

subroutine mtx_check(Mtx)
real*4, intent(in) :: Mtx(n,n)
real*4 ::  diag
diag=0.0e0
do i=1,n
  diag = diag + Mtx(i,i)
enddo
write(*,*)  'aver. sum of diag:',diag/float(n),'1,1:',Mtx(1,1),'1,n:',Mtx(1,n),'n,1:',Mtx(n,1),'n,n',Mtx(n,n)
end subroutine mtx_check

subroutine mtx_check_d(Mtx)
real*8, intent(in) :: Mtx(n,n)
real*8 ::  diag
diag=0.0d0
do i=1,n
  diag = diag + Mtx(i,i)
enddo
write(*,*)  'aver. sum of diag:',diag/dfloat(n),'1,1:',Mtx(1,1),'1,n:',Mtx(1,n),'n,1:',Mtx(n,1),'n,n',Mtx(n,n)
end subroutine mtx_check_d

end module gpu_stuff

!------------------------------------------------------
Program GPU_math_tests
  use gpu_stuff
  use fortran_timing

  call fortran_compiler_info
  call c_compiler_info

  call CudaAbouts
  call CulaAbouts

  call print_time_info

  call read_arguments
  select case (method)
    case("sgemm_cuda")
      call sgemm_cuda_test
    case("sgemm_cuda_c")
      call  c_gpu_math_tests(n,method,verb)
    case("sgemm_cula")
      call sgemm_cula_test
    case("sgemm_cula_c")
      call c_gpu_math_tests(n,method,verb)
    case("dgemm_cuda")
      call dgemm_cuda_test
    case("dgemm_cula")
      call dgemm_cula_test
    case("sgemv_cuda")
      call sgemv_cuda_test
    case("sgemv_cuda_c")
      call  c_gpu_math_tests(n,method,verb)
    case("sgemv_cula")
      call sgemv_cula_test
    case("dgemv_cula")
      call dgemv_cula_test
    case("sgemv_cula_c")
      call  c_gpu_math_tests(n,method,verb)
    case("sgeqrf_cula")
      call sgeqrf_test
    case("sgeqrf_cula_c")
      call  c_gpu_math_tests(n,method,verb)
    case("sgesv_cula")
      call sgesv_test
    case("sgesv_cula_c")
      call  c_gpu_math_tests(n,method,verb)
    case default
      write(*,*) 'no method for GPU math testing selected'
      call usage
      stop
  end select
End Program GPU_math_tests
