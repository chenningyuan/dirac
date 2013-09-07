program test_xcopy_with_big_array

use memory_allocator
 implicit none
 integer, parameter :: kreal    = 8
 real (kind=kreal), pointer ::     a_v4(:)
 real (kind=kreal), allocatable, target :: a_vt(:)
!---------------------
! COMMON /COMPI/ RCW
! INTEGER RCW
!---------------------
#include "../src/relccsd/complex.inc"
 integer :: NVTMIN 

 RCW = 1

 NVTMIN = 2147483948 ! ! larger than 2^31-1=2147483648
 print *,'NVTMIN=',NVTMIN,' RCW=',RCW
 print *,'a_VT of size NVTMIN*RCW=',NVTMIN*RCW

 call allocator_init

 !allocate(a_VT(NVTMIN*RCW))
 call alloc(a_VT ,NVTMIN*RCW, id="vt" )
 
 a_V4  => a_VT
 call my_rdints(a_V4)

 end program test_xcopy_with_big_array

 subroutine my_rdints(VVVVV)
implicit none
 real*8 :: A0r=0.0d0
 integer :: N = 2147483650 ! larger than 2^31-1
!---------------------
! COMMON /COMPI/ RCW
! INTEGER RCW
!---------------------
#include "../src/relccsd/complex.inc"
 REAL*8 :: VVVVV(RCW,*)
 print *,'rdints: befor call dcopy RCW=',RCW
 print *,' big N=',N
 print *,' A0r=',A0r
 call xcopy(N,A0r,0,VVVVV,1)
 print *,' after dcopy with big N=',N
 return
 end subroutine my_rdints
