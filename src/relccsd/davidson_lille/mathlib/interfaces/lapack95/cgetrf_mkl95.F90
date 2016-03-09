subroutine cgetrf_mkl95(a,ipiv,info)
!
!***********************************************************************
! fortran77 call:
! call cgetrf(m,n,a,lda,ipiv,info)
!***********************************************************************
!
!   use mkl77_lapack, only: mkl77_getrf, mkl77_xerbla
   implicit none
!     ----- parameter -----
   integer, parameter :: WP=kind(1.0E+00)
!     ----- arguments -----
   integer, intent(out), optional :: info
   complex(WP), intent(inout), dimension (:,:) :: a
   integer, intent(out), optional, target, dimension (:) :: ipiv
!     ----- local variables -----
   character(len=5), parameter :: srname='getrf'
   integer :: o_info,m,n,lda,l_stat_alloc,l_stat_dealloc
   integer, pointer :: o_ipiv(:)
!     ----- intrinsic functions -----
!  intrinsic max, present, size
!
!     +++++
!     start
!     +++++
!
!     ----- init optional and skipped scalars -----
   lda = max(1,size(a,1))
   m = size(a,1)
   n = size(a,2)
!     ----- init allocate status -----
   l_stat_alloc = 0
!     ----- allocate local and work arrays -----
   if(present(ipiv)) then
      o_ipiv => ipiv
   else
      allocate(o_ipiv(min(m,n)),stat=l_stat_alloc)
   endif
!     ----- call lapack77 routine -----
   if(l_stat_alloc.eq.0) then
      call cgetrf(m,n,a,lda,o_ipiv,o_info)
   else
      o_info = -1000
   endif
!     ----- deallocate local and work arrays -----
   if(.not. present(ipiv)) then
      deallocate(o_ipiv,stat=l_stat_dealloc)
   endif
!     ----- error handler -----
   if(present(info)) then
      info = o_info
   else if(o_info.le.-1000) then
      call xerbla(srname,-o_info)
   endif
   return
end subroutine cgetrf_mkl95
