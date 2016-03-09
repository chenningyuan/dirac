subroutine zgetri_mkl95(a,ipiv,info)
!
!***********************************************************************
! fortran77 call:
! call zgetri(n,a,lda,ipiv,work,lwork,info)
!***********************************************************************
!
!  use mkl77_lapack, only: mkl77_getri, mkl77_xerbla
   implicit none
!     ----- parameter -----
   integer, parameter :: WP=kind(1.0D+00)
!     ----- arguments -----
   integer, intent(out), optional :: info
   integer, intent(in), dimension (:) :: ipiv
   complex(WP), intent(inout), dimension (:,:) :: a
!     ----- local variables -----
   character(len=5), parameter :: srname='getri'
   integer :: o_info,n,lda,lwork,l_stat_alloc,l_stat_dealloc
   complex(WP), pointer, dimension (:) :: work
   complex(WP) :: s_work(1)
!     ----- intrinsic functions -----
!  intrinsic max, present, size
!
!     +++++
!     start
!     +++++
!
!     ----- init optional and skipped scalars -----
   lda = max(1,size(a,1))
   n = size(a,2)
!     ----- request work array(s) size -----
   lwork = -1
   call zgetri(n,a,lda,ipiv,s_work,lwork,o_info)
!     ----- exit if error: bad parameters -----
   if(o_info.eq.0) then
      lwork = s_work(1)
!     ----- init allocate status -----
      l_stat_alloc = 0
!     ----- allocate local and work arrays -----
      allocate(work(lwork),stat=l_stat_alloc)
!     ----- call lapack77 routine -----
      if(l_stat_alloc.eq.0) then
         call zgetri(n,a,lda,ipiv,work,lwork,o_info)
      else
          o_info = -1000
      endif
!     ----- deallocate local and work arrays -----
      deallocate(work,stat=l_stat_dealloc)
   endif
!     ----- error handler -----
   if(present(info)) then
      info = o_info
   else if(o_info.le.-1000) then
      call xerbla(srname,-o_info)
   endif
   return
end subroutine zgetri_mkl95
