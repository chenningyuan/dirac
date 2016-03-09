subroutine sspev_mkl95(a,w,uplo,z,info)
!
!***********************************************************************
! mkl fortran77 call:
! sspev(jobz,uplo,n,a,w,z,ldz,work,info)
!***********************************************************************
!
!  use dlapack, only: mkl77_spev, mkl77_xerbla
   implicit none
!     ----- parameter -----
   integer, parameter :: WP=kind(1.0E+00)
!     ----- arguments -----
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   real(WP), intent(out), dimension (:) :: w
   real(WP), intent(inout), dimension (:) :: a
   real(WP), intent(out), optional, target, dimension (:,:) :: z
!     ----- local variables -----
   character(len=4), parameter :: srname='spev'
   character(len=1) :: jobz,o_uplo
   integer :: o_info,n,ldz,l_stat_alloc,l_stat_dealloc,l_nn
#ifdef _ESSL_
   integer :: jobess
#endif
   real(WP), pointer, dimension (:,:) :: o_z
   real(WP), pointer, dimension (:) :: work
   real(WP), target :: l_a2_real(1,1)
!  intrinsic int, max, present, real, size, sqrt
!
!     +++++
!     start
!     +++++
!
!     ----- init optional and skipped scalars -----
   l_nn = size(a)
   if(l_nn.le.0) then
      n = l_nn
   else
! packed matrix "a(n*(n+1)/2)", therefore: n=(-1+8*(size(a)))/2
      n = int((-1.0E+00+sqrt(1.0E+00+8.0E+00*real(l_nn,WP)))*0.5E+00)
   endif
#ifdef _LAPACK_
   if(present(uplo)) then
      o_uplo = uplo
   else
      o_uplo = 'u'
   endif
   if(present(z)) then
      jobz = 'v'
   else
      jobz = 'n'
   endif
   if(present(z)) then
       ldz = max(1,size(z,1))
   else
      ldz = 1
   endif
!     ----- init allocate status -----
   l_stat_alloc = 0
!     ----- allocate local and work arrays -----
   if(present(z)) then
      o_z => z
   else
      o_z => l_a2_real
   endif
   allocate(work(3*n),stat=l_stat_alloc)
!     ----- call lapack77 routine -----
   if(l_stat_alloc==0) then
      call sspev(jobz,o_uplo,n,a,w,o_z,ldz,work,o_info)
   else
      o_info = -1000
   endif
#endif
#ifdef _ESSL_
   o_info = 0
   l_stat_alloc = 0
   if(present(uplo)) then
      jobess = 20
   else
      jobess = 0
   endif
   if(present(z)) then
      jobess = jobess+1
      ldz = max(1,size(z,1))
      o_z => z
      allocate(work(2*n),stat=l_stat_alloc)
   else
      ldz = 1
      o_z => l_a2_real
      allocate(work(n),stat=l_stat_alloc)
   endif
!     ----- init allocate status -----
!     ----- allocate local and work arrays -----
   if(l_stat_alloc==0) then
      call sspev(jobess,a,w,o_z,ldz,n,work,size(work))
   else
      o_info = -1000
   endif
#endif
!     ----- deallocate local and work arrays -----
   deallocate(work,stat=l_stat_dealloc)
!     ----- error handler -----
   if(present(info)) then
      info = o_info
   else if(o_info.le.-1000) then
      call xerbla(srname,-o_info)
   endif
end subroutine sspev_mkl95
