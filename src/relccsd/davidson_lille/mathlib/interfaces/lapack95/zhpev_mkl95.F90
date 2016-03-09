subroutine zhpev_mkl95(a,w,uplo,z,info)
!
!***********************************************************************
! fortran77 call:
! call zhpev(jobz,uplo,n,a,w,z,ldz,work,rwork,info)
!***********************************************************************
!
!  use mkl77_lapack, only: mkl77_hpev, mkl77_xerbla
   implicit none
!     ----- parameter -----
   integer, parameter :: WP=kind(1.0D+00)
!     ----- arguments -----
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   real(WP), intent(out), dimension (:) :: w
   complex(WP), intent(inout), dimension (:) :: a
   complex(WP), intent(out), optional, target, dimension (:,:) :: z
!     ----- local variables -----
   character(len=4), parameter :: srname='hpev'
   character(len=1) :: o_uplo,jobz
   integer :: o_info,n,ldz,l_stat_alloc,l_stat_dealloc,l_nn
#ifdef _ESSL_
   integer :: jobess
#endif
!     ----- allocatable arrays -----
   complex(WP), pointer, dimension (:,:) :: o_z
   complex(WP), pointer, dimension (:) :: work
   real(WP), pointer, dimension (:) :: rwork
   complex(WP), target :: l_a2_comp(1,1)
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
!     packed matrix "a(n*(n+1)/2)", therefore: n=(-1+sqrt(1+8*(size(a))))/2 -----
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
      o_z => z
   else
      ldz = 1
      o_z => l_a2_comp
   endif
!     ----- init allocate status -----
   l_stat_alloc = 0
!     ----- allocate local and work arrays -----
   if(present(z)) then
      o_z => z
   else
      o_z => l_a2_comp
   endif
   allocate(rwork(max(1,3*n-2)),stat=l_stat_alloc)
   if(l_stat_alloc.eq.0) then
      allocate(work(max(1,2*n-1)),stat=l_stat_alloc)
   endif
!     ----- call lapack77 routine -----
   if(l_stat_alloc.eq.0) then
      call zhpev(jobz,o_uplo,n,a,w,o_z,ldz,work,rwork,o_info)
   else
      o_info = -1000
   endif
!     ----- deallocate local and work arrays -----
   deallocate(rwork,stat=l_stat_dealloc)
   deallocate(work,stat=l_stat_dealloc)
#endif
#ifdef _ESSL_
   l_stat_alloc = 0
   o_info = 0
   if(present(uplo)) then
      jobess = 20
   else
      jobess = 0
   endif
   if(present(z)) then
      jobess = jobess+1
      ldz = max(1,size(z,1))
      o_z => z
      allocate(rwork(4*n),stat=l_stat_alloc)
   else
      ldz = 1
      o_z => l_a2_comp
      allocate(rwork(3*n),stat=l_stat_alloc)
   endif
!     ----- init allocate status -----
!     ----- allocate local and work arrays -----
   if(l_stat_alloc.eq.0) then
      call zhpev(jobess,a,w,o_z,ldz,n,work,size(work))
      deallocate(rwork)
   else
      o_info = -1000
   endif
#endif
!     ----- error handler -----
   if(present(info)) then
      info = o_info
   else if(o_info.le.-1000) then
      call xerbla(srname,-o_info)
   endif
end subroutine zhpev_mkl95
