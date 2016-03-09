subroutine cgeev_mkl95(a,w,vl,vr,info)
!
!***********************************************************************
! fortran77 call:
! call cgeev(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
!***********************************************************************
!
!  use mkl77_lapack, only: mkl77_geev, mkl77_xerbla
   implicit none
!     ----- parameter -----
   integer, parameter :: WP=kind(1.0E+00)
   double precision, parameter :: ZERO=0.0E+00,ONE=1.0E+00
!     ----- arguments -----
   integer, intent(out), optional :: info
   complex(WP), intent(out), dimension (:) :: w
   complex(WP), intent(inout), dimension (:,:) :: a
   complex(WP), intent(out), optional, target, dimension (:,:) :: vl,vr
!     ----- local variables -----
   character(len=4), parameter :: srname='geev'
   character(len=1) :: jobvl,jobvr
   integer :: n,lda,ldvl,ldvr,lwork,l_stat_alloc,l_stat_dealloc,o_info
#ifdef _ESSL_
   integer :: jobess,jc,k,i
   real (WP) :: dumr,dumi
   complex (WP) :: dum
#endif
!     ----- local arrays -----
   real(WP), pointer, dimension (:) :: rwork
   complex(WP), pointer, dimension (:) :: work
   complex(WP), pointer, dimension (:,:) :: o_vl,o_vr
   complex(WP) :: s_work(1)
   complex(WP), target :: l_a2_comp(1,1)
#ifdef _ESSL_
   logical, allocatable, dimension (:) :: selec
#endif
!  intrinsic max, present, size
!
!     +++++
!     start
!     +++++
!
!     ----- init optional and skipped scalars -----
   n = size(a,2)
#ifdef _LAPACK_
   if(present(vl)) then
       jobvl = 'v'
   else
       jobvl = 'n'
   endif
   if(present(vr)) then
       jobvr = 'v'
   else
       jobvr = 'n'
   endif
   lda = max(1,size(a,1))
   if(present(vl)) then
       ldvl = max(1,size(vl,1))
   else
       ldvl = 1
   endif
   if(present(vr)) then
       ldvr = max(1,size(vr,1))
   else
       ldvr = 1
   endif
!     ----- init allocate status -----
   l_stat_alloc = 0
!     ----- allocate local and work arrays -----
   if(present(vl)) then
      o_vl => vl
   else
      o_vl => l_a2_comp
   endif
   if(present(vr)) then
      o_vr => vr
   else
      o_vr => l_a2_comp
   endif
!     ----- request work array(s) size -----
   allocate(rwork(2*n), stat=l_stat_alloc)
   lwork = -1
   call cgeev(jobvl,jobvr,n,a,lda,w,o_vl,ldvl,o_vr,ldvr,s_work,lwork,rwork,  &
              o_info)
!     ----- exit if error: bad parameters -----
   if(o_info.ne.0) then
      if(present(info)) then
         info = o_info
      else if(o_info.le.-1000) then
         call xerbla(srname,-o_info)
         return
      endif
   endif
   lwork = s_work(1)
   if(l_stat_alloc.eq.0) then
       allocate(work(lwork),stat=l_stat_alloc)
   endif
!     ----- call lapack77 routine -----
   if(l_stat_alloc.eq.0) then
      call cgeev(jobvl,jobvr,n,a,lda,w,o_vl,ldvl,o_vr,ldvr,work,lwork,rwork,  &
                 o_info)
   else
      o_info = -1000
   endif
#endif
#ifdef _ESSL_
   if(present(vr)) then
      ldvr = max(1,size(vr,1))
      o_vr => vr
      jobess = 1
      allocate(rwork(3*n),work(2),selec(n),stat=l_stat_alloc)
   else
      ldvr = 1
      o_vr => l_a2_comp
      jobess = 0
      allocate(rwork(2*n),work(2),selec(2),stat=l_stat_alloc)
   endif
   if(l_stat_alloc.eq.0) then
      lda = max(1,size(a,1))
      o_info = 0
      call cgeev(jobess,a,lda,w,o_vr,ldvr,selec,n,rwork,size(rwork))
!
!     ----- let's stick to the MKL normalization     -----
!     ----- The computed eigenvectors are normalized -----
!     ----- to have Euclidean norm equal to 1        -----
!     ----- and largest component real               -----
!
      do jc=1,n
         k = 1
         dumr = abs(vr(1,jc))
         do i=2,n
            dumi = abs(vr(i,jc))
            if(dumi.gt.dumr) then
               dumr = dumi
               k = i
            endif
         enddo
         dum = conjg(vr(k,jc))/dumr
         vr(1:n,jc) = dum*vr(1:n,jc)
         dumr = sqrt(sum(vr(1:n,jc)*conjg(vr(1:n,jc))))
         vr(1:n,jc) = dcmplx((ONE/dumr),ZERO)*vr(1:n,jc)
      enddo
      deallocate(selec)
   else
      o_info = -1000
   endif
#endif
!     ----- deallocate local and work arrays -----
   deallocate(rwork,work,stat=l_stat_dealloc)
!     ----- error handler -----
   if(present(info)) then
      info = o_info
   else if(o_info.le.-1000) then
      call xerbla(srname,-o_info)
   endif
   return
end subroutine cgeev_mkl95
