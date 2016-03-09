subroutine sgeev_mkl95(a,wr,wi,vl,vr,info)
!
!***********************************************************************
! fortran77 call:
! call sgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
!***********************************************************************
!
!  use dlapack, only: mkl77_geev, mkl77_xerbla
   implicit none
!     ----- parameter -----
   integer, parameter :: WP=kind(1.0E+00)
   double precision, parameter :: ZERO=0.0E+00,ONE=1.0E+00
!     ----- arguments -----
   integer, intent(out), optional :: info
   real(WP), intent(inout), dimension (:,:) :: a
   real(WP), intent(out), dimension (:) :: wr,wi
   real(WP), intent(out), optional, target, dimension (:,:) :: vl,vr
!     ----- local variables -----
   character(len=4), parameter :: srname='geev'
   character(len=1) :: jobvl,jobvr
   integer :: n,lda,ldvl,ldvr,lwork,l_stat_alloc, l_stat_dealloc,o_info
#ifdef _ESSL_
   integer :: jobess,jc,i,k
   real (WP) :: dumr,dumi,ar,ai,eps
   complex (WP) :: dum
#endif
   real(WP), pointer, dimension (:) :: work
   real(WP), pointer, dimension (:,:) :: o_vl,o_vr
   real(WP) :: s_work(1)
   real(WP), target :: l_a2_real(1,1)
#ifdef _ESSL_
   complex(WP), allocatable, dimension (:) :: w
   complex(WP), allocatable, dimension (:,:) :: vvr
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
   if(present(vl)) then
      o_vl => vl
   else
      o_vl => l_a2_real
   endif
   if(present(vr)) then
      o_vr => vr
   else
      o_vr => l_a2_real
   endif
!     ----- request work array(s) size -----
   lwork = -1
   call sgeev(jobvl,jobvr,n,a,lda,wr,wi,o_vl,ldvl,o_vr,ldvr,s_work,lwork,o_info)
!     ----- exit if error: bad parameters -----
   if(o_info.ne.0) then
      if(present(info)) then
         info = o_info
      else if(o_info.le.-1000) then
         call xerbla(srname,-o_info)
      endif
   endif
   lwork = s_work(1)
!     ----- allocate local and work arrays -----
   allocate(work(lwork),stat=l_stat_alloc)
!     ----- call lapack77 routine -----
   if(l_stat_alloc.eq.0) then
      call sgeev(jobvl,jobvr,n,a,lda,wr,wi,o_vl,ldvl,o_vr,ldvr,work,lwork,  &
                 o_info)
   else
      o_info = -1000
   endif
#endif
#ifdef _ESSL_
   if(present(vr)) then
      ldvr = max(1,size(vr,1))
      jobess = 1
      allocate(vvr(n,n),w(n),work(2*n),selec(n),stat=l_stat_alloc)
   else
      ldvr = 1
      jobess = 0
      allocate(vvr(2,2),w(n),work(n),selec(2),stat=l_stat_alloc)
   endif
   if(l_stat_alloc.eq.0) then
      lda = max(1,size(a,1))
      o_info = 0
      eps = epsilon(ONE)
      call sgeev(jobess,a,lda,w,o_vr,ldvr,selec,n,work,size(work))
      wr(1:n) = real(w,WP)
      wi(1:n) = aimag(w)
!
!     ----- let's stick to the MKL normalization     -----
!     ----- The computed eigenvectors are normalized -----
!     ----- to have Euclidean norm equal to 1        -----
!     ----- and largest component real               -----
!
      o_info = 0
      jc = 1
      do while(jc.le.n)
         if(abs(wi(jc)).lt.eps) then
            vr(1:n,jc) = dble(vvr(1:n,jc))
            dumr = sqrt(sum(vr(1:n,jc)**2))
            vr(1:n,jc) = (ONE/dumr)*vr(1:n,jc)
            jc = jc+1
         else if(abs(wi(jc)+wi(jc+1)).lt.eps) then
            k = 1
            dumr = abs(vvr(1,jc))
            do i=2,n
               dumi = abs(vvr(i,jc))
               if(dumi.gt.dumr) then
                  dumr = dumi
                  k = i
               endif
            enddo
            dum = conjg(vvr(k,jc))/dumr
            vvr(1:n,jc) = dum*vvr(1:n,jc)
            dumr = sqrt(sum(vvr(1:n,jc)*conjg(vvr(1:n,jc))))
            vvr(1:n,jc) = dcmplx((ONE/dumr),ZERO)*vvr(1:n,jc)
            vr(1:n,jc) = real(vvr(1:n,jc))
            vr(1:n,jc+1) = aimag(vvr(1:n,jc))
            jc = jc+2
         else
            print *,' strange. pas de complexe conjugee'
         endif
      enddo
      deallocate(w,selec)
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
   return
end subroutine sgeev_mkl95
