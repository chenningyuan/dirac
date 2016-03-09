subroutine cgesvd_mkl95(a,s,u,vt,ww,job,info)
!
!*******************************************************************************
! fortran77 call:
!   call cgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
!*******************************************************************************
!
!   use mkl77_lapack, only: mkl77_gesvd, mkl77_xerbla
    implicit none
!
!     ----- parameter -----
!
    character(len=5), parameter :: srname='GESVD'
    integer, parameter :: WP=kind(1.0E+00)
!
!     ----- arguments -----
!
    character(len=1), intent(in), optional :: job
    integer, intent(out), optional :: info
    complex(WP), intent(inout), dimension (:,:) :: a
    complex(WP), intent(out), dimension (:) :: s
    complex(WP), intent(out), optional, target, dimension (:,:) :: u,vt
    complex(WP), intent(out), optional, dimension (:) :: ww
!
!     ----- local variables -----
!
    character(len=1) :: o_job,jobu,jobvt
    character(len=8) :: sv2107
    integer :: o_info,m,n,lda,ldu,ldvt,lwork,l_stat_alloc,l_stat_dealloc
#ifdef _ESSL_
    character(len=8) :: sv2017
    integer :: iopt,ierr1,ierr2,nb,i
#endif
    complex(WP) :: s_work(1)
!
!     ----- alloctable variables -----
!
    real(WP), pointer, dimension (:) :: rwork
    complex(WP), pointer, dimension (:) :: work
    complex(WP), pointer, dimension (:,:) :: o_u,o_vt
!
!     ----- "allocate" optional arrays -----
!
    complex(WP), target :: l_a2_complex(1,1)
    intrinsic max, present, size
!
!     +++++
!     start
!     +++++
!
!     ----- executable statements -----
!     ----- init optional and skipped scalars -----
!
    if(present(job)) then
        o_job = job
    else
        o_job = 'N'
    endif
    lda = max(1,size(a,1))
    if(present(u)) then
        ldu = max(1,size(u,1))
    else
        ldu = 1
    endif
    if(present(vt)) then
        ldvt = max(1,size(vt,1))
    else
        ldvt = 1
    endif
    m = size(a,1)
    n = size(a,2)
    if(present(u)) then
        if(size(u,2).eq.m) then
            jobu = 'A'
        else
            jobu = 'S'
        endif
    else
        if((o_job.eq.'U'.or.o_job.eq.'u')) then
            jobu = 'O'
        else
            jobu = 'N'
        endif
    endif
    if(present(vt)) then
        if(size(vt,1).eq.n) then
            jobvt = 'A'
        else
            jobvt = 'S'
        endif
    else
        if((o_job.eq.'V'.or.o_job.eq.'v')) then
            jobvt = 'O'
        else
            jobvt = 'N'
        endif
    endif
!
!     ----- init allocate status -----
!
    l_stat_alloc = 0
!
!     ----- allocate local and work arrays -----
!
    if(present(u)) then
        o_u => u
    else
        o_u => l_a2_complex
    endif
    if(present(vt)) then
        o_vt => vt
    else
        o_vt => l_a2_complex
    endif
#ifdef _LAPACK_
!
!     ----- request work array(s) size -----
!
    lwork = -1
    allocate(rwork(5*min(m,n)),stat=l_stat_alloc)
    if(l_stat_alloc.ne.0) then
       o_info = -1000
      call xerbla(srname,-o_info)
      return
    endif
    l_stat_alloc = 0
    call cgesvd(jobu,jobvt,m,n,a,lda,s,o_u,ldu,o_vt,ldvt,s_work,lwork,rwork,  &
                o_info)
!
!     ----- exit if error: bad parameters -----
!
    if(o_info.ne.0) then
!
!     ----- error handler -----
!
       if(present(info)) then
           info = o_info
       else if(o_info.le.-1000) then
          call xerbla(srname,-o_info)
       endif
       return
    endif
    lwork = s_work(1)
!
!     ----- allocate work arrays with requested sizes -----
!
    allocate(work(lwork),stat=l_stat_alloc)
!
!     ----- call lapack routine -----
!
    if(l_stat_alloc.eq.0) then
        call cgesvd(jobu,jobvt,m,n,a,lda,s,o_u,ldu,o_vt,ldvt,work,lwork,  &
                    rwork,o_info)
    else
       o_info = -1000
    endif
#endif
#ifdef _ESSL_
     write(6,"(5x,a)") 'sorry: interface for singular value decomposition&
                        & for complex matrices not implemented'
!#   if(jobvt.ne.'N') iopt = 11
!#   if(jobu.ne.'N') iopt = 12
!#   if(mod(iopt,10).ne.2) then
!#      nb = 0
!#   else
!#      nb = m
!#      o_u = 0.0E+00
!#      do i=1,m
!#         o_u(i,i) = 1.0D+00
!#      enddo
!#   endif
!#   if(iopt.eq.0.or.iopt.eq.10) then
!#      lwork = n+max(m,n)
!#   else if(iopt.eq.1.or.iopt.eq.11) then
!#      lwork = 2*n+max(m,n)
!#   else if(iopt.eq.2.or.iopt.eq.12) then
!#      lwork = 2*n+max(m,n,nb)
!#   endif
!#   allocate(work(lwork),stat=l_stat_alloc)
!#   call einfo(0)
!#!  call errsav(2015, sv2015)
!#!  call errset(2015,1,-1,0,0,2015)
!#   call errsav(2107, sv2107)
!#   call errset(2107,1,-1,0,0,2107)
!#   call cgesvf(iopt,a,lda,o_u,ldu,nb,s,m,n,work,lwork)
!#!
!#!     ----- error handler -----
!#!
!#   call einfo(2107,ierr1,ierr2)
!#   if(ierr1.gt.0) then
!#      o_info = min(m,n)-ierr1
!#   else
!#      o_info = 0
!#   endif
!#   if(present(info)) then
!#      info = o_info
!#   else if(o_info.le.-1000) then
!#      write(6,"(5x,a)") 'DGESVF terminated with error ',ierr1
!#   endif
!#   if(jobvt.eq.'A') then
!#      vt(1:n,1:n) = a(1:n,1:n)
!#   else if(jobvt.eq.'S') then
!#      vt(1:min(m,n),1:n) = a(1:min(m,n),1:n)
!#   endif
!#
!#   call errstr(2107,sv2107)
!#!  call errstr(2015,sv2015)
#endif
!
!     ----- set output optional array parameters -----
!
    if(present(ww)) ww = work(2:min(m,n))
!
!     ----- deallocate work arrays with requested sizes -----
!
    deallocate(work,stat=l_stat_dealloc)
!
!     ----- error handler -----
!
#ifdef _LAPACK_
   if(present(info)) then
       info = o_info
   else if(o_info.le.-1000) then
      call xerbla(srname,-o_info)
   endif
#endif
   return
end subroutine cgesvd_mkl95
