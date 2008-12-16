!-----------------------------------------------------------------------------------
!     1. diagonalize a matrix (symmetric or non-symmetric)
!     2. check for presence of complex eigenvalues (if any)
!     3. deal with according to abort_if_complex flag
!     4. if necessary, sort remaining eigenvalues and
!        eigenvectors in order of increasing eigenvalues
!     5. orthonormalize/biorthogonalize eigenvectors
!
!     sym              => symmetric or non-symmetric
!     n                => dimension of the matrix
!     ldmat            => leading dimension of the matrix
!     matrix(ldmat,n)  => matrix to be diagoanlized
!     nev              => number of real eigenvalues found
!     ev(n)            => first nev entries are eigenvalues
!     rvec(ldmat,n)    => first nev columns are right eigenvectors
!     lvec(ldmat,n)    => first nev columns are left  eigenvectors
!
!     abort_if_complex => if .true. then run is aborted if complex
!                         eigenvalues are found.
!                         if .false. then run continues if complex
!                         eigenvalues are found, they are eliminated
!                         and number of real eigenvalues is
!                         returned in nev
!
!     sort             => sort eigenvalues in the order
!                         of increasing eigenvalues
!
!     biortho          => applies only for non-symmetric case
!
!                         if true then left and right
!                         eigenvectors should be biorthonormal.
!
!                         if false then left and right
!                         eigenvectors should be biorthogonal,
!                         but each left and right vector should
!                         be seperately normalized.
!-----------------------------------------------------------------------------------
      subroutine diagonalize_hscc       &
                 (                      &
                    sym,                &
                    n, ldmat,           &
                    matrix,             &
                    nev,                &
                    ev,                 &
                    rvec,               &
                    lvec,               &
                    abort_if_complex,   &
                    sort,               &
                    biortho             &
                 )
!-----------------------------------------------------------------------------------
!
!          Called from DiagHEFF routine
!
!-----------------------------------------------------------------------------------

      use IOModule
      implicit none
      logical sym
      integer n
      integer ldmat
      double precision  matrix(ldmat,n)
      integer nev
      double precision  ev(n)
      double precision  rvec(ldmat,n)  ! right vector
      double precision  lvec(ldmat,n)  ! left vector
      logical abort_if_complex
      logical sort
      logical biortho

!     ... for dsyev and dgeev ...

      integer info
      double precision  a(ldmat,n)
      double precision  wr(n), wi(n)
      double precision  vr(ldmat,n)
      double precision  vl(ldmat,n)
      integer lwork
      double precision  work(8*n)
      double precision  ovl(n,n)
      integer ipiv(n)

!     ... local variables ...

      integer i, j
      integer ir
      double precision  eps
      double precision  DLAMCH
      double precision  swptmp
      double precision  DDOT
      double precision  DNRM2
      double precision  norm
      double precision  dtmp

      info = 0
      lwork = 8*n

      if (sym) then

         a = matrix
         wr = dble(0)
         call DSYEV('v','l',n,a,ldmat,wr,work,lwork,info)
         if (info /= 0) then
            write(UOUT,*) "diagonalize: Failure in LAPACK routine DSYEV"
            write(UOUT,*) "INFO  = ",info
            call flushout
            call quit('diagonalize: Failure in LAPACK routine DSYEV')
         endif

         ev = wr
         rvec = a
         lvec = a
         nev = n

      else

         a = matrix
         wr = dble(0)
         wi = dble(0)
         vr = dble(0)
         vl = dble(0)

! MI, Nov.2007/TAU
! ... this strange error occurs (opteron amd64/pgf90 5.2.4) when some dirac object
!       files are present in compilation: 
!  ** On entry to DORGQR parameter number  5 had an illegal value
! ... problem is deep in the linked lapack library
!
! .. from library dgeev it calls the dirac routine DORGQR with HURT(i8 problem) parameters
!
!     #1  dorghr_ address: 0x4DCA15 
!=> #0  dorgqr line: "gp_dlapack.F"@7341 address: 0x4AC862 
!     m = 21474836483 , n = 21474836483 , k = 21474836483 , a = 0x6662948, lda = 4 , tau = 0x66629F0, work = 0x6662A10, lwork = 584115552280 , info = 16106122464 
!
!
! system library routine DGEEV calls dirac -i8 compiled routine DORGQR with wrong parameters
!
         call DGEEV('v','v',n,a,ldmat,wr,wi,vl,ldmat,vr,ldmat,work,lwork,info)
         if (info.ne.0) then
            write(UOUT,*)  &
            "diagonalize: Failure in LAPACK routine DGEEV "
            write(UOUT,*) "INFO  = ",info
            call flush(UOUT)
            !call abort
            call Quit('diagonalize: Failure in LAPACK routine DGEEV !')
         endif

!        deal with complex eigenvalues

         eps = DLAMCH('e')
         ir = 0
         do i = 1, n
            if (abs(wi(i)).lt.eps) then
               ir = ir + 1
               ev(ir) = wr(i)
               if (ir.ne.i) then
                  vr(:,ir) = vr(:,i)
                  vl(:,ir) = vl(:,i)
               endif
            endif
         enddo

         nev = ir

         if (abort_if_complex.and.(nev.ne.n)) then
            write(UOUT,*) "Complex eigenvalues during diagonalization"
            write(UOUT,'(a,i3,a,i3)')               &
            '# of real eigevalues found, nev=',nev, &
            ' matrix dimension, n=',n 
            !write(UOUT,*) ( ( wr(i),wi(i) ),i=1,n)
            do i=1,n
              write(UOUT,*)  wr(i), wi(i) 
            enddo
            write(UOUT,*) "...aborting..."
            call flushout
            ! call abort
            call QUIT( 'diagonalize :'//  &
           'Complex eigenvalues during diagonalization...aborting')
         endif

!        sort real eigenvalues in increasing order
!        corresponding eigenvectors are also sorted.

         if (sort) then

            do ir = 1, nev
               do i = ir+1, nev
                  if (ev(i).lt.ev(ir)) then
                     swptmp =  ev(ir)
                     ev(ir) = ev(i)
                     ev(i) = swptmp
                     call DSWAP(n,vr(1,ir),1,vr(1,i),1)
                     call DSWAP(n,vl(1,ir),1,vl(1,i),1)
                  endif
               enddo
            enddo
         endif

!        make sure that left and right eigenvectors form
!        biorthogonal set and are normalized to unity.
!
!        multiply right eigenvectors by the inverse
!        of overlap matrix of left and right eigenvectors.

         call DGEMM('t','n',nev,nev,n,1.d00,vl,ldmat,vr,ldmat,0.0d00,ovl,n)

         call DGETRF(nev,nev,ovl,n,ipiv,info)
         if (info.ne.0) then
            write(UOUT,*) "Failure in LAPACK routine DGETRF"
            write(UOUT,*) "INFO  = ",info
            write(UOUT,*) 'nev=',nev,' n=',n
            write(UOUT,*) 'ipiv=',(ipiv(i),i=1,n)
            do i=1,n
            do j=1,n
              write(UOUT,'(2x,a,i3,a,i3,a,f12.6)') ' ovl(',i,';',j,')=',ovl(i,j)
            enddo           
            enddo           
            call flushout
            call quit('Failure in LAPACK routine DGETRF')
         endif

         call DGETRI(nev,ovl,n,ipiv,work,lwork,info)
         if (info.ne.0) then
            write(UOUT,*) "Failure in LAPACK routine DGETRI"
            write(UOUT,*) "INFO  = ",info
            call flush(UOUT)
            call abort
         endif

!        build new set of left and right eigenvectors

         call DGEMM('n','n',n,nev,nev,1.0d00,vr,ldmat,ovl,n,0.0d00,rvec,ldmat)
         call DCOPY(ldmat*nev,vl,1,lvec,1)

!        normalize right eigenvectors seperately.
!        left eigenvectors are already normalized.

         if (.not.biortho) then
            do i = 1, nev
               norm = DNRM2(n,rvec(1,i),1)
               call DSCAL(n,1.0d00/norm,rvec(1,i),1)
            enddo
         endif

      endif

      return
      end subroutine diagonalize_hscc
