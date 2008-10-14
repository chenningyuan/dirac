
      subroutine selftest_lapack
c     Call Lapack with some small tests to see if things look ok
      implicit none
      integer N, N2, N3, lwork, liwork, info, idum, M,
     &     idum2,idum3,idum4
      double precision A,DDUM,DDUM2,EIG,VEC,workdum
      N = 11
      N2 = N
      N3 = N
c     Query dsyevr work workspace size
C     If we use integer*8 for info and dsyevr expects integer*4
C     only parts of info will changed, and we can detect this later.
      INFO = -1
      liwork = -1
      CALL DSYEVR('V','A','U',N,A,N2,DDUM,DDUM2,IDUM,IDUM2,0.0D0,
     &     IDUM3,EIG,VEC,N3,IDUM4,WORKDUM,-1,liwork,-1,INFO)
      if (info.ne.0.or.liwork.lt.0) then
         print *,'ERRROR: Lapack test failed with info =',INFO,
     &        'LIWORK=',liwork
#ifdef INT_STAR8
         print *,'You are using 64 bit integers in Dirac'//
     &        ', make sure your Lapack library does the same!'
#endif
         print *,'Cannot continue without a working Lapack'
         call quit('Lapack test failed')
      endif
      end
