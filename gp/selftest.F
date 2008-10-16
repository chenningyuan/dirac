
      subroutine selftest_lapack
c     Call Lapack with some small tests to see if things look ok
      implicit none
      integer N, N2, N3, lwork, liwork, info, idum, M,
     &     idum2,idum3,idum4
#ifdef INT_STAR8
      integer*4 N_4, N2_4, N3_4, lwork_4, liwork_4, info_4, idum_4, M_4,
     &     idum2_4,idum3_4,idum4_4
#endif
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
#ifdef INT_STAR8
C     integer*8 call failed, try with integer*4
      N_4 = 11
      N2 = N_4
      N3 = N_4
      INFO_4 = -1
      liwork_4 = -1
      CALL DSYEVR('V','A','U',N_4,A,N2_4,DDUM,DDUM2,IDUM_4,
     &     IDUM2_4,0.0D0, IDUM3_4,EIG,VEC,N3_4,IDUM4_4,WORKDUM,
     &     -1,liwork_4,-1,INFO_4)
      if (info_4.eq.0.and.liwork_4.gt.0) then
         print *,'Call to LAPACK(DSYEVR) failed with integer*8,'//
     &        'but worked with integer*4 arguments.'
         print *,'However, DIRAC is configured to use integer*8.'
         print *,'You can fix this mismatch by linking to a integer*8'
         print *,'LAPACK library, or recompile DIRAC using integer*4.'
         print *
      else
         print *,'ERRROR: Lapack integer*8 test failed with info =',
     &        INFO,'LIWORK=',liwork
         print *,'ERRROR: Lapack integer*4 test failed with info =',
     &        INFO,'LIWORK=',liwork
      endif
#else
         print *,'ERRROR: Lapack test failed with info =',INFO,
     &        'LIWORK=',liwork
#endif

         print *,'This test can be disabled in configure.'
         print *,'Cannot continue without a working LAPACK, quitting.'
         call quit('LAPACK test failed')
      endif
      end
