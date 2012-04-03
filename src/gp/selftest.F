!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

      subroutine selftest_lapack
c     Call Lapack with some small tests to see if things look ok
      implicit none
      integer N, N2, N3, lwork, liwork, info, idum, M,
     &     idum2,idum3,idum4
#ifdef INT_STAR8
      integer*4 N_4, N2_4, N3_4, lwork_4, liwork_4, info_4, idum_4, M_4,
     &     idum2_4,idum3_4,idum4_4, mone_4, mone2_4
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
      N2_4 = N_4
      N3_4 = N_4
      INFO_4 = -1
      liwork_4 = -1
      mone_4 = -1
      mone2_4 = -1
      CALL DSYEVR('V','A','U',N_4,A,N2_4,DDUM,DDUM2,IDUM_4,
     &     IDUM2_4,0.0D0, IDUM3_4,EIG,VEC,N3_4,IDUM4_4,WORKDUM,
     &     mone_4,liwork_4,mone2_4,INFO_4)
      if (info_4.eq.0.and.liwork_4.gt.0) then
         print *,'WARNING: Self test failed!'
         print *,'Call to LAPACK(DSYEVR) failed with integer*8, but'
         print *,'worked with integer*4 arguments. However, DIRAC is'
         print *,'configured to use integer*8 (64 bit integers).'
         print *,' '
         print *,'You can fix this mismatch by linking to a integer*8'
         print *,'LAPACK library, or recompile DIRAC using integer*4.'
      else
         print *,'ERRROR: Lapack integer*8 test failed with info =',
     &        INFO,'LIWORK=',liwork
         print *,'ERRROR: Lapack integer*4 test failed with info =',
     &        INFO_4,'LIWORK=',liwork_4
      endif
#else
         print *,'ERRROR: Lapack test failed with info =',INFO,
     &        'LIWORK=',liwork
#endif
         print *,'Cannot continue without a working LAPACK, quitting.'
         call quit('LAPACK integer*4/8 test failed')
      endif
      end
C
C     MPI compatibility test
C
#ifdef VAR_MPI
      subroutine selftest_mpi
c     Call MPI routines with some small tests to see if things look ok
      use interface_to_mpi
      implicit real*8 (A-H,O-Z)
#include "infpar.h"
      INTEGER MPIERR, IREQ
      DIMENSION ISTAT(df_MPI_STATUS_SIZE)
C
C     send dummy message (from MASTER to MASTER...)
C
      CALL interface_MPI_ISEND(MYTID,1,0,25,
     &                         global_communicator,IREQ)
      CALL interface_MPI_PROBE(df_MPI_ANY_SOURCE,df_MPI_ANY_TAG,
     &                         global_communicator,ISTAT)
      IF( ISTAT(df_MPI_SOURCE) .ne. MPARID )THEN
         print *,' *** ERROR: MPI self test failed! *** '
         print *,' sender of message should have ID 0'
         print *,' but has ID ',ISTAT(df_MPI_SOURCE)
         print *,' '
         print *,' Make sure that you are using the right mpirun for'
         print *,' your particular MPI library version.'
         print *,' '
#ifdef INT_STAR8
         print *,' DIRAC is configured to use integer*8 '
     &        //'(64 bit integers).'
         print *,' You need to link a MPI library that is compiled'
         print *,' using integer*8, make sure that this is the case.'
#else
         print *,' DIRAC is configured to use integer*4 '
     &        //' (32 bit integers).'
         print *,' You need to link a MPI library that is compiled'
         print *,' using integer*4, make sure that this is the case.'
#endif
         print *,' '
         print *,' Cannot continue without a working MPI library.'
         print *,' Quitting.'
         CALL QUIT(' *** ERROR: MPI self test failed ***')
      END IF
C
C     ... it looks ok, receive dummy message
C
      CALL interface_MPI_IRECV(IMESS,1,0,df_MPI_ANY_TAG,
     &                         global_communicator,IREQ)
      CALL interface_MPI_WAIT(IREQ,ISTAT)
      END 
#endif
