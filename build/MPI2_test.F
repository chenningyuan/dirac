!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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

      program MPI2_TEST
      implicit real*8 (A-H,O-Z)
      include 'mpif.h'
      INTEGER MPIERR
      INTEGER ISTATUS(MPI_STATUS_SIZE)

      INTEGER MYID, NELEMENT
      INTEGER NPROC, NAMELENGTH
      INTEGER AMATL, XMATL, MY_X_WIN
      PARAMETER (mx_size_m=100 000 000)
#if ! defined (SYS_AIX) && ! defined (SYS_IBM)
      INTEGER(KIND=MPI_ADDRESS_KIND) ptr1, ptr2
#endif
      POINTER (ptr1, xmat(mx_size_m))
      POINTER (ptr2, ymat(mx_size_m))
      INTEGER(kind=MPI_ADDRESS_KIND) MY_MM_NEEDSSCR, MY_yy_NEEDSSCR
      INTEGER*8 NELEMENT8, MY_yy_NEEDS, MY_MM_NEEDS
      INTEGER ISSIZE_R8, I_ONE, I_FILE_LEN, MY_MPI_FH
      CHARACTER*72 MACHINENAME
      CHARACTER*9 MY_FILE_NM_L
      CHARACTER*11 MY_FILE_NM
C
      NAMELENGTH = 0
      NELEMENT   = 0
      MYID       = 0
      NPROC      = 0
      ISSIZE_R8  = 0
C
      CALL MPI_INIT(MPIERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID,  MPIERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, MPIERR)
C
      CALL MPI_GET_PROCESSOR_NAME(MACHINENAME,NAMELENGTH,MPIERR)
C
C     write my hostname
C
      write(*,*) ' Processor ',MYID,' is running on host ',
     & MACHINENAME(1:NAMELENGTH)
C
      MY_MM_NEEDS    = 0
      MY_MM_NEEDSSCR = 0
C
      NELEMENT = 10 000 000
C
C     determine number of elements NELEMENT REAL*8 (8 byte long)
C
      NELEMENT8 = NELEMENT
C
C     real*8 in bytes
      CALL MPI_TYPE_EXTENT(MPI_REAL8,ISSIZE_R8,MPIERR)
C
C     memory allocation (in bytes): NELEMENT8 * byte-length
C
      MY_MM_NEEDSSCR = ISSIZE_R8 * NELEMENT8
C
      MY_MM_NEEDS    = MY_MM_NEEDSSCR
      MY_yy_NEEDS    = MY_MM_NEEDS
C
C     allocate memory
C
      CALL MPI_ALLOC_MEM(MY_MM_NEEDS,MPI_INFO_NULL,ptr1,MPIerr)
      CALL MPI_ALLOC_MEM(MY_yy_NEEDS,MPI_INFO_NULL,ptr2,MPIerr)
C
C     fill matrices
C
      DO JI = 1, NELEMENT
        xmat(JI) = 1.0D0
      END DO
C
C     open memory window
C
      CALL MPI_WIN_CREATE(xmat,MY_MM_NEEDS,8,MPI_INFO_NULL,
     &                    MPI_COMM_WORLD,MY_X_WIN,MPIERR)
C
C     lock window (MPI_LOCK_SHARED mode)
C
C     select target ==> if itarget == myid: no 1-sided communication 
      ITARGET = MYID + 1
      IF( ITARGET .ge. NPROC ) ITARGET = 0
C
      CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,ITARGET,MPI_MODE_NOCHECK,
     &                  MY_X_WIN,MPIERR)
C
C     transfer data from target ITARGET
C
      JCOUNT_T = NELEMENT8
C
      JCOUNT   = JCOUNT_T
      CALL MPI_GET(YMAT,JCOUNT,MPI_REAL8,ITARGET,0,JCOUNT_T,
     &              MPI_REAL8, MY_X_WIN, MPIERR)
C
C     unlock
C
      CALL MPI_WIN_UNLOCK(ITARGET,MY_X_WIN,MPIERR)
C
C     free memory window
C
      CALL MPI_WIN_FREE(MY_X_WIN,MPIERR)
C
C     de-allocate
C
      CALL MPI_FREE_MEM(xmat,MPIERR)
      CALL MPI_FREE_MEM(ymat,MPIERR)
C
C     test MPI file handling
C
      MY_FILE_NM_L = 'localfile'
      I_ONE        = 1
C
      WRITE (MY_FILE_NM,'(A9,A1,I1)') MY_FILE_NM_L,'.',I_ONE
C
      I_FILE_LEN = 11
C
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD, MY_FILE_NM(1:I_FILE_LEN),
     &                   MPI_MODE_WRONLY + MPI_MODE_CREATE,
     &                   MPI_INFO_NULL, MY_MPI_FH, MPIERR)
C
      CALL MPI_FILE_CLOSE(MY_MPI_FH,MPIERR)
C
      call MPI_FINALIZE(MPIERR)

      end program
