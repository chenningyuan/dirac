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

C aspg: at this point (2007-05-08) i have not yet setup a way yo 
C       automatically enable VAR_PFS in the make system, so if
C       you need this you have to add "-DVAR_PFS" to the
C
C       DEPEND_FLAGS=
C       CPP_FLAGS=
C       CCP_FLAGS=
C
C       variables in Makefile.config, after it has been generated
C

#if defined (VAR_MPI2) 
C
C Allow for the use of a global file (or parallel file system) in
C moltra/tradr6.F and friends.
C This eliminates the message passing for the matrix transposition
C during the MOLTRA step completely. More importantly, it removes
C the synchronization wait time (most probably due to the limited
C number of buffers).
C
C The current implementation uses (partly asynchronous) MPI I/O
C (part of MPI-2), which is at least available on IBM and SGI Altix
C machines.
C
C ***IMPORTANT*** ***IMPORTANT*** ***IMPORTANT*** ***IMPORTANT***
C
C In order to use this modification you have to define VAR_PFS
C (furthermore, we assume that you already have VAR_MPI2 defined).
C
C ***IMPORTANT*** ***IMPORTANT*** ***IMPORTANT*** ***IMPORTANT***
C
C Note (WML TODO?): it is possible to make a more portable implementation
C using Fortran direct access and multiple files (possibly local) and do
C an intermediate more-or-less perfectly load-balanced message passing
C and re-ordering of local files. (Comparable with the current
C SELECT_LOCAL_FROM_GLOBAL routine.)
C
C Implemented as "enabling work" for the TRIPOD project
C (a DECI / DEISA project).
C Written by Walter.Lioen@sara.nl, 2007.
C
      SUBROUTINE GLOBAL_FILE_OPEN_WRITE
      IMPLICIT NONE

      INCLUDE "mpif.h"
#include "bucket_sort.h"
#include "infpar.h"

      INTEGER IERR
      INTEGER GF_INFO
      INTEGER I, J, K

      CALL MPI_INFO_CREATE(GF_INFO, IERR)
#if defined (SYS_AIX) || defined (SYS_IBM)
C     ... IBMs GPFS specific file hint - SK Aug 2008
      CALL MPI_INFO_SET(GF_INFO, "IBM_largeblock_io", "true",IERR)
#else
      CALL MPI_INFO_SET(GF_INFO, "direct_write", "true", IERR)
      CALL MPI_INFO_SET(GF_INFO, "romio_cb_write", "disable", IERR)
      CALL MPI_INFO_SET(GF_INFO, "romio_ds_write", "disable", IERR)
#endif

      CALL MPI_FILE_OPEN(
     +   MPI_COMM_WORLD, GF_NAME,
     +   MPI_MODE_WRONLY + MPI_MODE_CREATE + MPI_MODE_EXCL,
     +   GF_INFO, GF_HANDLE, IERR
     +)
      CALL MPI_INFO_FREE(GF_INFO, IERR)

      DO K = 0, NUMNOD
         DO J = 0, NUMNOD
            GF_NREC(J,K) = 0
            DO I = 1, N_BUFFER
               GF_BREC(I,J,K) = 0
            ENDDO
         ENDDO
      ENDDO

      WRITE_REQUEST = MPI_REQUEST_NULL

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_OPEN_READ
      IMPLICIT NONE

      INCLUDE "mpif.h"
#include "bucket_sort.h"
#include "infpar.h"

      INTEGER IERR
      INTEGER GF_INFO
      INTEGER I, NMPROC, NMBUFF

      NMPROC = NUMNOD + 1
      NMBUFF = N_BUFFER * NMPROC

      CALL MPI_INFO_CREATE(GF_INFO, IERR)
#if defined (SYS_AIX) || defined (SYS_IBM)
C     ... IBMs GPFS specific file hint - SK Aug 2008
      CALL MPI_INFO_SET(GF_INFO, "IBM_largeblock_io", "true",IERR)
#else
      CALL MPI_INFO_SET(GF_INFO, "ind_rd_buffer_size", "16777216", IERR)
      CALL MPI_INFO_SET(GF_INFO, "romio_cb_read", "disable", IERR)
      CALL MPI_INFO_SET(GF_INFO, "romio_ds_read", "enable", IERR)
#endif

      CALL MPI_FILE_OPEN(
     +   MPI_COMM_WORLD, GF_NAME,
     +   MPI_MODE_RDONLY + MPI_MODE_DELETE_ON_CLOSE,
     +   GF_INFO, GF_HANDLE, IERR
     +)

      CALL MPI_INFO_FREE(GF_INFO, IERR)

      DO I=0,NUMNOD
         CALL MPI_BCAST(GF_NREC(0,I),NMPROC,my_MPI_INTEGER,I,
     +                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(GF_BREC(1,0,I),NMBUFF,my_MPI_INTEGER,I,
     +                  MPI_COMM_WORLD,IERR)
      ENDDO

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_CLOSE_WRITE
      IMPLICIT NONE

      INCLUDE "mpif.h"
#include "bucket_sort.h"

      INTEGER IERR

      CALL MPI_WAIT(WRITE_REQUEST,MPI_STATUS_IGNORE,IERR)
      CALL MPI_FILE_CLOSE(GF_HANDLE, IERR)

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_CLOSE_READ
      IMPLICIT NONE

      INCLUDE "mpif.h"
#include "bucket_sort.h"

      INTEGER IERR

      CALL MPI_FILE_CLOSE(GF_HANDLE, IERR)

      RETURN
      END

      INTEGER FUNCTION GLOBAL_FILE_RECORD(
     +   TARGET_NODE, SOURCE_NODE, LOCAL_RECORD)
      IMPLICIT NONE
      INTEGER TARGET_NODE, SOURCE_NODE, LOCAL_RECORD

      INCLUDE "mpif.h"
#include "bucket_sort.h"
#include "infpar.h"
      INTEGER IERR

      IF (TARGET_NODE .LT. 0 .OR. TARGET_NODE .GT. NUMNOD) THEN
         PRINT *, 'ERROR: GLOBAL_FILE_RECORD: TARGET_NODE = ',
     +      TARGET_NODE
         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
      ENDIF
      IF (SOURCE_NODE .LT. 1 .OR. SOURCE_NODE .GT. NUMNOD) THEN
         PRINT *, 'ERROR: GLOBAL_FILE_RECORD: SOURCE_NODE = ',
     +      SOURCE_NODE
         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
      ENDIF
      IF (LOCAL_RECORD .LT. 1 .OR. LOCAL_RECORD .GT. GF_MAXREC) THEN
C increase GF_MAXREC in bucket_sort.h
         PRINT *, 'ERROR: GLOBAL_FILE_RECORD: LOCAL_RECORD = ',
     +      LOCAL_RECORD
         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
      ENDIF

      GLOBAL_FILE_RECORD =
     +   (
     +      (
     +         (
     +            TARGET_NODE * NUMNOD
     +         ) + SOURCE_NODE - 1
     +      ) * GF_MAXREC
     +   ) + LOCAL_RECORD

      RETURN
      END

      FUNCTION GLOBAL_FILE_OFFSET(
     +   TARGET_NODE, SOURCE_NODE, LOCAL_RECORD)
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INTEGER(KIND=MPI_OFFSET_KIND) GLOBAL_FILE_OFFSET
      INTEGER TARGET_NODE, SOURCE_NODE, LOCAL_RECORD

#include "bucket_sort.h"

      INTEGER GLOBAL_FILE_RECORD
      EXTERNAL GLOBAL_FILE_RECORD

C     Prevent INTEGER overflow by splitting the following computation
C     in two parts (the first part performs an implicit cast).
      GLOBAL_FILE_OFFSET =
     +   GLOBAL_FILE_RECORD(TARGET_NODE,SOURCE_NODE,LOCAL_RECORD) - 1
      GLOBAL_FILE_OFFSET = GLOBAL_FILE_OFFSET * GF_RECLEN

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_SEEK(TARGET_NODE,SOURCE_NODE,LOCAL_RECORD)
      IMPLICIT NONE
      INTEGER TARGET_NODE, SOURCE_NODE, LOCAL_RECORD

      INCLUDE "mpif.h"
#include "bucket_sort.h"
      INTEGER(KIND=MPI_OFFSET_KIND) RECORD_OFFSET
      INTEGER IERR

      INTEGER(KIND=MPI_OFFSET_KIND) GLOBAL_FILE_OFFSET
      EXTERNAL GLOBAL_FILE_OFFSET

      RECORD_OFFSET =
     +   GLOBAL_FILE_OFFSET(TARGET_NODE,SOURCE_NODE,LOCAL_RECORD)

      CALL MPI_FILE_SEEK(
     +   GF_HANDLE,RECORD_OFFSET,MPI_SEEK_SET,IERR)

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_WRITE(
     +   RGBUF,IGBUF,N,SOURCE_NODE,TARGET_NODE,IBUF)
      IMPLICIT NONE
      INTEGER N, SOURCE_NODE, TARGET_NODE, IBUF
      REAL*8 RGBUF
      INTEGER IGBUF

      INCLUDE "mpif.h"
#include "bucket_sort.h"

      DIMENSION RGBUF(NGBFSZ)
      DIMENSION IGBUF(2,NGBFSZ)
      INTEGER PREV_REC, THIS_REC, INCOUNT, POSITION, IERR
      INTEGER(KIND=MPI_OFFSET_KIND) RECORD_OFFSET

      INTEGER(KIND=MPI_OFFSET_KIND) GLOBAL_FILE_OFFSET
      EXTERNAL GLOBAL_FILE_OFFSET

      CALL MPI_WAIT(WRITE_REQUEST,MPI_STATUS_IGNORE,IERR)

      PREV_REC = GF_BREC(IBUF,TARGET_NODE,SOURCE_NODE)
      THIS_REC = GF_NREC(TARGET_NODE,SOURCE_NODE) + 1

      POSITION = 0

      INCOUNT = 1
      CALL MPI_PACK(PREV_REC,INCOUNT,my_MPI_INTEGER,
     +  GF_RECBUF,GF_RECLEN,POSITION,MPI_COMM_SELF,IERR)

      INCOUNT = 1
      CALL MPI_PACK(N,INCOUNT,my_MPI_INTEGER,
     +  GF_RECBUF,GF_RECLEN,POSITION,MPI_COMM_SELF,IERR)

      INCOUNT = NGBFSZ
      CALL MPI_PACK(RGBUF,INCOUNT,MPI_DOUBLE_PRECISION,
     +  GF_RECBUF,GF_RECLEN,POSITION,MPI_COMM_SELF,IERR)

      INCOUNT = 2*NGBFSZ
      CALL MPI_PACK(IGBUF,INCOUNT,my_MPI_INTEGER,
     +  GF_RECBUF,GF_RECLEN,POSITION,MPI_COMM_SELF,IERR)

      RECORD_OFFSET =
     +  GLOBAL_FILE_OFFSET(TARGET_NODE,SOURCE_NODE,THIS_REC)

      IERR = 0

      CALL MPI_FILE_IWRITE_AT(
     +  GF_HANDLE,RECORD_OFFSET,GF_RECBUF,GF_RECLEN,MPI_PACKED,
     +  WRITE_REQUEST,IERR)

!     check for successful writing. SK - Aug 2010
!     if it fails you need to change the general record offset 
!     for a given node; for more information consult the 
!     comments in bucket_sort.h.
      IF(IERR .ne. 0)THEN
        PRINT *, 'ERROR: GLOBAL_FILE_WRITE: writing failed at
     +  RECORD_OFFSET (mytid) =',RECORD_OFFSET,SOURCE_NODE
!       errorcode 112
        CALL MPI_ABORT(MPI_COMM_WORLD,112,IERR)
      END IF

      GF_BREC(IBUF,TARGET_NODE,SOURCE_NODE) = THIS_REC
      GF_NREC(TARGET_NODE,SOURCE_NODE) = THIS_REC

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_GET_RECORD(
     +   TARGET_NODE,SOURCE_NODE,LOCAL_RECORD)
      IMPLICIT NONE
      INTEGER TARGET_NODE, SOURCE_NODE, LOCAL_RECORD

      INCLUDE "mpif.h"
#include "bucket_sort.h"
      INTEGER   STATUS(MPI_STATUS_SIZE)
      INTEGER NUMREAD
      INTEGER IERR
      INTEGER(KIND=MPI_OFFSET_KIND) RECORD_OFFSET

      INTEGER(KIND=MPI_OFFSET_KIND) GLOBAL_FILE_OFFSET
      EXTERNAL GLOBAL_FILE_OFFSET

      RECORD_OFFSET = GLOBAL_FILE_OFFSET(
     +  TARGET_NODE,SOURCE_NODE,LOCAL_RECORD)

      CALL MPI_FILE_READ_AT(
     +  GF_HANDLE,RECORD_OFFSET,GF_RECBUF,GF_RECLEN,MPI_PACKED,
     +  STATUS,IERR)
CSK   ... zero NUMREAD, otherwise it may fail
      NUMREAD = 0
      CALL MPI_GET_COUNT(STATUS,MPI_PACKED,NUMREAD,IERR)
      WRITE(6,*) ' GLOBAL_FILE_GET... NUMREAD',NUMREAD
      IF (NUMREAD .NE. GF_RECLEN) THEN
         PRINT *, 'ERROR: GLOBAL_FILE_GET_RECORD: NUMREAD = ', NUMREAD
         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
      ENDIF

      RETURN
      END

      SUBROUTINE GLOBAL_FILE_READ(A,NROW,ICOL,IBUF)
      IMPLICIT NONE
      INTEGER NROW, ICOL,IBUF
      REAL*8 A(NROW)
C
C     Read back the sorted integrals and put them at the right position
C
      INCLUDE "mpif.h"
#include "bucket_sort.h"
#include "infpar.h"
      INTEGER SOURCE_NODE
      REAL*8 RGBUF(NGBFSZ)
      INTEGER IGBUF(2,NGBFSZ)
      INTEGER I, THIS_REC, PREV_REC, N, POSITION, OUTCOUNT, IERR

      DO I=1,NROW
         A(I) = 0D0
      ENDDO

      DO SOURCE_NODE = 1, NUMNOD

         THIS_REC = GF_BREC(IBUF,MYTID,SOURCE_NODE)

         IF (THIS_REC .EQ. 0) GOTO 20

   10    CALL LOCAL_FILE_GET_RECORD(SOURCE_NODE,THIS_REC)

            POSITION = 0

            OUTCOUNT = 1
            CALL MPI_UNPACK(GF_RECBUF,GF_RECLEN,POSITION,
     +         PREV_REC,OUTCOUNT,my_MPI_INTEGER,MPI_COMM_WORLD,IERR)

            OUTCOUNT = 1
            CALL MPI_UNPACK(GF_RECBUF,GF_RECLEN,POSITION,
     +         N,OUTCOUNT,my_MPI_INTEGER,MPI_COMM_WORLD,IERR)

            OUTCOUNT = NGBFSZ
            CALL MPI_UNPACK(GF_RECBUF,GF_RECLEN,POSITION,
     +         RGBUF,OUTCOUNT,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)

            OUTCOUNT = 2*NGBFSZ
            CALL MPI_UNPACK(GF_RECBUF,GF_RECLEN,POSITION,
     +         IGBUF,OUTCOUNT,my_MPI_INTEGER,MPI_COMM_WORLD,IERR)

C           Pick out the integrals that we need.

            IF (N.GT.NGBFSZ) THEN
               PRINT *, 'ERROR: GLOBAL_FILE_READ: N, NGBFSZ =',
     +                  N, NGBFSZ
               CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
            ENDIF
            DO I = 1, N
               IF (IGBUF(2,I).EQ.ICOL) THEN
                  IF (IGBUF(1,I) .LT. 1 .OR. IGBUF(1,I) .GT. NROW) THEN
                     PRINT *, 'ERROR: GLOBAL_FILE_READ: IGBUF(1,I) =',
     +                        IGBUF(1,I),'I =',I
                     CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
                  ENDIF
                  A(IGBUF(1,I)) = RGBUF(I)
               ENDIF
            ENDDO

            THIS_REC = PREV_REC

         IF (THIS_REC .NE. 0) GOTO 10
   20    CONTINUE 

      ENDDO

      RETURN
      END

      SUBROUTINE LOCAL_FILE_OPEN_WRITE
      IMPLICIT NONE

      INCLUDE "mpif.h"
#include "bucket_sort.h"
#include "infpar.h"

      INTEGER IERR
      INTEGER GF_INFO

      WRITE (LF_NAME,'(a,i5.5)') 'local', MYTID

      CALL MPI_INFO_CREATE(GF_INFO, IERR)
#if defined (SYS_AIX) || defined (SYS_IBM)
C     ... IBMs GPFS specific file hint - SK Aug 2008
      CALL MPI_INFO_SET(GF_INFO, "IBM_largeblock_io", "true",IERR)
#else
      CALL MPI_INFO_SET(GF_INFO, "direct_write", "true", IERR)
      CALL MPI_INFO_SET(GF_INFO, "romio_cb_write", "disable", IERR)
      CALL MPI_INFO_SET(GF_INFO, "romio_ds_write", "disable", IERR)
#endif

      CALL MPI_FILE_OPEN(
     +   MPI_COMM_SELF, LF_NAME,
     +   MPI_MODE_WRONLY + MPI_MODE_CREATE,
     +   GF_INFO, LF_HANDLE, IERR
     +)

      CALL MPI_INFO_FREE(GF_INFO, IERR)

      WRITE_REQUEST = MPI_REQUEST_NULL

      RETURN
      END

      SUBROUTINE LOCAL_FILE_CLOSE_WRITE
      IMPLICIT NONE

      INCLUDE "mpif.h"
#include "bucket_sort.h"

      INTEGER IERR

      CALL MPI_WAIT(WRITE_REQUEST,MPI_STATUS_IGNORE,IERR)
      CALL MPI_FILE_CLOSE(LF_HANDLE, IERR)

      RETURN
      END

      SUBROUTINE SELECT_LOCAL_FROM_GLOBAL
      IMPLICIT NONE

      INCLUDE "mpif.h"
#include "bucket_sort.h"
#include "infpar.h"
      CHARACTER LF_RECBUF(GF_RECLEN)
      INTEGER SOURCE_NODE,LOCAL_RECORD
      INTEGER I
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      INTEGER NUMREAD
      INTEGER IERR

      CALL GLOBAL_FILE_OPEN_READ
      CALL LOCAL_FILE_OPEN_WRITE

      DO SOURCE_NODE = 1, NUMNOD

         CALL GLOBAL_FILE_SEEK(MYTID,SOURCE_NODE,1)

         DO LOCAL_RECORD = 1, GF_NREC(MYTID,SOURCE_NODE)

            CALL MPI_FILE_READ(
     +        GF_HANDLE,GF_RECBUF,GF_RECLEN,MPI_PACKED,
     +        ISTAT,IERR)
CSK         ... zero NUMREAD, otherwise it may fail
            NUMREAD = 0
            CALL MPI_GET_COUNT(ISTAT,MPI_PACKED,NUMREAD,IERR)
            IF (NUMREAD .NE. GF_RECLEN) THEN
               PRINT *, 'ERROR: SELECT_LOCAL_FROM_GLOBAL: NUMREAD = ',
     +         NUMREAD
               CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
            ENDIF

            CALL MPI_WAIT(WRITE_REQUEST,MPI_STATUS_IGNORE,IERR)

            DO I = 1, GF_RECLEN
              LF_RECBUF(I) = GF_RECBUF(I)
            ENDDO

            CALL MPI_FILE_IWRITE(
     +        LF_HANDLE,LF_RECBUF,GF_RECLEN,MPI_PACKED,
     +        WRITE_REQUEST,IERR)

         ENDDO

      ENDDO

      CALL LOCAL_FILE_CLOSE_WRITE
      CALL GLOBAL_FILE_CLOSE_READ

      RETURN
      END

      SUBROUTINE LOCAL_FILE_OPEN_READ
      IMPLICIT NONE

      INCLUDE "mpif.h"
#include "bucket_sort.h"
#include "infpar.h"

      INTEGER SOURCE_NODE

      OPEN(UNIT=LF_UNIT,FILE=LF_NAME,
     +     STATUS='OLD',ACTION='READ',
     +     ACCESS='DIRECT',RECL=GF_RECLEN)

      LF_OFFSET(0) = 0
      LF_OFFSET(1) = 0
      DO SOURCE_NODE = 2, NUMNOD
         LF_OFFSET(SOURCE_NODE) = LF_OFFSET(SOURCE_NODE-1) +
     +      GF_NREC(MYTID,SOURCE_NODE-1)
      ENDDO

      RETURN
      END

      INTEGER FUNCTION LOCAL_FILE_RECORD(SOURCE_NODE, LOCAL_RECORD)
      IMPLICIT NONE
      INTEGER SOURCE_NODE, LOCAL_RECORD

      INCLUDE "mpif.h"
#include "bucket_sort.h"
#include "infpar.h"
      INTEGER IERR

      IF (SOURCE_NODE .LT. 1 .OR. SOURCE_NODE .GT. NUMNOD) THEN
         PRINT *, 'ERROR: LOCAL_FILE_RECORD: SOURCE_NODE = ',
     +      SOURCE_NODE
         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
      ENDIF
      IF (LOCAL_RECORD .LT. 1 .OR. LOCAL_RECORD .GT. GF_MAXREC) THEN
C increase GF_MAXREC in bucket_sort.h
         PRINT *, 'ERROR: LOCAL_FILE_RECORD: LOCAL_RECORD = ',
     +      LOCAL_RECORD
         CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
      ENDIF

      LOCAL_FILE_RECORD = LF_OFFSET(SOURCE_NODE) + LOCAL_RECORD

      RETURN
      END

      SUBROUTINE LOCAL_FILE_GET_RECORD(SOURCE_NODE,LOCAL_RECORD)
      IMPLICIT NONE
      INTEGER SOURCE_NODE, LOCAL_RECORD

#include "bucket_sort.h"

      INTEGER IREC
      INTEGER LOCAL_FILE_RECORD
      EXTERNAL LOCAL_FILE_RECORD

      IREC = LOCAL_FILE_RECORD(SOURCE_NODE,LOCAL_RECORD)
      READ(UNIT=LF_UNIT,REC=IREC) GF_RECBUF

      RETURN
      END

      SUBROUTINE LOCAL_FILE_CLOSE_READ
      IMPLICIT NONE

#include "bucket_sort.h"

      CLOSE(UNIT=LF_UNIT,STATUS='DELETE')

      RETURN
      END
C
#else /* VAR_MPI2 */ 
C aspg: dummy subroutine to avoid compilation errors in case mpi is 
C       not enabled and this file somehow ends up in a makefile
      SUBROUTINE GLOBAL_FILE_STUB
      END
#endif /* VAR_MPI2 */
