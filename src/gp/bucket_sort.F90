!-----------------------------------------
! Module with variables previously defined
! at src/gp/bucket_sort.h include file (COMMON /BF_FLNAME/, /BF_RECS/, /BF_SENDI/,
! /BF_SENDI8/, /KL_BUF_INFO/, /GF_STATE/, /GF_BUFFER/, /LF_NAME/)
!-----------------------------------------
!     Define common blocks and other data for bucket sort
   module bucket_sort_mod
   implicit none
!  Offset for the file unit
   integer, parameter :: I_BF_UNIT = 70
!  The maximum number of files that can be opened.
   integer, parameter :: MAX_BUFF_FILES = 5
!  The length of records, measured in integer words.
!  Performance tip : Make this an multiple of 1024 MINUS 1 (to leave space for the 2 additional administrative ints)
   integer, parameter :: NADMIN = 1
!  !!! CAUTION!!! the current value for this parameter leads to
!  crashes on various machines BRUTUS/ETH Zuerich, Horseshoe/SDU Odense, my laptop, ...
!  (I could not test IBM machines since I do not have access anymore).
!  SK - Aug 2010.
!
!  example: node 1 will always try to read/write the first record at:
!  (( NKWORDS * 1024 *4)**2) * 4 * 64 > 4.38 * 10^12.
   integer, parameter :: NKWORDS = 32
   integer, parameter :: NKBUFFS = 1
   integer, parameter :: NGBFSZ = NKWORDS*1024-NADMIN
!  The maximum number of buckets per file.
   integer, parameter :: N_BUFFER = NKBUFFS*1024
!  The memory needed to store all buckets for one sort (need MAX_MEM
!  times one real and two integer numbers)
   integer, parameter :: MAX_MEM = NGBFSZ*N_BUFFER
!  The file name for each file. This name is made unique
!  by appending the ID of the node in parallel runs, in
!  in serial runs we append 0.
   character(len=10) :: BF_NAME(MAX_BUFF_FILES)
!  Information about the sort and the files :
!  - The highest record written
!  - For each buffer
!     - the last record written
!     - the amount of unflushed data in the buffer
   integer :: I_BF_REC(MAX_BUFF_FILES), LGREC(N_BUFFER,MAX_BUFF_FILES), &
              LGBUF(N_BUFFER,MAX_BUFF_FILES)
!  Information necessary for distributed write and associated
!  communication
!  The maximum number of send buffers that is used.
   integer, parameter :: MAX_SEND_BUFFERS = 5
!  Buffers for sending integrals
   integer :: ISGBUF(4+4*NGBFSZ,MAX_SEND_BUFFERS), IREQ(MAX_SEND_BUFFERS)
!  Statistics
   integer*8 :: N_WRITTEN, N_DELETED
!  aspg, 2009-04-24
!  variables used in the fine-sorting of the 1HT of scheme6
   integer, parameter :: MAX_BUF_files = 1
   integer, parameter :: I_SORTED_UNIT = 400
   integer :: KLpairsInBuffer, max_records_per_buffer, records_in_buffer(N_BUFFER), &
              node_for_buffer(N_BUFFER), irec_md_lastwritten(N_BUFFER), &
              irec_md_global_lastwritten
   integer, parameter :: i_md_unit = 931
!
   integer, parameter :: GF_UNIT = 91
   character(len=*), parameter :: GF_NAME = 'gf-matrix'
   integer, parameter :: LF_UNIT = 92
   character(len=10) :: LF_NAME
   integer, parameter :: BLOCKSIZE = (NGBFSZ+1)*4
   integer, parameter :: MAX_MPI_SIZE = 32
#if defined (INT_STAR8)
   integer, parameter :: GF_RECLEN = 8 * BLOCKSIZE
#else
   integer, parameter :: GF_RECLEN = 4 * BLOCKSIZE
#endif
   integer, parameter :: GF_MAXREC = 64 * BLOCKSIZE
   integer :: GF_BREC(N_BUFFER, 0:MAX_MPI_SIZE-1, 0:MAX_MPI_SIZE-1), &
              GF_NREC(0:MAX_MPI_SIZE-1, 0:MAX_MPI_SIZE-1), LF_OFFSET(0:MAX_MPI_SIZE-1), &
              GF_HANDLE, LF_HANDLE, WRITE_REQUEST
   character :: GF_RECBUF(GF_RECLEN)
   end module bucket_sort_mod
