C
C     Define common blocks and other data for bucket sort
C
C     Offset for the file unit
      PARAMETER (I_BF_UNIT=70)
C
C     The maximum number of files that can be opened.
      PARAMETER (MAX_BUFF_FILES=5)
C
C     The length of records, measured in integer words.
C     Performance tip : Make this an multiple of 1024 MINUS 1 (to leave space for the 2 additional administrative ints)
C
      PARAMETER (NGBFSZ=2*1024-1)
C
C     The maximum number of buckets per file.
      PARAMETER (N_BUFFER=1024)
C
C     The memory needed to store all buckets for one sort (need MAX_MEM
C     times one real and two integer numbers)
C
      PARAMETER (MAX_MEM=NGBFSZ*N_BUFFER)
C
C     The file name for each file. This name is made unique
C     by appending the ID of the node in parallel runs, in
C     in serial runs we append 0.
C
      COMMON/ BF_FLNAME/ BF_NAME(MAX_BUFF_FILES)
      CHARACTER*10 BF_NAME
C
C     Information about the sort and the files :
C
C     - The highest record written
C     - For each buffer
C        - the last record written
C        - the amount of unflushed data in the buffer 
C
      COMMON /BF_RECS/ I_BF_REC(MAX_BUFF_FILES),
     &                 LGREC(N_BUFFER,MAX_BUFF_FILES),
     &                 LGBUF(N_BUFFER,MAX_BUFF_FILES)
C
C     Information necessary for distributed write and associated
C     communication
C
C     The maximum number of send buffers that is used.
      PARAMETER (MAX_SEND_BUFFERS=5)
C     Buffers for sending integrals
      COMMON /BF_SENDI/ ISGBUF(4+4*NGBFSZ,MAX_SEND_BUFFERS),
     &                  IREQ(MAX_SEND_BUFFERS)
C
C     Statistics
      INTEGER*8 N_WRITTEN, N_DELETED
      COMMON /BF_SENDI8/ N_WRITTEN, N_DELETED
