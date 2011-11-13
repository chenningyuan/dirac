C
C VERSION : $Revision: 5599 $
C DATE    : $Date: 2007-03-06 15:22:35 +0100 (Tue, 06 Mar 2007) $
C FILE    : inpt.inc
C
C MPI information for RELCCSD
C
C     General MPI information
C
      COMMON/CCPAR/ MASTER,MYPROC,NMPROC
C
C     MASTER: contains node number of master (0 in general)
C     MYPROC:  -"-     own process number. According to MPI
C             counting starts at zero.....# proc -1 !
C     NMPROC:  -"-     total number of processes (including master !)
C

C Information for integral distribution
C
C     Distribution of two-electron integrals
C     First index : start, stop of block on node
C     Second index : integral type (VOVV,VVVV)
C     Third index : irrep
C
      INTEGER*8 NV6PAR,NV5PAR,IDIST
      COMMON/CCPARV/ NV6PAR,NV5PAR,IDIST(6,10,32)
C
C DIRAC-specific preprocessor statement to solve 32/64 bit compatibility problems 
C
#if defined (INT_STAR8)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

