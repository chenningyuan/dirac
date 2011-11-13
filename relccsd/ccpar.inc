C
C VERSION : $Revision: 1.7 $
C DATE    : $Date: 2003/02/06 09:16:54 $
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
