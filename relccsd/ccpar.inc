C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : inpt.inc
C SOURCE  : $Source$
C
C MPI information for RELCCSD
C
C     General MPI information
C
      COMMON/CCPAR/ MASTER,MYPROC,NMPROC
C
C     Distribution of two-electron integrals
C     First index : start, stop of block on node
C     Second index : integral type (VOVV,VVVV)
C     Third index : irrep
C
      COMMON/CCPARV/ NV6PAR,IDIST(6,2,8)
