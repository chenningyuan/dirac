C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : lucipar_rel.inc
C
C MPI information for LUCIAREL
C
C     General MPI information
C
      COMMON/LUCIPARREL/ MASTER,MYPROC,NMPROC
C
C     MASTER: contains node number of master (0 in general)
C
C     MYPROC:  -"-     own process number. According to MPI
C             counting starts at zero.....# proc -1 !
C
C     NMPROC:  -"-     total number of processes (including master !)
C
