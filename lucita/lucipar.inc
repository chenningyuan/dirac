C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : inpt.inc
C
C MPI information for LUCITA
C
C     General MPI information
C
      INTEGER IRUNPA,INI_ITS,IRDSCF,ICOMPUTE,
     &        NPARWAY,IRUNPCG,MYNEW_COMM,
     &        MYNEW_ID,NEWCOMM_PROC,LUSCR_MPI
      COMMON/LUCIPAR/ MASTER,MYPROC,NMPROC,IRUNPA,INI_ITS,
     &                IRDSCF,ICOMPUTE,
     &                NPARWAY,IRUNPCG,MYNEW_COMM,MYNEW_ID,
     &                NEWCOMM_PROC,LUSCR_MPI
C
C     MASTER: contains node number of master (0 in general)
C     MYPROC:  -"-     own process number. According to MPI
C             counting starts at zero.....# proc -1 !
C     NMPROC:  -"-     total number of processes (including master !)
C
C     IRUNPA: IRUNPA = 0 --> serial run; IRUNPA = 1 --> parallel run
C     INI_ITS: sending specified c-coeff. in the first loop over S batches
C     IRDSCF: read in SCALFAC in MV7 and NOT in SBLOCKS (redundant reading)
C               1 --> MV7; 0 --> SBLOCKS (allocation in SBLOCK!)
C     ICOMPUTE: controls all computation after the call to sblock
C               1 --> computation is done; 0 --> just checking block distribution
C     NPARWAY: controls the enabling/disabling of the fg-version
C              0 --> no fg! 1 --> fg incl. MASTER (see IUSE_VER = 1) 
C              2 --> fg wo MASTER (see IUSE_VER = 3)
C
C     IRUNPCG: controls the enabling/disabling of the cg-version
C              0 --> no cg! 1 --> cg incl. MASTER 2 --> cg wo MASTER
C              3 --> cg wo MASTER in combination with fg
C
C     IRTNDRAT: ratio between numnod and nroot
C
