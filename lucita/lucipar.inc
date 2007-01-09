C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : lucipar.inc
C
C MPI information for LUCITA
C
C     General MPI information
C
      CHARACTER*12 CVECTOTN,SCLFTOTN
      CHARACTER*10 SVECNAME
      CHARACTER*6 PLU1,PLU2,PLU3,PLU4,PLU5,PLU6,PLU7,PLU8
      INTEGER IRUNPA,INI_ITS,IRDSCF,ICOMPUTE,
     &        NPARWAY,IRUNPCG,MYNEW_COMM,
     &        MYNEW_ID,NEWCOMM_PROC,LUSCR_MPI,CVECFLL,SCLFFLL,
     &        IRDWR,ICVECFL,ISCLFCFL,ISVECFL,ILU1,ILU2,ILU3,ILU4,
     &        ILU5,ILU6,ILU7,ILU8,JCOMCOM
      COMMON/LUCIPAR/ MASTER,MYPROC,NMPROC,IRUNPA,INI_ITS,
     &                IRDSCF,ICOMPUTE,
     &                NPARWAY,IRUNPCG,MYNEW_COMM,MYNEW_ID,
     &                NEWCOMM_PROC,LUSCR_MPI,CVECFLL,SCLFFLL,
     &                IRDWR,ICVECFL,ISCLFCFL,ISVECFL,
     &                ILU1,ILU2,ILU3,ILU4,ILU5,ILU6,ILU7,ILU8,
     &                JCOMCOM,
     &                PLU1,PLU2,PLU3,PLU4,PLU5,PLU6,PLU7,PLU8,
     &                CVECTOTN,SCLFTOTN,SVECNAME
C
C     MASTER: contains node number of master (0 in general)
C
C     MYPROC:  -"-     own process number. According to MPI
C             counting starts at zero.....# proc -1 !
C
C     NMPROC:  -"-     total number of processes (including master !)
C
C     IRUNPA: IRUNPA = 0 --> serial run; IRUNPA = 1 --> parallel run
C
C     INI_ITS: enables sending of c-coeff. in the first loop over s-batches;
C              only for 0 < NPARWAY < 4 !
C
C     IRDSCF: read in SCLFAC in RASSG3 (removes redundant reading in SBLOCKS)
C               1 --> enabled; allocation in MV7! 
C
C     ICOMPUTE: controls all computation after the call to sblock
C               1 --> computation is done; 0 --> just checking block distribution
C               only for 3 < NPARWAY <= 6 !
C
C     NPARWAY: controls the enabling/disabling of the fg-version
C              0 --> no fg! 1 - 3 --> fg incl. master 
C              4 - 6 --> fg wo master 
C
C     IRUNPCG: controls the enabling/disabling of the cg-version
C              0 --> no cg! 1 --> cg incl. MASTER 2 --> cg wo MASTER
C              3 --> cg wo MASTER in combination with fg
C
C     MYNEW_COMM: new communicator used in parallel calculations 
C                 (split of MPI_COMM_WORLD)
C
C     MYNEW_ID: contains own process number in MYNEW_COMM. According to MPI 
C               counting starts at zero !
C
C     NEWCOMM_PROC: contains total number of processes (including MYNEW_ID = 0 !) 
C                   in the subgroup MYNEW_COMM!
C
C     LUSCR_MPI: contains file handle for parallel communication file SIGMA.MYPROC 
C
