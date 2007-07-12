C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : lucipar_rel.inc
C
C MPI information for LUCIAREL
C
C     General MPI information
C
      INTEGER MSLVOUT_REL, LUWRT, NPTEST_VAR, NUM_BLOCKS,
     &        IIOMOD_REL, NFLGRPS_REL, NROOT_SAVE, NVEC_SAVE,
     &        IRC_SAVE, JVEC_SF, NUM_BLOCKS2,L_COMBI,IAM_NOT_INV,
     &        I_NZERO_LEN, IBLOCKCOMM, ISHARE_DENSM, IPARALLELIO,
     &        LF2_ZERO, ISMEMFAC, LEN_ALL_INT
	PARAMETER (MSLVOUT_REL = 7)
      COMMON/LUCIPARREL/ MASTER, MYPROC, NMPROC, LUWRT, 
     &                   NPTEST_VAR, NUM_BLOCKS, IIOMOD_REL,
     &                   NFLGRPS_REL, NROOT_SAVE, NVEC_SAVE,
     &                   IRC_SAVE, JVEC_SF, NUM_BLOCKS2,L_COMBI,
     &                   IAM_NOT_INV, I_NZERO_LEN, IBLOCKCOMM, 
     &                   ISHARE_DENSM, IPARALLELIO, LF2_ZERO,
     &                   ISMEMFAC, LEN_ALL_INT
C
C     MASTER: contains node number of master (0 in general)
C
C     MYPROC:  -"-     own process number. According to MPI
C             counting starts at zero.....# proc -1 !
C
C     NMPROC:  -"-     total number of processes (including master !)
C
C
C     IIOMOD_REL: determines integral I/O mode in parallel calculations
C
      INTEGER MYNEW_ID, MYNEW_COMM, NEWCOMM_PROC, ICOMM, 
     &        ICOMM_ID, ICOMM_SIZE, MY_GROUPN
      INTEGER IALL_LU1, IALL_LU2, IALL_LU3, IALL_LU4, IALL_LU5, 
     &        IALL_LUC, LUDUMPF
C
      COMMON/LUPARGROUP/ NEWCOMM_PROC, MYNEW_ID, MYNEW_COMM, ICOMM,
     &	             ICOMM_ID, ICOMM_SIZE, MY_GROUPN
      COMMON/LUCIPFILE/ IDIA, ILU1, ILU2, ILU3, ILU4, ILU5, ILU6, ILU7,
     &                  ILUC,
     &                  IALL_LU1, IALL_LU2, IALL_LU3, IALL_LU4,
     &                  IALL_LU5, IALL_LUC, LUDUMPF
#if defined (VAR_MPI)
      INTEGER MY_ACT_BLK1, MY_ACT_BLK2, MY_ACT_BLK_ALL
      INTEGER*8 MY_VEC1_IOFF, MY_VEC2_IOFF, MY_DIA_OFF, MY_LU1_OFF,
     &          MY_LU2_OFF, MY_LU3_OFF, MY_LU4_OFF, MY_LU5_OFF,
     &          MY_LUC_OFF
	COMMON/LUCIFILEOFF/ MY_VEC1_IOFF, MY_VEC2_IOFF, 
     &                    MY_DIA_OFF, MY_LU1_OFF,
     &                    MY_LU2_OFF, MY_LU3_OFF, MY_LU4_OFF, 
     &                    MY_LU5_OFF, MY_LUC_OFF,
     &                    MY_ACT_BLK1, MY_ACT_BLK2, MY_ACT_BLK_ALL
#endif
