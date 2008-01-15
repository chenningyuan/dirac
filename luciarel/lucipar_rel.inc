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
     &        LF2_ZERO, ISMEMFAC, LMEMFREE_PTR, IIJKL_ROD_LL
	PARAMETER (MSLVOUT_REL = 7)
CSK   COMMON/LUCIPARREL/         MYPROC, NMPROC, LUWRT, 
      COMMON/LUCIPARREL/ MASTER, MYPROC, NMPROC, LUWRT, 
     &                   NPTEST_VAR, NUM_BLOCKS, IIOMOD_REL,
     &                   NFLGRPS_REL, NROOT_SAVE, NVEC_SAVE,
     &                   IRC_SAVE, JVEC_SF, NUM_BLOCKS2,L_COMBI,
     &                   IAM_NOT_INV, I_NZERO_LEN, IBLOCKCOMM, 
     &                   ISHARE_DENSM, IPARALLELIO, LF2_ZERO,
     &                   ISMEMFAC, LMEMFREE_PTR, IIJKL_ROD_LL
       INTEGER*8 LEN_ALL_INT, IS_LENGTH_TT, NINT_TOTAL, LEN_T_BUFF_NZ
       COMMON/LUCIPAR_I8/ LEN_ALL_INT, IS_LENGTH_TT, 
     &                    NINT_TOTAL, LEN_T_BUFF_NZ
       REAL*8 mem_cp_time
       COMMON/LUCI_REAL/ mem_cp_time
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
     &        ICOMM_ID, ICOMM_SIZE, MY_GROUPN, N_MASTER,
     &        MYNEW_ID_SM, MYNEW_COMM_SM, NEWCOMM_PROC_SM,
     &        N_MASTER_SM
      INTEGER IALL_LU1, IALL_LU2, IALL_LU3, IALL_LU4, IALL_LU5, 
     &        IALL_LUC, IIJKL_ROD
C
      COMMON/LUPARGROUP/ NEWCOMM_PROC, MYNEW_ID, MYNEW_COMM, ICOMM,
     &	                 ICOMM_ID, ICOMM_SIZE, MY_GROUPN, N_MASTER,
     &                   MYNEW_ID_SM, MYNEW_COMM_SM, NEWCOMM_PROC_SM,
     &                   N_MASTER_SM
      COMMON/LUCIPFILE/ IDIA, ILU1, ILU2, ILU3, ILU4, ILU5, ILU6, ILU7,
     &                  ILUC,
     &                  IALL_LU1, IALL_LU2, IALL_LU3, IALL_LU4,
     &                  IALL_LU5, IALL_LUC, IIJKL_ROD
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
      LOGICAL SHARED_M, REORD_IJKL
      COMMON/LUCI_LOG/ SHARED_M, REORD_IJKL
!     character block
      CHARACTER*9 NIIJKL_ROD
      COMMON/LUCI_CHAR/ NIIJKL_ROD
