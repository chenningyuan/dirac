!
! VERSION : $Revision$
! DATE    : $Date$
! FILE    : krmc_shmem.h
!
!                MPI memory pointer 
!     ======================================
!
#if defined (VAR_MPI2)
      INTEGER(KIND=MPI_ADDRESS_KIND) MY_T_PTR
!
      COMMON/MPI_MEMPOINTER/ MY_T_PTR
!
!             MPI memory window handler 
!     ======================================
!
      INTEGER MY_T_WIN
!
      COMMON/MPI_WIN_HAND/ MY_T_WIN
!
!             array length
!     ======================================
!
      INTEGER(KIND=MPI_OFFSET_KIND) LEN_T_BUFF
!
      COMMON/IBUFF_LENGTH/ LEN_T_BUFF
!
#endif
#endif
