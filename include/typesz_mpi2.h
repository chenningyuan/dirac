!
! VERSION : $Revision: 6563 $
! DATE    : $Date: 2008-03-03 14:23:33 +0100 (Mon, 03 Mar 2008) $
! FILE    : typesz_mpi2.h
!
!     MPI-2 type size for memory allocation
!     ======================================
!
#if defined (VAR_MPI2)
!     INTEGER(KIND=MPI_ADDRESS_KIND) ISIZE_dp, ISIZE_my_int
      INTEGER ISIZE_dp, ISIZE_my_int, ISIZE_intdef, ISIZE_int8,
     &                  ISIZE_log
!
      COMMON/TYPE_SZMPI2/ ISIZE_dp, ISIZE_my_int, 
     &                    ISIZE_intdef, ISIZE_int8, ISIZE_log
!
#endif
