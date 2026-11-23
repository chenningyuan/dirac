!-----------------------------------------
!  Module with variables previously defined
!  at src/relccsd/ccpar.inc include file
!-----------------------------------------
!  MPI information for RELCCSD
!      General MPI information
!  (the header's #define my_MPI_INTEGER was used by none of its users: dropped)
   module ccpar_mod
   implicit none
   integer :: MASTER, MYPROC, NMPROC
!  MASTER: contains node number of master (0 in general)
!  MYPROC:  -"-     own process number. According to MPI
!          counting starts at zero.....# proc -1 !
!  NMPROC:  -"-     total number of processes (including master !)

!  Information for integral distribution
!     Distribution of two-electron integrals
!     First index : start, stop of block on node
!     Second index : integral type (VOVV,VVVV)
!     Third index : irrep
   integer*8 :: NV6PAR, NV5PAR, IDIST(6,10,32)
   end module ccpar_mod
