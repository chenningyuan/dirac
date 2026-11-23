!-----------------------------------------
! Module with variables previously defined
! at src/relccsd/files.inc include file (COMMON /FILES/, /RSFILE/, /PARNAM/)
!-----------------------------------------
!  Files to be used in RELCCSD
   module files_relccsd_mod
   implicit none
   integer, parameter :: NCCFIL = 11    ! Number of word addressable files
!  Define the unit numbers
   integer, parameter :: MSLVOUT = 12   ! Output from slaves
   integer, parameter :: MRCONEE = 13   ! Effective one-electron integrals
   integer, parameter :: MDINT   = 14   ! Two-electron integrals
   integer, parameter :: MDPROP  = 15   ! Property integrals
   integer, parameter :: MCCRES  = 16   ! Restart information
   integer, parameter :: MCCDENS = 17   ! Effective 1-electron density matrix
   integer, parameter :: IUHEFF  = 18   ! Effective Hamiltonian (Fock space)
   integer, parameter :: ITAMPL  = 26   ! T1 and T2 amplitudes (used in DIIS only)
   integer, parameter :: ITAMPD  = 27   ! Error vector for T1 and T2 amplitudes
   integer, parameter :: IOAMPL  = 28   ! Orbital gradients
   integer, parameter :: IOAMPS  = 29   ! Sigma vector for orbital gradients
   integer, parameter :: ITAPT   = 40   ! Sorted integral files
   integer, parameter :: ITAPL   = 50   ! DZ intermediate for Lambda equations
   integer, parameter :: ITAPADC = 60   ! Base for ADC-files
!  The standard input and output units may be redefined in the code
   integer :: IN, IW
!  Length of the records (integer words) on the restart file
   integer, parameter :: RSRCLEN = 1024
   integer :: IRECIN, IRECST, IRECEN, IRECT1, NRECT1, IRECT2, NRECT2, LRSFIL, LCCFIL
   character(len=6)  :: RSBASF   ! generic name for RESTART file
   character(len=11) :: RSFILN   ! node-unique RESTART file name of length LRSFIL
   character(len=11) :: CCBASF   ! generic name for RELCCSD output file
   character(len=15) :: CCFILN   ! node-unique RELCCSD output file name
   end module files_relccsd_mod
