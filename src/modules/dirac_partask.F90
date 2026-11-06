!-----------------------------------------
! Module with the constants previously defined
! as macros at dirac_partask.h include file
! (task codes the master sends to the co-workers)
!-----------------------------------------
   module dirac_partask_mod
   implicit none
   integer, parameter :: FINISH_PAR_SIGNAL = -2
   integer, parameter :: EXIT_NODEMENU     = -1
   integer, parameter :: RELEASE_NODES     = 0
   integer, parameter :: HERFCK_PAR        = 1
   integer, parameter :: MOLTRA_PAR        = 2
   integer, parameter :: RELCC_PAR         = 3
   integer, parameter :: XCINT_PAR         = 4
   integer, parameter :: LUCITA_PAR        = 5
   integer, parameter :: LUCIAREL_PAR      = 6
   integer, parameter :: KRMCSCF_PAR       = 7
   integer, parameter :: INTEREST_INIT_PAR = 8
   integer, parameter :: INTEREST_4C_PAR   = 9
   integer, parameter :: INTEREST_2C_PAR   = 10
   integer, parameter :: GASCIP_PAR        = 11
   integer, parameter :: PELIB_PAR         = 12
   integer, parameter :: POLPRP_PAR        = 13
   integer, parameter :: GEN1INT_INIT      = 14
   integer, parameter :: GEN1INT_GET_INT   = 15
   integer, parameter :: GEN1INT_GET_EXPT  = 16
   integer, parameter :: GEN1INT_GET_CUBE  = 17
   integer, parameter :: FDE_PAR           = 18
   integer, parameter :: VISUAL_PAR        = 24
   end module dirac_partask_mod
