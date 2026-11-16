!-----------------------------------------
! Module with variables previously defined
! at files_r.inc include file (COMMON /IPARNAMREL/, /PARNAMREL/)
!-----------------------------------------
   module files_r_mod
   implicit none
   integer :: LUSC929_MPIL, LUSC94_MPIL, LUSC95_MPIL, LUMOUTRI_MPIL, LUSC344_MPIL
   character(len=10) :: LUC_MPIN, LUHC_MPIN, LUSC34_MPIN, LUSC35_MPIN, LUSC36_MPIN, &
                        LUSC37_MPIN, LUSC38_MPIN, LUSC39_MPIN, LUSC61_MPIN, LUSC62_MPIN, &
                        LUSC93_MPIN, LUSC92_MPIN
   character(len=11) :: NLUSC95_MPI, NLUMOUTR_MPI, NLUMOUTI_MPI
   character(len=13) :: NLUDIA_MPI, NLUSC1_MPI, NLUSC2_MPI, NLUSC3_MPI, NLUSC94_MPI
   character(len=14) :: NLUC_MPI, NLUHC_MPI, NLUSC34_MPI, NLUSC35_MPI, NLUSC36_MPI, &
                        NLUSC37_MPI, NLUSC38_MPI, NLUSC39_MPI, NLUSC41_MPI, NLUSC61_MPI, &
                        NLUSC62_MPI, NLUSC93_MPI, NLUSC92_MPI
   character(len=7) :: LUSC95_MPIN, LUMOUTR_MPIN, LUMOUTI_MPIN
   character(len=9) :: LUDIA_MPIN, LUSC1_MPIN, LUSC2_MPIN, LUSC3_MPIN, LUSC94_MPIN
   end module files_r_mod
