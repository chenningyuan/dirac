!-----------------------------------------
! Module with variables previously defined
! at dcbana.h include file (COMMON /DCLPRP/)
!-----------------------------------------
!     wave function analysis module
!     =============================
!     DOPOP  - perform Mulliken population analysis
!     DOPRJ  - perform projection onto another solution
!     DOVEC  - print MO-coefficients
!     DORHO1 - density plot in one dimension
!              for visualization
!     DOLOC  - perform orbital localization
   module dcbana_mod
   implicit none
   logical :: DOPOP, DOPRJ, DOVEC, DO1RHO, DO3RHO, DO1WT, DOLOC
   end module dcbana_mod
