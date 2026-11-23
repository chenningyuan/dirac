!-----------------------------------------
! Module with variables previously defined
! at src/dirrci/files.inc include file (COMMON /FNUMBS/, /FNAMES/)
!-----------------------------------------
!  SOURCE  : $Source$
   module files_dirrci_mod
   implicit none
!  file numbers and names (initialized in BLOCK DATA DIRDAT before)
   integer :: MRCONEE = 8, MDCINT = 9, MDBINT = 10, MRCTRIV = 11, MRCFINV = 12, &
              MRCEXCS = 13, MRCVECS = 14, MRCTWOE = 15, MDPROP = 16
   character(len=7) :: FNAM(9) = (/ 'MRCONEE', 'MDCINT ', 'MDBINT ', 'MDTRIV ', 'MRCFINV', &
                                    'MRCEXCS', 'MRCVECS', 'MRCTWOE', 'MDPROP ' /)
   end module files_dirrci_mod
