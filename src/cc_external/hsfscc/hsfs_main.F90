Program HIGH_SECTORS_FSCC_METHODS
! standalone program, which also reads command line arguments
use high_sect_execs
 call print_hsfscc_header
 call get_line_arguments
 call HIGH_SECTOR_FSCC_METHODS
End Program

subroutine print_hsfscc_header
implicit none
#include "priunit.h"
  write(LUPRI,'(2x,a)') "---------------------------------------------------------------------------------"
  write(LUPRI,*)
  write(LUPRI,'(2x,a)')  'Welcome to the standalone High Sectors Fock Space Coupled Cluster methods program'
  write(LUPRI,*)
  write(LUPRI,'(2x,a)') ' Written by Uzi Kaldor, S.R. Hughes  '
  write(LUPRI,'(2x,a)') ' Later modied and interfaced by Ephraim Eliav, Miroslav Ilias'
  write(LUPRI,*)
  write(LUPRI,'(2x,a)')  '---------------------------------------------------------------------------------'
end subroutine print_hsfscc_header
