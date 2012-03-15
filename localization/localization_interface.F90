  subroutine localization_interface
  !==============================================================================!
  !           Interface of the minimization module with Dirac program            !
  !==============================================================================!
  implicit none

#include "dcbloc.h"
#include "priunit.h"

  ! Internal variables
  integer :: ier, iunit

  !------------------------------------------------------------------------------!
  iunit = 60

  ! open file
  open(iunit,file='OPTIM.INP',status='new',form='formatted',access='sequential',iostat=ier)

  if(ier /= 0) call quit('init_minimization: ERROR while creating new OPTIM.INP file')
  rewind 60

  write(iunit,*) iprloc, itrloc, lupri
  write(iunit,*) thfull, thdiag, thgrad
  write(iunit,*) lgfull, lgdiag, lggrad, lgchck
  write(iunit,*) hesloc

  ! close file
  close(iunit,status='keep')

  end subroutine localization_interface
