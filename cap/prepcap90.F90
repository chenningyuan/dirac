!=======================================================================
!
! file contains pure Fortran90 classes related to CAP project
!
!=======================================================================

  module cap

! private ! must be public that subroutines can be used in F77 parts !
#include "dcbcap.h"

!  logical, private :: cap_ovlp
!   integer(kind=4), private :: iprint
!   integer(kind=4), private :: iprcap

   type cap_int_param
   end type cap_int_param

  contains

  subroutine cap_init
!#include "implicit.h"
#include "priunit.h"
!#include "dcbcap.h"

! iprcap = IPRCAP
! cap_ovlp = CAP_OVLP

  if (IPRCAP >= 2) then
   write(LUPRI,*) 'in init routine !'
   write(LUPRI,*) 'module iprcap =',IPRCAP
   write(LUPRI,*) 'module cap_ovlp =',CAPOVLP
  endif

! call quit('in cap_init...')

   end subroutine cap_init

  end module cap

 subroutine cap_ovlp(WORK,LWORK)
use Memory
! check if not make check of DIRAC/CAP overlap integrals

 end subroutine cap_ovlp

