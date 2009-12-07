module interface_for_pammag
! Radovan's dft_cb_interface.F90 used as a model

!  use memory_allocator
!  use dirac_interface

  implicit none

  public common_blocks_vars

  save

!new names for test routine:
  integer, public :: fcklon_unit = 70
  integer, public :: nr_ao_total
  integer, public :: nr_ao_large
  integer, public :: nr_qu
  integer, public :: print_prp_level

  private

contains

  subroutine common_blocks_vars()
!-----------------------------------------------------------------------
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "dummy.h"
#include "dcbfir.h"
#include "dcbham.h"
!-----------------------------------------------------------------------

    print_prp_level = iprprp 
    fcklon_unit     = lufckl
    nr_ao_total     = ntbas(0)
    nr_ao_large     = ntbas(1)
    nr_qu           = nz

  end subroutine

end module
