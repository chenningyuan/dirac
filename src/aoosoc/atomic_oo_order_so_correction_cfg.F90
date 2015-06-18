!dirac_copyright_start
!      Copyright (c) 2014 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

module atomic_oo_order_so_correction_cfg

  implicit none

  save

! parameters

  integer, parameter, public :: aoo_fsoc                    =  26    ! file unit for aoo2esoc 

! character block
  character (len=6),  public :: file_name_mo_coefficients   = 'undef '
  character (len=8),  public :: file_name_aoo_contributions = 'aoo2esoc'

! logical block
  logical, public :: aoomod                                 =  .false. ! atomic-oo-order module activated
  logical, public :: aoomod_debug                           =  .false. ! debug mode (set print level high)
  logical, public :: aoo_2c_mmf_mos                         =  .false. ! .true. if we are running a 4-component SCF

! integer block
  integer, public :: aoo_prt                                =  0 ! print level
  integer, public :: aoo_soc_order                          =  2 ! default: spin-same orbit contributions

end module
