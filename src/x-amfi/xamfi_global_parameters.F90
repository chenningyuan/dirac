!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

module xamfi_global_parameters

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
  logical, public :: aoofocc                                =  .true.  ! X-AMFI with fractional-occupation 2e-Fock-matrix
  logical, public :: aoo2ediag                              =  .false. ! amfX2C:  X-AMFI with atomic 2e-Fock contributions created/added to diagonal of defining X2C h1
  logical, public :: aooeamf                                =  .false. ! eamfX2C: X-AMFI with first-order molecular 2e-Fock contributions created/added to defining X2C h1 

! integer block
  integer, public :: aoo_prt                                =  0 ! print level
  integer, public :: aoo_soc_order                          =  2 ! default: spin-same orbit contributions

! real*8 block
  real*8, public :: xamfi_energy_contributions              = 0.0d0

end module
