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

module srdft_cfg


! Manu: this module contains info about MCSCF-srDFT 

  implicit none

  save

!               compute long-range C integrals (needed for MCSCF-srDFT)
  logical, public :: srdft_cfg_lrcoulomb_int = .false.
!               compute short-range C integrals (needed for MCSCF-srDFT)
  logical, public :: srdft_cfg_srcoulomb_int = .false.

  private

end module
