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

module dirac_cfg

!radovan: this module can contain info about the chosen method and hamiltonian
!         the individual methods and hamiltonians can store their configuration
!         in individual modules (like for instance dft_cfg)
!         if slaves need to access the content make sure it is broadcast to them
!         (now it's not) this should be done by a method
!         which should be part of this module

  implicit none

  save

!                    kohn-sham dft
  logical, public :: dirac_cfg_dft_calculation    = .false.
!                    hartree-fock
  logical, public :: dirac_cfg_hf_calculation     = .false.
!                    hf or ks
  logical, public :: dirac_cfg_scf_calculation    = .false.
!                    MCSCF-srDFT
  logical, public :: dirac_cfg_mcsrdft_calculation = .false.
!                    openrsp module
  logical, public :: dirac_cfg_openrsp            = .false.
!                    visualization module
  logical, public :: dirac_cfg_visual             = .false.
!                    exacc module
  logical, public :: dirac_cfg_exacc              = .false.
!                    relcc module
  logical, public :: dirac_cfg_relcc              = .false.
!                    embedding module
  logical, public :: dirac_cfg_fde                = .false.
  logical, public :: dirac_cfg_fde_response       = .false.
  logical, public :: dirac_cfg_fde_export         = .false.
!                    Polarizable Continuum Model module
  logical, public :: dirac_cfg_pcm                = .false.
  private

end module
