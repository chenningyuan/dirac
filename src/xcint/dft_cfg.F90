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

module dft_cfg

   use xc_mpi

   implicit none

   public sync_dft_cfg

   logical, public :: dft_cfg_alda_hs                 = .false.
   logical, public :: dft_cfg_alda_ha                 = .false.
   logical, public :: dft_cfg_xalda                   = .false.
   logical, public :: dft_cfg_no_sdft                 = .false.
   logical, public :: dft_cfg_sdft_collinear          = .false.
   logical, public :: dft_cfg_grac                    = .false.
   logical, public :: dft_cfg_grac_is_active          = .false.
   logical, public :: dft_cfg_saop                    = .false.
   logical, public :: dft_cfg_saop_is_active          = .false.
   logical, public :: dft_cfg_saop_with_response_part = .false.
   logical, public :: dft_cfg_asymptote_is_lb94       = .false.
   logical, public :: dft_cfg_asymptote_is_lbalpha    = .false.
   logical, public :: dft_cfg_blocked                 = .false.
   logical, public :: dft_cfg_pointwise               = .false.
   logical, public :: dft_cfg_overlap_diagnostic      = .false.

   real(8), public :: dft_cfg_tinydens                = 1.0d-14
   real(8), public :: dft_cfg_screening               = 1.0d-18
   real(8), public :: dft_cfg_ac_ip                   = 0.0d0
   real(8), public :: dft_cfg_ac_threshold            = 0.0d0
   real(8), public :: dft_cfg_grac_alpha              = 0.0d0
   real(8), public :: dft_cfg_grac_beta               = 0.0d0

   private

contains

   subroutine sync_dft_cfg()

#ifdef VAR_MPI
      call xc_mpi_bcast(dft_cfg_alda_hs                 )
      call xc_mpi_bcast(dft_cfg_alda_ha                 )
      call xc_mpi_bcast(dft_cfg_xalda                   )
      call xc_mpi_bcast(dft_cfg_no_sdft                 )
      call xc_mpi_bcast(dft_cfg_sdft_collinear          )
      call xc_mpi_bcast(dft_cfg_grac                    )
      call xc_mpi_bcast(dft_cfg_grac_is_active          )
      call xc_mpi_bcast(dft_cfg_saop                    )
      call xc_mpi_bcast(dft_cfg_saop_is_active          )
      call xc_mpi_bcast(dft_cfg_saop_with_response_part )
      call xc_mpi_bcast(dft_cfg_asymptote_is_lb94       )
      call xc_mpi_bcast(dft_cfg_asymptote_is_lbalpha    )
      call xc_mpi_bcast(dft_cfg_blocked                 )
      call xc_mpi_bcast(dft_cfg_pointwise               )
      call xc_mpi_bcast(dft_cfg_overlap_diagnostic      )
      call xc_mpi_bcast(dft_cfg_tinydens                )
      call xc_mpi_bcast(dft_cfg_screening               )
      call xc_mpi_bcast(dft_cfg_ac_ip                   )
      call xc_mpi_bcast(dft_cfg_ac_threshold            )
      call xc_mpi_bcast(dft_cfg_grac_alpha              )
      call xc_mpi_bcast(dft_cfg_grac_beta               )
#endif /* ifdef VAR_MPI */

   end subroutine

end module
