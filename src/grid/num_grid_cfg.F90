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

module num_grid_cfg

   implicit none

   public report_num_grid

   character(30), public :: num_grid_cfg_gridfile
   integer,       public :: num_grid_cfg_angint                  = 41
   integer,       public :: num_grid_cfg_angmin                  = 15
   integer,       public :: num_grid_cfg_integration_check_level = 0
   logical,       public :: num_grid_cfg_estimate_radii          = .false.
   logical,       public :: num_grid_cfg_force_4c_grid           = .false.
   logical,       public :: num_grid_cfg_import_grid             = .false.
   logical,       public :: num_grid_cfg_multigrid               = .false.
   logical,       public :: num_grid_cfg_no_pruning              = .false.
   logical,       public :: num_grid_cfg_zipgrid                 = .true.
   real(8),       public :: num_grid_cfg_radint                  = 1.0d-13
!*effQED
   logical,       public :: num_grid_cfg_numint_print            = .true.   
   logical,       public :: num_grid_cfg_grasp                   = .false.   
   logical,       public :: num_grid_cfg_uehcut                  = .false.   
   logical,       public :: num_grid_cfg_allatom                 = .false.   
   logical,       public :: num_grid_cfg_printpot                = .false.   
   integer,       public :: num_grid_cfg_angqed                  = 15
   integer,       public :: num_grid_cfg_minz                    = 19
   real(8),       public :: num_grid_cfg_screen                  = 1.0D-10

   private

contains

   subroutine report_num_grid()

      write(*, *) ' ===== Numerical integration grid ====='

      if (num_grid_cfg_import_grid) then

         write(*, '(a)')        '   - integration grid imported from file '//num_grid_cfg_gridfile

      else

         write(*, '(a)')        '   - radial quadrature according to '// &
                                'R. Lindh, P.-Aa. Malmqvist, and L. Gagliardi'
         write(*, '(a, e12.5)') '     precision of radial quadrature set to: ', num_grid_cfg_radint
         write(*, '(a, i2)')    '   - angular quadrature using the Lebedev scheme, '// &
                                'exact up to order L = ', num_grid_cfg_angint

         if (num_grid_cfg_no_pruning) then
            write(*, '(a)')      '   - pruning of angular grid turned off'
         end if

         if (num_grid_cfg_integration_check_level > 0) then
            write(*, '(a)')      '   - numerical integration will be tested'
         end if

         if (num_grid_cfg_estimate_radii) then
            write(*, '(a)')      '   - estimate relative atomic sizes for use in the Becke'
            write(*, '(a)')      '     partitioning scheme from atomic contributions '// &
                                 'to the density'
         end if

      end if

   end subroutine

end module
