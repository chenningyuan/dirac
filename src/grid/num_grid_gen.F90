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

module num_grid_gen

   implicit none

   public delete_num_grid
   public reset_num_grid
   public generate_num_grid
   public stash_zipgrid
   public unstash_zipgrid

   private

   logical :: grid_is_done     = .false.
   logical :: remember_zipgrid = .false.

contains

   subroutine delete_num_grid()

   logical :: tobe
   integer :: ludft

   !> delete DFT-grid
   inquire(file='numerical_grid',exist=tobe)
   if(tobe)then
     inquire(file='numerical_grid',opened=tobe)
     if(tobe)then
       inquire(file='numerical_grid',number=ludft)
       CLOSE(ludft,STATUS='DELETE')
     else
       ludft = 109
       CALL OPNFIL(ludft,'numerical_grid','OLD','xxxxxx')
       CLOSE(ludft,STATUS='DELETE')
     end if
   end if

   end subroutine

   subroutine reset_num_grid()

      grid_is_done = .false.

   end subroutine

   subroutine stash_zipgrid(zipgrid)

      logical, intent(inout) :: zipgrid

      remember_zipgrid = zipgrid
      zipgrid          = .false.

   end subroutine

   subroutine unstash_zipgrid(zipgrid)

      logical, intent(out) :: zipgrid

      zipgrid = remember_zipgrid

   end subroutine

   subroutine generate_num_grid(dmat)

      use dirac_cfg
!     --------------------------------------------------------------------------
      real(8), intent(in) :: dmat(*)
!     --------------------------------------------------------------------------
! dcbham.h : iprham
#include "dcbham.h"

      if (.not. grid_is_done) then
!       calculate grid
        call dft_grid(dmat, iprham)
        grid_is_done = .true.
      end if

   end subroutine

end module
