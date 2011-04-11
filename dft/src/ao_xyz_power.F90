!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

module ao_xyz_power

   implicit none

   public ao_xyz_power_init
   public ao_xyz_power_delete

   private

   integer, allocatable, public :: x_power(:)
   integer, allocatable, public :: y_power(:)
   integer, allocatable, public :: z_power(:)
   integer, allocatable, public :: power_offset(:)

contains

   subroutine ao_xyz_power_init(l_max)

!     --------------------------------------------------------------------------
      integer, intent(in) :: l_max
!     --------------------------------------------------------------------------
      integer             :: l, m, n
      integer             :: length
!     --------------------------------------------------------------------------

!>    reset arrays
      call ao_xyz_power_delete()

!>    get length of arrays for subsequent allocation
      length = 0
      do l = 1, l_max
         do m = 1, l
            do n = 1, m
               length = length + 1
            end do
         end do
      end do

!>    allocate arrays
      allocate(x_power(length))
      allocate(y_power(length))
      allocate(z_power(length))
      allocate(power_offset(l_max))

!>    fill arrays
      length = 0
      do l = 1, l_max
         power_offset(l) = length + 1
         do m = 1, l
            do n = 1, m
               length = length + 1
               x_power(length) = l - m
               y_power(length) = m - n
               z_power(length) = n - 1
            end do
         end do
      end do

   end subroutine

   subroutine ao_xyz_power_delete()

      if (allocated(x_power     )) deallocate(x_power     )
      if (allocated(y_power     )) deallocate(y_power     )
      if (allocated(z_power     )) deallocate(z_power     )
      if (allocated(power_offset)) deallocate(power_offset)

   end subroutine

end module
