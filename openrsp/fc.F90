module fc

   use interface_ao
   use interface_mo
   use ao_eval

   implicit none

   public get_fc_integrals

   private

contains

   subroutine get_fc_integrals(           &
                               mat_dim,   &
                               mat,       &
                               center,    &
                               take_derv, &
                               center_d,  &
                               ixyz_d     &
                              )

!     --------------------------------------------------------------------------
      integer, intent(in)  :: mat_dim
      real(8)              :: mat(*)
      logical, intent(in)  :: take_derv
      integer, intent(in)  :: center
      integer, intent(in)  :: center_d
      integer, intent(in)  :: ixyz_d
!     --------------------------------------------------------------------------
      real(8), allocatable :: ao(:, :)
      real(8), allocatable :: buffer(:, :)
      integer              :: i, j, iblock
      integer              :: nr1, nr2, st1, st2
      real(8)              :: d
!     --------------------------------------------------------------------------

      call interface_ao_read(.false.)
      call interface_mo_read()
      call ao_eval_init(1, 0, .true.)

      allocate(ao(1,     nr_ao_slices*nr_ao_cartesian))
      ao = 0.0d0
      allocate(buffer(1, nr_ao_slices*nr_ao_cartesian))
      buffer = 0.0d0

      call get_ao(1, center_xyz(center, 1), center_xyz(center, 2), center_xyz(center, 3), ao, buffer)

      do iblock = 1, nr_ao_blocks
         nr1 = ao_block_nr(iblock)
         st1 = ao_block_start(iblock)
         nr2 = ao_block_nr(lssl_block_partner(iblock, 0, 0))
         st2 = ao_block_start(lssl_block_partner(iblock, 0, 0))
         do i = st1, st1 + nr1 - 1
            do j = st2, st2 + nr2 - 1
               if (take_derv) then
                  d = 0.0d0
                  if (ao_center(i) == center_d) then
                     d = d - ao(1, ao_off_g1_m0(ixyz_d, 0) + i)*ao(1, j)
                  end if
                  if (ao_center(j) == center_d) then
                     d = d - ao(1, ao_off_g1_m0(ixyz_d, 0) + j)*ao(1, i)
                  end if
               else
                  d = ao(1, i)*ao(1, j)
               end if
               mat((i-1)*mat_dim + j) = d
            end do
         end do
      end do

      deallocate(ao)
      deallocate(buffer)

   end subroutine

end module
