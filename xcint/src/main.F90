program main

   use interface_ao
   use interface_file_io
   use xc_driver
   use xc_derv

   implicit none

   real(8), allocatable :: dmat(:, :)
   real(8), allocatable :: fmat(:, :, :, :)
   real(8), allocatable :: fmat_reference(:, :)
   integer, parameter   :: unit_grid = 44

   call read_input()
   call interface_ao_read(.false.)

   fun_is_automatic = .true.

   allocate(dmat(nr_ao, nr_ao))
   dmat = 0.0d0

   allocate(fmat(nr_ao, nr_ao, nr_quaternion_blocks, 1))
   fmat = 0.0d0

   call matrix_from_file('xc_dmat0', nr_ao, dmat)

   open(unit_grid,                 &
        file   = 'numerical_grid', &
        status = 'unknown',        &
        form   = 'formatted',      &
        access = 'sequential')
   rewind(unit_grid)

   call xc_integrator(                                        &
                      xc_mat_dim      = nr_ao,                &
                      xc_nz           = nr_quaternion_blocks, &
                      xc_dmat_0       = dmat,                 &
                      xc_nr_dmat      = 0,                    &
                      xc_nr_fmat      = 1,                    &
                      xc_fmat         = fmat,                 &
                      xc_do_potential = .true.                &
                     )

   deallocate(dmat)

   allocate(fmat_reference(nr_ao, nr_ao))
   call matrix_from_file('xc_fmat', nr_ao, fmat_reference)
   call compare_matrices(nr_ao, fmat, fmat_reference)

   deallocate(fmat)
   deallocate(fmat_reference)

contains

   subroutine compare_matrices(nr_ao, mat, mat_reference)

!     --------------------------------------------------------------------------
      integer, intent(in) :: nr_ao
      real(8), intent(in) :: mat(*)
      real(8), intent(in) :: mat_reference(*)
!     --------------------------------------------------------------------------
      integer             :: i
      real(8)             :: diff
      real(8)             :: max_diff
!     --------------------------------------------------------------------------

      max_diff = 0.0d0
      do i = 1, nr_ao*nr_ao
         diff = abs(mat(i) - mat_reference(i))
         if (diff > max_diff) then
            max_diff = diff
         end if
      end do
      print *, 'compare matrices: largest deviation =', max_diff

   end subroutine

   subroutine read_input()

      use keyword
      use character_processing
      use dirac_cfg
      use interface_input
      use interface_functional_read

!     --------------------------------------------------------------------------
      character(kw_length) :: word
      character(kw_length) :: kw_section
      character(80)        :: line
      real(8)              :: hfx_out, mu_out, beta_out
!     --------------------------------------------------------------------------

!     this is to catch keywords that appear before some section starts
      kw_section = '       '

      open(unit_in, file='xc_input')

      rewind unit_in
      do while (.true.)
         read(unit_in, '(a7)', end=1) word

         if (word == '       ') then
!           blank line
         else
            select case (word(1:1))
               case ('!', '#')
!                 comment

               case ('*')
!                 section
                  kw_section = uppercase(word)

               case default
!                 keyword
                  select case (kw_section)
                     case ('*DFT   ')
                        call read_input_dft(word, kw_section)
                     case ('*END OF')
                        go to 1
                     case default
                  end select
            end select
         end if
      end do

1     rewind unit_in
      do while (.true.)
         read(unit_in, '(a80)', end=2) line

         if (line(1:7) == '.DFTAUT') then
            read(unit_in, '(a80)') line
            call parse_functional(line, xc_fun, hfx_out, mu_out, beta_out, .true.)
         end if
      end do

2     close(unit_in)

   end subroutine

end program
