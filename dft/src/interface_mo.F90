
module interface_mo

   implicit none

#ifdef PRG_DIRAC
   public interface_mo_write
#endif
   public interface_mo_read
   public interface_mo_clear

   logical, public :: interface_mo_set

   integer, public :: nr_open_shells
   integer, public :: nr_electrons_inactive
   integer, public :: nr_electrons_active
   integer, public :: nr_electrons_total
   integer, public :: nr_mo
   integer, public :: nr_mo_gerade
   integer, public :: nr_mo_gerade_negative_secondary
   integer, public :: nr_mo_gerade_positive_inactive
   integer, public :: nr_mo_gerade_positive_secondary
   integer, public :: nr_mo_ungerade
   integer, public :: nr_mo_ungerade_negative_secondary
   integer, public :: nr_mo_ungerade_positive_inactive
   integer, public :: nr_mo_ungerade_positive_secondary
   integer, public :: nr_df_open
   logical, public :: mo_eigenvalues_available
   logical, public :: mo_coef_available

   real(8), allocatable, public :: mo_eigenvalues(:)
   real(8), allocatable, public :: mo_coef(:)
   real(8), allocatable, public :: df_open(:)

   private

   character(*), parameter :: file_name_mo      = 'interface_mo'
   character(*), parameter :: file_name_mo_coef = 'interface_mo_coef'
   integer,      parameter :: file_unit         = 137
   logical                 :: file_exists

contains

   subroutine interface_mo_clear()

      interface_mo_set                  = .false.

      nr_open_shells                    = 0
      nr_electrons_inactive             = 0
      nr_electrons_active               = 0
      nr_electrons_total                = 0
      nr_mo                             = 0
      nr_mo_gerade                      = 0
      nr_mo_gerade_negative_secondary   = 0
      nr_mo_gerade_positive_inactive    = 0
      nr_mo_gerade_positive_secondary   = 0
      nr_mo_ungerade                    = 0
      nr_mo_ungerade_negative_secondary = 0
      nr_mo_ungerade_positive_inactive  = 0
      nr_mo_ungerade_positive_secondary = 0
      nr_df_open                        = 0
      mo_eigenvalues_available          = .false.
      mo_coef_available                 = .false.

      if (allocated(df_open        )) deallocate(df_open        )
      if (allocated(mo_coef        )) deallocate(mo_coef        )
      if (allocated(mo_eigenvalues )) deallocate(mo_eigenvalues )

   end subroutine

#ifdef PRG_DIRAC
   subroutine interface_mo_write()

!     --------------------------------------------------------------------------
      real(8) :: dummy
      integer :: i
      logical :: lucoef_opened_by_dirac
!     --------------------------------------------------------------------------

#include "dcbbas.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbdhf.h"

      call interface_mo_clear()

      inquire(file = 'DFCOEF', exist = file_exists)
      if (file_exists) then
         inquire(unit=lucoef, opened=lucoef_opened_by_dirac)
         if (lucoef_opened_by_dirac) then
            close(lucoef, status = 'keep')
         end if

         open(file_unit,              &
              file   = 'DFCOEF',      &
              status = 'unknown',     &
              form   = 'unformatted', &
              access = 'sequential')
         rewind(file_unit)
         allocate(mo_eigenvalues(norbt))
         mo_eigenvalues = 0.0d0
         allocate(mo_coef(n2bbasxq))
         mo_coef = 0.0d0

         call reacmo_no_work(file_unit,      &
                             'DFCOEF',       &
                             mo_coef,        &
                             mo_eigenvalues, &
                             (/0/),          &
                             dummy,          &
                             6)
         mo_eigenvalues_available = .true.
         mo_coef_available        = .true.

         close(file_unit, status = 'keep')
         if (lucoef_opened_by_dirac) then
            open(lucoef,                 &
                 file   = 'DFCOEF',      &
                 status = 'unknown',     &
                 form   = 'unformatted', &
                 access = 'sequential')
         end if
      end if

!     always regenerate this file
      inquire(file = file_name_mo, exist = file_exists)
      if (file_exists) then
         open(file_unit,             &
              file   = file_name_mo, &
              status = 'unknown',    &
              form   = 'formatted',  &
              access = 'sequential')
         close(file_unit, status = 'delete')
      end if
      open(file_unit,             &
           file   = file_name_mo, &
           status = 'new',        &
           form   = 'formatted',  &
           access = 'sequential')
      rewind(file_unit)

      nr_open_shells        = nasht
      nr_electrons_inactive = nelect
      nr_electrons_active   = naelec
      nr_electrons_total    = nr_electrons_inactive + nr_electrons_active
      nr_mo                 = norbt

      nr_mo_gerade_negative_secondary = npsh(1)
      nr_mo_gerade_positive_inactive  = nish(1)
      nr_mo_gerade_positive_secondary = nesh(1) - nr_mo_gerade_positive_inactive
      nr_mo_gerade = nr_mo_gerade_negative_secondary &
                   + nr_mo_gerade_positive_inactive  &
                   + nr_mo_gerade_positive_secondary

      nr_mo_ungerade_negative_secondary = npsh(2)
      nr_mo_ungerade_positive_inactive  = nish(2)
      nr_mo_ungerade_positive_secondary = nesh(2) - nr_mo_ungerade_positive_inactive
      nr_mo_ungerade = nr_mo_ungerade_negative_secondary &
                     + nr_mo_ungerade_positive_inactive  &
                     + nr_mo_ungerade_positive_secondary

      if (aoc) then
         nr_df_open = nopen
      else
         nr_df_open = 0
      end if

      write(file_unit, '(4x, a)') 'nr_open_shells'
      write(file_unit, '(4x, i6)') nr_open_shells
      write(file_unit, '(4x, a)') 'nr_electrons_inactive'
      write(file_unit, '(4x, i6)') nr_electrons_inactive
      write(file_unit, '(4x, a)') 'nr_electrons_active'
      write(file_unit, '(4x, i6)') nr_electrons_active
      write(file_unit, '(4x, a)') 'nr_electrons_total'
      write(file_unit, '(4x, i6)') nr_electrons_total
      write(file_unit, '(4x, a)') 'nr_mo'
      write(file_unit, '(4x, i6)') nr_mo
      write(file_unit, '(4x, a)') 'nr_mo_gerade'
      write(file_unit, '(4x, i6)') nr_mo_gerade
      write(file_unit, '(4x, a)') 'nr_mo_gerade_negative_secondary'
      write(file_unit, '(4x, i6)') nr_mo_gerade_negative_secondary
      write(file_unit, '(4x, a)') 'nr_mo_gerade_positive_inactive'
      write(file_unit, '(4x, i6)') nr_mo_gerade_positive_inactive
      write(file_unit, '(4x, a)') 'nr_mo_gerade_positive_secondary'
      write(file_unit, '(4x, i6)') nr_mo_gerade_positive_secondary
      write(file_unit, '(4x, a)') 'nr_mo_ungerade'
      write(file_unit, '(4x, i6)') nr_mo_ungerade
      write(file_unit, '(4x, a)') 'nr_mo_ungerade_negative_secondary'
      write(file_unit, '(4x, i6)') nr_mo_ungerade_negative_secondary
      write(file_unit, '(4x, a)') 'nr_mo_ungerade_positive_inactive'
      write(file_unit, '(4x, i6)') nr_mo_ungerade_positive_inactive
      write(file_unit, '(4x, a)') 'nr_mo_ungerade_positive_secondary'
      write(file_unit, '(4x, i6)') nr_mo_ungerade_positive_secondary
      write(file_unit, '(4x, a)') 'nr_df_open'
      write(file_unit, '(4x, i6)') nr_df_open
      if (nr_df_open > 0) then
         do i = 1, nr_df_open
            write(file_unit, '(4x, i6, e26.16)') i, df_open(i)
         end do
      end if

      if (mo_eigenvalues_available) then
         write(file_unit, '(4x, a)') 'mo_eigenvalues'
         do i = 1, nr_mo
            write(file_unit, '(4x, i6, e26.16)') i, mo_eigenvalues(i)
         end do
      end if

      close(file_unit, status = 'keep')

      if (mo_coef_available) then
!        always regenerate this file
         inquire(file = file_name_mo_coef, exist = file_exists)
         if (file_exists) then
            open(file_unit,                  &
                 file   = file_name_mo_coef, &
                 status = 'unknown',         &
                 form   = 'formatted',       &
                 access = 'sequential')
            close(file_unit, status = 'delete')
         end if
         open(file_unit,                  &
              file   = file_name_mo_coef, &
              status = 'new',             &
              form   = 'formatted',       &
              access = 'sequential')
         rewind(file_unit)
         write(file_unit, *) size(mo_coef)
         do i = 1, size(mo_coef)
            write(file_unit, *) mo_coef(i)
         end do
         close(file_unit, status = 'keep')
      end if

      call interface_mo_clear()

   end subroutine
#endif

   subroutine interface_mo_read()

      use xc_mpi

!     --------------------------------------------------------------------------
      character(120) :: line
      integer        :: i, idummy, n, size_mo_coef
!     --------------------------------------------------------------------------

      call interface_mo_clear()

      if (xc_mpi_is_master()) then
         inquire(file = file_name_mo, exist = file_exists)
         if (.not. file_exists) then
            print *, 'interface_mo file not found'
            stop
         end if
         open(file_unit,             &
              file   = file_name_mo, &
              status = 'unknown',    &
              form   = 'formatted',  &
              access = 'sequential')
         rewind(file_unit)

         read(file_unit, *)
         read(file_unit, *) nr_open_shells
         read(file_unit, *)
         read(file_unit, *) nr_electrons_inactive
         read(file_unit, *)
         read(file_unit, *) nr_electrons_active
         read(file_unit, *)
         read(file_unit, *) nr_electrons_total
         read(file_unit, *)
         read(file_unit, *) nr_mo
         read(file_unit, *)
         read(file_unit, *) nr_mo_gerade
         read(file_unit, *)
         read(file_unit, *) nr_mo_gerade_negative_secondary
         read(file_unit, *)
         read(file_unit, *) nr_mo_gerade_positive_inactive
         read(file_unit, *)
         read(file_unit, *) nr_mo_gerade_positive_secondary
         read(file_unit, *)
         read(file_unit, *) nr_mo_ungerade
         read(file_unit, *)
         read(file_unit, *) nr_mo_ungerade_negative_secondary
         read(file_unit, *)
         read(file_unit, *) nr_mo_ungerade_positive_inactive
         read(file_unit, *)
         read(file_unit, *) nr_mo_ungerade_positive_secondary
         read(file_unit, *)
         read(file_unit, *) nr_df_open
         if (nr_df_open > 0) then
            allocate(df_open(nr_df_open))
            df_open = 0.0d0
            do i = 1, nr_df_open
               read(file_unit, *) idummy, df_open(i)
            end do
         end if

         rewind(file_unit)
         do while (.true.)

            read(file_unit, '(a120)', end=1) line

            if (word_contains(line, 'mo_eigenvalues')) then
               allocate(mo_eigenvalues(nr_mo))
               do i = 1, nr_mo
                  read(file_unit, *) idummy, mo_eigenvalues(i)
               end do
               mo_eigenvalues_available = .true.
            end if

         end do
1        continue

         inquire(file = file_name_mo_coef, exist = file_exists)
         if (file_exists) then
            open(file_unit,                  &
                 file   = file_name_mo_coef, &
                 status = 'unknown',         &
                 form   = 'formatted',       &
                 access = 'sequential')
            rewind(file_unit)
            read(file_unit, *) n
            allocate(mo_coef(n))
            do i = 1, n
               read(file_unit, *) mo_coef(i)
            end do
            close(file_unit, status = 'keep')
            mo_coef_available = .true.
         end if

         interface_mo_set = .true.
      end if

#ifdef VAR_MPI
      if (xc_mpi_get_nr_proc() > 1) then

         call xc_mpi_bcast(interface_mo_set                  )

         call xc_mpi_bcast(nr_open_shells                    )
         call xc_mpi_bcast(nr_electrons_inactive             )
         call xc_mpi_bcast(nr_electrons_active               )
         call xc_mpi_bcast(nr_electrons_total                )
         call xc_mpi_bcast(nr_mo                             )
         call xc_mpi_bcast(nr_mo_gerade                      )
         call xc_mpi_bcast(nr_mo_gerade_negative_secondary   )
         call xc_mpi_bcast(nr_mo_gerade_positive_inactive    )
         call xc_mpi_bcast(nr_mo_gerade_positive_secondary   )
         call xc_mpi_bcast(nr_mo_ungerade                    )
         call xc_mpi_bcast(nr_mo_ungerade_negative_secondary )
         call xc_mpi_bcast(nr_mo_ungerade_positive_inactive  )
         call xc_mpi_bcast(nr_mo_ungerade_positive_secondary )
         call xc_mpi_bcast(nr_df_open                        )
         call xc_mpi_bcast(mo_eigenvalues_available          )
         call xc_mpi_bcast(mo_coef_available                 )

         if (nr_df_open > 0) then
            if (.not. xc_mpi_is_master()) then
               allocate(df_open(nr_df_open))
            end if
            call xc_mpi_bcast(df_open)
         end if

         if (mo_coef_available) then
            if (xc_mpi_is_master()) then
               size_mo_coef = size(mo_coef)
            end if
            call xc_mpi_bcast(size_mo_coef)
            if (.not. xc_mpi_is_master()) then
               allocate(mo_coef(size_mo_coef))
            end if
            call xc_mpi_bcast(mo_coef)
         end if

         if (mo_eigenvalues_available) then
            if (.not. xc_mpi_is_master()) then
               allocate(mo_eigenvalues(nr_mo))
            end if
            call xc_mpi_bcast(mo_eigenvalues)
         end if
      end if
#endif

   end subroutine

   function word_contains(word, substring)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: word
      character(*), intent(in) :: substring
!     --------------------------------------------------------------------------
      logical                  :: word_contains
!     --------------------------------------------------------------------------

      if (index(word, substring) > 0) then
         word_contains = .true.
      else
         word_contains = .false.
      end if

   end function

end module
