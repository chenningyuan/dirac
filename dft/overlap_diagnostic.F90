module overlap_diagnostic
! pbht mo overlap diagnostic

   use memory_allocator

   implicit none

   public integrate_lambda
   public overlap_diagnostic_init
   public overlap_diagnostic_write
   public overlap_diagnostic_print

   private

#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h" 

   integer                 :: nr_mo
   real(8), allocatable    :: mo_coef(:)
   real(8), allocatable    :: lambda_matrix(:, :)

   character(*), parameter :: file_name = 'lambda_matrix'
   integer,      parameter :: file_unit = 137
   logical                 :: file_exists

contains

   subroutine overlap_diagnostic_init(nz, nr_mo_in, nr_ao)

!     --------------------------------------------------------------------------
      integer, intent(in)    :: nz
      integer, intent(in)    :: nr_mo_in
      integer, intent(in)    :: nr_ao
!     --------------------------------------------------------------------------

      nr_mo = nr_mo_in

    ! call alloc(mo_coef, nr_ao*nr_mo*nz)
      call alloc(mo_coef, nr_ao*nr_ao*nz)
      call alloc(lambda_matrix, nr_mo, nr_mo)
      call read_mo_coef(mo_coef)

   end subroutine

   subroutine overlap_diagnostic_write()

      inquire(file = file_name, exist = file_exists)

      if (file_exists) then
!        panic
!        close(file_unit, status = 'keep')
!        return
      else
         open(file_unit,              &
              file   = file_name,     &
              status = 'new',         &
              form   = 'unformatted', &
              access = 'sequential')
         rewind(file_unit)
      end if

      write(file_unit) nr_mo
      write(file_unit) lambda_matrix
      close(file_unit, status = 'keep')

      if (allocated(mo_coef))       call dealloc(mo_coef)
      if (allocated(lambda_matrix)) call dealloc(lambda_matrix)

   end subroutine

   subroutine overlap_diagnostic_print()

      inquire(file = file_name, exist = file_exists)

      if (file_exists) then
         open(file_unit,              &
              file   = file_name,     &
              status = 'old',         &
              form   = 'unformatted', &
              access = 'sequential')
         rewind(file_unit)
      else
!        panic
!        close(file_unit, status = 'keep')
!        return
      end if

      read(file_unit) nr_mo

      call alloc(lambda_matrix, nr_mo, nr_mo)

      read(file_unit) lambda_matrix

      close(file_unit, status = 'keep')

      write(*, *) 'blub lambda'
      call prqmat(lambda_matrix,  &
                  norbt,          &
                  norbt,          &
                  norbt,          &
                  norbt,          &
                  1,              &
                  (/1, 2, 3, 4/), &
                  6)

      call dealloc(lambda_matrix)

   end subroutine

   subroutine integrate_lambda(w, nz, nr_mo, nr_ao, ao, buffer)

!     --------------------------------------------------------------------------
      real(8), intent(in)    :: w
      integer, intent(in)    :: nz
      integer, intent(in)    :: nr_mo
      integer, intent(in)    :: nr_ao
      real(8), intent(in)    :: ao(*)
      real(8), intent(inout) :: buffer(*)
!     --------------------------------------------------------------------------
      integer                :: p, q
      integer                :: i, k
      integer                :: iz
!     --------------------------------------------------------------------------

      i = 0
      do iz = 1, nz
         do p = 1, nr_mo
            buffer(p) = 0.0d0
            do k = 1, nr_ao
               i = i + 1
               buffer(p) = buffer(p) + mo_coef(i)*ao(k)
            end do
         end do

         do p = 1, nr_mo
            do q = 1, nr_mo
              lambda_matrix(p, q) = lambda_matrix(p, q) + w*dabs(buffer(p)*buffer(q))
            end do
         end do
      end do

   end subroutine

   subroutine read_mo_coef(mo_coef)

!     --------------------------------------------------------------------------
      real(8), intent(out) :: mo_coef(*)
!     --------------------------------------------------------------------------
      real(8), allocatable :: eigenvalues(:)
      integer, parameter   :: file_unit = 66
!     --------------------------------------------------------------------------

      call alloc(eigenvalues, norbt)

      call opnfil(file_unit, 'DFCOEF', 'OLD', 'read_mo_coef')

      call reacmo(file_unit,   &
                  'DFCOEF',    &
                  mo_coef,     &
                  eigenvalues, &
                  (/0/),       &
                  (/0.0d0/),   &
                  6,           &
                  (/0.0d0/),   &
                  (/0/),       &
                  (/0/))

      close(file_unit, status = 'keep')

      call dealloc(eigenvalues)

   end subroutine

end module
