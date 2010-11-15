!dirac_copyright_start
!dirac_copyright_end

module basis_info_ng

   implicit none

   public basis_info_write_ng

   private

   integer              :: nr_shells(0:2)
   integer, parameter   :: nr_types = 2
   integer              :: max_nr_primitives(0:2)
   integer              :: first_shell(0:2)
   integer              :: last_shell(0:2)
   logical              :: is_spherical(0:2)
   integer, allocatable :: shell_type(:)
   integer, allocatable :: shell_center(:)
   integer, allocatable :: l_quantum_nr(:)
   integer, allocatable :: nr_primitives(:)
   real(8), allocatable :: primitive_exp(:, :)
   real(8), allocatable :: contraction_coef(:, :)

   character(*), parameter :: file_name = 'basis_info_ng'
   integer,      parameter :: file_unit = 137
   logical                 :: file_exists

   character(*), parameter :: type_label(0:2) = (/'total', 'large', 'small'/)

contains

   subroutine basis_info_write_ng()

!     --------------------------------------------------------------------------
      integer :: i, j, i_c, i_s, ic
      integer :: ishell, icount, isym, iprim, irep
!     --------------------------------------------------------------------------

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "shells.h"
#include "symmet.h"
#include "aovec.h"
#include "primit.h"
#include "sphtrm.h"
#include "dcbbas.h"
#include "pincom.h"
#include "dgroup.h"
#include "dcbham.h"

!     call basis_info_delete()

      inquire(file = file_name, exist = file_exists)

      if (file_exists) then
         close(file_unit, status = 'keep')
         return
      else
         open(file_unit,            &
              file   = file_name,   &
              status = 'new',       &
              form   = 'formatted', &
              access = 'sequential')
         rewind(file_unit)
      end if

      nr_shells(1) = nlrgsh
      nr_shells(2) = kmax - nr_shells(1)
      nr_shells(0) = nr_shells(1) + nr_shells(2)

      first_shell(1) = 1
      last_shell(1)  = nr_shells(1)
      first_shell(2) = nr_shells(1)  + 1
      last_shell(2)  = last_shell(1) + nr_shells(2)

      allocate(l_quantum_nr(nr_shells(0)))
      allocate(shell_center(nr_shells(0)))
      do ishell = 1, nr_shells(0)
         l_quantum_nr(ishell) = nhkt(ishell)
         shell_center(ishell) = ncent(ishell)
      end do

      allocate(shell_type(nr_shells(0)))
      do ishell = 1, nr_shells(0)
         do ic = 1, nr_types
            if (ishell >= first_shell(ic) .and. ishell <= last_shell(ic)) then
               shell_type(ishell) = ic
            end if
         end do
      end do

      max_nr_primitives = 0
      do ishell = 1, nr_shells(0)
         i = jstrt(ishell) + 1
         j = jstrt(ishell) + nuco(ishell)
         icount = 0
         do iprim = i, j
            if (dabs(priccf(iprim, numcf(ishell))) > tiny(0.0d0)) then
               icount = icount + 1
               if (icount > max_nr_primitives(shell_type(ishell))) then
                  max_nr_primitives(shell_type(ishell)) = icount
               end if
            end if
         end do
      end do
      max_nr_primitives(0) = maxval(max_nr_primitives)

      allocate(nr_primitives(nr_shells(0)))
      allocate(primitive_exp(nr_shells(0), max_nr_primitives(0)))
      allocate(contraction_coef(nr_shells(0), max_nr_primitives(0)))

      do ishell = 1, nr_shells(0)
         i = jstrt(ishell) + 1
         j = jstrt(ishell) + nuco(ishell)
         icount = 0
         do iprim = i, j
            if (dabs(priccf(iprim, numcf(ishell))) > tiny(0.0d0)) then
               icount = icount + 1
               nr_primitives(ishell)            = icount
               contraction_coef(ishell, icount) = priccf(iprim, numcf(ishell))
               primitive_exp(ishell, icount)    = priexp(iprim)
            end if
         end do
      end do

      do ic = 1, nr_types

         is_spherical(ic) = .false.
         do ishell = first_shell(ic), last_shell(ic)
            if (sphr(ishell)) then
               is_spherical(ic) = .true.
            end if
         end do
         
         write(file_unit, *) '*** basis_'//type_label(ic)
         
         write(file_unit, *) 'is_spherical'
         write(file_unit, *)  is_spherical(ic)
         
         write(file_unit, *) 'nr_shells'
         write(file_unit, *)  nr_shells(ic)
         
         write(file_unit, *) 'l_quantum_nr'
         do ishell = first_shell(ic), last_shell(ic)
            write(file_unit, *) l_quantum_nr(ishell)
         end do
         
         write(file_unit, *) 'shell_center'
         do ishell = first_shell(ic), last_shell(ic)
            write(file_unit, *) shell_center(ishell)
         end do
         
         write(file_unit, *) 'max_nr_primitives'
         write(file_unit, *)  max_nr_primitives(ic)
         
         write(file_unit, *) 'nr_primitives'
         do ishell = first_shell(ic), last_shell(ic)
            write(file_unit, *) nr_primitives(ishell)
         end do
         
         write(file_unit, *) 'primitive_exp_and_contraction_coef'
         do ishell = first_shell(ic), last_shell(ic)
            do i = 1, nr_primitives(ishell)
               write(file_unit, '(2f26.16)') primitive_exp(ishell, i), contraction_coef(ishell, i)
            end do
         end do

      end do

      deallocate(l_quantum_nr)
      deallocate(shell_center)
      deallocate(shell_type)
      deallocate(nr_primitives)
      deallocate(primitive_exp)
      deallocate(contraction_coef)

      close(file_unit, status = 'keep')

   end subroutine

end module
