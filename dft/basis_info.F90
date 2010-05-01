
!> \brief   write/read basis set information to/from file
!> \author  radovan bast
!> \date    2010-05-01
!> \todo    mpi bcasting
!> \todo    complete this module so that dft_ao gets all info from file
!> \warning incomplete, work in progress

module basis_info

   implicit none

   public basis_info_write
   public basis_info_read
   public basis_info_delete

   private

!  -----------------------------------------------------------------------------
!  public non-allocatables

   integer,              public :: nr_shells
   integer,              public :: nr_cartesian_ao

!  -----------------------------------------------------------------------------
!  public allocatables

   integer, allocatable, public :: cartesian_deg(:)
   integer, allocatable, public :: cartesian_off(:)
   integer, allocatable, public :: spherical_deg(:)
   integer, allocatable, public :: spherical_off(:)
   integer, allocatable, public :: n_quantum_nr(:)
   integer, allocatable, public :: shell_center(:)
   logical, allocatable, public :: shell_is_spherical(:)
   real(8), allocatable, public :: shell_extent(:)
   integer, allocatable, public :: nr_nonzero_contraction_coefs(:)
   real(8), allocatable, public :: contraction_coef(:, :)
   real(8), allocatable, public :: primitive_exp(:, :)
   integer, allocatable, public :: l_exp(:)
   integer, allocatable, public :: m_exp(:)
   integer, allocatable, public :: n_exp(:)

!  -----------------------------------------------------------------------------
!  private

   character(*), parameter :: file_name = 'basis_info'
   integer,      parameter :: file_unit = 137
   logical                 :: file_exists

   real(8), parameter :: cell_size          = 2.25d0
!                                             cell_size*0.5d0*sqrt(3.0d0)
   real(8), parameter :: half_cell_diagonal = cell_size*0.866025404d0
   real(8), parameter :: ao_threshold       = 1.0d-15

!  this is probably eq. 8.2.25 in the pink book
   real(8), parameter :: g(10) = (/1.0d0, 1.3333d0, 1.6d0, 1.83d0, 2.03d0, 2.22d0, 2.39d0, 2.55d0, 2.70d0, 2.84d0/)

contains

!> \brief deallocate arrays
   subroutine basis_info_delete()

      nr_shells       = 0
      nr_cartesian_ao = 0

      if (allocated(cartesian_deg               )) deallocate(cartesian_deg               )
      if (allocated(cartesian_off               )) deallocate(cartesian_off               )
      if (allocated(spherical_deg               )) deallocate(spherical_deg               )
      if (allocated(spherical_off               )) deallocate(spherical_off               )
      if (allocated(n_quantum_nr                )) deallocate(n_quantum_nr                )
      if (allocated(shell_center                )) deallocate(shell_center                )
      if (allocated(shell_is_spherical          )) deallocate(shell_is_spherical          )
      if (allocated(nr_nonzero_contraction_coefs)) deallocate(nr_nonzero_contraction_coefs)
      if (allocated(contraction_coef            )) deallocate(contraction_coef            )
      if (allocated(primitive_exp               )) deallocate(primitive_exp               )
      if (allocated(shell_extent                )) deallocate(shell_extent                )
      if (allocated(l_exp                       )) deallocate(l_exp                       )
      if (allocated(m_exp                       )) deallocate(m_exp                       )
      if (allocated(n_exp                       )) deallocate(n_exp                       )

   end subroutine

!> \brief write basis set info to file
   subroutine basis_info_write()

!     --------------------------------------------------------------------------
      integer :: i, j, i_c, i_s
      integer :: ishell, icount, isym, iprim
      integer :: contracted_shell_nr
      real(8) :: r
      integer :: maxs
      integer :: maxc
      integer :: maxp
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

      inquire(file  = file_name, exist = file_exists)

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

      maxs = mxshel
      maxc = mxcont
      maxp = mxprim

      call allocate_arrays(maxs, &
                           maxc, &
                           maxp)

      nr_shells = kmax

      nr_cartesian_ao = 0
      shell_extent    = 0.0d0

      i_c = 1
      i_s = 1

      do ishell = 1, nr_shells

         cartesian_deg(ishell)      = kckt(ishell)
         cartesian_off(ishell)      = i_c - 1
         spherical_deg(ishell)      = khkt(ishell)
         spherical_off(ishell)      = i_s - 1
         n_quantum_nr(ishell)       = nhkt(ishell)
         shell_center(ishell)       = ncent(ishell)
         shell_is_spherical(ishell) = sphr(ishell)

         contracted_shell_nr = numcf(ishell)

         i = jstrt(ishell) + 1
         j = jstrt(ishell) + nuco(ishell)

         icount = 0
         do iprim = i, j
            if (dabs(priccf(iprim, contracted_shell_nr)) > 1.0d-30) then
               icount = icount + 1
               nr_nonzero_contraction_coefs(ishell) = icount
               contraction_coef(ishell, icount)     = priccf(iprim, contracted_shell_nr)
               primitive_exp(ishell, icount)        = priexp(iprim)

!              get shell extents
               r = (log(abs(priccf(iprim, contracted_shell_nr))) &
                  - log(ao_threshold))/priexp(iprim)

               if (shell_extent(ishell) < r) shell_extent(ishell) = r
            end if
         end do

         shell_extent(ishell) = sqrt(shell_extent(ishell))*g(n_quantum_nr(ishell)) &
                              + half_cell_diagonal

!        better square this once than sqrt p2 in each point
         shell_extent(ishell) = shell_extent(ishell)**2.0d0

         do isym = 0, maxopr
            if (iand(isym, istbao(ishell)) == 0) then
               nr_cartesian_ao = nr_cartesian_ao + cartesian_deg(ishell)
            end if
         end do

         call lmnval(n_quantum_nr(ishell), cartesian_deg(ishell), l_exp(i_c), m_exp(i_c), n_exp(i_c))

         i_c = i_c + cartesian_deg(ishell)
         i_s = i_s + spherical_deg(ishell)

      end do

      write(file_unit, *) maxs
      write(file_unit, *) maxc
      write(file_unit, *) maxp

      write(file_unit, *) nr_shells

      do ishell = 1, nr_shells
         write(file_unit, *) cartesian_deg(ishell),      &
                             cartesian_off(ishell),      &
                             spherical_deg(ishell),      &
                             spherical_off(ishell),      &
                             n_quantum_nr(ishell),       &
                             shell_center(ishell),       &
                             shell_is_spherical(ishell), &
                             shell_extent(ishell),       &
                             nr_nonzero_contraction_coefs(ishell)
         do i = 1, nr_nonzero_contraction_coefs(ishell)
            write(file_unit, *) contraction_coef(ishell, i), &
                                primitive_exp(ishell, i)
         end do
      end do

      write(file_unit, *) nr_cartesian_ao
      do i = 1, nr_cartesian_ao
         write(file_unit, *) l_exp(i), &
                             m_exp(i), &
                             n_exp(i)
      end do

      close(file_unit, status = 'keep')

   end subroutine

!> \brief read basis set info from file
   subroutine basis_info_read()

!     --------------------------------------------------------------------------
      integer :: i, ishell
      integer :: maxs
      integer :: maxc
      integer :: maxp
!     --------------------------------------------------------------------------

      inquire(file  = file_name, exist = file_exists)

      if (file_exists) then
         open(file_unit,            &
              file   = file_name,   &
              status = 'unknown',   &
              form   = 'formatted', &
              access = 'sequential')
         rewind(file_unit)
      else
         call quit('cannot read from file basis_info')
      end if

      read(file_unit, *) maxs
      read(file_unit, *) maxc
      read(file_unit, *) maxp

      call allocate_arrays(maxs, &
                           maxc, &
                           maxp)

      read(file_unit, *) nr_shells

      do ishell = 1, nr_shells
         read(file_unit, *) cartesian_deg(ishell),      &
                            cartesian_off(ishell),      &
                            spherical_deg(ishell),      &
                            spherical_off(ishell),      &
                            n_quantum_nr(ishell),       &
                            shell_center(ishell),       &
                            shell_is_spherical(ishell), &
                            shell_extent(ishell),       &
                            nr_nonzero_contraction_coefs(ishell)
         do i = 1, nr_nonzero_contraction_coefs(ishell)
            read(file_unit, *) contraction_coef(ishell, i), &
                               primitive_exp(ishell, i)
         end do
      end do

      read(file_unit, *) nr_cartesian_ao
      do i = 1, nr_cartesian_ao
         read(file_unit, *) l_exp(i), &
                            m_exp(i), &
                            n_exp(i)
      end do

      close(file_unit, status = 'keep')

   end subroutine

   subroutine allocate_arrays(maxs, &
                              maxc, &
                              maxp)

!     --------------------------------------------------------------------------
      integer, intent(in) :: maxs
      integer, intent(in) :: maxc
      integer, intent(in) :: maxp
!     --------------------------------------------------------------------------

      allocate(cartesian_deg(maxs))
      allocate(cartesian_off(maxs))
      allocate(spherical_deg(maxs))
      allocate(spherical_off(maxs))
      allocate(n_quantum_nr(maxs))
      allocate(shell_center(maxs))
      allocate(shell_is_spherical(maxs))
      allocate(nr_nonzero_contraction_coefs(maxs))
      allocate(contraction_coef(maxs, maxc))
      allocate(primitive_exp(maxs, maxc))
      allocate(shell_extent(maxs))
      allocate(l_exp(maxp))
      allocate(m_exp(maxp))
      allocate(n_exp(maxp))

   end subroutine

end module
