
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
   integer,              public :: nr_boson_sym
   integer,              public :: nr_cartesian_ao
   integer,              public :: nr_so_ao

!  -----------------------------------------------------------------------------
!  public allocatables

   integer, allocatable, public :: cartesian_deg(:)
   integer, allocatable, public :: cartesian_off(:)
   integer, allocatable, public :: spherical_deg(:)
   integer, allocatable, public :: spherical_off(:)
   integer, allocatable, public :: l_quantum_nr(:)
   integer, allocatable, public :: shell_center(:)
   real(8), allocatable, public :: shell_center_xyz(:, :)
!>                                 shell stabiliser: sym that does not move center
   integer, allocatable, public :: shell_stabilizer(:)
   logical, allocatable, public :: shell_is_spherical(:)
   real(8), allocatable, public :: shell_extent(:)
   integer, allocatable, public :: nr_nonzero_contraction_coefs(:)
   real(8), allocatable, public :: contraction_coef(:, :)
   real(8), allocatable, public :: primitive_exp(:, :)
   integer, allocatable, public :: l_exp(:)
   integer, allocatable, public :: m_exp(:)
   integer, allocatable, public :: n_exp(:)
   integer, allocatable, public :: i_so_ao(:, :)
   real(8), allocatable, public :: f_so_ao(:)
   real(8), allocatable, public :: bitstring_parity(:)
   real(8), allocatable, public :: cartesian_to_spherical(:)
   integer, allocatable, public :: spherical_address(:)

!  -----------------------------------------------------------------------------
!  private

   character(*), parameter :: file_name = 'basis_info'
   integer,      parameter :: file_unit = 137
   logical                 :: file_exists

   logical                 :: only_large
   logical,      parameter :: force_small_ao = .false.

   real(8),      parameter :: cell_size          = 2.25d0
!>                                                 cell_size*0.5d0*sqrt(3.0d0)
   real(8),      parameter :: half_cell_diagonal = cell_size*0.866025404d0
   real(8),      parameter :: ao_threshold       = 1.0d-15

!>                                            this is probably eq. 8.2.25 in the pink book
   real(8),      parameter :: g(10) = (/1.0d0, 1.3333d0, 1.6d0, 1.83d0, 2.03d0, 2.22d0, 2.39d0, 2.55d0, 2.70d0, 2.84d0/)

   integer,      parameter :: max_nr_so_ao = 10000

contains

!> \brief deallocate arrays
   subroutine basis_info_delete()

      nr_shells       = 0
      nr_boson_sym    = 0
      nr_cartesian_ao = 0
      nr_so_ao        = 0

      if (allocated(cartesian_deg               )) deallocate(cartesian_deg               )
      if (allocated(cartesian_off               )) deallocate(cartesian_off               )
      if (allocated(spherical_deg               )) deallocate(spherical_deg               )
      if (allocated(spherical_off               )) deallocate(spherical_off               )
      if (allocated(l_quantum_nr                )) deallocate(l_quantum_nr                )
      if (allocated(shell_center                )) deallocate(shell_center                )
      if (allocated(shell_center_xyz            )) deallocate(shell_center_xyz            )
      if (allocated(shell_stabilizer            )) deallocate(shell_stabilizer            )
      if (allocated(shell_is_spherical          )) deallocate(shell_is_spherical          )
      if (allocated(shell_extent                )) deallocate(shell_extent                )
      if (allocated(nr_nonzero_contraction_coefs)) deallocate(nr_nonzero_contraction_coefs)
      if (allocated(contraction_coef            )) deallocate(contraction_coef            )
      if (allocated(primitive_exp               )) deallocate(primitive_exp               )
      if (allocated(l_exp                       )) deallocate(l_exp                       )
      if (allocated(m_exp                       )) deallocate(m_exp                       )
      if (allocated(n_exp                       )) deallocate(n_exp                       )
      if (allocated(i_so_ao                     )) deallocate(i_so_ao                     )
      if (allocated(f_so_ao                     )) deallocate(f_so_ao                     )
      if (allocated(bitstring_parity            )) deallocate(bitstring_parity            )
      if (allocated(cartesian_to_spherical      )) deallocate(cartesian_to_spherical      )
      if (allocated(spherical_address           )) deallocate(spherical_address           )

   end subroutine

!> \brief write basis set info to file
   subroutine basis_info_write()

!     --------------------------------------------------------------------------
      integer :: i, j, i_c, i_s, ic
      integer :: ishell, icount, isym, iprim, irep
      integer :: iso, iso_bs, iao
      integer :: contracted_shell_nr, nr_ao_irep, nshell, mul, n
      real(8) :: r, f
      integer :: ns
      integer :: nc
      integer :: np
      integer :: nsoao
      integer :: ncs
      integer :: nq
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

      only_large = .false.
      if (mc == 1) then
         only_large = .true.
      end if
      if (only_large .and. force_small_ao) then
         only_large = .false.
      end if

      if (only_large) then
         nr_shells = nlrgsh
      else
         nr_shells = kmax
      end if

      nr_boson_sym = nbsym

      ns    = nr_shells
      nc    = mxcont
      np    = mxprim
      nsoao = max_nr_so_ao
      ncs   = ncsp
      nq    = mxqn

      call allocate_arrays(ns,    &
                           nc,    &
                           np,    &
                           nsoao, &
                           ncs,   &
                           nq)

      bitstring_parity       = pt
      cartesian_to_spherical = csp
      spherical_address      = ispadr

      nr_cartesian_ao = 0
      shell_extent    = 0.0d0

      i_c = 1
      i_s = 1

      do ishell = 1, nr_shells

         cartesian_deg(ishell)       = kckt(ishell)
         cartesian_off(ishell)       = i_c - 1
         spherical_deg(ishell)       = khkt(ishell)
         spherical_off(ishell)       = i_s - 1
         l_quantum_nr(ishell)        = nhkt(ishell)
         shell_center(ishell)        = ncent(ishell)
         shell_center_xyz(1, ishell) = cent(ishell, 1, 1)
         shell_center_xyz(2, ishell) = cent(ishell, 2, 1)
         shell_center_xyz(3, ishell) = cent(ishell, 3, 1)
         shell_stabilizer(ishell)    = istbao(ishell)
         shell_is_spherical(ishell)  = sphr(ishell)

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

         shell_extent(ishell) = sqrt(shell_extent(ishell))*g(l_quantum_nr(ishell)) &
                              + half_cell_diagonal

!        better square this once than sqrt p2 in each point
         shell_extent(ishell) = shell_extent(ishell)**2.0d0

         do isym = 0, maxopr
            if (iand(isym, istbao(ishell)) == 0) then
               nr_cartesian_ao = nr_cartesian_ao + cartesian_deg(ishell)
            end if
         end do

         call lmnval(l_quantum_nr(ishell), cartesian_deg(ishell), l_exp(i_c), m_exp(i_c), n_exp(i_c))

         i_c = i_c + cartesian_deg(ishell)
         i_s = i_s + spherical_deg(ishell)

      end do

!     ao to so business

      i = 1
      nr_so_ao = 0
      do irep = 0, maxrep
         icount = 0
         nr_ao_irep = nbbas(irep, 1) + nbbas(irep, 2)
         do iso = i, i + nr_ao_irep - 1

            ishell = iand(ishft(ipind(iso), -16), 65535)
            nshell = iand(ishft(ipind(iso),  -8),   255)

            mul = istbao(ishell)
            n   = l_quantum_nr(ishell)

            iao = kstrt(ishell) + nshell

            icount = icount + 1
            if (icount > nbbas(irep, 1)) then
!              is small
               ic = 2
            else
!              is large
               ic = 1
            end if

!           do not plan to AO-SO transform small component functions
!           in case of zero small-small metric (.LEVY-LEBLOND and .X2C4)
            if (.not. (only_large .and. (ic == 2))) then

!              shift by difference boson sorted minus boson unsorted
               iso_bs = iso + ibbas(irep, ic) - icos(irep + 1, ic)

               do isym = 0, maxopr
                  if (iand(isym, mul) == 0) then

                     f = pt(iand(isym, ieor(irep, isymao(n, nshell))))

                     nr_so_ao = nr_so_ao + 1

                     if (nr_so_ao > max_nr_so_ao) then
                        call quit('increase max_nr_so_ao')
                     end if

                     i_so_ao(1, nr_so_ao) = iso_bs
                     i_so_ao(2, nr_so_ao) = iao
                     f_so_ao(nr_so_ao)    = f

                     iao = iao + spherical_deg(ishell)

                  endif
               end do

            end if

         end do
         i = i + nr_ao_irep
      end do
      nsoao = nr_so_ao

      write(file_unit, *) ns   
      write(file_unit, *) nc   
      write(file_unit, *) np   
      write(file_unit, *) nsoao
      write(file_unit, *) ncs  
      write(file_unit, *) nq   

      write(file_unit, *) nr_shells      
      write(file_unit, *) nr_boson_sym   
      write(file_unit, *) nr_cartesian_ao
      write(file_unit, *) nr_so_ao       

      do ishell = 1, nr_shells
         write(file_unit, *) cartesian_deg(ishell),       &
                             cartesian_off(ishell),       &
                             spherical_deg(ishell),       &
                             spherical_off(ishell),       &
                             l_quantum_nr(ishell),        &
                             shell_center(ishell),        &
                             shell_center_xyz(1, ishell), &
                             shell_center_xyz(2, ishell), &
                             shell_center_xyz(3, ishell), &
                             shell_stabilizer(ishell),    &
                             shell_is_spherical(ishell),  &
                             shell_extent(ishell),        &
                             nr_nonzero_contraction_coefs(ishell)


         do i = 1, nr_nonzero_contraction_coefs(ishell)
            write(file_unit, *) contraction_coef(ishell, i), &
                                primitive_exp(ishell, i)
         end do
      end do
      do i = 1, nr_cartesian_ao
         write(file_unit, *) l_exp(i), &
                             m_exp(i), &
                             n_exp(i)
      end do
      do i = 1, nr_so_ao
         write(file_unit, *) i_so_ao(1, i), &
                             i_so_ao(2, i), &
                             f_so_ao(i)
      end do
      do i = 0, 7
         write(file_unit, *) bitstring_parity(i)
      end do
      do i = 1, ncs
         write(file_unit, *) cartesian_to_spherical(i)
      end do
      do i = 1, nq
         write(file_unit, *) spherical_address(i)
      end do

      close(file_unit, status = 'keep')

   end subroutine

!> \brief read basis set info from file
   subroutine basis_info_read()

!     --------------------------------------------------------------------------
      integer :: i, ishell
      integer :: ns
      integer :: nc
      integer :: np
      integer :: nsoao
      integer :: ncs
      integer :: nq
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

      read(file_unit, *) ns   
      read(file_unit, *) nc   
      read(file_unit, *) np   
      read(file_unit, *) nsoao
      read(file_unit, *) ncs  
      read(file_unit, *) nq   

      call allocate_arrays(ns,    &
                           nc,    &
                           np,    &
                           nsoao, &
                           ncs,   &
                           nq)

      read(file_unit, *) nr_shells      
      read(file_unit, *) nr_boson_sym   
      read(file_unit, *) nr_cartesian_ao
      read(file_unit, *) nr_so_ao       

      do ishell = 1, nr_shells
         read(file_unit, *) cartesian_deg(ishell),       &
                            cartesian_off(ishell),       &
                            spherical_deg(ishell),       &
                            spherical_off(ishell),       &
                            l_quantum_nr(ishell),        &
                            shell_center(ishell),        &
                            shell_center_xyz(1, ishell), &
                            shell_center_xyz(2, ishell), &
                            shell_center_xyz(3, ishell), &
                            shell_stabilizer(ishell),    &
                            shell_is_spherical(ishell),  &
                            shell_extent(ishell),        &
                            nr_nonzero_contraction_coefs(ishell)


         do i = 1, nr_nonzero_contraction_coefs(ishell)
            read(file_unit, *) contraction_coef(ishell, i), &
                               primitive_exp(ishell, i)
         end do
      end do
      do i = 1, nr_cartesian_ao
         read(file_unit, *) l_exp(i), &
                             m_exp(i), &
                             n_exp(i)
      end do
      do i = 1, nr_so_ao
         read(file_unit, *) i_so_ao(1, i), &
                            i_so_ao(2, i), &
                            f_so_ao(i)
      end do
      do i = 0, 7
         read(file_unit, *) bitstring_parity(i)
      end do
      do i = 1, ncs
         read(file_unit, *) cartesian_to_spherical(i)
      end do
      do i = 1, nq
         read(file_unit, *) spherical_address(i)
      end do

      close(file_unit, status = 'keep')

   end subroutine

   subroutine allocate_arrays(ns,    &
                              nc,    &
                              np,    &
                              nsoao, &
                              ncs,   &
                              nq)

!     --------------------------------------------------------------------------
      integer, intent(in) :: ns
      integer, intent(in) :: nc
      integer, intent(in) :: np
      integer, intent(in) :: nsoao
      integer, intent(in) :: ncs
      integer, intent(in) :: nq
!     --------------------------------------------------------------------------

      allocate(cartesian_deg(ns))
      allocate(cartesian_off(ns))
      allocate(spherical_deg(ns))
      allocate(spherical_off(ns))
      allocate(l_quantum_nr(ns))
      allocate(shell_center(ns))
      allocate(shell_center_xyz(3, ns))
      allocate(shell_stabilizer(ns))
      allocate(shell_is_spherical(ns))
      allocate(shell_extent(ns))
      allocate(nr_nonzero_contraction_coefs(ns))
      allocate(contraction_coef(ns, nc))
      allocate(primitive_exp(ns, nc))
      allocate(l_exp(np))
      allocate(m_exp(np))
      allocate(n_exp(np))
      allocate(i_so_ao(2, nsoao))
      allocate(f_so_ao(nsoao))
      allocate(bitstring_parity(0:7))
      allocate(cartesian_to_spherical(ncs))
      allocate(spherical_address(nq))

   end subroutine

end module
