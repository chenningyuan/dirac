
!    radovan bast <radovan.bast@uit.no> 2011-04-20
!    module for the export/import of basis set information

!    general design
!    --------------
!  - file is formatted
!  - file is read in free format
!  - sections have labels
!  - order of sections does not matter
!  - section symmetry can be omitted if irrelevant (e.g. C1)
!  - normalization is absorbed in contraction coefficients
!    (we can discuss this point)
!    AOs are normalized based on angular momentum and exponent
!    such that each shell has a common normalization
!    with this < AO | AO > =   1 for s, px, py, pz, dxy, dxz, dyz, fxyz
!                              3 for dxx, dyy, dzz, fxxy, ...
!                             15     fxxx, fyyy, fzzz, gxxxy, ...
!                            105     gxxxx, ...
!                            ...     ...
!                              9 for gxxyy, ...
!                            etc     ...
!  - ireps go from 0 to nr_ireps-1

!    todo and known issues
!    ---------------------
!  - remove shell stabilizers and refer to center stabilizers
!  - move all symmetry related info to *** symmetry, also center degeneracy and
!    sym independent nuclei
!  - make it possible to read from basis_info without *** symmetry assuming C1
!  - remove dependency on dft_cfg_betasigma

module basis_info

   implicit none

   public basis_info_write
   public basis_info_read

   integer, public, parameter :: l_max = 20

   integer, public :: ao_nr(0:7, 2)
   integer, public :: ao_off(0:7, 2)
   integer, public :: ao_off_hermit(8, 2)
   integer, public :: c2s_index(l_max)
   integer, public :: iq(0:3)
   integer, public :: irep_of_axes(3)
   integer, public :: irep_of_rotations(3)
   integer, public :: llss_nr_ao(2, 16, 0:3, 0:7)
   integer, public :: llss_nr_density_blocks(0:3, 0:7)
   integer, public :: llss_st_ao(2, 16, 0:3, 0:7)
   integer, public :: lssl_nr_ao(2, 16, 0:3, 0:7)
   integer, public :: lssl_nr_density_blocks(0:3, 0:7)
   integer, public :: lssl_st_ao(2, 16, 0:3, 0:7)
   integer, public :: nr_ao
   integer, public :: nr_ao_cartesian
   integer, public :: nr_ao_spherical
   integer, public :: nr_ireps
   integer, public :: nr_shells
   integer, public :: nr_shells_large
   integer, public :: nr_shells_small
   integer, public :: nr_so_ao
   integer, public :: nr_sym_dependent_centers
   integer, public :: nr_sym_independent_centers
   integer, public :: pq_to_uq(4, 0:7)
   integer, public :: uq_to_pq(4, 0:7)
   logical, public :: basis_is_spherical
   logical, public :: use_only_large
   real(8), public :: llss_prefactor(16, 0:3, 0:7)

   integer, allocatable, public :: cartesian_deg(:)
   integer, allocatable, public :: cartesian_off(:)
   integer, allocatable, public :: i_so_ao(:, :)
   integer, allocatable, public :: l_quantum_nr(:)
   integer, allocatable, public :: nr_nonzero_contraction_coefs(:)
   integer, allocatable, public :: shell_center(:)
   integer, allocatable, public :: shell_stabilizing_irep(:)
   integer, allocatable, public :: spherical_deg(:)
   integer, allocatable, public :: spherical_off(:)
   integer, allocatable, public :: sym_independent_center_deg(:)
   real(8), allocatable, public :: c2s_mat(:)
   real(8), allocatable, public :: contraction_coef(:, :)
   real(8), allocatable, public :: f_so_ao(:)
   real(8), allocatable, public :: primitive_exp(:, :)
   real(8), allocatable, public :: sym_dependent_center_charge(:)
   real(8), allocatable, public :: sym_dependent_center_xyz(:, :)
   real(8), allocatable, public :: sym_independent_center_charge(:)
   real(8), allocatable, public :: sym_independent_center_xyz(:, :)

   private

   character(*), parameter :: file_name = 'basis_info'
   integer,      parameter :: file_unit = 137
   logical                 :: file_exists

contains

   subroutine basis_info_clear()

      ao_nr                      = 0
      ao_off                     = 0
      ao_off_hermit              = 0
      c2s_index                  = 0
      iq                         = 0
      irep_of_axes               = 0
      irep_of_rotations          = 0
      llss_nr_ao                 = 0
      llss_nr_density_blocks     = 0
      llss_st_ao                 = 0
      lssl_nr_ao                 = 0
      lssl_nr_density_blocks     = 0
      lssl_st_ao                 = 0
      nr_ao                      = 0
      nr_ao_cartesian            = 0
      nr_ao_spherical            = 0
      nr_ireps                   = 0
      nr_shells                  = 0
      nr_shells_large            = 0
      nr_shells_small            = 0
      nr_so_ao                   = 0
      nr_sym_dependent_centers   = 0
      nr_sym_independent_centers = 0
      pq_to_uq                   = 0
      uq_to_pq                   = 0
      basis_is_spherical         = .false.
      use_only_large             = .false.
      llss_prefactor             = 0.0d0

      if (allocated(cartesian_deg                )) deallocate(cartesian_deg                )
      if (allocated(cartesian_off                )) deallocate(cartesian_off                )
      if (allocated(i_so_ao                      )) deallocate(i_so_ao                      )
      if (allocated(l_quantum_nr                 )) deallocate(l_quantum_nr                 )
      if (allocated(nr_nonzero_contraction_coefs )) deallocate(nr_nonzero_contraction_coefs )
      if (allocated(shell_center                 )) deallocate(shell_center                 )
      if (allocated(shell_stabilizing_irep       )) deallocate(shell_stabilizing_irep       )
      if (allocated(spherical_deg                )) deallocate(spherical_deg                )
      if (allocated(spherical_off                )) deallocate(spherical_off                )
      if (allocated(sym_independent_center_deg   )) deallocate(sym_independent_center_deg   )
      if (allocated(c2s_mat                      )) deallocate(c2s_mat                      )
      if (allocated(contraction_coef             )) deallocate(contraction_coef             )
      if (allocated(f_so_ao                      )) deallocate(f_so_ao                      )
      if (allocated(primitive_exp                )) deallocate(primitive_exp                )
      if (allocated(sym_dependent_center_charge  )) deallocate(sym_dependent_center_charge  )
      if (allocated(sym_dependent_center_xyz     )) deallocate(sym_dependent_center_xyz     )
      if (allocated(sym_independent_center_charge)) deallocate(sym_independent_center_charge)
      if (allocated(sym_independent_center_xyz   )) deallocate(sym_independent_center_xyz   )

   end subroutine

   subroutine basis_info_write()

!     --------------------------------------------------------------------------
      integer                 :: i, j, ic
      integer                 :: ishell, icount, iprim, irep, icent
      integer, allocatable    :: nr_primitives(:)
      character(1), parameter :: ls(2) = (/'L', 'S'/)
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
#include "nuclei.h"

      call basis_info_clear()

!     always regenerate this file
      inquire(file = file_name, exist = file_exists)
      if (file_exists) then
         open(file_unit,            &
              file   = file_name,   &
              status = 'unknown',   &
              form   = 'formatted', &
              access = 'sequential')
         close(file_unit, status = 'delete')
      end if

      open(file_unit,            &
           file   = file_name,   &
           status = 'new',       &
           form   = 'formatted', &
           access = 'sequential')
      rewind(file_unit)

      nr_ireps  = nbsym
      nr_shells = kmax

      allocate(l_quantum_nr(nr_shells))
      allocate(shell_center(nr_shells))
      allocate(shell_stabilizing_irep(nr_shells))
      do ishell = 1, nr_shells
         l_quantum_nr(ishell)           = nhkt(ishell) - 1
         shell_center(ishell)           = ncent(ishell)
         shell_stabilizing_irep(ishell) = istbao(ishell)
      end do

      allocate(nr_primitives(nr_shells))
      do ishell = 1, nr_shells
         i = jstrt(ishell) + 1
         j = jstrt(ishell) + nuco(ishell)
         icount = 0
         do iprim = i, j
            if (dabs(priccf(iprim, numcf(ishell))) > tiny(0.0d0)) then
               icount = icount + 1
               nr_primitives(ishell) = icount
            end if
         end do
      end do

      allocate(primitive_exp(nr_shells,    sum(nr_primitives)))
      allocate(contraction_coef(nr_shells, sum(nr_primitives)))
      do ishell = 1, nr_shells
         i = jstrt(ishell) + 1
         j = jstrt(ishell) + nuco(ishell)
         icount = 0
         do iprim = i, j
            if (dabs(priccf(iprim, numcf(ishell))) > tiny(0.0d0)) then
               icount = icount + 1
               contraction_coef(ishell, icount) = priccf(iprim, numcf(ishell))
               primitive_exp(ishell, icount)    = priexp(iprim)
            end if
         end do
      end do

      write(file_unit, '(a)') '*** geometry'

      write(file_unit, '(4x, a)') 'nr_sym_independent_centers'
      write(file_unit, '(4x, i6)') nucind

      write(file_unit, '(4x, a)') 'charges_and_coordinates_independent'
      do icent = 1, nucind
         write(file_unit, '(4f27.16)') charge(icent),  &
                                       cord(1, icent), &
                                       cord(2, icent), &
                                       cord(3, icent)
      end do

      write(file_unit, '(4x, a)') 'degeneracy_of_center'
      do icent = 1, nucind
         write(file_unit, '(4x, 2i6)') icent, nucdeg(icent)
      end do

      write(file_unit, '(4x, a)') 'nr_sym_dependent_centers'
      write(file_unit, '(4x, i6)') nucdep

      write(file_unit, '(4x, a)') 'charges_and_coordinates_dependent'
      do icent = 1, nucind
         if (mult(istbnu(icent)) == 1) then
               write(file_unit, '(4f27.16)') charge(icent),  &
                                             cord(1, icent), &
                                             cord(2, icent), &
                                             cord(3, icent)
         else
            do irep = 0, nr_ireps - 1
               if (iand(irep, istbnu(icent)) == 0) then
                  write(file_unit, '(4f27.16)') charge(icent),                               &
                                                cord(1, icent)*pt(iand(isymax(1, 1), irep)), &
                                                cord(2, icent)*pt(iand(isymax(2, 1), irep)), &
                                                cord(3, icent)*pt(iand(isymax(3, 1), irep))
               end if
            end do
         end if
      end do

      write(file_unit, '(a)') '*** basis'

      basis_is_spherical = .false.
      do ishell = 1, kmax
         if (sphr(ishell)) then
            basis_is_spherical = .true.
         end if
      end do
      write(file_unit, '(4x, a)') 'is_spherical'
      write(file_unit, '(4x, l6)') basis_is_spherical

      use_only_large = .false.
      if (mc == 1) then
         use_only_large = .true.
      end if
      write(file_unit, '(4x, a)') 'use_only_large'
      write(file_unit, '(4x, l6)') use_only_large

      write(file_unit, '(4x, a)') 'nr_primitive_exp'
      write(file_unit, '(4x, i6)') sum(nr_primitives)

      write(file_unit, '(4x, 5a6, 2a26)') '# cent', &
                                          'L/S',    &
                                          'shell',  &
                                          'l',      &
                                          'irep',   &
                                          'prim',   &
                                          'coef'
      do ishell = 1, kmax
         if (ishell > nlrgsh) then
            ic = 2
         else
            ic = 1
         end if
         do i = 1, nr_primitives(ishell)
            write(file_unit, '(4x, i6, a6, 3i6, 2e26.16)') shell_center(ishell),           &
                                                           ls(ic),                         &
                                                           ishell,                         &
                                                           l_quantum_nr(ishell),           &
                                                           shell_stabilizing_irep(ishell), &
                                                           primitive_exp(ishell, i),       &
                                                           contraction_coef(ishell, i)
         end do
      end do
      deallocate(nr_primitives)

      write(file_unit, '(a)') '*** symmetry'

      write(file_unit, '(4x, a)') 'nr_ireps'
      write(file_unit, '(4x, i6)') nr_ireps

      write(file_unit, '(4x, a)') 'irep_of_axes'
      write(file_unit, '(4x, i6)') isymax(1, 1)
      write(file_unit, '(4x, i6)') isymax(2, 1)
      write(file_unit, '(4x, i6)') isymax(3, 1)

      write(file_unit, '(4x, a)') 'irep_of_rotations'
      write(file_unit, '(4x, i6)') isymax(1, 2)
      write(file_unit, '(4x, i6)') isymax(2, 2)
      write(file_unit, '(4x, i6)') isymax(3, 2)

      write(file_unit, '(4x, 7a8)') '# irep',  &
                                    'nr-L',    &
                                    'nr-S',    &
                                    'D-off-L', &
                                    'D-off-S', &
                                    'H-off-L', &
                                    'H-off-S'
      do irep = 0, nr_ireps - 1
         write(file_unit, '(4x, 7i8)') irep,            &
                                       nbbas(irep, 1),  &
                                       nbbas(irep, 2),  &
                                       ibbas(irep, 1),  &
                                       ibbas(irep, 2),  &
                                       icos(irep+1, 1), &
                                       icos(irep+1, 2)
      end do

      pq_to_uq = ipqtoq
      uq_to_pq = iqtopq
      write(file_unit, '(4x, a)') 'quaternion_packing'
      do i = 1, 4
         write(file_unit, '(8x, 8i2)') pq_to_uq(i, 0:7)
      end do
      do i = 1, 4
         write(file_unit, '(8x, 8i2)') uq_to_pq(i, 0:7)
      end do

      iq(0) = 1
      do i = 1, 3
         iq(i) = jm4pos(4 - i)
         iq(i) = iqtopq(iq(i), isymax(i, 2))
      end do
      write(file_unit, '(4x, a)') 'quaternion_units'
      write(file_unit, '(8x, 4i2)') iq(0:3)

      close(file_unit, status = 'keep')

      call basis_info_clear()

   end subroutine

   subroutine basis_info_read(use_mpi_if_possible)

      use dft_cfg
      use xc_mpi

!     --------------------------------------------------------------------------
      logical, intent(in)  :: use_mpi_if_possible
!     --------------------------------------------------------------------------
      character(120)       :: line
      integer              :: i, j, k, l, m, n
      integer              :: sum_nr_primitives
      integer              :: ioff_c, ioff_s
      integer              :: idummy, irep, cent
      logical              :: symmetry_info_found
      real(8)              :: prim, coef
      character(1)         :: ch
      integer              :: inter_irep(0:3), jxyz, irep1, irep2, ic, nr1, nr2, mc, icount
      integer              :: nr_ao_irep, isym, deg
      integer              :: lvar, mvar, nvar
      integer              :: iao, iso, iso_dirac_sorting
      integer              :: ishell, nshell
      integer              :: nround
      integer              :: size_c2s_mat
      integer              :: lssl_ic(2) = (/2, 1/)
      real(8)              :: f
      integer, allocatable :: isymao(:, :)
      integer, allocatable :: irepao(:)
      integer, allocatable :: ipind(:)
      integer, allocatable :: kstart(:)
      real(8), parameter   :: bitstring_parity(0:7) = (/1.0d0, -1.0d0, -1.0d0, 1.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0/)
!     --------------------------------------------------------------------------

      call basis_info_clear()

      if (xc_mpi_is_master()) then

      symmetry_info_found = .false.

      inquire(file = file_name, exist = file_exists)
      if (.not. file_exists) then
         print *, 'basis_info file not found'
         stop
      end if

      open(file_unit,            &
           file   = file_name,   &
           status = 'unknown',   &
           form   = 'formatted', &
           access = 'sequential')
      rewind(file_unit)

      do while (.true.)

         read(file_unit, '(a120)', end=1) line

         if (word_contains(line, '*** geometry')) then
            read(file_unit, *)
            read(file_unit, *) nr_sym_independent_centers
            read(file_unit, *)
            allocate(sym_independent_center_charge(nr_sym_independent_centers))
            allocate(sym_independent_center_xyz(nr_sym_independent_centers, 3))
            do i = 1, nr_sym_independent_centers
               read(file_unit, *) sym_independent_center_charge(i), &
                                  sym_independent_center_xyz(i, 1:3)
            end do
            read(file_unit, *)
            allocate(sym_independent_center_deg(nr_sym_independent_centers))
            do i = 1, nr_sym_independent_centers
               read(file_unit, *) idummy, sym_independent_center_deg(i)
            end do
            read(file_unit, *)
            read(file_unit, *) nr_sym_dependent_centers
            read(file_unit, *)
            allocate(sym_dependent_center_charge(nr_sym_dependent_centers))
            allocate(sym_dependent_center_xyz(nr_sym_dependent_centers, 3))
            do i = 1, nr_sym_dependent_centers
               read(file_unit, *) sym_dependent_center_charge(i), &
                                  sym_dependent_center_xyz(i, 1:3)
            end do
         end if

         if (word_contains(line, '*** symmetry')) then
            symmetry_info_found = .true.
            read(file_unit, *)
            read(file_unit, *) nr_ireps
            read(file_unit, *)
            read(file_unit, *) irep_of_axes(1:3)
            read(file_unit, *)
            read(file_unit, *) irep_of_rotations(1:3)
            read(file_unit, *)
            do irep = 0, nr_ireps - 1
               read(file_unit, *) idummy,                   &
                                  ao_nr(irep, 1),           &
                                  ao_nr(irep, 2),           &
                                  ao_off(irep, 1),          &
                                  ao_off(irep, 2),          &
                                  ao_off_hermit(irep+1, 1), &
                                  ao_off_hermit(irep+1, 2)
            end do
            read(file_unit, *)
            do i = 1, 4
               read(file_unit, *) pq_to_uq(i, 0:7)
            end do
            do i = 1, 4
               read(file_unit, *) uq_to_pq(i, 0:7)
            end do
            read(file_unit, *)
            read(file_unit, *) iq(0:3)
         end if

         if (word_contains(line, '*** basis')) then
            read(file_unit, *)
            read(file_unit, *) basis_is_spherical
            read(file_unit, *)
            read(file_unit, *) use_only_large
            read(file_unit, *)
            read(file_unit, *) sum_nr_primitives
            read(file_unit, *)

!           overkill allocation for convenience
!           it will "never become a bottleneck" (2011-04-29)
            allocate(shell_center(sum_nr_primitives))
            allocate(l_quantum_nr(sum_nr_primitives))
            allocate(nr_nonzero_contraction_coefs(sum_nr_primitives))
            allocate(shell_stabilizing_irep(sum_nr_primitives))
            allocate(cartesian_deg(sum_nr_primitives))
            allocate(cartesian_off(sum_nr_primitives))
            allocate(spherical_deg(sum_nr_primitives))
            allocate(spherical_off(sum_nr_primitives))
            allocate(primitive_exp(sum_nr_primitives, sum_nr_primitives))
            allocate(contraction_coef(sum_nr_primitives, sum_nr_primitives))

            j = 0
            k = 0
            ioff_c = 0
            ioff_s = 0
            nr_shells_large = 0
            nr_shells_small = 0
            do i = 1, sum_nr_primitives
               read(file_unit, *) cent,   &
                                  ch,     &
                                  ishell, &
                                  l,      &
                                  irep,   &
                                  prim,   &
                                  coef
               if (ishell > k) then
                  j = 1
                  shell_center(ishell)           = cent
                  l_quantum_nr(ishell)           = l
                  shell_stabilizing_irep(ishell) = irep
                  cartesian_deg(ishell) = (l+1)*(l+2)/2
                  cartesian_off(ishell) = ioff_c
                  ioff_c = ioff_c + cartesian_deg(ishell)
                  spherical_deg(ishell) = 2*l + 1
                  spherical_off(ishell) = ioff_s
                  ioff_s = ioff_s + spherical_deg(ishell)
                  if (ch == 'L') then
                     nr_shells_large = nr_shells_large + 1
                  end if
                  if (ch == 'S') then
                     nr_shells_small = nr_shells_small + 1
                  end if
               else
                  j = j + 1
               end if
               primitive_exp(ishell, j)    = prim
               contraction_coef(ishell, j) = coef
               nr_nonzero_contraction_coefs(ishell) = j
               k = ishell
            end do
         end if

      end do
1     continue

      close(file_unit, status = 'keep')

      if (use_only_large) then
         mc = 1
      else
         mc = 2
      end if

!     transform as rotation
      llss_nr_density_blocks = 0

      inter_irep(0) = 0
      do jxyz = 1, 3
         inter_irep(jxyz) = irep_of_rotations(jxyz)
      end do

      do irep = 0, nr_ireps - 1
         do jxyz = 0, 3
            i = 0
            do ic = 1, mc
               do irep1 = 0, nr_ireps - 1
                  irep2 = ieor(irep1, ieor(inter_irep(jxyz), irep))

                  nr1 = ao_nr(irep1, ic)
                  nr2 = ao_nr(irep2, ic)

                  if ((nr1 > 0) .and. (nr2 > 0)) then
                     i = i + 1
                     llss_nr_ao(1, i, jxyz, irep) = nr1
                     llss_st_ao(1, i, jxyz, irep) = ao_off(irep1, ic) + 1
                     llss_nr_ao(2, i, jxyz, irep) = nr2
                     llss_st_ao(2, i, jxyz, irep) = ao_off(irep2, ic) + 1
                     llss_prefactor(i, jxyz, irep) = 1.0d0
!                                        fixme ugly dependency
                     if ((ic == 2) .and. dft_cfg_betasigma) then
                        llss_prefactor(i, jxyz, irep) = -1.0d0
                     end if
                  end if
               end do
            end do
            llss_nr_density_blocks(jxyz, irep) = i
         end do
      end do

!     transform as translation
      lssl_nr_density_blocks = 0

      inter_irep(0) = 0
      do jxyz = 1, 3
         inter_irep(jxyz) = irep_of_axes(jxyz)
      end do

      do irep = 0, nr_ireps - 1
         do jxyz = 0, 3
            i = 0
            do ic = 1, mc
               do irep1 = 0, nr_ireps - 1
                  irep2 = ieor(irep1, ieor(inter_irep(jxyz), irep))

                  nr1 = ao_nr(irep1, ic)
                  nr2 = ao_nr(irep2, lssl_ic(ic))

                  if ((nr1 > 0) .and. (nr2 > 0)) then
                     i = i + 1
                     lssl_nr_ao(1, i, jxyz, irep) = nr1
                     lssl_st_ao(1, i, jxyz, irep) = ao_off(irep1, ic) + 1
                     lssl_nr_ao(2, i, jxyz, irep) = nr2
                     lssl_st_ao(2, i, jxyz, irep) = ao_off(irep2, lssl_ic(ic)) + 1
                  end if
               end do
            end do
            lssl_nr_density_blocks(jxyz, irep) = i
         end do
      end do

      if (basis_is_spherical) then
         allocate(isymao(l_max, 2*l_max + 1))
         isymao = 0
         do l = 0, l_max
            do i = 1, 2*l + 1
               isymao(l + 1, i) = ireplm(l, mdef(l, i))
            end do
         end do
      else
         allocate(isymao(l_max, (l_max+1)*(l_max+2)/2))
         isymao = 0
         do l = 1, l_max
            k = 0
            do m = 1, l
               do n = 1, m
                  k = k + 1
                  lvar = mod(l-m, 2)*irep_of_axes(1)
                  mvar = mod(m-n, 2)*irep_of_axes(2)
                  nvar = mod(n-1, 2)*irep_of_axes(3)
                  isymao(l, k) = ieor(lvar, ieor(mvar, nvar))
               end do
            end do
         end do
      end if

      nr_shells = nr_shells_large + nr_shells_small

      nr_ao_cartesian = 0
      nr_ao_spherical = 0
      do ishell = 1, nr_shells
         do irep = 0, nr_ireps - 1
            if (iand(irep, shell_stabilizing_irep(ishell)) == 0) then
               nr_ao_cartesian = nr_ao_cartesian + cartesian_deg(ishell)
               nr_ao_spherical = nr_ao_spherical + spherical_deg(ishell)
            end if
         end do
      end do

      if (basis_is_spherical) then
         nr_ao = nr_ao_spherical
      else
         nr_ao = nr_ao_cartesian
      end if

      allocate(irepao(nr_ao))
      allocate(ipind(nr_ao))
      irepao = 0
      ipind = 0
      iso = 0
      do irep = 0, nr_ireps - 1
         iao = 0
         do ishell = 1, nr_shells
            if (basis_is_spherical) then
               deg = spherical_deg(ishell)
            else
               deg = cartesian_deg(ishell)
            end if
            do k = 1, deg
               iao = iao + 1
               if (iand(shell_stabilizing_irep(ishell), ieor(irep, isymao(l_quantum_nr(ishell) + 1, k))) == 0) then
                  iso = iso + 1
                  ipind(iso) = ishell*2**16 + k*2**8 + irepao(iao)
                  irepao(iao) = irepao(iao) + 1
               end if
            end do
         end do
      end do
      deallocate(irepao)

      allocate(kstart(nr_shells))
      kstart = 0
      do ishell = 2, nr_shells
         if (basis_is_spherical) then
            deg = spherical_deg(ishell-1)
         else
            deg = cartesian_deg(ishell-1)
         end if
         kstart(ishell) = kstart(ishell - 1) + sym_independent_center_deg(shell_center(ishell-1))*deg
      end do

      nr_so_ao = 0
      do nround = 1, 2
!        two rounds: first round gets dimensions
!                    second round allocates and sets arrays

         if (nround == 2) then
            allocate(i_so_ao(2, nr_so_ao))
            allocate(f_so_ao(nr_so_ao))
            nr_so_ao = 0
         end if

         i = 1
         do irep = 0, nr_ireps - 1
            icount = 0
            nr_ao_irep = ao_nr(irep, 1) + ao_nr(irep, 2)
            do iso = i, i + nr_ao_irep - 1

               ishell = iand(ishft(ipind(iso), -16), 65535)
               nshell = iand(ishft(ipind(iso),  -8),   255)

               iao = kstart(ishell) + nshell

               icount = icount + 1
               if (icount > ao_nr(irep, 1)) then
!                 is small
                  ic = 2
               else
!                 is large
                  ic = 1
               end if

!              do not plan to AO-SO transform small component functions
!              in case of zero small-small metric (.LEVY-LEBLOND and .X2C4)
               if (.not. (use_only_large .and. (ic == 2))) then

!                 shift by difference dirac sorting minus hermit sorting
                  iso_dirac_sorting = iso + ao_off(irep, ic) - ao_off_hermit(irep + 1, ic)

                  do isym = 0, nr_ireps - 1
                     if (iand(isym, shell_stabilizing_irep(ishell)) == 0) then

                        f = bitstring_parity(iand(isym, ieor(irep, isymao(l_quantum_nr(ishell) + 1, nshell))))

                        nr_so_ao = nr_so_ao + 1
                        if (nround == 2) then
                           i_so_ao(1, nr_so_ao) = iso_dirac_sorting
                           i_so_ao(2, nr_so_ao) = iao
                           f_so_ao(nr_so_ao)    = f
                        end if

                        if (basis_is_spherical) then
                           iao = iao + spherical_deg(ishell)
                        else
                           iao = iao + cartesian_deg(ishell)
                        end if

                     endif
                  end do
               end if
            end do
            i = i + nr_ao_irep
         end do
      end do

      deallocate(isymao)
      deallocate(ipind)
      deallocate(kstart)

      if (.not. symmetry_info_found) then
!        fixme
      end if

!     limit active nr of shells
      if (use_only_large) then
         nr_shells = nr_shells_large
      end if

      call setup_c2s_mat()
      size_c2s_mat = size(c2s_mat)

      end if

#ifdef VAR_MPI
      if (use_mpi_if_possible .and. (xc_mpi_get_nr_proc() > 1)) then

         call xc_mpi_bcast(nr_so_ao                  )
         call xc_mpi_bcast(nr_sym_dependent_centers  )
         call xc_mpi_bcast(nr_sym_independent_centers)
         call xc_mpi_bcast(sum_nr_primitives         )
         call xc_mpi_bcast(size_c2s_mat              )

         if (.not. xc_mpi_is_master()) then
            allocate(sym_independent_center_charge(nr_sym_independent_centers))
            allocate(sym_independent_center_xyz(nr_sym_independent_centers, 3))
            allocate(sym_independent_center_deg(nr_sym_independent_centers))

            allocate(sym_dependent_center_charge(nr_sym_dependent_centers))
            allocate(sym_dependent_center_xyz(nr_sym_dependent_centers, 3))

            allocate(shell_center(sum_nr_primitives))
            allocate(l_quantum_nr(sum_nr_primitives))
            allocate(nr_nonzero_contraction_coefs(sum_nr_primitives))
            allocate(shell_stabilizing_irep(sum_nr_primitives))
            allocate(cartesian_deg(sum_nr_primitives))
            allocate(cartesian_off(sum_nr_primitives))
            allocate(spherical_deg(sum_nr_primitives))
            allocate(spherical_off(sum_nr_primitives))
            allocate(primitive_exp(sum_nr_primitives, sum_nr_primitives))
            allocate(contraction_coef(sum_nr_primitives, sum_nr_primitives))
            allocate(c2s_mat(size_c2s_mat))

            allocate(i_so_ao(2, nr_so_ao))
            allocate(f_so_ao(nr_so_ao))
         end if

         call xc_mpi_bcast(ao_nr                     )
         call xc_mpi_bcast(ao_off                    )
         call xc_mpi_bcast(ao_off_hermit             )
         call xc_mpi_bcast(c2s_index                 )
         call xc_mpi_bcast(iq                        )
         call xc_mpi_bcast(irep_of_axes              )
         call xc_mpi_bcast(irep_of_rotations         )
         call xc_mpi_bcast(llss_nr_ao                )
         call xc_mpi_bcast(llss_nr_density_blocks    )
         call xc_mpi_bcast(llss_st_ao                )
         call xc_mpi_bcast(lssl_nr_ao                )
         call xc_mpi_bcast(lssl_nr_density_blocks    )
         call xc_mpi_bcast(lssl_st_ao                )
         call xc_mpi_bcast(nr_ao                     )
         call xc_mpi_bcast(nr_ao_cartesian           )
         call xc_mpi_bcast(nr_ao_spherical           )
         call xc_mpi_bcast(nr_ireps                  )
         call xc_mpi_bcast(nr_shells                 )
         call xc_mpi_bcast(nr_shells_large           )
         call xc_mpi_bcast(nr_shells_small           )
         call xc_mpi_bcast(pq_to_uq                  )
         call xc_mpi_bcast(uq_to_pq                  )
         call xc_mpi_bcast(basis_is_spherical        )
         call xc_mpi_bcast(use_only_large            )
         call xc_mpi_bcast(llss_prefactor            )

         call xc_mpi_bcast(cartesian_deg                )
         call xc_mpi_bcast(cartesian_off                )
         call xc_mpi_bcast(i_so_ao                      )
         call xc_mpi_bcast(l_quantum_nr                 )
         call xc_mpi_bcast(nr_nonzero_contraction_coefs )
         call xc_mpi_bcast(shell_center                 )
         call xc_mpi_bcast(shell_stabilizing_irep       )
         call xc_mpi_bcast(spherical_deg                )
         call xc_mpi_bcast(spherical_off                )
         call xc_mpi_bcast(sym_independent_center_deg   )
         call xc_mpi_bcast(c2s_mat                      )
         call xc_mpi_bcast(contraction_coef             )
         call xc_mpi_bcast(f_so_ao                      )
         call xc_mpi_bcast(primitive_exp                )
         call xc_mpi_bcast(sym_dependent_center_charge  )
         call xc_mpi_bcast(sym_dependent_center_xyz     )
         call xc_mpi_bcast(sym_independent_center_charge)
         call xc_mpi_bcast(sym_independent_center_xyz   )
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

   function fac(n)

!     --------------------------------------------------------------------------
      integer, intent(in) :: n
      real(8)             :: fac
!     it is real(8) so that expressions don't overflow!
!     --------------------------------------------------------------------------
      integer             :: i
!     --------------------------------------------------------------------------

      fac = 1.0d0
      do i = 1, n
         fac = fac*dble(i)
      end do

   end function

   function fac2(n)

!     --------------------------------------------------------------------------
      integer, intent(in) :: n
      real(8)             :: fac2
!     it is real(8) so that expressions don't overflow!
!     --------------------------------------------------------------------------
      integer             :: i
!     --------------------------------------------------------------------------

      if (n < 0) then
         fac2 = dble(n + 2)
         do i = n + 4, 1, 2
            fac2 = fac2*dble(i)
         end do
         fac2 = 1.0d0/fac2
      else if (n == 0) then
         fac2 = 1.0d0
      else
         fac2 = dble(n)
         do i = n - 2, 1, -2
            fac2 = fac2*dble(i)
         end do
      end if

   end function

   function binom(i, j)

!     --------------------------------------------------------------------------
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer             :: binom
      real(8)             :: f
!     --------------------------------------------------------------------------

      f = fac(i)/(fac(i - j)*fac(j))
!     f and fac are real(8) so that binom doesn't overflow!
      binom = int(f)

   end function

   subroutine setup_c2s_mat()

!     --------------------------------------------------------------------------
      integer              :: nc, ns
      integer              :: i, j, k, l, m, n
      integer              :: ix, ilm
      real(8)              :: cm, cmk, cmki, cmkij
      real(8), allocatable :: tmat(:, :)
      real(8), allocatable :: cossin_coef(:, :)
      real(8), allocatable :: legendre_coef(:)
!     --------------------------------------------------------------------------


      c2s_index = 0
      n = 0
      do l = 2, l_max
         c2s_index(l) = n + 1
         n = n + ((l+1)*(l+2)/2)**2
      end do
      allocate(c2s_mat(n))
      c2s_mat = 0.0d0

      do l = 2, l_max

         nc = (l+1)*(l+2)/2
         ns = 2*l + 1

         allocate(legendre_coef(0:l))
         legendre_coef = 0.0d0
         do k = 0, l/2
            legendre_coef(l - 2*k) = (dble((-1)**k)/dble(2**l))*binom(l, k)*binom(2*(l - k), l)
         end do

         allocate(cossin_coef(0:l, 0:l))
         cossin_coef = 0.0d0
         do m = 0, l
            cossin_coef(m, 0) = 1.0d0
            do k = 1, m
               cossin_coef(m, k) = cossin_coef(m - 1, k - 1)*dble((-1)**(k - 1))
               if (m > k) cossin_coef(m, k) = cossin_coef(m, k) + cossin_coef(m - 1, k)
            end do
         end do

         allocate(tmat(nc, ns))
         tmat = 0.0d0
         do m = 0, l
            cm = sqrt(2.0d0*fac(l - m)/fac(l + m))
            if (m  ==  0) cm = 1.0d0
            cm = cm/dsqrt(dble(fac2(2*l - 1)))
            do k = mod(l - m, 2), l - m, 2
               if (m > 0) legendre_coef(k) = dble((k + 1))*legendre_coef(k + 1)
               cmk = cm*legendre_coef(k)
               do i = 0, (l - k - m)/2
                  cmki = cmk*binom((l - k - m)/2,i)
                  do j = 0, i
                     cmkij = cmki*binom(i, j)
                     do n = 0, m
                        ix = l - 2*j - m + n
                        ix = ix*(ix + 1)/2 + l + 1 - m - 2*i
                        if (mod(n, 2) == 1) then
                           ilm = 1 + l - m
                        else
                           ilm = 1 + l + m
                        end if
                        tmat(ix, ilm) = tmat(ix, ilm) + cmkij*cossin_coef(m, n)
                     end do
                  end do
               end do
            end do
         end do

         k = c2s_index(l)
         do i = 1, nc
            do j = 1, ns
               c2s_mat(k) = tmat(i, j)
               k = k + 1
            end do
         end do

         deallocate(tmat)
         deallocate(cossin_coef)
         deallocate(legendre_coef)

      end do

   end subroutine

   function mdef(l, i)

!     --------------------------------------------------------------------------
      integer, intent(in) :: l
      integer, intent(in) :: i
      integer             :: mdef
!     --------------------------------------------------------------------------

      if (l == 1) then
          select case (i)
             case (1)
                mdef = 1
             case (2)
                mdef = -1
             case (3)
                mdef = 0
          end select
      else
         mdef = i - l - 1
      endif

   end function

   function ireplm(l, m)

!     --------------------------------------------------------------------------
      integer, intent(in) :: l
      integer, intent(in) :: m
      integer             :: ireplm
!     --------------------------------------------------------------------------

      ireplm = 0
      if (mod(l + m,  2) == 1) then
         ireplm = irep_of_axes(3)
      end if
      if (mod(abs(m), 2) == 1) then
         ireplm = ieor(ireplm, irep_of_axes(1))
      end if
      if (m < 0) then
         ireplm = ieor(ireplm, ieor(irep_of_axes(1), irep_of_axes(2)))
      end if

   end function

end module
