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

module overlap_diagnostic

!  MO overlap diagnostic according to Peach, Benfield, Helgaker, and Tozer:
!  M. J. G. Peach, P. Benfield, T. U. Helgaker, and D. J. Tozer
!  J. Chem. Phys. 128, 044118 (2008)

!  implemented 2010/06/30 by Radovan Bast <radovan.bast@uit.no>

!  known problems: so far no mpi, some point groups make problems with lvcorr

   use memory_allocator

   implicit none

   public sum_lambda
   public integrate_lambda_matrix
   public overlap_diagnostic_init
   public overlap_diagnostic_delete
   public lambda_matrix_to_disc
   public lambda_matrix_from_disc

   private

#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"

   integer                 :: nr_g_mo_ns
   integer                 :: nr_g_mo_pi
   integer                 :: nr_g_mo_ps
   integer                 :: nr_u_mo_ns
   integer                 :: nr_u_mo_pi
   integer                 :: nr_u_mo_ps
   integer                 :: nr_g_mo
   integer                 :: nr_u_mo
   integer                 :: nr_g_ao
   integer                 :: nr_u_ao

   integer                 :: nr_g_ao_blocks
   integer                 :: nr_u_ao_blocks
   integer                 :: g_ao_start(8)
   integer                 :: u_ao_start(8)
   integer                 :: g_ao_end(8)
   integer                 :: u_ao_end(8)

   real(8), allocatable    :: mo_coef(:)
   real(8), allocatable    :: lambda_matrix(:, :)

   character(*), parameter :: file_name = 'lambda_matrix'
   integer,      parameter :: file_unit = 137
   logical                 :: file_exists

contains

   subroutine sum_lambda(lambda_numerator,   &
                         lambda_denominator, &
                         i,                  &
                         s,                  &
                         kappa)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: lambda_numerator
      real(8), intent(inout) :: lambda_denominator
      integer, intent(in)    :: i
      integer, intent(in)    :: s
      real(8), intent(in)    :: kappa
!     --------------------------------------------------------------------------

      lambda_numerator   = lambda_numerator   + kappa*lambda_matrix(npsh(1) + i, npsh(1) + s)
      lambda_denominator = lambda_denominator + kappa

   end subroutine

   subroutine overlap_diagnostic_init()

!     --------------------------------------------------------------------------
      integer :: k, irep, ic, nc
!     --------------------------------------------------------------------------

      nr_g_mo_ns = npsh(1)
      nr_g_mo_pi = nish(1)
      nr_g_mo_ps = nesh(1) - nr_g_mo_pi
      nr_g_mo    = nr_g_mo_ns + nr_g_mo_pi + nr_g_mo_ps

      nr_u_mo_ns = npsh(2)
      nr_u_mo_pi = nish(2)
      nr_u_mo_ps = nesh(2) - nr_u_mo_pi
      nr_u_mo    = nr_u_mo_ns + nr_u_mo_pi + nr_u_mo_ps

      nr_g_ao    = nfbas(1, 0)
      nr_u_ao    = nfbas(2, 0)

      nr_g_ao_blocks = 0
      nr_u_ao_blocks = 0
      g_ao_start     = 0
      u_ao_start     = 0
      g_ao_end       = 0
      u_ao_end       = 0

      if (twocomp) then
         nc = 1
      else
         nc = 2
      end if

      do ic = 1, nc
         do irep = 0, maxrep
            if (nbbas(irep, ic) > 0) then
               if (jbtof(irep, ic) == 1) then
                  nr_g_ao_blocks = nr_g_ao_blocks + 1
                  g_ao_start(nr_g_ao_blocks) = ibbas(irep, ic) + 1
                  g_ao_end(nr_g_ao_blocks)   = ibbas(irep, ic) + nbbas(irep, ic)
               else
                  nr_u_ao_blocks = nr_u_ao_blocks + 1
                  u_ao_start(nr_u_ao_blocks) = ibbas(irep, ic) + 1
                  u_ao_end(nr_u_ao_blocks)   = ibbas(irep, ic) + nbbas(irep, ic)
               end if
            end if
         end do
      end do

      call alloc(mo_coef, ncmotq)
      call alloc(lambda_matrix, nr_g_mo + nr_u_mo, nr_g_mo + nr_u_mo)
      mo_coef       = 0.0d0
      lambda_matrix = 0.0d0
      call read_mo_coef(mo_coef)

   end subroutine

   subroutine overlap_diagnostic_delete()

      nr_g_mo_ns = 0
      nr_g_mo_pi = 0
      nr_g_mo_ps = 0
      nr_u_mo_ns = 0
      nr_u_mo_pi = 0
      nr_u_mo_ps = 0
      nr_g_mo    = 0
      nr_u_mo    = 0
      nr_g_ao    = 0
      nr_u_ao    = 0

      nr_g_ao_blocks = 0
      nr_u_ao_blocks = 0
      g_ao_start     = 0
      u_ao_start     = 0
      g_ao_end       = 0
      u_ao_end       = 0

      if (allocated(mo_coef))       call dealloc(mo_coef)
      if (allocated(lambda_matrix)) call dealloc(lambda_matrix)

   end subroutine

   subroutine lambda_matrix_to_disc()

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

      write(file_unit) nr_g_mo_ns
      write(file_unit) nr_g_mo_pi
      write(file_unit) nr_g_mo_ps
      write(file_unit) nr_u_mo_ns
      write(file_unit) nr_u_mo_pi
      write(file_unit) nr_u_mo_ps
      write(file_unit) nr_g_mo
      write(file_unit) nr_u_mo
      write(file_unit) nr_g_ao
      write(file_unit) nr_u_ao
      write(file_unit) lambda_matrix

      close(file_unit, status = 'keep')

      if (allocated(mo_coef))       call dealloc(mo_coef)
      if (allocated(lambda_matrix)) call dealloc(lambda_matrix)

   end subroutine

   subroutine lambda_matrix_from_disc()

!     --------------------------------------------------------------------------
      integer :: p
      real(8) :: trace
!     --------------------------------------------------------------------------

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

      read(file_unit) nr_g_mo_ns
      read(file_unit) nr_g_mo_pi
      read(file_unit) nr_g_mo_ps
      read(file_unit) nr_u_mo_ns
      read(file_unit) nr_u_mo_pi
      read(file_unit) nr_u_mo_ps
      read(file_unit) nr_g_mo
      read(file_unit) nr_u_mo
      read(file_unit) nr_g_ao
      read(file_unit) nr_u_ao

      call alloc(lambda_matrix, nr_g_mo + nr_u_mo, nr_g_mo + nr_u_mo)
      read(file_unit) lambda_matrix

      close(file_unit, status = 'keep')

      write(*, *) 'debug lambda matrix'
      call prqmat(lambda_matrix,  &
                  norbt,          &
                  norbt,          &
                  norbt,          &
                  norbt,          &
                  1,              &
                  (/1, 2, 3, 4/), &
                  6)

     trace = 0.0d0
     do p = 1, norbt
        trace = trace + lambda_matrix(p, p)
     end do
     write(*, *) 'lambda trace (should be nr of orbitals):', trace
     write(*, *) 'lambda trace divided by nr of orbitals:', trace/norbt

   end subroutine

   subroutine integrate_lambda_matrix(w,  &
                                      ao, &
                                      buffer)

!     --------------------------------------------------------------------------
      real(8), intent(in)    :: w
      real(8), intent(in)    :: ao(*)
      real(8), intent(inout) :: buffer(:, :)
!     --------------------------------------------------------------------------
      integer                :: s, t
      integer                :: i, k, n
      integer                :: iz
      real(8)                :: q1, q2, q3, q4, q_norm
!     --------------------------------------------------------------------------
!
!     structure of mo_coef (assuming no active orbitals):
!
!      _____  1 ... nz
!     |  ___
!     | |  _  1 ... nr gerade mo negative energy secondary
!     | | |_  1 ... nr gerade ao (order according to ibbas)
!     | |___
!     |  ___
!     | |  _  1 ... nr gerade mo positive energy inactive
!     | | |_  1 ... nr gerade ao (order according to ibbas)
!     | |___
!     |  ___
!     | |  _  1 ... nr gerade mo positive energy secondary
!     | | |_  1 ... nr gerade ao (order according to ibbas)
!     | |___
!     |_____
!
!      _____  1 ... nz
!     |  ___
!     | |  _  1 ... nr ungerade mo negative energy secondary
!     | | |_  1 ... nr ungerade ao (order according to ibbas)
!     | |___
!     |  ___
!     | |  _  1 ... nr ungerade mo positive energy inactive
!     | | |_  1 ... nr ungerade ao (order according to ibbas)
!     | |___
!     |  ___
!     | |  _  1 ... nr ungerade mo positive energy secondary
!     | | |_  1 ... nr ungerade ao (order according to ibbas)
!     | |___
!     |_____

      i = 0
      do iz = 1, nz
!        gerade
         do s = 1, nr_g_mo
            buffer(s, iz) = 0.0d0
            do n = 1, nr_g_ao_blocks
               do k = g_ao_start(n), g_ao_end(n)
                  i = i + 1
                  buffer(s, iz) = buffer(s, iz) + ao(k)*mo_coef(i)
               end do
            end do
         end do
      end do
      do iz = 1, nz
         if (nfsym > 1) then
!           ungerade
            do s = nr_g_mo + 1, nr_g_mo + nr_u_mo
               buffer(s, iz) = 0.0d0
               do n = 1, nr_u_ao_blocks
                  do k = u_ao_start(n), u_ao_end(n)
                     i = i + 1
                     buffer(s, iz) = buffer(s, iz) + ao(k)*mo_coef(i)
                  end do
               end do
            end do
         end if
      end do

      do s = 1, nr_g_mo + nr_u_mo
         do t = s, nr_g_mo + nr_u_mo
            select case (nz)
               case (1)
                  q1 = buffer(s, 1)*buffer(t, 1)
                  q_norm = dsqrt(q1*q1)
               case (2)
                  q1 = buffer(s, 1)*buffer(t, 1) &
                     - buffer(s, 2)*buffer(t, 2)
                  q2 = buffer(s, 1)*buffer(t, 2) &
                     + buffer(s, 2)*buffer(t, 1)
                  q_norm = dsqrt(q1*q1 + q2*q2)
               case (4)
                  q1 = buffer(s, 1)*buffer(t, 1) &
                     - buffer(s, 2)*buffer(t, 2) &
                     - buffer(s, 3)*buffer(t, 3) &
                     - buffer(s, 4)*buffer(t, 4)
                  q2 = buffer(s, 1)*buffer(t, 2) &
                     + buffer(s, 2)*buffer(t, 1) &
                     + buffer(s, 3)*buffer(t, 4) &
                     - buffer(s, 4)*buffer(t, 3)
                  q3 = buffer(s, 1)*buffer(t, 3) &
                     + buffer(s, 3)*buffer(t, 1) &
                     - buffer(s, 2)*buffer(t, 4) &
                     + buffer(s, 4)*buffer(t, 2)
                  q4 = buffer(s, 1)*buffer(t, 4) &
                     + buffer(s, 4)*buffer(t, 1) &
                     + buffer(s, 2)*buffer(t, 3) &
                     - buffer(s, 3)*buffer(t, 2)
                  q_norm = dsqrt(q1*q1 + q2*q2 + q3*q3 + q4*q4)
            end select
            lambda_matrix(s, t) = lambda_matrix(s, t) + w*q_norm
         end do
      end do

   end subroutine

   subroutine read_mo_coef(mo_coef)

!     --------------------------------------------------------------------------
      real(8), intent(out) :: mo_coef(*)
!     --------------------------------------------------------------------------
      integer, parameter   :: file_unit = 66
      real(8)              :: energy
!     --------------------------------------------------------------------------

      call opnfil(file_unit, 'DFCOEF', 'OLD', 'read_mo_coef')

      call reacmo_no_work(file_unit, &
                          'DFCOEF',  &
                          mo_coef,   &
                          (/0.0d0/), &
                          (/0/),     &
                          energy,    &
                          2)

      close(file_unit, status = 'keep')

   end subroutine

end module
