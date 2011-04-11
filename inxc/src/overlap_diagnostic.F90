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

!  implemented 2011/01/11 by Radovan Bast <radovan.bast@uit.no>

   use memory_allocator

   implicit none

   public sum_lambda
   public integrate_o_matrix
   public overlap_diagnostic_init
   public overlap_diagnostic_collect
   public overlap_diagnostic_delete
   public o_matrix_to_disc
   public o_matrix_from_disc

   private

   integer                 :: nr_qu
   integer                 :: nr_fermi

   integer                 :: nr_mo
   integer                 :: nr_g_mo
   integer                 :: nr_g_mo_ns
   integer                 :: nr_g_mo_pi
   integer                 :: nr_g_mo_ps
   integer                 :: nr_u_mo
   integer                 :: nr_u_mo_ns
   integer                 :: nr_u_mo_pi
   integer                 :: nr_u_mo_ps

   integer                 :: nr_g_ao
   integer                 :: nr_u_ao
   integer                 :: nr_g_ao_blocks
   integer                 :: nr_u_ao_blocks
   integer                 :: g_ao_start(8)
   integer                 :: u_ao_start(8)
   integer                 :: g_ao_end(8)
   integer                 :: u_ao_end(8)

   real(8), allocatable    :: mo_coef(:)
   real(8), allocatable    :: o_matrix(:, :)
   real(8), allocatable    :: buffer(:, :)

   character(*), parameter :: file_name = 'o_matrix'
   integer,      parameter :: file_unit = 137
   logical                 :: file_exists

contains

   subroutine sum_lambda(lambda_numerator,   &
                         lambda_denominator, &
                         i,                  &
                         s,                  &
                         ifermi,             &
                         sfermi,             &
                         kappa)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: lambda_numerator
      real(8), intent(inout) :: lambda_denominator
      integer, intent(in)    :: i
      integer, intent(in)    :: s
      integer, intent(in)    :: ifermi
      integer, intent(in)    :: sfermi
      real(8), intent(in)    :: kappa
!     --------------------------------------------------------------------------
      integer                :: ioff, soff
!     --------------------------------------------------------------------------

      select case (ifermi)
         case (1)
            ioff = nr_g_mo_ns
         case (2)
            ioff = nr_g_mo + nr_u_mo_ns
      end select

      select case (sfermi)
         case (1)
            soff = nr_g_mo_ns
         case (2)
            soff = nr_g_mo + nr_u_mo_ns
      end select

      lambda_numerator   = lambda_numerator   + kappa*o_matrix(ioff + i, soff + s)
      lambda_denominator = lambda_denominator + kappa

   end subroutine

   subroutine overlap_diagnostic_init(parallel_xc, i_am_master)

      use xc_mpi

!     --------------------------------------------------------------------------
      logical, intent(in) :: parallel_xc, i_am_master
!     --------------------------------------------------------------------------
      integer             :: k, irep, ic, nc
!     --------------------------------------------------------------------------

#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"

      if (i_am_master) then
         nr_qu    = nz
         nr_fermi = nfsym
         
         nr_g_mo_ns = npsh(1)
         nr_g_mo_pi = nish(1)
         nr_g_mo_ps = nesh(1) - nr_g_mo_pi
         nr_g_mo    = nr_g_mo_ns + nr_g_mo_pi + nr_g_mo_ps
         
         nr_u_mo_ns = npsh(2)
         nr_u_mo_pi = nish(2)
         nr_u_mo_ps = nesh(2) - nr_u_mo_pi
         nr_u_mo    = nr_u_mo_ns + nr_u_mo_pi + nr_u_mo_ps
         
         nr_mo      = nr_g_mo + nr_u_mo
         
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
      end if

#ifdef VAR_MPI
      if (parallel_xc) then
         call xc_mpi_bcast(nr_qu         )
         call xc_mpi_bcast(nr_fermi      )
         call xc_mpi_bcast(nr_mo         )
         call xc_mpi_bcast(nr_g_mo       )
         call xc_mpi_bcast(nr_g_mo_ns    )
         call xc_mpi_bcast(nr_g_mo_pi    )
         call xc_mpi_bcast(nr_g_mo_ps    )
         call xc_mpi_bcast(nr_u_mo       )
         call xc_mpi_bcast(nr_u_mo_ns    )
         call xc_mpi_bcast(nr_u_mo_pi    )
         call xc_mpi_bcast(nr_u_mo_ps    )
         call xc_mpi_bcast(nr_g_ao       )
         call xc_mpi_bcast(nr_u_ao       )
         call xc_mpi_bcast(nr_g_ao_blocks)
         call xc_mpi_bcast(nr_u_ao_blocks)
         call xc_mpi_bcast(g_ao_start    )
         call xc_mpi_bcast(u_ao_start    )
         call xc_mpi_bcast(g_ao_end      )
         call xc_mpi_bcast(u_ao_end      )
      end if
#endif

      call alloc(mo_coef, nr_qu*(nr_g_mo*nr_g_ao + nr_u_mo*nr_u_ao))
      call alloc(o_matrix, nr_mo, nr_mo)
      call alloc(buffer, nr_mo, nr_qu)

      mo_coef  = 0.0d0
      o_matrix = 0.0d0

      if (i_am_master) then
         call read_mo_coef(mo_coef)
      end if

#ifdef VAR_MPI
      if (parallel_xc) then
         call xc_mpi_bcast(mo_coef)
      end if
#endif

   end subroutine

   subroutine overlap_diagnostic_collect(i_am_master)

      use xc_mpi

!     --------------------------------------------------------------------------
      logical, intent(in) :: i_am_master
!     --------------------------------------------------------------------------
      integer             :: ierr
      real(8)             :: dummy
!     --------------------------------------------------------------------------

#ifdef VAR_MPI
#include "mpif.h"

      if (i_am_master) then
         call mpi_reduce(MPI_IN_PLACE,         &
                         o_matrix,             &
                         nr_mo*nr_mo,          &
                         mpi_double_precision, &
                         mpi_sum,              &
                         0,                    &
                         mpi_comm_world,       &
                         ierr)
      else
         call mpi_reduce(o_matrix,             &
                         dummy,                &
                         nr_mo*nr_mo,          &
                         mpi_double_precision, &
                         mpi_sum,              &
                         0,                    &
                         mpi_comm_world,       &
                         ierr)
      end if
#endif

   end subroutine

   subroutine overlap_diagnostic_delete()

      nr_qu          = 0
      nr_fermi       = 0

      nr_mo          = 0
      nr_g_mo        = 0
      nr_g_mo_ns     = 0
      nr_g_mo_pi     = 0
      nr_g_mo_ps     = 0
      nr_u_mo        = 0
      nr_u_mo_ns     = 0
      nr_u_mo_pi     = 0
      nr_u_mo_ps     = 0

      nr_g_ao        = 0
      nr_u_ao        = 0
      nr_g_ao_blocks = 0
      nr_u_ao_blocks = 0
      g_ao_start     = 0
      u_ao_start     = 0
      g_ao_end       = 0
      u_ao_end       = 0

      if (allocated(mo_coef))  call dealloc(mo_coef)
      if (allocated(o_matrix)) call dealloc(o_matrix)
      if (allocated(buffer))   call dealloc(buffer)

   end subroutine

   subroutine o_matrix_to_disc()

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
      write(file_unit) o_matrix

      close(file_unit, status = 'keep')

      if (allocated(mo_coef))       call dealloc(mo_coef)
      if (allocated(o_matrix)) call dealloc(o_matrix)

   end subroutine

   subroutine o_matrix_from_disc()

!     --------------------------------------------------------------------------
      integer :: p
      real(8) :: trace, diff, diff_tolerance
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

      call alloc(o_matrix, nr_g_mo + nr_u_mo, nr_g_mo + nr_u_mo)
      read(file_unit) o_matrix

      close(file_unit, status = 'keep')

      write(*, *) 'O_ia matrix (for PBHT MO Overlap Diagnostic)'
      call prqmat(o_matrix,       &
                  nr_mo,          &
                  nr_mo,          &
                  nr_mo,          &
                  nr_mo,          &
                  1,              &
                  (/1, 2, 3, 4/), &
                  6)

     trace = 0.0d0
     do p = 1, nr_mo
        trace = trace + o_matrix(p, p)
     end do
     write(*, *) 'trace (should be nr of orbitals):', trace
     write(*, *) 'trace/nr of orbitals (should be very close to 1):', trace/nr_mo
     diff = dabs((trace/nr_mo) - 1.0d0)
     diff_tolerance = 1.0d-4
     if (diff > diff_tolerance) then
        call quit('overlap_diagnostic: avg O_ia trace significantly different from 1.0')
     end if

   end subroutine

   subroutine integrate_o_matrix(w, ao)

!     --------------------------------------------------------------------------
      real(8), intent(in)    :: w
      real(8), intent(in)    :: ao(*)
!     --------------------------------------------------------------------------
      integer                :: s, t
      integer                :: i, k, n
      integer                :: iz
      real(8)                :: q1, q2, q3, q4, q_norm
!     --------------------------------------------------------------------------
!
!     structure of mo_coef (assuming no active orbitals):
!
!      _____  1 ... nr_qu
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
!      _____  1 ... nr_qu
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
      do iz = 1, nr_qu
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
      do iz = 1, nr_qu
         if (nr_fermi > 1) then
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
            select case (nr_qu)
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
            o_matrix(s, t) = o_matrix(s, t) + w*q_norm
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

end module overlap_diagnostic
