!      Copyright (c) 2011 by the authors of DIRAC.
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

module xc_geometric

   use interface_ao
   use ao_eval

   implicit none

   public xc_geo_0_lda_c1
   public xc_geo_0_gga_c1

   private

   integer :: pos(9) = (/1, 2, 3, 2, 4, 5, 3, 5, 6/)

contains

   subroutine xc_geo_0_lda_c1(g,       &
                              d_n,     &
                              mat_dim, &
                              dmat_0,  &
                              ao)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: g(*)
      real(8), intent(in)    :: d_n
      integer, intent(in)    :: mat_dim
      real(8), intent(in)    :: dmat_0(mat_dim, mat_dim)
      real(8), intent(in)    :: ao(*)
!     --------------------------------------------------------------------------
      integer                :: k, ks, kc
      integer                :: l, ls, lc
      integer                :: ixyz, ic
      real(8)                :: u
!     --------------------------------------------------------------------------

      do k = 1, nr_ao
         u = 0.0d0
         do l = 1, nr_ao
            u = u + 2.0d0*ao(l)*dmat_0(l, k)
         end do
         do ixyz = 1, 3
            ic = 3*(ao_center(k) - 1) + ixyz
            g(ic) = g(ic) - d_n*u*ao(ao_off_g1_m0(ixyz, 0) + k)
         end do
      end do

   end subroutine

   subroutine xc_geo_0_gga_c1(g,       &
                              d_n,     &
                              d_z,     &
                              gn,      &
                              mat_dim, &
                              dmat_0,  &
                              ao)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: g(*)
      real(8), intent(in)    :: d_n
      real(8), intent(in)    :: d_z
      real(8), intent(in)    :: gn(3)
      integer, intent(in)    :: mat_dim
      real(8), intent(in)    :: dmat_0(mat_dim, mat_dim)
      real(8), intent(in)    :: ao(*)
!     --------------------------------------------------------------------------
      integer                :: k, ks, kc
      integer                :: l, ls, lc
      integer                :: ixyz, jxyz, ic
      real(8)                :: t1, t2, t3, t4(3), t5
!     --------------------------------------------------------------------------

      do k = 1, nr_ao
         t1 = 0.0d0
         do l = 1, nr_ao
            t1 = t1 + ao(l)*dmat_0(l, k)
         end do
         t1 = 2.0d0*t1

         t2 = 0.0d0
         do jxyz = 1, 3
            t3 = 0.0d0
            do l = 1, nr_ao
               t3 = t3 + ao(ao_off_g1_m0(jxyz, 0) + l)*dmat_0(l, k)
            end do
            t2 = t2 + t3*gn(jxyz)
         end do
         t2 = 4.0d0*t2

         t4 = 0.0d0
         do ixyz = 1, 3
            do jxyz = 1, 3
               t5 = 0
               do l = 1, nr_ao
                  t5 = t5 + ao(ao_off_g2_m0(pos((jxyz-1)*3 + ixyz), 0) + k)*ao(l)*dmat_0(l, k)
               end do
               t4(ixyz) = t4(ixyz) + t5*gn(jxyz)
            end do
         end do
         t4 = 4.0d0*t4

         do ixyz = 1, 3
            ic = 3*(ao_center(k) - 1) + ixyz
            g(ic) = g(ic) - d_n*t1*ao(ao_off_g1_m0(ixyz, 0) + k) - d_z*t2*ao(ao_off_g1_m0(ixyz, 0) + k) - d_z*t4(ixyz)
         end do
      end do

   end subroutine


end module
