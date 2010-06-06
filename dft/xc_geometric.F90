module xc_geometric

   use basis_info

   implicit none 

   public xc_geo_0_lda_c1
   public xc_geo_0_gga_c1
   public xc_geo_1_lda_c1
   public xc_geo_2_lda_c1

   private

   integer :: pos(9) = (/1, 2, 3, 2, 4, 5, 3, 5, 6/)

contains

   subroutine xc_geo_0_lda_c1(g,        &
                              d_n,      &
                              mat_dim,  &
                              dmat_0,   &
                              ao,       &
                              ao1)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: g(*)
      real(8), intent(in)    :: d_n
      integer, intent(in)    :: mat_dim
      real(8), intent(in)    :: dmat_0(mat_dim, mat_dim)
      real(8), intent(in)    :: ao(mat_dim)
      real(8), intent(in)    :: ao1(mat_dim, 3)
!     --------------------------------------------------------------------------
      integer                :: k, ks, kc
      integer                :: l, ls, lc
      integer                :: ixyz, ic, nr_center
      real(8)                :: t0
!     --------------------------------------------------------------------------

      k = 0
      do ks = 1, nr_shells
         nr_center = shell_center(ks)
         do kc = 1, cartesian_deg(ks)
            k = k + 1
      
            t0 = 0.0d0
            l = 0
            do ls = 1, nr_shells
               do lc = 1, cartesian_deg(ls)
                  l = l + 1
                  t0 = t0 + ao(l)*dmat_0(l, k)
               end do
            end do
      
            do ixyz = 1, 3
               ic = 3*(nr_center - 1) + ixyz
               g(ic) = g(ic) - 2.0d0*d_n*ao1(k, ixyz)*t0
            end do
         end do
      end do
   
   end subroutine

   subroutine xc_geo_0_gga_c1(g,        &
                              d_n,      &
                              d_z,      &
                              gn,       &
                              mat_dim,  &
                              dmat_0,   &
                              ao,       &
                              ao1,      &
                              ao2)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: g(*)
      real(8), intent(in)    :: d_n
      real(8), intent(in)    :: d_z
      real(8), intent(in)    :: gn(3)
      integer, intent(in)    :: mat_dim
      real(8), intent(in)    :: dmat_0(mat_dim, mat_dim)
      real(8), intent(in)    :: ao(mat_dim)
      real(8), intent(in)    :: ao1(mat_dim, 3)
      real(8), intent(in)    :: ao2(mat_dim, 6)
!     --------------------------------------------------------------------------
      integer                :: k, ks, kc
      integer                :: l, ls, lc
      integer                :: ixyz, jxyz, ic, nr_center
      real(8)                :: t0, t1, f0, f1
!     --------------------------------------------------------------------------

      k = 0
      do ks = 1, nr_shells
         nr_center = shell_center(ks)
         do kc = 1, cartesian_deg(ks)
            k = k + 1
      
            do ixyz = 1, 3

               t0 = 0.0d0
               t1 = 0.0d0
               l = 0
               do ls = 1, nr_shells
                  do lc = 1, cartesian_deg(ls)
                     l = l + 1
                     f0 = ao1(k, ixyz)*ao(l)
                     f1 = 0.0
                     do jxyz = 1, 3
                        f1 = f1 + gn(jxyz)*ao1(k,      ixyz             )*ao1(l, jxyz)
                        f1 = f1 + gn(jxyz)*ao2(k, pos((jxyz-1)*3 + ixyz))*ao (l)
                     end do
                     t0 = t0 + dmat_0(l, k)*f0
                     t1 = t1 + dmat_0(l, k)*f1
                  end do
               end do

               ic = 3*(nr_center - 1) + ixyz
               g(ic) = g(ic) - 2.0d0*t0*d_n - 2.0d0*t1*d_z
            end do
         end do
      end do
   
   end subroutine

   subroutine xc_geo_1_lda_c1(g,        &
                              d_n,      &
                              d_nn,     &
                              n_t1,     &
                              mat_dim,  &
                              dmat_0,   &
                              tmat_1,   &
                              ao,       &
                              ao1)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: g(*)
      real(8), intent(in)    :: d_n
      real(8), intent(in)    :: d_nn
      real(8), intent(in)    :: n_t1
      integer, intent(in)    :: mat_dim
      real(8), intent(in)    :: dmat_0(mat_dim, mat_dim)
      real(8), intent(in)    :: tmat_1(mat_dim, mat_dim)
      real(8), intent(in)    :: ao(mat_dim)
      real(8), intent(in)    :: ao1(mat_dim, 3)
!     --------------------------------------------------------------------------
      integer                :: k, ks, kc
      integer                :: l, ls, lc
      integer                :: ixyz, ic, nr_center
      real(8)                :: t0, t1
!     --------------------------------------------------------------------------

      k = 0
      do ks = 1, nr_shells
         nr_center = shell_center(ks)
         do kc = 1, cartesian_deg(ks)
            k = k + 1
      
            t0 = 0.0d0
            t1 = 0.0d0
            l = 0
            do ls = 1, nr_shells
               do lc = 1, cartesian_deg(ls)
                  l = l + 1
                  t0 = t0 + ao(l)*(dmat_0(l, k) + dmat_0(k, l))
                  t1 = t1 + ao(l)*(tmat_1(l, k) + tmat_1(k, l))
               end do
            end do
      
            do ixyz = 1, 3
               ic = 3*(nr_center - 1) + ixyz
               g(ic) = g(ic) - 2.0d0* d_n*ao1(k, ixyz)*t1 &
                             - 2.0d0*d_nn*ao1(k, ixyz)*t0*n_t1
            end do
         end do
      end do
   
   end subroutine

   subroutine xc_geo_2_lda_c1(g,        &
                              d_nn,     &
                              d_nnn,    &
                              n_t1,     &
                              n_t2,     &
                              mat_dim,  &
                              dmat_0,   &
                              tmat_1,   &
                              tmat_2,   &
                              ao,       &
                              ao1)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: g(*)
      real(8), intent(in)    :: d_nn
      real(8), intent(in)    :: d_nnn
      real(8), intent(in)    :: n_t1
      real(8), intent(in)    :: n_t2
      integer, intent(in)    :: mat_dim
      real(8), intent(in)    :: dmat_0(mat_dim, mat_dim)
      real(8), intent(in)    :: tmat_1(mat_dim, mat_dim)
      real(8), intent(in)    :: tmat_2(mat_dim, mat_dim)
      real(8), intent(in)    :: ao(mat_dim)
      real(8), intent(in)    :: ao1(mat_dim, 3)
!     --------------------------------------------------------------------------
      integer                :: k, ks, kc
      integer                :: l, ls, lc
      integer                :: ixyz, ic, nr_center
      real(8)                :: t0, t1, t2
!     --------------------------------------------------------------------------

      k = 0
      do ks = 1, nr_shells
         nr_center = shell_center(ks)
         do kc = 1, cartesian_deg(ks)
            k = k + 1
      
            t0 = 0.0d0
            t1 = 0.0d0
            t2 = 0.0d0
            l = 0
            do ls = 1, nr_shells
               do lc = 1, cartesian_deg(ls)
                  l = l + 1
                  t0 = t0 + ao(l)*(dmat_0(l, k) + dmat_0(k, l))
                  t1 = t1 + ao(l)*(tmat_1(l, k) + tmat_1(k, l))
                  t2 = t2 + ao(l)*(tmat_2(l, k) + tmat_2(k, l))
               end do
            end do
      
            do ixyz = 1, 3
               ic = 3*(nr_center - 1) + ixyz
               g(ic) = g(ic) - 4.0d0*d_nn*ao1(k, ixyz)*t1*n_t2 &
                             - 4.0d0*d_nn*ao1(k, ixyz)*t2*n_t1 &
                             - 4.0d0*d_nnn*ao1(k, ixyz)*t0*n_t1*n_t2
            end do
         end do
      end do
   
   end subroutine

end module
