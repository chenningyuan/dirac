!dirac_copyright_start
!dirac_copyright_end

module xc_london_c1

! radovan: this is the "equation on paper approach" to LAO XC
!
!      (+) the code is very close to equations
!      (+) code does not need "magnetic derivatives" of AOs
!          and works with undifferentiated pointwise AO overlap
!      (+) i think it should be very easy to do higher order LAO
!          contributions to KS matrix and to KS response
!          using this approach
!          it should boil down to additional factorization of what we already have
!      (-) it does not work in symmetry adapted (SA) AO basis
!          because we cannot calculate SA AOs and "SA" R_kl in an uncoupled way
!          one could still use R_kl but introduce additional loops
!          in distribution routines and work with non-SA AOs
!          and transform to SA in each matrix element
!          but i think this would be very expensive
!          anyway this code can be useful for development in C1

! todo (*) work only on matrix blocks that couple AOs of different centers
!      (*) figure out centers of AOs in dft_ao

  use memory_allocator
  use basis_info

  implicit none

  public xc_london_c1_init
  public xc_london_c1_lda
  public xc_london_c1_gga

  private

!radovan: these common blocks will go away
!         p array should be set up in dft_ao
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

!                         ideally calculated only once
  real(8), allocatable :: p(:, :)
!                         ideally calculated only once
  real(8), allocatable :: r_kl(:, :, :)
  logical              :: all_set = .false.

  integer              :: b(3) = (/2, 3, 1/)
  integer              :: c(3) = (/3, 1, 2/)

contains

  subroutine xc_london_c1_init(mat_dim)

!   ----------------------------------------------------------------------------
    integer, intent(in) :: mat_dim
!   ----------------------------------------------------------------------------
    integer             :: i, j, k, l, ixyz, ishell, isym
    real(8)             :: px, py, pz
!   ----------------------------------------------------------------------------

    if (.not. all_set) then

      call alloc(r_kl, mat_dim, mat_dim, 3)
      call alloc(p, mat_dim, 3)

      i = 1
      do ishell = 1, kmax
        do isym = 0, maxopr
          if (iand(isym, istbao(ishell)) == 0) then

            px = pt(iand(isymax(1, 1), isym))*cent(ishell, 1, 1)
            py = pt(iand(isymax(2, 1), isym))*cent(ishell, 2, 1)
            pz = pt(iand(isymax(3, 1), isym))*cent(ishell, 3, 1)

            do j = i, i + spherical_deg(ishell) - 1
              p(j, 1) = px
              p(j, 2) = py
              p(j, 3) = pz
            end do

            i = i + spherical_deg(ishell)
          end if
        end do
      end do

      do ixyz = 1, 3
        do k = 1, mat_dim
          do l = 1, mat_dim
            r_kl(k, l, ixyz) = p(k, ixyz) - p(l, ixyz)
          end do
        end do
      end do

      call dealloc(p)

      all_set = .true.
    end if

  end subroutine

  subroutine xc_london_c1_lda(d_r,     &
                              mat_dim, &
                              nz,      &
                              fmat,    &
                              ao,      &
                              r)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: d_r
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(inout) :: fmat(mat_dim, mat_dim, *)
    real(8), intent(in)    :: ao(mat_dim)
    real(8), intent(in)    :: r(3)
!   ----------------------------------------------------------------------------
    integer                :: k, l, bxyz
    real(8)                :: f
    real(8)                :: rkl_cross_r
    real(8)                :: omega
!   ----------------------------------------------------------------------------

    do k = 1, mat_dim
!     work only on the triangle without diagonal elements
      do l = 1, k - 1

        omega = ao(k)*ao(l)

        do bxyz = 1, 3
          rkl_cross_r = r_kl(k, l, b(bxyz))*r(c(bxyz)) &
                      - r_kl(k, l, c(bxyz))*r(b(bxyz))
         
          f = 0.5d0*d_r*rkl_cross_r*omega
         
          fmat(k, l, (bxyz-1)*nz + 1) = fmat(k, l, (bxyz-1)*nz + 1) + f
          fmat(l, k, (bxyz-1)*nz + 1) = fmat(l, k, (bxyz-1)*nz + 1) - f
        end do
      end do
    end do

  end subroutine

  subroutine xc_london_c1_gga(d_r,     &
                              d_z,     &
                              g,       &
                              mat_dim, &
                              nz,      &
                              fmat,    &
                              ao,      &
                              ao_g1,   &
                              r)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: d_r
    real(8), intent(in)    :: d_z
    real(8), intent(in)    :: g(3)
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(inout) :: fmat(mat_dim, mat_dim, *)
    real(8), intent(in)    :: ao(mat_dim)
    real(8), intent(in)    :: ao_g1(mat_dim, 3)
    real(8), intent(in)    :: r(3)
!   ----------------------------------------------------------------------------
    integer                :: k, l, bxyz, gxyz
    real(8)                :: f
    real(8)                :: rkl_cross_r, rkl_cross_g
    real(8)                :: omega, nabla_omega(3)
    real(8), external      :: ddot
!   ----------------------------------------------------------------------------

    do k = 1, mat_dim
!     work only on the triangle without diagonal elements
      do l = 1, k - 1

        omega = ao(k)*ao(l)
        do gxyz = 1, 3
          nabla_omega(gxyz) = ao(k)*ao_g1(l, gxyz) + ao_g1(k, gxyz)*ao(l)
        end do

        do bxyz = 1, 3
          rkl_cross_r = r_kl(k, l, b(bxyz))*r(c(bxyz)) &
                      - r_kl(k, l, c(bxyz))*r(b(bxyz))
          rkl_cross_g = r_kl(k, l, b(bxyz))*g(c(bxyz)) &
                      - r_kl(k, l, c(bxyz))*g(b(bxyz))
         
          f = 0.5d0*d_r*rkl_cross_r*omega &
            +       d_z*rkl_cross_g*omega &
            +       d_z*rkl_cross_r*ddot(3, g, 1, nabla_omega, 1)
         
          fmat(k, l, (bxyz-1)*nz + 1) = fmat(k, l, (bxyz-1)*nz + 1) + f
          fmat(l, k, (bxyz-1)*nz + 1) = fmat(l, k, (bxyz-1)*nz + 1) - f
        end do
      end do
    end do

  end subroutine

end module
