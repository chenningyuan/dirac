module xc_london_c1

! radovan: this is the "equation on paper approach" to LAO XC
!
!      (+) the code is very close to equations
!      (+) code does not need "magnetic derivatives" of AOs
!      (+) i think it should be very easy to do higher order LAO
!          contributions to KS matrix and to KS response
!          using this approach
!      (-) it does not work in symmetry adapted (SA) AO basis
!          because we cannot calculate SA AOs and "SA" R_kl in an uncoupled way
!          one could still use R_kl but introduce additional loops
!          in distribution routines and work with non-SA AOs
!          and transform to SA in each matrix element
!          but i think this would be very expensive
!          anyway this code can be useful for development in C1

  use memory_allocator
  use dft_ao

  implicit none

  public xc_london_c1_init
  public xc_london_c1_lda

  private

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
  real(8), allocatable :: R_kl(:, :, :)
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

      call alloc(R_kl, mat_dim, mat_dim, 3)
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
            R_kl(k, l, ixyz) = p(k, ixyz) - p(l, ixyz)
          end do
        end do
      end do

      call dealloc(p)

      all_set = .true.
    end if

  end subroutine

  subroutine xc_london_c1_lda(s, mat_dim, nz, fmat, ao, cx, cy, cz)

    use dft_cb_interface

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(inout) :: fmat(mat_dim, mat_dim, *)
    real(8), intent(in)    :: ao(mat_dim)
    real(8), intent(in)    :: cx, cy, cz
!   ----------------------------------------------------------------------------
    integer                :: k, l, ixyz
    real(8)                :: f, r(3)
!   ----------------------------------------------------------------------------

    r(1) = cx
    r(2) = cy
    r(3) = cz

    do ixyz = 1, 3
      do k = 1, mat_dim
!       work only on the triangle without diagonal elements
        do l = 1, k - 1
          f = (R_kl(k, l, b(ixyz))*r(c(ixyz)) &
             - R_kl(k, l, c(ixyz))*r(b(ixyz)))

          f = 0.5d0*s*f*ao(k)*ao(l)

          fmat(k, l, (ixyz-1)*nz + 1) = fmat(k, l, (ixyz-1)*nz + 1) + f
          fmat(l, k, (ixyz-1)*nz + 1) = fmat(l, k, (ixyz-1)*nz + 1) - f
        end do
      end do

    end do

  end subroutine

end module
