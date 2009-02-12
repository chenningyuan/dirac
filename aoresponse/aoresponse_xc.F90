module aoresponse_xc_m

  use kinds_m
  use slater_m
  use scf_density_m
  use grid_m
  use memory_allocator
! use matrix_oper
  use matrix_defop

  implicit none

  real(kind_dp), parameter :: n_0_threshold             = 1.0d-14
  real(kind_dp), parameter :: xc_contribution_threshold = 1.0d-14

contains

! subroutine xc_contribution(e, order, D_0, D, F)
  subroutine xc_contribution(order, D_0, D, F)

!   real(kind_dp)    :: e(*)
    integer          :: order
    type(matrix)     :: D_0, D(:), F(:)
!   ----------------------------------------------------------------------------
!   type(matrix)     :: M(9) !fixme later, hardcoded for 2n+1 QR
    type(grid_t)     :: grid
    type(gao_data_t) :: gao_data
    integer          :: i, j, k, ij, ijk, iz, iq

#include "dgroup.h"


!   scale matrices
!   ==============

    D_0 = 4.0d0*D_0
    do i = 1, size(D)
      D(i) = 4.0d0*D(i)
    end do


!   half-phase
!   ==========

    if (nz < 4) then
      do iz = 1, nz
        do i = 1, size(D)
          iq = ipqtoq(iz, D(i)%pg_sym - 1)
          call q2bphase('D', iq, 1, D(i)%elms((iz-1)*D(i)%nrow*D(i)%ncol + 1))
        end do
        do i = 1, size(F)
          iq = ipqtoq(iz, F(i)%pg_sym - 1)
          call q2bphase('F', iq, 1, F(i)%elms((iz-1)*F(i)%nrow*F(i)%ncol + 1))
        end do
      end do
    end if


!   allocate "mix" matrices
!   =======================

!   if (order == 3) then
!     do i = 1, 3
!       do j = 1, 3
!         ij = i + (j-1)*3
!         M(ij) = 1.0d-100*D_0
!         M(ij)%pg_sym = 1 + ieor(D(i)%pg_sym - 1, D(j+3)%pg_sym - 1)
!       end do
!     end do
!   end if


!   integrate
!   =========

    call init_dft_grid(D_0, grid)
    call init_gao_data(gao_data, 0)

    call xc_contribution_backend(order, D_0, D, F, grid, gao_data)

    call del_grid(grid)
    call del_gao_data(gao_data)


!   get contribution from "mix" matrices and deallocate them
!   ========================================================

!   if (order == 3) then
!     do i = 1, 3
!       do j = 1, 3
!         ij = i + (j-1)*3
!         do k = 1, 3
!           ijk = ij + (k-1)*9

!           e(ijk) = e(ijk) + 0.5d0*(M(ij)**D(k+6))
!         end do
!         M(ij) = 0
!       end do
!     end do
!   end if


!   scale matrices
!   ==============

    D_0 = 0.25d0*D_0
    do i = 1, size(D)
      D(i) = 0.25d0*D(i)
    end do


!   half-phase
!   ==========

    if (nz < 4) then
      do iz = 1, nz
        do i = 1, size(D)
          iq = ipqtoq(iz, D(i)%pg_sym - 1)
          call q2bphase('D', iq, 1, D(i)%elms((iz-1)*D(i)%nrow*D(i)%ncol + 1))
        end do
        do i = 1, size(F)
          iq = ipqtoq(iz, F(i)%pg_sym - 1)
          call q2bphase('F', iq, 1, F(i)%elms((iz-1)*F(i)%nrow*F(i)%ncol + 1))
        end do
      end do
    end if

  end subroutine

  subroutine xc_contribution_backend(order, D_0, D, F, grid, gao_data)
  
    integer          :: order
    type(matrix)     :: D_0, D(:), F(:)
    type(grid_t)     :: grid 
    type(gao_data_t) :: gao_data
!   ----------------------------------------------------------------------------
    type(density_t)  :: n_0
    type(density_t)  :: n(100)
    type(density_t)  :: n_xc(100)
    type(density_t)  :: n_xc_2(100)
    integer          :: i, j, ipoint
    real(kind_dp)    :: derivative
    real(kind_dp)    :: worth_it
    real(kind_dp)    :: dftenergy

    call init_density(n_0, '1')

    do i = 1, size(D)
      call init_density(n(i),    '1')
      call init_density(n_xc(i), '1')
    end do 
    
!   do i = 1, size(M)
!     call init_density(n_xc_2(i), '1')
!   end do 
    
    do ipoint = 1, grid%nr_points

      call get_gao_data(gao_data,       &
                        grid%x(ipoint), &
                        grid%y(ipoint), &
                        grid%z(ipoint)) 

      call get_or_put_density(n_0, D_0, .true., gao_data)

      if (n_0%scalar > n_0_threshold) then

!       e(1) = e(1) + dftenergy(n_0%scalar, 0.0d0)*grid%w(ipoint)


!       evaluate perturbed densities
!       ============================

        do i = 1, size(D)
          call get_or_put_density(n(i), D(i), .true., gao_data)
        end do


!       evaluate prefactors
!       ===================

        derivative     = slater_derv(order, n_0%scalar, grid%w(ipoint))

        select case (order)

          case (2)
            n_xc(1)%scalar = 2.0d0*derivative*n(1)%scalar

          case (3)
            n_xc(1)%scalar = 2.0d0*derivative*n(1)%scalar*n(2)%scalar

        end select


!       fill matrices
!       =============

        do i = 1, size(F)
          if (abs(n_xc(i)%scalar) > xc_contribution_threshold) then
            call get_or_put_density(n_xc(i), F(i), .false., gao_data)
          end if
        end do
!       do i = 1, size(M)
!         if (abs(n_xc_2(i)%scalar) > xc_contribution_threshold) then
!           call get_or_put_density(n_xc_2(i), M(i), .false., gao_data)
!         end if
!       end do
      end if
    end do

  end subroutine

end module
