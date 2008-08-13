module aoresponse_xc_m

  use kinds_m
  use svwn_m
  use memory_allocator
  use matrix_oper
  use matrix_mod

  implicit none

  type dft_grid_t
    integer                    :: nr_points
    real(kind_dp), allocatable :: x(:)
    real(kind_dp), allocatable :: y(:)
    real(kind_dp), allocatable :: z(:)
    real(kind_dp), allocatable :: w(:)
  end type

  type gao_data_t
    integer                    :: max_derivative
    integer,       allocatable :: gao_center(:)
    real(kind_dp), allocatable :: gao(:)
    real(kind_dp), allocatable :: gao1(:)
    real(kind_dp), allocatable :: gao2(:)
    real(kind_dp), allocatable :: buffer(:)
  end type

  type density_t
    real(kind_dp) :: scalar
    logical       :: is_vector
    real(kind_dp) :: vector(3)
    integer       :: tr_sym
    integer       :: transforms_as
    logical       :: diagonal
    logical       :: times_minus_im
    real(kind_dp) :: factor_lx
    real(kind_dp) :: factor_xs
  end type

  real(kind_dp), parameter :: n_0_threshold             = 1.0d-14
  real(kind_dp), parameter :: xc_contribution_threshold = 1.0d-14

contains

  subroutine xc_contribution(e, order, D_0, D, F)

    real(kind_dp)    :: e(*)
    integer          :: order
    type(matrix)     :: D_0, D(:), F(:)
!   ----------------------------------------------------------------------------
    type(matrix)     :: M(9) !fixme later, hardcoded for 2n+1 QR
    type(dft_grid_t) :: grid
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

    if (order == 3) then
      do i = 1, 3
        do j = 1, 3
          ij = i + (j-1)*3
          M(ij) = 1.0d-100*D_0
          M(ij)%pg_sym = 1 + ieor(D(i)%pg_sym - 1, D(j+3)%pg_sym - 1)
        end do
      end do
    end if


!   integrate
!   =========

    call get_grid(D_0, grid)
    call init_gao_data(gao_data, 0)
    call xc_contribution_backend(order, D_0, D, F, M, grid, gao_data)
    call del_gao_data(gao_data)


!   get contribution from "mix" matrices and deallocate them
!   ========================================================

    if (order == 3) then
      do i = 1, 3
        do j = 1, 3
          ij = i + (j-1)*3
          do k = 1, 3
            ijk = ij + (k-1)*9

            e(ijk) = e(ijk) + 0.5d0*(M(ij)**D(k+6))
          end do
          M(ij) = 0
        end do
      end do
    end if


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

  subroutine xc_contribution_backend(order, D_0, D, F, M, grid, gao_data)

    integer          :: order
    type(matrix)     :: D_0, D(:), F(:), M(:)
    type(dft_grid_t) :: grid
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

    do i = 1, size(M)
      call init_density(n_xc_2(i), '1')
    end do

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

        derivative     = slater_recursive(order, n_0%scalar, grid%w(ipoint))
        n_xc(1)%scalar = 2.0d0*derivative*n(1)%scalar

        if (order > 2) then
          derivative       = slater_recursive(order, n_0%scalar, grid%w(ipoint))
          n_xc_2(1)%scalar = derivative*n(1)%scalar*n(4)%scalar
          n_xc_2(2)%scalar = derivative*n(1)%scalar*n(5)%scalar
          n_xc_2(3)%scalar = derivative*n(1)%scalar*n(6)%scalar
          n_xc_2(4)%scalar = derivative*n(2)%scalar*n(4)%scalar
          n_xc_2(5)%scalar = derivative*n(2)%scalar*n(5)%scalar
          n_xc_2(6)%scalar = derivative*n(2)%scalar*n(6)%scalar
          n_xc_2(7)%scalar = derivative*n(3)%scalar*n(4)%scalar
          n_xc_2(8)%scalar = derivative*n(3)%scalar*n(5)%scalar
          n_xc_2(9)%scalar = derivative*n(3)%scalar*n(6)%scalar

          derivative       = slater_recursive(order-1, n_0%scalar, grid%w(ipoint))
          do i = 1, size(D)
            n_xc(i)%scalar = 2.0d0*derivative*n(i)%scalar
          end do
        end if


!       fill matrices
!       =============

        do i = 1, size(F)
          if (abs(n_xc(i)%scalar) > xc_contribution_threshold) then
            call get_or_put_density(n_xc(i), F(i), .false., gao_data)
          end if
        end do
        do i = 1, size(M)
          if (abs(n_xc_2(i)%scalar) > xc_contribution_threshold) then
            call get_or_put_density(n_xc_2(i), M(i), .false., gao_data)
          end if
        end do
      end if
    end do

  end subroutine

  subroutine init_density(self, kind)

    type(density_t), intent(inout) :: self
    character,       intent(in)    :: kind

#include "dcbgen.h"

    self%scalar = 0.0d0
    self%vector = 0.0d0

    select case (kind)
      case ('s', 'a', 'iba', 'bs')
        self%is_vector = .true.
      case ('1', 'g', 'ibg', 'b')
        self%is_vector = .false.
    end select

    select case (kind)
      case ('1', 'g', 'iba', 'b')
        self%tr_sym =  1
      case ('s', 'a', 'ibg', 'bs')
        self%tr_sym = -1
    end select

    select case (kind)
      case ('1', 'g', 'ibg', 'b')
        self%transforms_as = 0
      case ('a', 'iba')
        self%transforms_as = 1
      case ('s', 'bs')
        self%transforms_as = 2
    end select

    select case (kind)
      case ('1', 's', 'b', 'bs')
        self%diagonal = .true.
      case ('g', 'a', 'ibg', 'iba')
        self%diagonal = .false.
    end select

    select case (kind)
      case ('ibg', 'iba')
        self%times_minus_im = .true.
      case ('1', 's', 'g', 'a', 'b', 'bs')
        self%times_minus_im = .false.
    end select

    select case (kind)
      case ('1', 's', 'g')
        self%factor_lx =  1.0d0
        self%factor_xs =  1.0d0
      case ('ibg', 'b', 'bs')
        self%factor_lx =  1.0d0
        self%factor_xs = -1.0d0
      case ('a')
        self%factor_lx = -cval
        self%factor_xs = -cval
      case ('iba')
        self%factor_lx = -cval
        self%factor_xs =  cval
    end select

  end subroutine

  subroutine init_gao_data(self, max_derivative)

    type(gao_data_t), intent(inout) :: self
    integer                         :: max_derivative

#include "dcbbas.h"

    self%max_derivative = max_derivative

    call alloc(self%gao_center, ntbas(0))
    call alloc(self%buffer, 4*ntbas(0))

    call alloc(self%gao, ntbas(0))
    if (self%max_derivative > 0) then
      call alloc(self%gao1, 3*ntbas(0))
    end if
    if (self%max_derivative > 1) then
      call alloc(self%gao2, 6*ntbas(0))
    end if

  end subroutine

  subroutine del_gao_data(self)

    type(gao_data_t), intent(inout) :: self

    call dealloc(self%gao_center)
    call dealloc(self%buffer)

    call dealloc(self%gao)
    if (self%max_derivative > 0) then
      call dealloc(self%gao1)
    end if
    if (self%max_derivative > 1) then
      call dealloc(self%gao2)
    end if

  end subroutine

  subroutine get_gao_data(self, x, y, z)

    type(gao_data_t), intent(inout) :: self
    real(kind_dp),    intent(in)    :: x
    real(kind_dp),    intent(in)    :: y
    real(kind_dp),    intent(in)    :: z
!   ----------------------------------------------------------------------------
    real(kind_dp)                   :: dummy_dp
    integer                         :: dummy_i

    self%buffer = 0.0d0

    call getsos(self%gao,            &
                self%gao1,           &
                self%gao2,           &
                dummy_dp,            &
                self%gao_center,     &
                x,                   &
                y,                   &
                z,                   &
                self%buffer,         &
                dummy_dp,            &
                dummy_i,             &
                self%gao_center,     &
                self%max_derivative, &
                .false., 0)

  end subroutine

  subroutine get_or_put_density(density, M, get_density, gao_data)

    type(density_t),  intent(inout), target :: density
    type(matrix),     intent(inout)         :: M
    logical,          intent(in)            :: get_density
    type(gao_data_t), intent(in),    target :: gao_data
!   ----------------------------------------------------------------------------
    real(kind_dp),    pointer               :: gao_slice_left(:)
    real(kind_dp),    pointer               :: gao_slice_right(:)
    real(kind_dp),    pointer               :: buffer(:)
    integer                                 :: component_lx
    integer                                 :: component_xs
    integer                                 :: ixyz, iz, offset

#include "dgroup.h"
#include "dcbbas.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "symmet.h"

    if (get_density) then
      density%scalar = 0.0d0
      density%vector = 0.0d0
    end if

    if (density%diagonal) then
      component_lx = 1
      component_xs = 2
    else
      component_lx = 2
      component_xs = 1
    end if

    gao_slice_left  => gao_data%gao
    gao_slice_right => gao_data%gao
    buffer          => gao_data%buffer

    if (density%is_vector) then

      do ixyz = 1, 3

        iz = jm4pos(4 - ixyz)
        iz = iqtopq(iz, isymax(ixyz, 2))

        offset = (iz-1)*ntbas(0)*ntbas(0)

        call general_density_block(density%vector(ixyz),                                  &
                                   M%elms(1+offset),                                      &
                                   get_density,                                           &
                                   ieor(isymax(ixyz, density%transforms_as), M%pg_sym-1), &
                                   gao_slice_left,                                        &
                                   gao_slice_right,                                       &
                                   buffer,                                                &
                                   1,                                                     &
                                   component_lx,                                          &
                                   density%factor_lx)

        call general_density_block(density%vector(ixyz),                                  &
                                   M%elms(1+offset),                                      &
                                   get_density,                                           &
                                   ieor(isymax(ixyz, density%transforms_as), M%pg_sym-1), &
                                   gao_slice_left,                                        &
                                   gao_slice_right,                                       &
                                   buffer,                                                &
                                   component_xs,                                          &
                                   2,                                                     &
                                   density%factor_xs)

      end do

    else

        call general_density_block(density%scalar,                                        &
                                   M%elms(1+offset),                                      &
                                   get_density,                                           &
                                   M%pg_sym - 1,                                          &
                                   gao_slice_left,                                        &
                                   gao_slice_right,                                       &
                                   buffer,                                                &
                                   1,                                                     &
                                   component_lx,                                          &
                                   density%factor_lx)

        call general_density_block(density%scalar,                                        &
                                   M%elms(1+offset),                                      &
                                   get_density,                                           &
                                   M%pg_sym - 1,                                          &
                                   gao_slice_left,                                        &
                                   gao_slice_right,                                       &
                                   buffer,                                                &
                                   component_xs,                                          &
                                   2,                                                     &
                                   density%factor_xs)

    end if

  end subroutine

  subroutine get_grid(D, grid)

    type(matrix),     intent(inout) :: D
    type(dft_grid_t), intent(inout) :: grid
!   ----------------------------------------------------------------------------
    real(kind_dp),    allocatable   :: work(:)
    integer                         :: kfree, lfree
    integer                         :: file_unit
    integer                         :: nr_points
    integer                         :: counter
    integer                         :: i

    file_unit = 44

    kfree = 1
    lfree = 10000000
    call alloc(work, lfree)

    call dftgrd(D%elms, work(kfree), lfree, 1)

    call dealloc(work)


!   first round: get the total nr_points
!   ====================================

    rewind file_unit
    nr_points = 137 ! just set it to something > 0

    counter = 0
    do while (nr_points > 0)

      read(file_unit) nr_points

      if (nr_points > 0) then

        read(file_unit)
        counter = counter + nr_points

      end if
    end do


!   second round: allocate and read
!   ===============================

    grid%nr_points = counter

    call alloc(grid%x, grid%nr_points)
    call alloc(grid%y, grid%nr_points)
    call alloc(grid%z, grid%nr_points)
    call alloc(grid%w, grid%nr_points)

    rewind file_unit
    nr_points = 137 ! just set it to something > 0

    counter = 0
    do while (nr_points > 0)

      read(file_unit) nr_points

      if (nr_points > 0) then

        call reaqua(grid%x(1+counter:1+counter+nr_points), &
                    grid%y(1+counter:1+counter+nr_points), &
                    grid%z(1+counter:1+counter+nr_points), &
                    grid%w(1+counter:1+counter+nr_points), &
                    file_unit,           &
                    nr_points)

        counter = counter + nr_points

      end if
    end do

    rewind file_unit
!   close(file_unit, status = 'keep')

  end subroutine

end module
