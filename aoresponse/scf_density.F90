module scf_density_m

  use kinds_m
  use memory_allocator
  use matrix_genop
  use matrix_defop

  implicit none

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

contains

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

        call general_density_block(density%scalar,    &
                                   M%elms(1),         &
                                   get_density,       &
                                   M%pg_sym - 1,      &
                                   gao_slice_left,    &
                                   gao_slice_right,   &
                                   buffer,            &
                                   1,                 &
                                   component_lx,      &
                                   density%factor_lx)

        call general_density_block(density%scalar,    &
                                   M%elms(1),         &
                                   get_density,       &
                                   M%pg_sym - 1,      &
                                   gao_slice_left,    &
                                   gao_slice_right,   &
                                   buffer,            &
                                   component_xs,      &
                                   2,                 &
                                   density%factor_xs)

    end if

  end subroutine

end module
