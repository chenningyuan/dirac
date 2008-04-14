
!-------------------------------------------------------------------------------
!      radovan bast (2008-04-06)
!-------------------------------------------------------------------------------

module scf_density_class

  use kinds_m
  use vector_algebra_m

  implicit none

#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"

  type gao_data_t
    integer(kind_i4)              :: max_derivative
    integer(kind_i4), allocatable :: gao_center(:)
    real(kind_dp),    allocatable :: gao(:)
    real(kind_dp),    allocatable :: gao1(:)
    real(kind_dp),    allocatable :: gao2(:)
    real(kind_dp),    allocatable :: buffer(:)
  end type gao_data_t

  type density_t
    real(kind_dp)     :: scalar
    logical           :: is_vector
    type(dp_vector_t) :: vector
    integer(kind_i4)  :: t_symmetry
    integer(kind_i4)  :: transforms_as
    logical           :: diagonal
    logical           :: times_minus_im
    real(kind_dp)     :: factor_lx
    real(kind_dp)     :: factor_xs
  end type density_t

  public gao_data_t, density_t

  public init_density
  public init_ao_dmat, del_ao_dmat, init_mo_coef, del_mo_coef, get_mo_coef, get_ao_dmat_zo
  public init_gao_data, del_gao_data, get_gao_data
  public get_or_put_density, general_density_block
  public print_square_q_matrix

  private

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
        self%t_symmetry =  1
      case ('s', 'a', 'ibg', 'bs')
        self%t_symmetry = -1
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
      case ('ibaa')
        self%factor_lx = -cval
        self%factor_xs =  cval
    end select
  end subroutine init_density

  subroutine init_ao_dmat(self)
    real(kind_dp), intent(inout), allocatable :: self(:, :, :)

    allocate(self(ntbas(0), ntbas(0), nz))
    self = 0.0d0
  end subroutine init_ao_dmat

  subroutine del_ao_dmat(self)
    real(kind_dp), intent(inout), allocatable :: self(:, :, :)

    deallocate(self)
  end subroutine del_ao_dmat

  subroutine init_mo_coef(self)
    real(kind_dp), intent(inout), allocatable :: self(:, :, :)

    allocate(self(ntbas(0), ntbas(0), nz))
    self = 0.0d0
  end subroutine init_mo_coef

  subroutine del_mo_coef(self)
    real(kind_dp), intent(inout), allocatable :: self(:, :, :)

    deallocate(self)
  end subroutine del_mo_coef

  subroutine get_mo_coef(self)
    real(kind_dp), intent(inout), target :: self(:, :, :)

    real(kind_dp), allocatable           :: eigenvalues(:)
    integer, parameter                   :: file_unit = 70
    logical                              :: file_exists
    real(kind_dp)                        :: dummy_dp
    integer(kind_i4)                     :: dummy_i4

    file_exists = .false.
    inquire(file = 'DFCOEF', exist = file_exists)
    if (.not. file_exists) then
      call quit('error in get_cmo, file DFCOEF not found')
    end if

    open(file_unit,            &
         file = 'DFCOEF',      &
         status = 'old',       &
         form = 'unformatted', &
         access = 'sequential')

    rewind file_unit

    allocate(eigenvalues(norbt))

    call reacmo(file_unit, 'DFCOEF', self, eigenvalues, dummy_i4, dummy_dp, 6, &
                dummy_dp, dummy_i4, dummy_i4)

    deallocate(eigenvalues)

    close(file_unit, status = 'keep')

!   call print_square_q_matrix(self, ntbas(0))
  end subroutine get_mo_coef

  subroutine get_ao_dmat_zo(ao_dmat_zo, mo_coef)
    real(kind_dp), intent(out) :: ao_dmat_zo(:, :, :)
    real(kind_dp), intent(in)  :: mo_coef(:, :, :)

    ao_dmat_zo = 0.0d0

    call genden(ao_dmat_zo, mo_coef, 1, 0)

    call dscal(n2bbasxq, 2.0d0, ao_dmat_zo, 1)
  end subroutine get_ao_dmat_zo

  subroutine init_gao_data(self)
    type(gao_data_t), intent(inout) :: self

    allocate(self%gao_center(ntbas(0)))
    allocate(self%buffer(4*ntbas(0)))

    allocate(self%gao(ntbas(0)))
    if (self%max_derivative > 0) then
      allocate(self%gao1(3*ntbas(0)))
    end if
    if (self%max_derivative > 1) then
      allocate(self%gao2(6*ntbas(0)))
    end if
  end subroutine init_gao_data

  subroutine del_gao_data(self)
    type(gao_data_t), intent(inout) :: self

    deallocate(self%gao_center)
    deallocate(self%buffer)

    deallocate(self%gao)
    if (self%max_derivative > 0) then
      deallocate(self%gao1)
    end if
    if (self%max_derivative > 1) then
      deallocate(self%gao2)
    end if
  end subroutine del_gao_data

  subroutine get_gao_data(self, p)
    type(gao_data_t),  intent(inout) :: self
    type(dp_vector_t), intent(in)    :: p

    real(kind_dp)                    :: dummy_dp
    integer(kind_i4)                 :: dummy_i4

    self%buffer = 0.0d0

    call getsos(self%gao, self%gao1, self%gao2,                   &
                dummy_dp, self%gao_center,                        &
                p%component(1),  p%component(2),  p%component(3), &
                self%buffer, dummy_dp, dummy_i4,                  &
                self%gao_center, self%max_derivative,             &
                .false., 0)
  end subroutine get_gao_data

  subroutine get_or_put_density(op, density, matrix, get_density, irepop, gao_data)
    character(*),     intent(in)            :: op
    type(density_t),  intent(inout), target :: density
    real(kind_dp),    intent(inout), target :: matrix(:, :, :)
    logical,          intent(in)            :: get_density
    integer(kind_i4), intent(in)            :: irepop
    type(gao_data_t), intent(in), target    :: gao_data

    real(kind_dp),    pointer               :: density_pointer
    real(kind_dp),    pointer               :: matrix_iz_block(:, :), gao_slice_left(:), gao_slice_right(:), buffer(:)
    integer(kind_i4)                        :: component_lx, component_xs, ixyz, iz

#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "symmet.h"

    if (density%diagonal) then
      component_lx = 1
      component_xs = 2
    else
      component_lx = 2
      component_xs = 1
    end if

    buffer => gao_data%buffer

    gao_slice_left  => gao_data%gao
    gao_slice_right => gao_data%gao

    if (density%is_vector) then
      do ixyz = 1, 3

        iz = jm4pos(4 - ixyz)
        iz = iqtopq(iz, isymax(ixyz, 2))

        matrix_iz_block => matrix(:, :, iz)

!       call general_density_block(density_pointer, matrix_iz_block,                               &
!                                  get_density, ieor(isymax(ixyz, density%transforms_as), irepop), &
!                                  gao_slice_left, gao_slice_right,                                &
!                                  buffer, 1, component_lx, density%factor_lx)
!       call general_density_block(density_pointer, matrix_iz_block,                               &
!                                  get_density, ieor(isymax(ixyz, density%transforms_as), irepop), &
!                                  gao_slice_left, gao_slice_right,                                &
!                                  buffer, component_xs, 2, density%factor_xs)

      end do
    else

      matrix_iz_block => matrix(:, :, 1)

!     call general_density_block(density_pointer, matrix_iz_block, &
!                                get_density, irepop,              &
!                                gao_slice_left, gao_slice_right,  &
!                                buffer, 1, component_lx, density%factor_lx)
!     call general_density_block(density_pointer, matrix_iz_block, &
!                                get_density, irepop,              &
!                                gao_slice_left, gao_slice_right,  &
!                                buffer, component_xs, 2, density%factor_xs)

    end if
  end subroutine get_or_put_density

  subroutine general_density_block(density, matrix, get_density, irepop,    &
                                   gao_slice_left, gao_slice_right, buffer, &
                                   component_left, component_right, scale_by)
    real(kind_dp),    intent(inout)      :: density
    real(kind_dp),    intent(inout)      :: matrix(:, :)
    logical,          intent(in)         :: get_density
    integer(kind_i4), intent(in)         :: irepop
    real(kind_dp),    pointer            :: gao_slice_left(:), gao_slice_right(:), buffer(:)
    integer(kind_i4), intent(in)         :: component_left, component_right
    real(kind_dp),    intent(in)         :: scale_by

    integer(kind_i4)                     :: i1, i2, n1, n2, j, irep1, irep2
    real(kind_dp)                        :: ddot, f

#include "dcbham.h"

    if (.not.((component_left + component_right > 2) .and. (mc == 1))) then
      do irep1 = 0, nbsym - 1
        n1 = nbbas(irep1, component_left)
        if (n1 > 0) then

          irep2 = ieor(irep1, irepop)
          n2 = nbbas(irep2, component_right)
          if (n2 > 0) then

            if (get_density) then
              i1 = ibbas(irep1, component_left)  + 1
              i2 = ibbas(irep2, component_right) + 1

              call dgemv('n', nbbas(irep1, component_left), nbbas(irep2, component_right), 1.0d0, &
                         matrix(i1, i2), ntbas(0), gao_slice_right(i2), 1, 0.0d0, buffer, 1)
              density = density + scale_by*ddot(nbbas(irep1, component_left), gao_slice_left(i1), 1, buffer, 1)
            else
              i1 = ibbas(irep1, component_left)
              i2 = ibbas(irep2, component_right)

              do j = 1, n2
                f = scale_by*density*gao_slice_right(i2 + j)
                call daxpy(n1, f, gao_slice_left(i1 + 1), 1, matrix(i1 + 1, i2 + j), 1)
              end do
            end if
          end if
        end if
      end do
    end if
  end subroutine general_density_block

  subroutine print_square_q_matrix(matrix, matrix_dimension)
    real(kind_dp),    intent(in) :: matrix(:, :, :)
    integer(kind_i4), intent(in) :: matrix_dimension

#include "priunit.h"

    call prqmat(matrix,           &
                matrix_dimension, &
                matrix_dimension, &
                matrix_dimension, &
                matrix_dimension, &
                nz,               &
                ipqtoq(1, 0),     &
                lupri)
  end subroutine print_square_q_matrix

end module scf_density_class
