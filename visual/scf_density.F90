
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

  public init_ao_dmat, del_ao_dmat
  public init_mo_coef, del_mo_coef
  public get_mo_coef, get_ao_dmat_zo

  public gao_data_t
  public init_gao_data, del_gao_data, get_gao_data

  public get_or_put_density, general_density_block

  public print_square_q_matrix

  private

! blub

contains

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

    call getsos(self%gao, self%gao1, self%gao2,           &
                dummy_dp, self%gao_center, p%x, p%y, p%z, &
                self%buffer, dummy_dp, dummy_i4,          &
                self%gao_center, self%max_derivative,     &
                .false., 0)
  end subroutine get_gao_data















  subroutine get_or_put_density(density, matrix, get_density, irepop, gao_data)
    real(kind_dp),    intent(inout)      :: density
    real(kind_dp),    intent(inout)      :: matrix(:, :, :)
    logical,          intent(in)         :: get_density
    integer(kind_i4), intent(in)         :: irepop
    type(gao_data_t), intent(in), target :: gao_data

    real(kind_dp),    pointer            :: gao_slice_left(:), gao_slice_right(:), buffer(:)
    logical                              :: diagonal

!   diagonal = .true.
    gao_slice_left  => gao_data%gao
    gao_slice_right => gao_data%gao
    buffer => gao_data%buffer

    call general_density_block(density, matrix, get_density, irepop, gao_slice_left, gao_slice_right, &
                               buffer, 1, 1)
    write(*, *) 'schaka', density
  end subroutine get_or_put_density

  subroutine general_density_block(density, matrix, get_density, irepop, gao_slice_left, gao_slice_right, &
                                   buffer, component_left, component_right)
    real(kind_dp),    intent(inout)      :: density
    real(kind_dp),    intent(inout)      :: matrix(:, :, :)
    logical,          intent(in)         :: get_density
    integer(kind_i4), intent(in)         :: irepop
    real(kind_dp),    pointer            :: gao_slice_left(:), gao_slice_right(:), buffer(:)
    integer(kind_i4), intent(in)         :: component_left, component_right

    integer(kind_i4)                :: i1, i2, n1, n2, j, irep1, irep2
    real(kind_dp)                   :: factor
    real(kind_dp)                   :: ddot
    real(kind_dp)                   :: f = 1.0d0


       do irep1 = 0, nbsym-1

         n1 = nbbas(irep1,component_left)
         if(n1 > 0) then

           irep2 = ieor(irep1,irepop)
           n2    = nbbas(irep2,component_right)
           if(n2 > 0) then

             if(get_density) then

               i1 = ibbas(irep1,component_left) + 1
               i2 = ibbas(irep2,component_right) + 1

               call dgemv('n',nbbas(irep1,component_left),nbbas(irep2,component_right),1.0d0, &
                          matrix(i1,i2,1),ntbas(0),gao_slice_right(i2),1,0.0d0,buffer,1)
               density = density + f*ddot(nbbas(irep1,component_left),gao_slice_left(i1),1,buffer,1)

             else

               i1 = ibbas(irep1,component_left)
               i2 = ibbas(irep2,component_right)

               do j = 1, n2
                 factor = f*density*gao_slice_right(i2+j)
                 call daxpy(n1,factor,gao_slice_left(i1+1),1,matrix(i1+1,i2+j,1),1)
               end do
             end if
           end if
         end if
       end do
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
