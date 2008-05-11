module dirac_interface_m

  use kinds_m
  use memory_allocator
  use dirac_matrix_operations_m
  use aoresponse_input_m

  implicit none

contains

  subroutine get_mo_coef(C)
    type(matrix), intent(inout) :: C

#include "dcborb.h"

    real(kind_dp), allocatable :: eigenvalues(:)
    integer, parameter         :: file_unit = 70
    logical                    :: file_exists
    real(kind_dp)              :: dummy_dp
    integer(kind_i4)           :: dummy_i4

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

    call alloc(eigenvalues, norbt)

    call reacmo(file_unit, 'DFCOEF', C%elements, eigenvalues, dummy_i4, dummy_dp, 6, &
                dummy_dp, dummy_i4, dummy_i4)

    call dealloc(eigenvalues)

    close(file_unit, &
          status = 'keep')

  end subroutine get_mo_coef

  subroutine get_zo_ao_dmat(C, D)
    type(matrix), intent(in)  :: C
    type(matrix), intent(inout) :: D

!   D = 0.0d0*D
    call genden(D%elements, C%elements, 1, 0)
    D = 2.0d0*D

  end subroutine get_zo_ao_dmat

  subroutine get_ao_pmat(op_name, A)
    character(*), intent(in)  :: op_name
    type(matrix), intent(out) :: A

#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbxpr.h"

    real(kind_dp), allocatable :: ptri(:)
    real(kind_dp), allocatable :: op1int(:)
    logical,       allocatable :: first(:)
    real(kind_dp), allocatable :: work(:)

    integer(kind_i4)           :: kfree, lfree
    logical                    :: file_exists, is_square
    real(kind_dp)              :: dummy_dp
    integer(kind_i4)           :: file_unit, dummy_i4, idimension

    is_square = .false.
    if (iprptim(find_op_index(op_name)) == 0) then
      is_square = .true.
    end if

    if (is_square) then
      idimension = n2bbasx
    else
      idimension = nnbbasx
    end if

    call alloc(ptri,   nz*idimension)
    call alloc(op1int, idimension)
    allocate(first(nz))

    kfree = 1
    lfree = 2*nz*n2bbasx
    call alloc(work, lfree)

    file_unit = lu1int

    file_exists = .false.
    inquire(file = 'AOPROPER', exist = file_exists)
    if (.not. file_exists) then
      call quit('file aoproper not found')
    end if

    open(file_unit,           &
         file   = 'AOPROPER', &
         status = 'unknown',  &
         form   = 'unformatted')

    rewind file_unit

    if (.not. isdef(A)) call mat_init(A, ntbas(0), ntbas(0))

    call prpex2(find_op_index(op_name), &
                dummy_dp,               &
                .false.,                &
                idimension,             &
                ptri,                   &
                A%elements,             &
                op1int,                 &
                first,                  &
                work, kfree, lfree, 0)

    close(file_unit, &
          status = 'keep')

    call dealloc(ptri)
    call dealloc(op1int)
    deallocate(first)
    call dealloc(work)

  end subroutine get_ao_pmat

end module dirac_interface_m
