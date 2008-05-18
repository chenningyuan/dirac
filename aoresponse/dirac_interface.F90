module dirac_interface_m

  use kinds_m
  use memory_allocator
  use dirac_matrix_operations_m
  use aoresponse_input_m

  implicit none

  public get_mo_coef
  public get_mo_fmat
  public get_zo_ao_dmat
  public get_ao_pmat
  public get_ao_smat

  private

contains

  subroutine get_mo_fmat(F)

    type(matrix), intent(inout) :: F

#include "dcbbas.h"

    type(matrix)                :: F1, F2

    integer, parameter          :: file_unit = 70
    logical                     :: file_exists

    file_exists = .false.
    inquire(file = 'DFFCK1', exist = file_exists)
    if (.not. file_exists) then
      call quit('error in get_mo_fmat: file DFFCK1 not found')
    end if

    file_exists = .false.
    inquire(file = 'DFFCK2', exist = file_exists)
    if (.not. file_exists) then
      call quit('error in get_mo_fmat: file DFFCK2 not found')
    end if

    if (.not. isdef(F)) call mat_init(F, ntbas(0), ntbas(0))

    F1 = 1.0d-100*F
    F2 = 1.0d-100*F

    open(file_unit,       &
         file = 'DFFCK1', &
         status = 'old',  &
         form = 'unformatted')

    rewind file_unit

    call reafck(file_unit, F1%elements, .true., 1)

    close(file_unit, status = 'keep')

    open(file_unit,       &
         file = 'DFFCK2', &
         status = 'old',  &
         form = 'unformatted')

    rewind file_unit

    call reafck(file_unit, F2%elements, .true., 1)

    close(file_unit, status = 'keep')

!   F = 2.0d0*F1 + 0.25d0*F2
    F = 2.0d0*F1 + F2

    F1 = 0
    F2 = 0

  end subroutine get_mo_fmat

  subroutine get_ao_smat(S)

    type(matrix), intent(inout) :: S

#include "dcbbas.h"

    real(kind_dp), allocatable :: work(:)
    integer(kind_i4)           :: kfree, lfree, kwork, lwork

    kfree = 1

#include "memint.h"

    if (.not. isdef(S)) call mat_init(S, ntbas(0), ntbas(0))

    lfree = nnbbasx + n2bbasx
    lfree = lfree*2 !fixme: factor 2 shouldn't be needed
    call alloc(work, lfree)

    call gtovlt(S%elements, 1.0d0, work, lfree)

    call dealloc(work)

  end subroutine get_ao_smat

  subroutine get_mo_coef(C)

    type(matrix), intent(inout) :: C

#include "dcborb.h"
#include "dcbbas.h"

    real(kind_dp), allocatable :: eigenvalues(:)
    integer, parameter         :: file_unit = 70
    logical                    :: file_exists
    real(kind_dp)              :: dummy_dp
    integer(kind_i4)           :: dummy_i4

    file_exists = .false.
    inquire(file = 'DFCOEF', exist = file_exists)
    if (.not. file_exists) then
      call quit('error in get_mo_coef: file DFCOEF not found')
    end if

    open(file_unit,            &
         file = 'DFCOEF',      &
         status = 'old',       &
         form = 'unformatted', &
         access = 'sequential')

    rewind file_unit

    call alloc(eigenvalues, norbt)

!   if (.not. isdef(C)) call mat_init(C, ntbas(0), norbt)
    if (.not. isdef(C)) call mat_init(C, ntbas(0), ntbas(0))

    call reacmo(file_unit,   &            
                'DFCOEF',    &           
                C%elements,  &             
                eigenvalues, &              
                dummy_i4,    &           
                dummy_dp,    &           
                6,           &     
                dummy_dp,    &           
                dummy_i4,    &           
                dummy_i4)

    call dealloc(eigenvalues)

    close(file_unit, status = 'keep')

  end subroutine get_mo_coef

  subroutine get_zo_ao_dmat(C, D)

    type(matrix), intent(in)    :: C
    type(matrix), intent(inout) :: D

#include "dcbbas.h"

    if (.not. isdef(D)) call mat_init(D, ntbas(0), ntbas(0))

    call genden(D%elements, C%elements, 1, 0)

    D = 2.0d0*D

  end subroutine get_zo_ao_dmat

  subroutine get_ao_pmat(op_name, A)

    character(*), intent(in)    :: op_name
    type(matrix), intent(inout) :: A

#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbxpr.h"

    real(kind_dp), allocatable :: ptri(:)
    real(kind_dp), allocatable :: op1int(:)
    logical,       allocatable :: first(:)

    real(kind_dp), allocatable :: work(:)
    integer(kind_i4)           :: kfree, lfree

    integer(kind_i4)           :: op_index
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
      call quit('error in get_ao_pmat: file AOPROPER not found')
    end if

    open(file_unit,           &
         file   = 'AOPROPER', &
         status = 'unknown',  &
         form   = 'unformatted')

    rewind file_unit

    if (.not. isdef(A)) call mat_init(A, ntbas(0), ntbas(0))

    op_index = find_op_index(op_name)

    call prpex2(op_index,   &
                dummy_dp,   &
                .false.,    &
                idimension, &
                ptri,       &
                A%elements, &
                op1int,     &
                first,      &
                work, kfree, lfree, 0)

    A%pg_sym        = iprpsym(op_index)
    A%tr_sym        = iprptim(op_index)
    A%fermi_pointer = jbtof(A%pg_sym - 1, 1)
!!!!!!!!!!!
    A%ih_sym        = 1
!!!!!!!!!!!

    close(file_unit, &
          status = 'keep')

    call dealloc(ptri)
    call dealloc(op1int)
    deallocate(first)
    call dealloc(work)

  end subroutine get_ao_pmat

end module dirac_interface_m
