module aoresponse_input_m

  use kinds_m
  use memory_allocator

  implicit none

  type perturbation_op_t
    integer(kind_i4) :: sym, tr_sym, nr_components, typ
    real(kind_dp)    :: prefactors(2)
    character(8)     :: op_name, labels(2)
  end type perturbation_op_t

  integer, parameter :: max_nr_op = 137

  type(perturbation_op_t), save :: perturbation_op(max_nr_op)

contains

  function find_op_index(op_name) result(op_index)
    character(*), intent(in) :: op_name
    integer(kind_i4)         :: op_index, i

    op_index = -1

    do i = 1, max_nr_op
      if (perturbation_op(i)%op_name == op_name) then
        op_index = i
      end if
    end do

    if (op_index == -1) then
      call quit('find_op_index: index not found')
    end if

  end function find_op_index

  subroutine add_op(op_name,            &
                    typ, nr_components, &
                    prefactors,         &
                    labels)
    character(*),     intent(in)  :: op_name
    integer(kind_i4), intent(in)  :: typ, nr_components
    real(kind_dp),    intent(in)  :: prefactors(nr_components)
    character(*),     intent(in)  :: labels(nr_components)

    integer(kind_i4)              :: indxpr, isyxpr, itrxpr, i

    real(kind_dp),    allocatable :: work(:)
    integer(kind_i4)              :: lwork

    lwork = 200000 !fixme
    call alloc(work, lwork)

    call xprind(op_name,                &
                typ, nr_components,     &
                prefactors,             &
                labels,                 &
                indxpr, isyxpr, itrxpr, &
                work, lwork, 0)

    call dealloc(work)

    perturbation_op(indxpr)%sym           = isyxpr
    perturbation_op(indxpr)%tr_sym        = itrxpr
    perturbation_op(indxpr)%nr_components = nr_components
    perturbation_op(indxpr)%typ           = typ
    perturbation_op(indxpr)%op_name       = op_name

    do i = 1, nr_components
      perturbation_op(indxpr)%prefactors(i) = prefactors(i)
      perturbation_op(indxpr)%labels(i)     = labels(i)
    end do

  end subroutine add_op

  subroutine aoresponse_input()

    integer(kind_i4) :: i

    do i = 1, max_nr_op
      perturbation_op(i)%op_name = ''
    end do

    call add_op('XDIPLEN', 1, 1, (/1.0d0/), (/'XDIPLEN '/))
    call add_op('YDIPLEN', 1, 1, (/1.0d0/), (/'YDIPLEN '/))
    call add_op('ZDIPLEN', 1, 1, (/1.0d0/), (/'ZDIPLEN '/))

  end subroutine aoresponse_input

end module aoresponse_input_m
