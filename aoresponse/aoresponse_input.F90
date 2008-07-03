module aoresponse_input_m

  use kinds_m
  use memory_allocator

  implicit none

  type perturbation_op_t
    integer       :: pg_sym, tr_sym, nr_components, typ
    real(kind_dp) :: prefactors(2)
    character(8)  :: op_name, labels(2)
  end type

  integer, parameter            :: max_nr_op = 137

  type(perturbation_op_t), save :: perturbation_op(max_nr_op)

  logical,                 save :: aoresponse_skip_pp
  logical,                 save :: aoresponse_skip_pn

  real(kind_dp),           save :: aoresponse_threshold_response
  real(kind_dp),           save :: aoresponse_threshold_rhs
  real(kind_dp),           save :: aoresponse_freq(4)
  real(kind_dp),           save :: aoresponse_freq_process

  logical,                 save :: aoresponse_alpha
  logical,                 save :: aoresponse_beta
  logical,                 save :: aoresponse_beta_2np1
  logical,                 save :: aoresponse_gamma

  logical,                 save :: aoresponse_gamma_static
  logical,                 save :: aoresponse_dc_shg
  logical,                 save :: aoresponse_thg
  logical,                 save :: aoresponse_idri
  logical,                 save :: aoresponse_dc_kerr

contains

  function find_op_index(op_name) result(op_index)

    character(*), intent(in) :: op_name
    integer                  :: op_index, i

    op_index = -1

    do i = 1, max_nr_op
      if (perturbation_op(i)%op_name == op_name) then
        op_index = i
      end if
    end do

    if (op_index == -1) then
      call quit('find_op_index: index not found')
    end if

  end function

  subroutine add_op(op_name,            &
                    typ, nr_components, &
                    prefactors,         &
                    labels)

    character(7),  intent(in)  :: op_name
    integer,       intent(in)  :: typ, nr_components
    real(kind_dp), intent(in)  :: prefactors(nr_components)
    character(8),  intent(in)  :: labels(nr_components)
!   ----------------------------------------------------------------------------
    integer                    :: indxpr, isyxpr, itrxpr, i
    real(kind_dp), allocatable :: work(:)
    integer                    :: lwork

    lwork = 200000
    call alloc(work, lwork)

    call xprind(op_name//'         ',   &
                typ, nr_components,     &
                prefactors,             &
                labels,                 &
                indxpr, isyxpr, itrxpr, &
                work, lwork, 0)

    call dealloc(work)

    perturbation_op(indxpr)%pg_sym        = isyxpr
    perturbation_op(indxpr)%tr_sym        = itrxpr
    perturbation_op(indxpr)%nr_components = nr_components
    perturbation_op(indxpr)%typ           = typ
    perturbation_op(indxpr)%op_name       = op_name

    do i = 1, nr_components
      perturbation_op(indxpr)%prefactors(i) = prefactors(i)
      perturbation_op(indxpr)%labels(i)     = labels(i)
    end do

  end subroutine

end module
