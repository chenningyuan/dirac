module aoresponse_input_m

  use kinds_m

  implicit none

  type perturbation_op_t
    integer(kind_i4) :: op_index, sym, tr_sym, nr_components, typ
    real(kind_dp)    :: prefactor(2)
    character(8)     :: op_name, label(2)
  end type perturbation_op_t
  public perturbation_op_t

  type perturbation_field_t
    type(perturbation_op_t) :: component(3)
  end type perturbation_field_t
  public perturbation_field_t

  type(perturbation_field_t), save :: perturbation_field(2)
  public perturbation_field

  public hash_field_to_index
  public aoresponse_input

  private

contains

  function hash_field_to_index(f) result(i)
    character(*), intent(in) :: f
    integer(kind_i4)         :: i

    select case (f)
      case ('EL')
        i = 1
      case ('NOLO')
        i = 2
      case default
        call quit('panic in hash_field_to_index')
    end select

  end function hash_field_to_index

  subroutine aoresponse_input(work, lwork)
    real(kind_dp),    intent(inout) :: work(lwork)
    integer(kind_i4), intent(inout) :: lwork

    integer(kind_i4)                :: indxpr, isyxpr, itrxpr, i, j

    perturbation_field(hash_field_to_index('EL'))%component(1)%op_name = 'ex'
    perturbation_field(hash_field_to_index('EL'))%component(2)%op_name = 'ey'
    perturbation_field(hash_field_to_index('EL'))%component(3)%op_name = 'ez'

    perturbation_field(hash_field_to_index('EL'))%component(1)%nr_components = 1
    perturbation_field(hash_field_to_index('EL'))%component(2)%nr_components = 1
    perturbation_field(hash_field_to_index('EL'))%component(3)%nr_components = 1

    perturbation_field(hash_field_to_index('EL'))%component(1)%typ = 1
    perturbation_field(hash_field_to_index('EL'))%component(2)%typ = 1
    perturbation_field(hash_field_to_index('EL'))%component(3)%typ = 1

    perturbation_field(hash_field_to_index('EL'))%component(1)%prefactor(1) = 1.0d0
    perturbation_field(hash_field_to_index('EL'))%component(2)%prefactor(1) = 1.0d0
    perturbation_field(hash_field_to_index('EL'))%component(3)%prefactor(1) = 1.0d0

    perturbation_field(hash_field_to_index('EL'))%component(1)%label(1) = 'XDIPLEN'
    perturbation_field(hash_field_to_index('EL'))%component(2)%label(1) = 'YDIPLEN'
    perturbation_field(hash_field_to_index('EL'))%component(3)%label(1) = 'ZDIPLEN'

    perturbation_field(hash_field_to_index('NOLO'))%component(1)%op_name = 'bx'
    perturbation_field(hash_field_to_index('NOLO'))%component(2)%op_name = 'by'
    perturbation_field(hash_field_to_index('NOLO'))%component(3)%op_name = 'bz'

    perturbation_field(hash_field_to_index('NOLO'))%component(1)%nr_components = 2
    perturbation_field(hash_field_to_index('NOLO'))%component(2)%nr_components = 2
    perturbation_field(hash_field_to_index('NOLO'))%component(3)%nr_components = 2

    perturbation_field(hash_field_to_index('NOLO'))%component(1)%typ = 5
    perturbation_field(hash_field_to_index('NOLO'))%component(2)%typ = 6
    perturbation_field(hash_field_to_index('NOLO'))%component(3)%typ = 7

    perturbation_field(hash_field_to_index('NOLO'))%component(1)%prefactor(1) = 68.5d0
    perturbation_field(hash_field_to_index('NOLO'))%component(1)%prefactor(2) = 68.5d0
    perturbation_field(hash_field_to_index('NOLO'))%component(2)%prefactor(1) = 68.5d0
    perturbation_field(hash_field_to_index('NOLO'))%component(2)%prefactor(2) = 68.5d0
    perturbation_field(hash_field_to_index('NOLO'))%component(3)%prefactor(1) = 68.5d0
    perturbation_field(hash_field_to_index('NOLO'))%component(3)%prefactor(2) = 68.5d0

    perturbation_field(hash_field_to_index('NOLO'))%component(1)%label(1) = 'ZDIPLEN'
    perturbation_field(hash_field_to_index('NOLO'))%component(1)%label(2) = 'YDIPLEN'
    perturbation_field(hash_field_to_index('NOLO'))%component(2)%label(1) = 'XDIPLEN'
    perturbation_field(hash_field_to_index('NOLO'))%component(2)%label(2) = 'ZDIPLEN'
    perturbation_field(hash_field_to_index('NOLO'))%component(3)%label(1) = 'YDIPLEN'
    perturbation_field(hash_field_to_index('NOLO'))%component(3)%label(2) = 'XDIPLEN'

    do i = 1, 2
      do j = 1, 3

        call xprind(perturbation_field(i)%component(j)%op_name,       &
                    perturbation_field(i)%component(j)%typ,           &
                    perturbation_field(i)%component(j)%nr_components, &
                    perturbation_field(i)%component(j)%prefactor,     &
                    perturbation_field(i)%component(j)%label,         &
                    indxpr, isyxpr, itrxpr,                           &
                    work, lwork, 0)

        perturbation_field(i)%component(j)%op_index = indxpr
        perturbation_field(i)%component(j)%sym      = isyxpr
        perturbation_field(i)%component(j)%tr_sym   = itrxpr

      end do
    end do

  end subroutine aoresponse_input

end module aoresponse_input_m
