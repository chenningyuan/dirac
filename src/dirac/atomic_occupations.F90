module atomic_occupations
  implicit none
  private

  integer, parameter :: max_open_shells = 2
  integer, parameter :: max_kappa = 7

  type :: atomic_occupation_t
    integer :: atomic_number = 0
    integer :: charge = 0
    character(len=3) :: symbol = '   '
    integer :: closed(2) = 0
    integer :: nopen = 0
    integer :: open_electrons(max_open_shells) = 0
    integer :: open_shells(2,max_open_shells) = 0
    integer :: nkappa = 0
    integer :: kappa(max_kappa) = 0
    integer :: closed_kappa(max_kappa) = 0
    integer :: open_kappa(max_kappa,max_open_shells) = 0
  end type atomic_occupation_t

  type(atomic_occupation_t), parameter :: neutral_atoms(120) = [ &
    atomic_occupation_t(  1, 0, 'H ', [ 0, 0], 1, [1,0], reshape([2,0,0,0], [2,2])), &
    atomic_occupation_t(  2, 0, 'He', [ 2, 0], 0), &
    atomic_occupation_t(  3, 0, 'Li', [ 2, 0], 1, [1,0], reshape([2,0,0,0], [2,2])), &
    atomic_occupation_t(  4, 0, 'Be', [ 4, 0], 0), &
    atomic_occupation_t(  5, 0, 'B ', [ 4, 0], 1, [1,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t(  6, 0, 'C ', [ 4, 0], 1, [2,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t(  7, 0, 'N ', [ 4, 0], 1, [3,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t(  8, 0, 'O ', [ 4, 0], 1, [4,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t(  9, 0, 'F ', [ 4, 0], 1, [5,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 10, 0, 'Ne', [ 4, 6], 0), &
    atomic_occupation_t( 11, 0, 'Na', [ 4, 6], 1, [1,0], reshape([2,0,0,0], [2,2])), &
    atomic_occupation_t( 12, 0, 'Mg', [ 6, 6], 0), &
    atomic_occupation_t( 13, 0, 'Al', [ 6, 6], 1, [1,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 14, 0, 'Si', [ 6, 6], 1, [2,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 15, 0, 'P ', [ 6, 6], 1, [3,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 16, 0, 'S ', [ 6, 6], 1, [4,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 17, 0, 'Cl', [ 6, 6], 1, [5,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 18, 0, 'Ar', [ 6,12], 0), &
    atomic_occupation_t( 19, 0, 'K ', [ 6,12], 1, [1,0], reshape([2,0,0,0], [2,2])), &
    atomic_occupation_t( 20, 0, 'Ca', [ 8,12], 0), &
    atomic_occupation_t( 21, 0, 'Sc', [ 8,12], 1, [1,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 22, 0, 'Ti', [ 8,12], 1, [2,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 23, 0, 'V ', [ 8,12], 1, [3,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 24, 0, 'Cr', [ 6,12], 2, [5,1], reshape([10,0,2,0], [2,2])), &
    atomic_occupation_t( 25, 0, 'Mn', [ 8,12], 1, [5,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 26, 0, 'Fe', [ 8,12], 1, [6,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 27, 0, 'Co', [ 8,12], 1, [7,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 28, 0, 'Ni', [ 8,12], 1, [8,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 29, 0, 'Cu', [ 8,12], 1, [9,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 30, 0, 'Zn', [18,12], 0), &
    atomic_occupation_t( 31, 0, 'Ga', [18,12], 1, [1,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 32, 0, 'Ge', [18,12], 1, [2,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 33, 0, 'As', [18,12], 1, [3,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 34, 0, 'Se', [18,12], 1, [4,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 35, 0, 'Br', [18,12], 1, [5,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 36, 0, 'Kr', [18,18], 0), &
    atomic_occupation_t( 37, 0, 'Rb', [18,18], 1, [1,0], reshape([2,0,0,0], [2,2])), &
    atomic_occupation_t( 38, 0, 'Sr', [20,18], 0), &
    atomic_occupation_t( 39, 0, 'Y ', [20,18], 1, [1,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 40, 0, 'Zr', [20,18], 1, [2,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 41, 0, 'Nb', [18,18], 2, [4,1], reshape([10,0,2,0], [2,2])), &
    atomic_occupation_t( 42, 0, 'Mo', [18,18], 2, [5,1], reshape([10,0,2,0], [2,2])), &
    atomic_occupation_t( 43, 0, 'Tc', [18,18], 2, [6,1], reshape([10,0,2,0], [2,2])), &
    atomic_occupation_t( 44, 0, 'Ru', [18,18], 2, [7,1], reshape([10,0,2,0], [2,2])), &
    atomic_occupation_t( 45, 0, 'Rh', [18,18], 2, [8,1], reshape([10,0,2,0], [2,2])), &
    atomic_occupation_t( 46, 0, 'Pd', [28,18], 0), &
    atomic_occupation_t( 47, 0, 'Ag', [28,18], 1, [1,0], reshape([2,0,0,0], [2,2])), &
    atomic_occupation_t( 48, 0, 'Cd', [30,18], 0), &
    atomic_occupation_t( 49, 0, 'In', [30,18], 1, [1,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 50, 0, 'Sn', [30,18], 1, [2,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 51, 0, 'Sb', [30,18], 1, [3,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 52, 0, 'Te', [30,18], 1, [4,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 53, 0, 'I ', [30,18], 1, [5,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 54, 0, 'Xe', [30,24], 0), &
    atomic_occupation_t( 55, 0, 'Cs', [30,24], 1, [1,0], reshape([2,0,0,0], [2,2])), &
    atomic_occupation_t( 56, 0, 'Ba', [32,24], 0), &
    atomic_occupation_t( 57, 0, 'La', [32,24], 1, [1,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 58, 0, 'Ce', [32,24], 2, [1,1], reshape([0,14,10,0], [2,2])), &
    atomic_occupation_t( 59, 0, 'Pr', [32,24], 1, [3,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 60, 0, 'Nd', [32,24], 1, [4,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 61, 0, 'Pm', [32,24], 1, [5,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 62, 0, 'Sm', [32,24], 1, [6,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 63, 0, 'Eu', [32,24], 1, [7,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 64, 0, 'Gd', [32,24], 1, [8,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 65, 0, 'Tb', [32,24], 1, [9,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 66, 0, 'Dy', [32,24], 1, [10,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 67, 0, 'Ho', [32,24], 1, [11,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 68, 0, 'Er', [32,24], 1, [12,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 69, 0, 'Tm', [32,24], 1, [13,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 70, 0, 'Yb', [32,38], 0), &
    atomic_occupation_t( 71, 0, 'Lu', [32,38], 1, [1,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 72, 0, 'Hf', [32,38], 1, [2,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 73, 0, 'Ta', [32,38], 1, [3,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 74, 0, 'W ', [32,38], 1, [4,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 75, 0, 'Re', [32,38], 1, [5,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 76, 0, 'Os', [32,38], 1, [6,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 77, 0, 'Ir', [32,38], 1, [7,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 78, 0, 'Pt', [32,38], 1, [8,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 79, 0, 'Au', [32,38], 1, [9,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 80, 0, 'Hg', [42,38], 0), &
    atomic_occupation_t( 81, 0, 'Tl', [42,38], 1, [1,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 82, 0, 'Pb', [42,38], 1, [2,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 83, 0, 'Bi', [42,38], 1, [3,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 84, 0, 'Po', [42,38], 1, [4,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 85, 0, 'At', [42,38], 1, [5,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t( 86, 0, 'Rn', [42,44], 0), &
    atomic_occupation_t( 87, 0, 'Fr', [42,44], 1, [1,0], reshape([2,0,0,0], [2,2])), &
    atomic_occupation_t( 88, 0, 'Ra', [44,44], 0), &
    atomic_occupation_t( 89, 0, 'Ac', [44,44], 1, [1,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 90, 0, 'Th', [44,44], 1, [2,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t( 91, 0, 'Pa', [44,44], 2, [2,1], reshape([0,14,10,0], [2,2])), &
    atomic_occupation_t( 92, 0, 'U ', [44,44], 2, [3,1], reshape([0,14,10,0], [2,2])), &
    atomic_occupation_t( 93, 0, 'Np', [44,44], 2, [4,1], reshape([0,14,10,0], [2,2])), &
    atomic_occupation_t( 94, 0, 'Pu', [44,44], 1, [6,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 95, 0, 'Am', [44,44], 1, [7,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 96, 0, 'Cm', [44,44], 1, [8,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 97, 0, 'Bk', [44,44], 1, [9,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 98, 0, 'Cf', [44,44], 1, [10,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t( 99, 0, 'Es', [44,44], 1, [11,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t(100, 0, 'Fm', [44,44], 1, [12,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t(101, 0, 'Md', [44,44], 1, [13,0], reshape([0,14,0,0], [2,2])), &
    atomic_occupation_t(102, 0, 'No', [44,58], 0), &
    atomic_occupation_t(103, 0, 'Lr', [44,58], 1, [1,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t(104, 0, 'Rf', [44,58], 1, [2,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t(105, 0, 'Db', [44,58], 1, [3,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t(106, 0, 'Sg', [44,58], 1, [4,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t(107, 0, 'Bh', [44,58], 1, [5,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t(108, 0, 'Hs', [44,58], 1, [6,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t(109, 0, 'Mt', [44,58], 1, [7,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t(110, 0, 'Ds', [44,58], 1, [8,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t(111, 0, 'Rg', [44,58], 1, [9,0], reshape([10,0,0,0], [2,2])), &
    atomic_occupation_t(112, 0, 'Cn', [54,58], 0), &
    atomic_occupation_t(113, 0, 'Nh', [54,58], 1, [1,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t(114, 0, 'Fl', [54,58], 1, [2,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t(115, 0, 'Mc', [54,58], 1, [3,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t(116, 0, 'Lv', [54,58], 1, [4,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t(117, 0, 'Ts', [54,58], 1, [5,0], reshape([0,6,0,0], [2,2])), &
    atomic_occupation_t(118, 0, 'Og ', [54,64], 0), &
    atomic_occupation_t(119, 0, 'Uue', [54,64], 1, [1,0], reshape([2,0,0,0], [2,2])), &
    atomic_occupation_t(120, 0, 'Ubn', [56,64], 0) &
  ]

  public :: get_atomic_occupation

contains

  subroutine get_atomic_occupation(atomic_number, nfsym, &
                                   use_kappa_selection, found, &
                                   closed, nopen, open_electrons, &
                                   open_shells, &
                                   nkappa, kappa, closed_kappa, open_kappa, &
                                   symbol)
    integer, intent(in) :: atomic_number
    integer, intent(in) :: nfsym
    logical, intent(in) :: use_kappa_selection
    logical, intent(out) :: found
    integer, allocatable, intent(out) :: closed(:)
    integer, intent(out) :: nopen
    integer, allocatable, intent(out) :: open_electrons(:)
    integer, allocatable, intent(out) :: open_shells(:,:)
    integer, optional, intent(out) :: nkappa
    integer, allocatable, optional, intent(out) :: kappa(:)
    integer, allocatable, optional, intent(out) :: closed_kappa(:)
    integer, allocatable, optional, intent(out) :: open_kappa(:,:)
    character(len=3), optional, intent(out) :: symbol

    type(atomic_occupation_t) :: atom

    found = atomic_number >= 1 .and. atomic_number <= size(neutral_atoms) &
            .and. (nfsym == 1 .or. nfsym == 2)
    if (.not. found) then
      allocate(closed(0), open_electrons(0), open_shells(0,0))
      nopen = 0
      if (present(nkappa)) nkappa = 0
      if (present(kappa)) allocate(kappa(0))
      if (present(closed_kappa)) allocate(closed_kappa(0))
      if (present(open_kappa)) allocate(open_kappa(0,0))
      if (present(symbol)) symbol = '   '
      return
    end if

    atom = neutral_atoms(atomic_number)
    if (use_kappa_selection) call set_kappa_occupation(atom)
    nopen = atom%nopen
    allocate(closed(nfsym), open_electrons(nopen))
    allocate(open_shells(nfsym,nopen))

    open_electrons = atom%open_electrons(:nopen)
    if (nfsym == 1) then
      closed(1) = sum(atom%closed)
      if (nopen > 0) &
        open_shells(1,:) = sum(atom%open_shells(:,:nopen), dim=1)
    else
      closed = atom%closed
      if (nopen > 0) open_shells = atom%open_shells(:,:nopen)
    end if

    if (present(nkappa)) nkappa = merge(atom%nkappa, 0, use_kappa_selection)
    if (present(kappa)) then
      allocate(kappa(merge(atom%nkappa, 0, use_kappa_selection)))
      if (use_kappa_selection) kappa = atom%kappa(:atom%nkappa)
    end if
    if (present(closed_kappa)) then
      allocate(closed_kappa(merge(atom%nkappa, 0, use_kappa_selection)))
      if (use_kappa_selection) &
        closed_kappa = atom%closed_kappa(:atom%nkappa)
    end if
    if (present(open_kappa)) then
      allocate(open_kappa(merge(atom%nkappa, 0, use_kappa_selection), &
                            merge(nopen, 0, use_kappa_selection)))
      if (use_kappa_selection .and. nopen > 0) &
        open_kappa = atom%open_kappa(:atom%nkappa,:nopen)
    end if

    if (present(symbol)) symbol = atom%symbol
  end subroutine get_atomic_occupation

  subroutine set_kappa_occupation(atom)
    type(atomic_occupation_t), intent(inout) :: atom

    integer :: iopen
    integer :: nd_closed
    integer :: nf_closed
    integer :: np_closed
    integer :: ns_closed

    atom%nkappa = max_kappa
    atom%kappa = [-1, 1, -2, 2, -3, 3, -4]

    ns_closed = 2*number_closed_s_shells(atom)
    nd_closed = (atom%closed(1)-ns_closed)/10
    np_closed = number_closed_p_shells(atom)
    nf_closed = (atom%closed(2)-6*np_closed)/14

    atom%closed_kappa = [ns_closed, 2*np_closed, 4*np_closed, &
                          4*nd_closed, 6*nd_closed, &
                          6*nf_closed, 8*nf_closed]

    do iopen = 1, atom%nopen
      select case (sum(atom%open_shells(:,iopen)))
      case (2)
        atom%open_kappa(:,iopen) = [2, 0, 0, 0, 0, 0, 0]
      case (6)
        atom%open_kappa(:,iopen) = [0, 2, 4, 0, 0, 0, 0]
      case (10)
        atom%open_kappa(:,iopen) = [0, 0, 0, 4, 6, 0, 0]
      case (14)
        atom%open_kappa(:,iopen) = [0, 0, 0, 0, 0, 6, 8]
      end select
    end do
  end subroutine set_kappa_occupation

  integer function number_closed_s_shells(atom)
    type(atomic_occupation_t), intent(in) :: atom

    integer :: period
    logical :: open_s_shell

    select case (atom%atomic_number)
    case (1:2)
      period = 1
    case (3:10)
      period = 2
    case (11:18)
      period = 3
    case (19:36)
      period = 4
    case (37:54)
      period = 5
    case (55:86)
      period = 6
    case (87:118)
      period = 7
    case default
      period = 8
    end select

    open_s_shell = any(sum(atom%open_shells, dim=1) == 2)
    number_closed_s_shells = period
    if (open_s_shell) number_closed_s_shells = period-1

    ! Palladium has an empty 5s shell in the configuration used here.
    if (atom%atomic_number == 46) number_closed_s_shells = 4
  end function number_closed_s_shells

  integer function number_closed_p_shells(atom)
    type(atomic_occupation_t), intent(in) :: atom

    integer :: np

    number_closed_p_shells = 0
    do np = 0, 6
      if (atom%closed(2)-6*np < 0) cycle
      if (mod(atom%closed(2)-6*np, 14) == 0) &
        number_closed_p_shells = np
    end do
  end function number_closed_p_shells

end module atomic_occupations
