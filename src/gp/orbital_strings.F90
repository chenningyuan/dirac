module orbital_strings

! Common parsing utilities for orbital selection strings.
!
! This module is intended to become the shared home for the modern
! string2mo implementation and, eventually, for a numlst-compatible wrapper.

  implicit none
  private

  public string2mo
  public read_orbital_energies

contains

  subroutine string2mo(inp_str, nmo, mo_list, eps)

    use fermi_level, only: fermi_energy, mo_position

    character(len=*), intent(in)       :: inp_str
    integer, intent(out), optional     :: nmo
    integer, intent(out), optional     :: mo_list(:)
    real(8), intent(in), optional      :: eps(:)

    integer :: str_len, sl, p_len, check, nmo_all
    integer :: infty_len
    integer :: mo1, mo2, i, j, k, i_mo, size_mo
    character(len=20) :: sav_str
    character(len=14) :: allowed_str
    character(len=1)  :: infty_str
    character(len=1)  :: p_str, c_str, e_str
    real(8)                 :: e_min, e_max, e_dist
    real(8), allocatable    :: dist_a(:), orbital_energies(:)
    logical                 :: do_nmo, do_mo

    do_nmo = .false.
    do_mo = .false.
    if (present(nmo)) then
       do_nmo = .true.
       if (present(mo_list)) stop 'use only one output in string2mo'
       nmo = 0
    else if (present(mo_list)) then
       do_mo = .true.
       mo_list = -9898
       size_mo = size(mo_list)
    else
       stop 'no output in string2mo'
    end if

    if (present(eps)) then
       nmo_all = size(eps)
       allocate(orbital_energies(nmo_all))
       orbital_energies = eps
    else
       call read_orbital_energies(orbital_energies, nmo=nmo_all)
    end if
    if (nmo_all < 1) stop 'invalid orbital count in string2mo'

    str_len = len(inp_str)
    allowed_str = '0123456789+-HL'
    infty_str = 'o'
    p_str = '.'
    c_str = ','
    e_str = ' '
    sl = 0
    p_len = 0
    infty_len = 0
    check = 0
    sav_str = ' '

    if (inp_str(1:min(6,str_len)) == 'energy') then
       i = 7
       do k = 1, 100
          if (inp_str(i:i) == ' ') then
             i = i + 1
          else
             exit
          end if
       end do
       j = index(inp_str(i:),' ')
       j = i + j
       if (to_upper(inp_str(i:j-1)) == 'FERMI') then
          e_min = fermi_energy()
       else
          read(inp_str(i:j-1), *) e_min
       end if

       i = j
       do k = 1, 100
          if (inp_str(i:i) == ' ') then
             i = i + 1
          else
             exit
          end if
       end do
       j = index(inp_str(i:),' ')
       j = i + j
       if (to_upper(inp_str(i:j-1)) == 'FERMI') then
          e_max = fermi_energy()
       else
          read(inp_str(i:j-1), *) e_max
       end if
       i = index(inp_str(j+1:),' ') + j
       read(inp_str(j-1:i), *) e_dist
       if (e_min >= e_max) stop 'string2mo use: energy e_min, e_max, e_dist'

       allocate(dist_a(nmo_all))

       dist_a = orbital_energies(1:nmo_all) - e_min
       where(dist_a < 0.0d0) dist_a = 1000000.0d0
       j = minloc(abs(dist_a),1)
       mo1 = 1
       do i = j - 1, 1, -1
          if (abs(orbital_energies(i+1) - orbital_energies(i)) > e_dist) then
             mo1 = i + 1
             exit
          end if
       end do

       dist_a = orbital_energies(1:nmo_all) - e_max
       where(dist_a > 0.0d0) dist_a = 1000000.0d0
       j = minloc(abs(dist_a),1)
       mo2 = nmo_all
       do i = j, nmo_all - 1
          if (abs(orbital_energies(i) - orbital_energies(i+1)) > e_dist) then
             mo2 = i
             exit
          end if
       end do

       if (do_nmo) nmo = mo2 - mo1 + 1

       if (do_mo) then
          i = mo2 - mo1 + 1
          if (i /= size_mo) stop 'mo_list wrong size in string2mo'
          i_mo = 0
          do k = mo1, mo2
             i_mo = i_mo + 1
             mo_list(i_mo) = k
          end do
       end if
       deallocate(dist_a)
    else
       i_mo = 0
       do i = 1, str_len
          if (inp_str(i:i) == e_str) then
             cycle
          else if (inp_str(i:i) == p_str) then
             p_len = p_len + 1
             if (p_len == 1) then
                mo1 = orbital_position(sav_str(1:sl))
                sl = 0
             else if (p_len > 2) then
                stop 'wrong number of points in string2mo'
             end if
          else if (inp_str(i:i) == infty_str) then
             infty_len = infty_len + 1
             if (infty_len == 1) then
                sl = 1
             else if (infty_len > 2) then
                stop 'infinity sign ill defined in string2mo'
             end if
          else if (inp_str(i:i) == c_str) then
             if (p_len == 0) then
                mo1 = orbital_position(sav_str(1:sl))
                call add_orbital(mo1, nmo_all, do_nmo, nmo, do_mo, i_mo, mo_list)
                sl = 0
             else if (p_len == 2) then
                if (infty_len == 2) then
                   mo2 = nmo_all
                else
                   mo2 = orbital_position(sav_str(1:sl))
                end if
                call add_orbital_range(mo1, mo2, nmo_all, do_nmo, nmo, do_mo, i_mo, mo_list)
                p_len = 0
                sl = 0
             else
                stop 'not enough points in string2mo'
             end if
             infty_len = 0
          else
             check = 0
             do j = 1, len(allowed_str)
                if (to_upper(inp_str(i:i)) == allowed_str(j:j)) then
                   sl = sl + 1
                   if (sl > len(sav_str)) stop 'orbital token too long in string2mo'
                   sav_str(sl:sl) = to_upper(inp_str(i:i))
                   check = check + 1
                end if
             end do
             if (check /= 1) stop 'wrong character in string2mo'
          end if
       end do

       if (sl > 0) then
          if (p_len == 0) then
             mo1 = orbital_position(sav_str(1:sl))
             call add_orbital(mo1, nmo_all, do_nmo, nmo, do_mo, i_mo, mo_list)
          else if (p_len == 2) then
             if (infty_len == 2) then
                mo2 = nmo_all
             else
                mo2 = orbital_position(sav_str(1:sl))
             end if
             call add_orbital_range(mo1, mo2, nmo_all, do_nmo, nmo, do_mo, i_mo, mo_list)
          else
             stop 'not enough points in string2mo'
          end if
       end if

       if (do_mo) then
          if (i_mo /= size_mo) stop 'mo_list wrong size in string2mo'
          j = 1
          call quicksort(mo_list, j, size_mo)
          do j = 1, size_mo - 1
             if (mo_list(j) == mo_list(j+1)) stop 'mo appearing twice in string2mo'
          end do
       end if
    end if

  end subroutine string2mo

  subroutine add_orbital(mo, nmo_all, do_nmo, nmo, do_mo, i_mo, mo_list)
    integer, intent(in)             :: mo, nmo_all
    logical, intent(in)             :: do_nmo, do_mo
    integer, intent(inout)          :: i_mo
    integer, intent(inout), optional :: nmo
    integer, intent(inout), optional :: mo_list(:)

    if (mo <= nmo_all) then
       if (do_nmo .and. present(nmo)) nmo = nmo + 1
       if (do_mo) then
          i_mo = i_mo + 1
          mo_list(i_mo) = mo
       end if
    end if
  end subroutine add_orbital

  subroutine add_orbital_range(mo1, mo2, nmo_all, do_nmo, nmo, do_mo, i_mo, mo_list)
    integer, intent(inout)           :: mo1, mo2
    integer, intent(in)              :: nmo_all
    logical, intent(in)              :: do_nmo, do_mo
    integer, intent(inout)           :: i_mo
    integer, intent(inout), optional :: nmo
    integer, intent(inout), optional :: mo_list(:)

    integer :: k, tmp

    if (mo1 > mo2) then
       tmp = mo1
       mo1 = mo2
       mo2 = tmp
    end if
    if (mo1 < nmo_all) then
       if (mo2 > nmo_all) mo2 = nmo_all
       if (do_nmo .and. present(nmo)) nmo = nmo + mo2 - mo1 + 1
       if (do_mo) then
          do k = mo1, mo2
             i_mo = i_mo + 1
             mo_list(i_mo) = k
          end do
       end if
    end if
  end subroutine add_orbital_range

  subroutine read_orbital_energies(eps_unbarred, eps_barred, nmo)
    use checkpoint
    use fermi_level, only: set_mosection

    real(8), allocatable, intent(out)            :: eps_unbarred(:)
    real(8), allocatable, intent(out), optional  :: eps_barred(:)
    integer, intent(out), optional               :: nmo

    character(len=:), allocatable        :: group
    integer                             :: i, idx, n_mo(2), n_po(2), nmo_all
    integer                             :: npo_all, nto
    real(8), allocatable                :: all_eps(:)
    logical                             :: mbpt_orbitals, data_present, c1_info

    group = '/result/wavefunctions/scf/'
    call checkpoint_query('/result/wavefunctions/mbpt/mobasis/n_mo', mbpt_orbitals)
    if (mbpt_orbitals) group = '/result/wavefunctions/mbpt/'
    call set_mosection(group)

    call checkpoint_query(group//'mobasis/n_mo', data_present)
    if (.not.data_present) then
       write(*,*) ' ERROR: no orbitals present on ', group
       stop ' Fatal error in string2mo'
    end if
    call checkpoint_read(group//'mobasis/n_mo', idata=n_mo)

    call checkpoint_query(group//'mobasis/n_po', data_present)
    if (data_present) then
       call checkpoint_read(group//'mobasis/n_po', idata=n_po)
    else
       n_po = 0
    end if

    call checkpoint_query(group//'mobasis/eigenvalues_C1', c1_info)
    if (c1_info) then
       nmo_all = sum(n_mo)
       npo_all = sum(n_po)
       nto = nmo_all + npo_all
       allocate(all_eps(nto))
       call checkpoint_read(group//'mobasis/eigenvalues_C1', rdata=all_eps)
    else
       if (n_mo(2) > 0 .or. n_po(2) > 0) then
          print*, ' Not yet possible to use inversion symmetry in string2mo'
          stop ' Fatal error in string2mo'
       end if
       nmo_all = n_mo(1)
       npo_all = n_po(1)
       nto = nmo_all + npo_all
       allocate(all_eps(nto))
       call checkpoint_read(group//'mobasis/eigenvalues', rdata=all_eps)
    end if

    if (present(nmo)) nmo = nmo_all

    allocate(eps_unbarred(nmo_all))
    if (present(eps_barred)) allocate(eps_barred(nmo_all))
    do i = 1, nmo_all
       idx = npo_all + i
       eps_unbarred(i) = all_eps(idx)
       if (present(eps_barred)) then
          eps_barred(i) = all_eps(idx)
       end if
    end do
    deallocate(all_eps)
  end subroutine read_orbital_energies

  integer function orbital_position(token) result(mo)
    use fermi_level, only: mo_position

    character(len=*), intent(in)     :: token

    if (token(1:1) == 'H') then
       mo = mo_position(token)
    else if (token(1:1) == 'L') then
       mo = mo_position(token)
    else
       read(token, *) mo
    end if
  end function orbital_position

  recursive subroutine quicksort(mo_list, i_beg, i_end)
    integer, intent(inout) :: mo_list(:)
    integer, intent(in)    :: i_beg, i_end

    integer :: i, j
    integer :: ref, tmp

    i = floor((real(i_end) + real(i_beg))/2.0)
    ref = mo_list(i)
    i = i_beg
    j = i_end

    do
       do while (mo_list(i) < ref)
          i = i + 1
       end do
       do while (mo_list(j) > ref)
          j = j - 1
       end do
       if (i < j) then
          tmp = mo_list(i)
          mo_list(i) = mo_list(j)
          mo_list(j) = tmp
          i = i + 1
          j = j - 1
       else if (i == j) then
          i = i + 1
          exit
       else
          exit
       end if
    end do

    if (i_beg < j) call quicksort(mo_list, i_beg, j)
    if (i < i_end) call quicksort(mo_list, i, i_end)
  end subroutine quicksort

  pure function to_upper(str) result(out)
    character(len=*), intent(in) :: str
    character(len=len(str))      :: out

    integer :: i, c

    out = str
    do i = 1, len(str)
       c = iachar(str(i:i))
       if (c >= iachar('a') .and. c <= iachar('z')) then
          out(i:i) = achar(c - 32)
       end if
    end do
  end function to_upper

end module orbital_strings
