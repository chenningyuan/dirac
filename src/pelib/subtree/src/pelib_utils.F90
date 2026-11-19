!
!   Polarizable Embedding (PE) library
!   Copyright (C) The PE library developers. See the CONTRIBUTORS file
!                 in the top-level directory of this distribution.
!
!   This file is part of the PE library.
!
!   The PE library is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as
!   published by the Free Software Foundation, either version 3 of the
!   License, or (at your option) any later version.
!
!   The PE library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the PE library. If not, see <http://www.gnu.org/licenses/>.
!
!   Contact information:
!
!   Jogvan Magnus Haugaard Olsen
!   E-mail: foeroyingur@gmail.com
!
module pelib_utils

    use pelib_precision

    implicit none

contains

!------------------------------------------------------------------------------

function elem2charge(elem) result(charge)

    use pelib_options, only: luout

    character(len=*), intent(in) :: elem

    integer(ip) :: i
    real(rp) :: charge
    character(len=2), dimension(109) :: elements

    elements = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
                & 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
                & 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',&
                & 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',&
                & 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',&
                & 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',&
                & 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',&
                & 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',&
                & 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',&
                & 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',&
                & 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt' /)

    if (elem == 'X') then
        charge = 0.0_rp
        return
    end if

    do i = 1, 109
        if (elem == trim(elements(i))) then
            charge = real(i, rp)
            exit
        else
            charge = 0.0_rp
        end if
    end do

    if (int(charge) == 0 .and. elem /= 'X') then
        write(luout, *) 'WARNING: charge not found for element: ', elem
    end if

end function elem2charge

!------------------------------------------------------------------------------

function charge2elem(charge) result(elem)

    real(rp), intent(in) :: charge

    character(len=2) :: elem
    character(len=2), dimension(109) :: elements

    elements = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
                & 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
                & 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',&
                & 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',&
                & 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',&
                & 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',&
                & 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',&
                & 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',&
                & 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',&
                & 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',&
                & 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt' /)

    if ((nint(charge) == 0) .or. (nint(charge) > 109)) then
        elem = 'X '
        return
    end if

    elem = elements(nint(charge))

end function charge2elem

!------------------------------------------------------------------------------

function charge2mass(charge) result(mass)

    use pelib_options, only: luout

    real(rp), intent(in) :: charge

    integer(ip) :: i
    real(rp) :: mass
    real(rp), dimension(109) :: masses

    masses = (/   1.007825_rp,   4.002603_rp,   7.016005_rp,   9.012183_rp,&
              &  11.009305_rp,  12.000000_rp,  14.003074_rp,  15.994915_rp,&
              &  18.998403_rp,  19.992439_rp,  22.989770_rp,  23.985045_rp,&
              &  26.981541_rp,  27.976928_rp,  30.973763_rp,  31.972072_rp,&
              &  34.968853_rp,  39.962383_rp,  38.963708_rp,  39.962591_rp,&
              &  44.955914_rp,  47.947947_rp,  50.943963_rp,  51.940510_rp,&
              &  54.938046_rp,  55.934939_rp,  58.933198_rp,  57.935347_rp,&
              &  62.929599_rp,  63.929145_rp,  68.925581_rp,  73.921179_rp,&
              &  74.921596_rp,  79.916521_rp,  78.918336_rp,  83.911506_rp,&
              &  84.911800_rp,  87.905625_rp,  88.905856_rp,  89.904708_rp,&
              &  92.906378_rp,  97.905405_rp,  97.907216_rp, 101.904348_rp,&
              & 102.905503_rp, 105.903475_rp, 106.905095_rp, 113.903361_rp,&
              & 114.903875_rp, 119.902199_rp, 120.903824_rp, 129.906229_rp,&
              & 126.904477_rp, 131.904148_rp, 132.905429_rp, 137.905232_rp,&
              & 138.906347_rp, 139.905433_rp, 140.907647_rp, 141.907719_rp,&
              & 144.912743_rp, 151.919728_rp, 152.921225_rp, 157.924019_rp,&
              & 158.925342_rp, 163.929171_rp, 164.930319_rp, 165.930290_rp,&
              & 168.934212_rp, 173.938859_rp, 174.940770_rp, 179.946546_rp,&
              & 180.947992_rp, 183.950928_rp, 186.955744_rp, 191.961467_rp,&
              & 192.962917_rp, 194.964766_rp, 196.966543_rp, 201.970617_rp,&
              & 204.974401_rp, 207.976627_rp, 208.980374_rp, 208.982404_rp,&
              & 209.987126_rp, 222.017571_rp, 223.019733_rp, 226.025403_rp,&
              & 227.027750_rp, 232.038051_rp, 231.035880_rp, 238.050785_rp,&
              & 237.048168_rp, 244.064199_rp, 243.061373_rp, 247.070347_rp,&
              & 247.070300_rp, 251.079580_rp, 252.082944_rp, 257.095099_rp,&
              & 258.098570_rp, 259.100931_rp, 260.105320_rp, 261.108690_rp,&
              & 262.113760_rp, 263.118220_rp, 262.122930_rp, 269.134100_rp,&
              & 267.138000_rp /)

    if (int(charge) == 0) then
        mass = 0.0_rp
        return
    end if

    do i = 1, 109
        if (int(charge) == i) then
            mass = masses(i)
            exit
        else
            mass = 0.0_rp
        end if
    end do

    if (int(mass) == 0 .and. int(charge) == 0) then
        write(luout, *) 'WARNING: mass not found for element with charge: ', charge
    end if

end function charge2mass

!------------------------------------------------------------------------------

subroutine chcase(string, uplo)

    character(len=*), intent(inout) :: string
    character(len=*), intent(in), optional :: uplo

    integer(ip) :: i, gap
    character(len=1) :: a, z, o_uplo

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'u'
    end if

    gap = iachar('a') - iachar('A')

    if (o_uplo == 'u' .or. o_uplo == 'U') then
        a = 'a'
        z = 'z'
    else if (o_uplo == 'l' .or. o_uplo == 'L') then
        a = 'A'
        z = 'Z'
        gap = - gap
    else
        error stop 'Unknown case specified'
    end if

    do i = 1, len_trim(string)
        if (lge(string(i:i), a) .and. lle(string(i:i), z)) then
            string(i:i) = achar(iachar(string(i:i)) - gap)
        end if
    end do

end subroutine chcase

!------------------------------------------------------------------------------

subroutine openfile(filename, lunit, stat, frmt)

    character(*), intent(in) :: filename, stat, frmt
    integer(ip), intent(out) :: lunit
    integer(ip) :: i
    logical(lp) :: lexist, lopen

    if (stat == 'old') then
        inquire(file=filename, exist=lexist)

        if (.not. lexist) then
            print *, filename, ' not found.'
            error stop 'file not found.'
        end if
    end if

    do i = 21, 99
        inquire(unit=i, opened=lopen)
        if (lopen) then
            cycle
        else
            lunit = i
            open(unit=lunit, file=filename, status=stat, form=frmt)
            exit
        end if
    end do

    return

end subroutine openfile

!------------------------------------------------------------------------------

end module pelib_utils
