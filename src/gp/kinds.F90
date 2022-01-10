!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

!      here we can collect variable type definitions
!      i'm not so sure about the decimals and exponent ranges so please edit
!
!      radovan bast (2008-04-05)
!      inspiration: discussions with jonas juselius
!-------------------------------------------------------------------------------

module kinds_m

  implicit none

  integer, public, parameter :: kind_i2 = selected_int_kind(4)
  integer, public, parameter :: kind_i4 = selected_int_kind(8)

  integer, public, parameter :: kind_sp = selected_real_kind(6, 30)
  integer, public, parameter :: kind_dp = selected_real_kind(14, 200)

end module kinds_m
