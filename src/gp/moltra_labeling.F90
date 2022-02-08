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

module moltra_labeling

  integer, parameter, private ::  MAX_KRAMERS_PAIRS = 3200

  save

! Arrays needed to toggle between different indexing (TODO: figure out how to use parameters below)
  integer, public  :: moltra_isp (-MAX_KRAMERS_PAIRS:MAX_KRAMERS_PAIRS)
  integer, public  :: moltra_ikr (2*MAX_KRAMERS_PAIRS)
  integer, public  :: moltra_nkr
 
contains

      subroutine Make_Kramer_to_SpinorIndex (nfsym,nstr)
      integer ::    nfsym, nstr(nfsym), ii, jj, ifsym, i
      ii = 0
      jj = 0
      do ifsym = 1, nfsym
         do i = 1, nstr(ifsym)
            ii = ii + 1
            jj = jj + 1
            moltra_isp( ii) = jj
            moltra_isp(-ii) = jj + nstr(ifsym)
            moltra_ikr(jj)              =  ii
            moltra_ikr(jj+nstr(ifsym))  = -ii
         end do
         jj = jj + nstr(ifsym)
         if (jj > MAX_KRAMERS_PAIRS) call quit ('Update moltra_labeling')
      end do
      moltra_nkr = ii
      end subroutine

      integer function Kramer_to_SpinorIndex (ikr)
      Kramer_to_SpinorIndex =  moltra_isp(ikr)
      end function

      integer function Spinor_to_KramerIndex (isp)
      Spinor_to_KramerIndex =  moltra_ikr(isp)
      end function

end module
