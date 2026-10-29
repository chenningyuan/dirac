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

module aoint_zero
!
! Blocks of AO and SO two-electron integrals that are zero by symmetry
! (parity of the Hermite expansion), so that they need not be zeroed,
! copied or scanned in direct SO Fock matrix construction
! (C2HIN1 -> SYM2 -> FCKOU1).
!
!   aoz_use    set by CAOINT for one batch of AO integrals when the only
!              consumer is SYM2 followed by FCKOUT; zero AO blocks are
!              then flagged in ao_zero and left undefined in AOINT
!   ao_zero(i) AO block i (column of AOINT) is zero, not stored
!   soz_ok     so_zero is valid for the SO integrals of this shell
!              quadruplet (set by SYM2, reset by SYMLOP)
!   so_zero(i) SO block i (column of SOINT) is zero, not stored
!
implicit none
logical, public :: aoz_use = .false.
logical, public :: soz_ok  = .false.
logical, allocatable, public :: ao_zero(:), so_zero(:)
public :: aoz_reset, soz_alloc
contains
  subroutine aoz_reset(nblk)
    ! all nblk AO blocks are stored (none flagged zero)
    integer, intent(in) :: nblk
    if (allocated(ao_zero)) then
      if (size(ao_zero) < nblk) deallocate(ao_zero)
    end if
    if (.not. allocated(ao_zero)) allocate(ao_zero(max(nblk,64)))
    ao_zero(1:nblk) = .false.
  end subroutine
  subroutine soz_alloc(nblk)
    integer, intent(in) :: nblk
    if (allocated(so_zero)) then
      if (size(so_zero) < nblk) deallocate(so_zero)
    end if
    if (.not. allocated(so_zero)) allocate(so_zero(max(nblk,64)))
  end subroutine
end module
