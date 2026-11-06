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

module gab_screening
!
! Cauchy-Schwarz factors on the level of pairs of (symmetry independent)
! contracted orbitals, used to prune orbital pairs of a shell quadruplet
! in direct Fock matrix construction before any integrals are computed
! (TWOOD1 -> ODCSEG).
!
!   qorb(ia,ib)  sqrt( max |(ij|ij)| ) over all AOs i of orbital ia and
!                j of orbital ib (all components and symmetry images)
!   ao2orb(i)    orbital (1..KMAX) of AO i
!   qorb_igtyp   integral classes (IGTYP of GABGEN) qorb was made for
!   qorb_use     set by the Fock matrix driver around TWOINT
!   prune_thr    orbital pairs (ab) with qorb(a,b) < prune_thr(1) and
!                (cd) with qorb(c,d) < prune_thr(2) are skipped;
!                set per shell quadruplet in TWOOD1, zero otherwise
!
implicit none
logical, public :: qorb_ok   = .false.
logical, public :: qorb_use  = .false.
integer, public :: qorb_igtyp = -1
real(8), public :: prune_thr(2) = 0.0d0
real(8), allocatable, public :: qorb(:,:)
integer, allocatable, public :: ao2orb(:)
end module
