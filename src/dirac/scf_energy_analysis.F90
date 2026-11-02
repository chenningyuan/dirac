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

subroutine scf_energy_analysis(work, lwork)

!      use dcbham_mod, only: INTV1C, ONECAP, SMLV1C
!      use dcbgen_mod, only: DOLVC
!      use dcbbas_mod, only: N2BBASX, N2BBASXQ
      use priunit_mod, only: LUPRI
!      use dgroup_mod, only: NZ
!      use dcbdhf_mod, only: AOC, DF, DHFERG, E1PART, E2PART, ELERGY, IFCKOP, IHRMOP, INTFLG, ISYMOP, LUDENS, LUFCK1, MFMAT, MXOPEN, NFMAT, NOPEN
!      use frame_mod, only: POTNUC
  use dirac_cfg

#include "implicit.h"

  real(8), parameter :: d0 = 0.0d0
  real(8), parameter :: d2 = 2.0d0

#include "frame.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbgen.h"

  logical :: density_exists
  logical :: lvcon
  character(len=6) :: fmt
  character(len=6) :: mxform
  integer, intent(in) :: lwork
  integer :: ifckop_x(mfmat)
  integer :: iopen
  integer :: ishell
  integer :: itop
  integer :: iz
  integer :: intflg_class, iclass
  integer :: jprint
  real(8) :: basis_elergy
  real(8) :: dhferg_out
  real(8) :: e(4, 0:mxopen), ec(4, 0:mxopen), ex(4, 0:mxopen)
  real(8) :: elergy_out
  real(8), intent(inout) :: work(lwork)
  real(8), allocatable, target :: dmat_buf(:)
  real(8), allocatable, target :: f1_buf(:)
  real(8), allocatable, target :: f2_buf(:)
  real(8), allocatable, target :: f2x_buf(:)
  real(8), pointer :: dmat(:, :, :)
  real(8), pointer :: f1(:, :)
  real(8), pointer :: f2(:, :, :)
  real(8), pointer :: f2x(:, :, :)
  character(len=2) :: class_label(4)

  inquire (file='DFDENS', exist=density_exists)
  if (.not. density_exists) return

  call qenter('SCF_ENERGY_ANALYSIS')

  class_label = (/'LL', 'SL', 'SS', 'GT'/)

  lvcon = dolvc .and. ((intflg == 3 .or. intflg == 11) .or. (onecap .and. intv1c == 2))
  jprint = -1
  call dzero(e, size(e))
  e1part = d0
  e2part = d0

  allocate (f1_buf(n2bbasxq), f2_buf(n2bbasxq * nfmat), f2x_buf(n2bbasxq * nfmat), dmat_buf(n2bbasxq * nfmat))
  f1(1:n2bbasx, 1:nz) => f1_buf
  f2(1:n2bbasx, 1:nz, 1:nfmat) => f2_buf
  f2x(1:n2bbasx, 1:nz, 1:nfmat) => f2x_buf
  dmat(1:n2bbasx, 1:nz, 1:nfmat) => dmat_buf

  call opnfil(ludens, 'DFDENS', 'OLD', 'SCF_ENERGY_ANALYSIS')
  call readns(ludens, dmat_buf)
  close (ludens, status='keep')

  call twofck(isymop, ihrmop, ifckop, f2_buf, dmat_buf, nfmat, intflg, jprint, work, lwork)

  call opnfil(lufck1, 'DFFCK1', 'OLD', 'SCF_ENERGY_ANALYSIS')
  call reafck(lufck1, f1_buf, .false., 1)
  close (lufck1, status='keep')
  if (onecap .and. .not. smlv1c) then
    call catcorr(f1_buf, f2_buf, dmat_buf, work, lwork, jprint)
  end if

  if (aoc) then
    itop = nopen
  else
    itop = 0
  end if

  ! Evaluate one-electron energy
  do ishell = 0, itop
    dfrac1 = d2 * df(ishell)
    do iz = 1, nz
      e(1, ishell) = e(1, ishell) + dfrac1 * ddot(n2bbasx, dmat(1, iz, ishell + 1), 1, f1(1, iz), 1)
    end do
  end do

  ! Evaluate two-electron energy with the full Fock matrix
  call scf_energy_evaluation(e(2, 0:mxopen), dmat, f2, itop)

  ! Recompute Fock matrix to analyse Coulomb energy
  ifckop_x(1:nfmat) = 2
  ec = 0.d0
  e(3,0:itop) = 0.d0
  intflg_class = 1
  do iclass = 1, 4
    if (iand(intflg, intflg_class) /= 0) then
      call twofck(isymop, ihrmop, ifckop_x, f2x_buf, dmat_buf, nfmat, intflg_class, jprint, work, lwork)
      call scf_energy_evaluation(ec(iclass, 0:mxopen), dmat, f2x, itop)
      e(3,0:itop) = e(3,0:itop) + ec(iclass, 0:itop)
    end if
    intflg_class = intflg_class * 2 ! follow the bitpack intflg structure
  end do

  ifckop_x(1:nfmat) = 3
  ex = 0.d0
  e(4,0:itop) = 0.d0
  intflg_class = 1
  do iclass = 1, 4
    if (iand(intflg, intflg_class) /= 0) then
      call twofck(isymop, ihrmop, ifckop_x, f2x_buf, dmat_buf, nfmat, intflg_class, jprint, work, lwork)
      call scf_energy_evaluation(ex(iclass, 0:mxopen), dmat, f2x, itop)
      e(4,0:itop) = e(4,0:itop) + ex(iclass, 0:itop)
    end if
    intflg_class = intflg_class * 2 ! follow the bitpack intflg structure
  end do

  do iopen = 0, itop
    e1part = e1part + e(1, iopen)
    e2part = e2part + e(2, iopen)
  end do

  basis_elergy = e1part + e2part

  call header('Detailed energy analysis', -1)
  ergmax = max(abs(e1part), abs(e2part)) * 10
  fmt = mxform(ergmax, 20)
  call scf_energy_extensions(dmat_buf, work, lwork, jprint, lvcon, fmt, basis_elergy, elergy_out, dhferg_out)
  elergy = elergy_out
  dhferg = dhferg_out

  write (lupri, '(A,' // fmt // ')') 'Nuclear potential energy          :  ', potnuc
  write (lupri, '(A,' // fmt // ')') '* 1-electron energy               :  ', e1part
  write (lupri, '(A,' // fmt // ')') '* 2-electron energy               :  ', e2part
  write (lupri, '(A,' // fmt // ')') '** Coulomb repulsion energy       :  ', sum(e(3, 0:itop))
  intflg_class = 1
  do iclass = 1, 4
   if (iand(intflg, intflg_class) /= 0) then
      write (lupri, '(A,' // fmt // ')') '*** '//class_label(iclass)//' contribution               :  ', &
      sum(ec(iclass, 0:itop))
   end if
   intflg_class = intflg_class * 2
  end do
  write (lupri, '(A,' // fmt // ')') '** Exchange energy                :  ', sum(e(4, 0:itop))
  intflg_class = 1
  do iclass = 1, 4
   if (iand(intflg, intflg_class) /= 0) then
      write (lupri, '(A,' // fmt // ')') '*** '//class_label(iclass)//' contribution               :  ', &
      sum(ex(iclass, 0:itop))
   end if
   intflg_class = intflg_class * 2
  end do
  write (lupri, '(A,' // fmt // ')') 'Sum of all electronic terms       :  ', elergy
  write (lupri, '(A,' // fmt // ')') 'Molecular energy                  :  ', dhferg

  if (itop > 0) then
    write (lupri, '(/A)') 'Breakdown of electronic energy in terms of shell contributions: '
    do iopen = 0, itop
      write (lupri, '(A,I2)') 'Shell ', iopen
      write (lupri, '(A,' // fmt // ')') '1e     :', e(1, iopen)
      write (lupri, '(A,' // fmt // ')') '2e(C+E):', e(2, iopen)
      write (lupri, '(A,' // fmt // ')') '2e(C)  :', e(3, iopen)
      write (lupri, '(A,' // fmt // ')') '2e(E)  :', e(4, iopen)
    end do
  end if

  deallocate (f1_buf, f2_buf, f2x_buf, dmat_buf)
  call qexit('SCF_ENERGY_ANALYSIS')

end subroutine scf_energy_analysis

subroutine scf_energy_evaluation(e_shell, dmat, fmat, itop)

!      use dcbdhf_mod, only: DA, DF, MXOPEN, NFMAT
!      use dcbbas_mod, only: N2BBASX
!      use dgroup_mod, only: NZ
#include "implicit.h"

  real(8), parameter :: d2 = 2.0d0
  real(8), parameter :: dp25 = 0.25d0

#include "dcbbas.h"
#include "dgroup.h"
#include "dcbdhf.h"

  integer, intent(in) :: itop
  integer :: ishell
  integer :: iz
  integer :: jshell
  real(8), intent(inout) :: e_shell(0:mxopen)
  real(8), intent(in) :: dmat(n2bbasx, nz, nfmat)
  real(8), intent(in) :: fmat(n2bbasx, nz, nfmat)
  real(8) :: dfrac2
  real(8) :: etemp

  do ishell = 0, itop
    dfrac2 = dp25 * (d2 * df(ishell)) * (d2 * df(ishell)) * da(ishell)
    do iz = 1, nz
      e_shell(ishell) = e_shell(ishell) + dfrac2 * ddot(n2bbasx, dmat(1, iz, ishell + 1), 1, &
                                                        fmat(1, iz, ishell + 1), 1)
    end do
    do jshell = 0, ishell - 1
      dfrac2 = df(ishell) * df(jshell)
      do iz = 1, nz
        etemp = dfrac2 * ddot(n2bbasx, dmat(1, iz, jshell + 1), 1, fmat(1, iz, ishell + 1), 1)
        e_shell(ishell) = e_shell(ishell) + etemp
        e_shell(jshell) = e_shell(jshell) + etemp
      end do
    end do
  end do

end subroutine scf_energy_evaluation

subroutine scf_energy_extensions(dmat_buf, work, lwork, jprint, lvcon, fmt, basis_elergy, elergy_out, dhferg_out)

!      use dcbham_mod, only: IPRSOL, SOLVEN
!      use dcbbas_mod, only: N2BBASX, N2BBASXQ, NTBAS
!      use dgroup_mod, only: NZ
      use priunit_mod, only: LUPRI
!      use frame_mod, only: CORRLV, POTNUC
!      use dcbdhf_mod, only: DA, DF, ESOLVE, ESOLVN, NFMAT
  use dirac_cfg
  use fde_mod
#ifdef MOD_XAMFI
  use xamfi_global_parameters, only: aoomod, aooeamf
  use xamfi_utils, only: xamfi_get_energy_correction
#endif
#ifdef HAS_PCMSOLVER
  use pcm_scf
#endif
#ifdef HAS_PELIB
  use pelib_dirac_interface, only: use_pelib, pelib_ifc_energy, pelib_ifc_fock
#endif

#include "implicit.h"

  real(8), parameter :: d0 = 0.0d0
  real(8), parameter :: d2 = 2.0d0
  real(8), parameter :: dp25 = 0.25d0

#include "frame.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbham.h"
!#include "dcbgen.h"

  integer, intent(in) :: jprint
  integer, intent(in) :: lwork
  character(len=*), intent(in) :: fmt
  logical, intent(in) :: lvcon
  real(8), intent(in) :: dmat_buf(*)
  real(8), intent(in) :: basis_elergy
  real(8), intent(inout) :: work(lwork)
  real(8), intent(out) :: dhferg_out
  real(8), intent(out) :: elergy_out
  integer :: nf
  logical :: xamfi_corrections
  real(8) :: fde_energy_shift
  real(8) :: lvcorr_energy
  real(8) :: pcm_energy
  real(8) :: pe_el_nrg
  real(8) :: pe_tot_nrg
  real(8) :: solvent_energy
  real(8) :: xamfi_energy
  real(8) :: xc_mat_energy
  real(8) :: xc_total_energy
  real(8), allocatable :: ext_fmat(:)
#ifdef HAS_PELIB
  real(8), allocatable :: pe_fmat(:)
#endif
#ifdef MOD_XAMFI
  logical :: xamfi_vxc_
  real(8) :: xamfi_vxc_delta
  real(8), allocatable :: dmat_xamfi(:, :, :)
  real(8), allocatable :: xamfi_fmat(:)
#endif
  type(fde_import) :: itmp

  solvent_energy = d0
  pcm_energy = d0
  pe_el_nrg = d0
  pe_tot_nrg = d0
  xc_mat_energy = d0
  xc_total_energy = d0
  fde_energy_shift = d0
  lvcorr_energy = d0
  xamfi_energy = d0
  xamfi_corrections = .false.

  if (solven) then
    nf = 1
    allocate (ext_fmat(n2bbasxq * nfmat))
    call dzero(ext_fmat, size(ext_fmat))
    call solfck(ext_fmat, dmat_buf, nf, esolve, esolvn, work, lwork, iprsol)
    solvent_energy = esolve + esolvn
    deallocate (ext_fmat)
  end if

#ifdef HAS_PELIB
  if (use_pelib()) then
    allocate (pe_fmat(n2bbasx))
    pe_fmat = d0
    call pelib_ifc_fock(dmat_buf, pe_fmat, pe_tot_nrg, pe_el_nrg)
    if (jprint >= 1) call pelib_ifc_energy(dmat_buf)
    deallocate (pe_fmat)
  end if
#endif

#ifdef MOD_XAMFI
  inquire (file='XAMFI-ss-soc-contributions', exist=xamfi_corrections)
  if (xamfi_corrections .and. aoomod) then
    allocate (dmat_xamfi(n2bbasx, nz, nfmat))
    dmat_xamfi = reshape(dmat_buf(1:n2bbasxq * nfmat), shape(dmat_xamfi))
    allocate (xamfi_fmat(ntbas(0) ** 2 * nz))
    open (99, file='XAMFI-ss-soc-contributions', status='old', form='unformatted', &
          access='sequential', action='readwrite', position='rewind')
    read (99) xamfi_fmat(1:ntbas(0) ** 2 * nz)
    close (99, status='keep')
    xamfi_energy = xamfi_get_energy_correction(dmat_xamfi, xamfi_fmat, dp25 * d2 * df(0) * d2 * df(0) * da(0), &
                                               0, 0, nz, n2bbasx, nfmat)
    deallocate (dmat_xamfi)
    deallocate (xamfi_fmat)
  end if
#endif

  if (dirac_cfg_dft_calculation) then
#ifdef MOD_XAMFI
    inquire (file='XAMFI-Vxc-contributions', exist=xamfi_vxc_)
    if (xamfi_vxc_) then
      allocate (xamfi_fmat(n2bbasxq))
      open (99, file='XAMFI-Vxc-contributions', status='old', form='unformatted', &
            access='sequential', action='readwrite', position='rewind')
      read (99) xamfi_fmat(1:n2bbasxq)
      close (99, status='keep')
      xamfi_vxc_delta = ddot(n2bbasxq, dmat_buf, 1, xamfi_fmat, 1)
      call xcint_set_xc_mat_energy_delta(xamfi_vxc_delta)
      deallocate (xamfi_fmat)
    end if
#endif
    xc_total_energy = xcint_get_xc_energy()
    xc_mat_energy = xcint_get_xc_mat_energy()
    if (jprint >= 5) print *, 'ERGCAL: xcint_get_xc_energy()=', xc_total_energy
  end if

#ifdef HAS_PCMSOLVER
  if (dirac_cfg_pcm) then
    call pcm_energy_driver(dmat_buf, pcm_energy, work, lwork)
  end if
#endif

  if (lvcon) then
    call lvcorr(dmat_buf, work, lwork, jprint)
    lvcorr_energy = corrlv
  end if

  if (dirac_cfg_fde) then
    call fde_get_import_info(itmp)
    if (itmp%im_update_vemb) then
      fde_energy_shift = fde_get_embpot_energy()
    else if (itmp%im_vemb .and. itmp%im_remove_vemb_e_dcount) then
      fde_energy_shift = fde_dirac_embpot_static_energy()
    end if
  end if

  elergy_out = basis_elergy + xamfi_energy - xc_mat_energy + xc_total_energy - 0.5d0 * pe_el_nrg - fde_energy_shift
  dhferg_out = potnuc + elergy_out + solvent_energy + pcm_energy + pe_tot_nrg - 0.5d0 * pe_el_nrg + lvcorr_energy

  if (lvcon) write (lupri, '(A,' // fmt // ')') 'SS Coulombic correction           :  ', lvcorr_energy
  if (solven) write (lupri, '(A,' // fmt // ')') 'Solvent energy                    :  ', solvent_energy
#ifdef HAS_PCMSOLVER
  if (dirac_cfg_pcm) write (lupri, '(A,' // fmt // ')') 'Polarization energy               :  ', pcm_energy
#endif
#ifdef HAS_PELIB
  if (pe_tot_nrg /= d0 .or. pe_el_nrg /= d0) write (lupri, '(A,' // fmt // ')') 'Embedding energy                  :  ', &
                                                                                  pe_tot_nrg
#endif
#ifdef MOD_XAMFI
  if (xamfi_corrections) then
    if (aooeamf) then
      write (lupri, '(A,' // fmt // ')') 'eamf contr. to 2e-energy :  ', xamfi_energy
    else
      write (lupri, '(A,' // fmt // ')') 'amf contr. to 2e-energy  :  ', xamfi_energy
    end if
  end if
#endif

end subroutine scf_energy_extensions
