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

  use dirac_cfg
  use fde_mod
  use interface_mo_specific, only: interface_mo_write
  use num_grid_gen, only: generate_num_grid
#ifdef MOD_XAMFI
  use xamfi_global_parameters, only: aoomod, aooeamf, xamfi_energy_contributions
  use xamfi_utils, only: xamfi_get_energy_correction
#endif
#ifdef HAS_PCMSOLVER
  use pcm_scf
#endif
#ifdef HAS_PELIB
  use pelib_dirac_interface, only: use_pelib, pelib_ifc_fock, pelib_ifc_energy
#endif

#include "implicit.h"
#include "priunit.h"

  real(8), parameter :: d0 = 0.0d0
  real(8), parameter :: d2 = 2.0d0
  real(8), parameter :: dp5 = 0.5d0
  real(8), parameter :: dp25 = 0.25d0

#include "frame.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbgen.h"

  logical :: lvcon
  logical :: xamfi_corrections
  logical :: xamfi_vxc_
  character(len=2) :: m(4)
  character(len=6) :: fmt
  character(len=6) :: mxform
  integer, intent(in) :: lwork
  logical :: tobe
  logical :: saveflags(4)
  integer :: ifmat
  integer :: iopen
  integer :: ishell
  integer :: itop
  integer :: iz
  integer :: jprint=-1
  integer :: jshell
  integer :: nf
  integer :: ntask
  integer :: ifckop_c(4),ifckop_e(4)
  dimension :: e(4, 0:mxopen), work(lwork)
  integer, save :: inum = 0
#ifdef HAS_PCMSOLVER
  real(8) :: pcm_energy
#endif
  real(8) :: pe_tot_nrg
  real(8) :: pe_el_nrg
  real(8) :: xamfi_vxc_delta
  real(8), allocatable :: aoo2esssoc(:)
  real(8), allocatable, target :: dmat_buf(:)
  real(8), allocatable, target :: f1_buf(:)
  real(8), allocatable, target :: f2_buf(:), f2c_buf(:), f2e_buf(:)
  integer, allocatable :: npos_local(:)
  real(8), dimension(:), allocatable :: pe_fmat
  real(8), dimension(:), allocatable :: xamfi_fmat
  real(8), pointer :: f1(:, :)
  real(8), pointer :: dmat(:, :, :)
  real(8), pointer :: f2(:, :, :), f2c(:, :, :), f2e(:, :, :)
  type(fde_import) :: itmp
  real(8), save :: besterg = 1.0d100

  data m /'AR', 'AI', 'BR', 'BI'/

  call qenter('SCF_ENERGY_ANALYSIS')
  inum = inum + 1

  lvcon = dolvc .and. ((intflg == 3 .or. intflg == 11) .or. (onecap .and. intv1c == 2))
  call dzero(e, 4 * (mxopen + 1))
  ergbuf = dhferg
  e1part = d0
  e2part = d0
  allocate (f1_buf(n2bbasxq), f2_buf(n2bbasxq * nfmat), dmat_buf(n2bbasxq * nfmat))
  allocate (f2c_buf(n2bbasxq), f2e_buf(n2bbasxq * nfmat))
  f1(1:n2bbasx, 1:nz) => f1_buf
  f2(1:n2bbasx, 1:nz, 1:nfmat) => f2_buf
  f2c(1:n2bbasx, 1:nz, 1:nfmat) => f2c_buf
  f2e(1:n2bbasx, 1:nz, 1:nfmat) => f2e_buf
  dmat(1:n2bbasx, 1:nz, 1:nfmat) => dmat_buf

  call opnfil(ludens, 'DFDENS', 'OLD', 'SCF_ENERGY_ANALYSIS')
  call readns(ludens, dmat_buf)
  close (ludens, status='keep')

  ! Compute Coulomb and exchange separately
  call SaveTaskDistribFlags(saveflags)
  call SetTaskDistribFlags((/ .true., .true., .true., .true. /))
  call SetIntTaskArrayDimension(ntask, parcal)
  allocate (npos_local(max(1, ntask)))
  ifckop_c = 2
  ifckop_e = 3
  call twofck(isymop, ihrmop, ifckop, f2_buf, dmat_buf, nfmat, npos_local, intflg, jprint, work, lwork)
  call twofck(isymop, ihrmop, ifckop_c, f2c_buf, dmat_buf, nfmat, npos_local, intflg, jprint, work, lwork)
  call twofck(isymop, ihrmop, ifckop_e, f2e_buf, dmat_buf, nfmat, npos_local, intflg, jprint, work, lwork)
  deallocate (npos_local)
  call SetTaskDistribFlags(saveflags)

  if (dirac_cfg_fde) then
    call fde_get_import_info(itmp)
    if (itmp%im_update_vemb) then
      call interface_mo_write()
      call fde_dirac_set_nz(nz)
#ifdef VAR_MPI
      if (parcal) call dirac_parctl(fde_par)
#endif
      call fde_calculate_emb_pot_mat(ntbas(0), dmat_buf, f2_buf)
    end if
  end if

  if (dirac_cfg_dft_calculation) then
    call interface_mo_write()
    call generate_num_grid(dmat_buf)
#ifdef VAR_MPI
      if (parcal) call dirac_parctl(xcint_par)
#endif
    call xcint_potential_rks(ntbas(0), dmat_buf, f2_buf)
  end if

  if (solven) then
    nf = 1
    call solfck(f2_buf, dmat_buf, nf, esolve, esolvn, work, lwork, iprsol)
  end if

#ifdef HAS_PELIB
  if (use_pelib()) then
    call pelib_ifc_fock(dmat_buf, f2_buf, pe_tot_nrg, pe_el_nrg)
  end if
#endif

  call opnfil(lufck1, 'DFFCK1', 'OLD', 'SCF_ENERGY_ANALYSIS')
  call reafck(lufck1, f1_buf, .false., 1)
  close (lufck1, status='keep')
  if (onecap .and. .not. smlv1c) then
    call catcorr(f1_buf, f2_buf, dmat_buf, work, lwork, jprint)
  end if

  inquire (file='XAMFI-ss-soc-contributions', exist=tobe)
  if (tobe) then
    allocate (aoo2esssoc(ntbas(0) ** 2 * nz))
    open (99, file='XAMFI-ss-soc-contributions', status='old', form='unformatted', &
          access='sequential', action='readwrite', position='rewind')
    read (99) aoo2esssoc(1:ntbas(0) ** 2 * nz)
    close (99, status='keep')
    call daxpy(ntbas(0) ** 2 * nz, -1.0d0, aoo2esssoc, 1, f1_buf, 1)
    call daxpy(ntbas(0) ** 2 * nz, 1.0d0, aoo2esssoc, 1, f2_buf, 1)
    deallocate (aoo2esssoc)
  end if

  if (dirac_cfg_fde) then
    call fde_get_import_info(itmp)
  end if

#ifdef MOD_XAMFI
  inquire (file='XAMFI-ss-soc-contributions', exist=xamfi_corrections)
  if (xamfi_corrections .and. aoomod) then
    open (99, file='XAMFI-ss-soc-contributions', status='old', form='unformatted', &
          access='sequential', action='readwrite', position='rewind')
    allocate (xamfi_fmat(ntbas(0) ** 2 * nz))
    read (99) xamfi_fmat(1:ntbas(0) ** 2 * nz)
    close (99, status='keep')
  end if
#endif

  if (aoc) then
    itop = nopen
  else
    itop = 0
  end if

  do ishell = 0, itop
    dfrac1 = d2 * df(ishell)
    dfrac2 = dp25 * dfrac1 * dfrac1 * da(ishell)

    do iz = 1, nz
      e(1, ishell) = e(1, ishell) + dfrac1 * ddot(n2bbasx, dmat(1, iz, ishell + 1), 1, f1(1, iz), 1)
      e(2, ishell) = e(2, ishell) + dfrac2 * ddot(n2bbasx, dmat(1, iz, ishell + 1), 1, &
                                       f2(1, iz, ishell + 1), 1) ! total
      e(3, ishell) = e(3, ishell) + dfrac2 * ddot(n2bbasx, dmat(1, iz, ishell + 1), 1, &
                                       f2c(1, iz, ishell + 1), 1) ! Coulomb-only

      e(4, ishell) = e(4, ishell) + dfrac2 * ddot(n2bbasx, dmat(1, iz, ishell + 1), 1, &
                                       f2e(1, iz, ishell + 1), 1) ! exchange-only
    end do

    do jshell = 0, ishell - 1
      dfrac2 = df(ishell) * df(jshell)
      do iz = 1, nz
        etemp = dfrac2 * ddot(n2bbasx, dmat(1, iz, jshell + 1), 1, f2(1, iz, ishell + 1), 1)
        e(2, ishell) = e(2, ishell) + etemp
        e(2, jshell) = e(2, jshell) + etemp
        ! Coulomb only
        etemp = dfrac2 * ddot(n2bbasx, dmat(1, iz, jshell + 1), 1, f2c(1, iz, ishell + 1), 1)
        e(3, ishell) = e(3, ishell) + etemp
        e(3, jshell) = e(3, jshell) + etemp
        ! Exchange only
       etemp = dfrac2 * ddot(n2bbasx, dmat(1, iz, jshell + 1), 1, f2e(1, iz, ishell + 1), 1)
        e(4, ishell) = e(4, ishell) + etemp
        e(4, jshell) = e(4, jshell) + etemp
      end do
    end do
  end do

#ifdef MOD_XAMFI
  if (xamfi_corrections .and. aoomod) then
    dfrac2 = dp25 * d2 * df(0) * d2 * df(0) * da(0)
    xamfi_energy_contributions = xamfi_get_energy_correction(dmat, xamfi_fmat, dfrac2, 0, 0, &
                                                             nz, n2bbasx, nfmat)
    deallocate (xamfi_fmat)
  end if
#endif

  if (dirac_cfg_dft_calculation) then
#ifdef MOD_XAMFI
    inquire (file='XAMFI-Vxc-contributions', exist=xamfi_vxc_)
    if (xamfi_vxc_) then
      open (99, file='XAMFI-Vxc-contributions', status='old', form='unformatted', &
            access='sequential', action='readwrite', position='rewind')
      allocate (xamfi_fmat(n2bbasxq))
      read (99) xamfi_fmat(1:n2bbasxq)
      close (99, status='keep')
      xamfi_vxc_delta = ddot(n2bbasxq, dmat_buf, 1, xamfi_fmat, 1)
      call xcint_set_xc_mat_energy_delta(xamfi_vxc_delta)
      deallocate (xamfi_fmat)
    end if
#endif
    xcint_xc_energy = xcint_get_xc_mat_energy()
    e(2, 0) = e(2, 0) - xcint_xc_energy
  end if

  
  do iopen = 0, itop
    e1part = e1part + e(1, iopen)
    e2part = e2part + e(2, iopen)
  end do
 

  dhferg = potnuc
  if (solven) then
    e2part = e2part - esolve
    ergsol = esolve + esolvn
    dhferg = dhferg + ergsol
  end if

#ifdef HAS_PCMSOLVER
  if (dirac_cfg_pcm) then
    call pcm_energy_driver(dmat_buf, pcm_energy, work, lwork)
    dhferg = dhferg + pcm_energy
  end if
#endif

#ifdef HAS_PELIB
  if (use_pelib()) then
    allocate (pe_fmat(n2bbasx))
    e_pe = 0.0d0
    call pelib_ifc_fock(dmat_buf, pe_fmat, pe_tot_nrg, pe_el_nrg)
    if (jprint >= 1) call pelib_ifc_energy(dmat_buf)
    dhferg = dhferg + pe_tot_nrg - 0.5d0 * pe_el_nrg
    e_pe = pe_tot_nrg
    deallocate (pe_fmat)
  end if
#endif

  if (lvcon) then
    call lvcorr(dmat_buf, work, lwork, jprint)
    dhferg = dhferg + corrlv
  end if

  elergy = e1part + e2part

  if (dirac_cfg_dft_calculation) then
    xcint_elergy = xcint_get_xc_energy()
    elergy = elergy + xcint_elergy
    if (jprint >= 5) then
      print *, 'ERGCAL: xcint_get_xc_energy()=', xcint_elergy
    end if
  end if

  if (dirac_cfg_fde) then
    call fde_get_import_info(itmp)
    if (itmp%im_update_vemb) then
      elergy = elergy - fde_get_embpot_energy()
    else if (itmp%im_vemb .and. itmp%im_remove_vemb_e_dcount) then
      elergy = elergy - fde_dirac_embpot_static_energy()
    end if
  end if

  dhferg = dhferg + elergy

#ifdef HAS_PELIB
  if (use_pelib()) then
    elergy = elergy - 0.5d0 * pe_el_nrg
  end if
#endif

  CALL HEADER('Detailed energy analysis',-1)
  ergmax = max(abs(e1part), abs(e2part)) * 10
  fmt = mxform(ergmax, 20)
  write (lupri, '(A,' // fmt // ')') 'Nuclear potential energy          :  ', potnuc
  if (lvcon) write (lupri, '(A,' // fmt // ')') 'SS Coulombic correction            :  ', corrlv
  if (solven) write (lupri, '(A,' // fmt // ')') 'Solvent energy                     :  ', ergsol
#ifdef HAS_PCMSOLVER
  if (dirac_cfg_pcm) write (lupri, '(A,' // fmt // ')') 'Polarization energy                :  ', &
                                                         get_pcm_energy()
#endif
#ifdef HAS_PELIB
  if (use_pelib()) write (lupri, '(A,' // fmt // ')') 'Embedding energy                   :  ', &
                                                       pe_tot_nrg
#endif


    write (lupri, '(A,' // fmt // ')') '* 1-electron energy               :  ', e1part
    write (lupri, '(A,' // fmt // ')') '* 2-electron energy               :  ', e2part
    write (lupri, '(A,' // fmt // ')') '** Coulomb repulsion energy       :  ',sum(e(3, 0:iopen))
    write (lupri, '(A,' // fmt // ')') '** Exchange energy                :  ',sum(e(4, 0:iopen))
    write (lupri, '(A,' // fmt // ')') 'Total electronic energy           :  ', elergy
#ifdef MOD_XAMFI
    if (xamfi_corrections) then
    if (aooeamf) then
      write (lupri, '(A,' // fmt // ')') '... with eamf contr. to 2e-energy :  ', &
                                          xamfi_energy_contributions
    else
      write (lupri, '(A,' // fmt // ')') '... with amf contr. to 2e-energy  :  ', &
                                          xamfi_energy_contributions
    end if
  end if
#endif
  write (lupri, '(A,' // fmt // ')') 'Total energy                      :  ', dhferg


  if (itop >0) then
    write (lupri,'(/A)') 'Breakdown of electronic energy in terms of shell contributions: '
    do iopen = 0, itop
        write (lupri,'(A,I2)') 'Shell ',iopen
        write (lupri,'(A,' // fmt // ')') '1e     :',e(1, iopen)
        write (lupri,'(A,' // fmt // ')') '2e(C+E):',e(2, iopen)
        write (lupri,'(A,' // fmt // ')') '2e(C)  :',e(3, iopen)
        write (lupri,'(A,' // fmt // ')') '2e(E)  :',e(4, iopen)
      end do
  end if


  deallocate (f1_buf, f2_buf, dmat_buf)
  call qexit('SCF_ENERGY_ANALYSIS')

end subroutine scf_energy_analysis
