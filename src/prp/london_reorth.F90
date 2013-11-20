!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

module london_reorth_hb

  use memory_allocator 
  use dirac_cfg
  use dft_cfg
  use num_grid_gen
  use dirac_interface
#ifdef MOD_LAO_REARRANGED
  use london_helper  
  use london_utils   
#endif

  implicit none

  public london_reort_hb
!  public check_symmetry 
  save

  private

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "dcbprp.h"
#include "aovec.h"
#include "blocks.h"
#include "dcbfir.h"
#include "dcbham.h"
#include "iratdef.h"
#include "cbihr2.h"
#include "dcbdhf.h"
#include "shells.h"
#include "dcbxlr.h"

#include "infpar.h"



  integer, parameter :: max_nr_xyz_comp = 3
  integer            ::     nr_xyz_comp

  integer, parameter :: file_unit = 137
  integer, parameter :: tbmo_unit = 71
  logical            :: debug_me

contains

  subroutine london_reort_hb(cmo, ibeig, work, lwork)

! calculates the reorthonormalization contribution to dH(0)/dB|B=0
! when London atomic orbitals are used:
! {T^B, H(0)} = {T^B, h(0)} + {T^B, L}
! please see two last terms in eq. (31) in JCP 131, 124119 (2009)

! gosia: in connection-independent formulation, we calculate only:
! -{S^B, L}

!-----------------------------------------------------------------------
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    integer, intent(in)    :: lwork
    real(8), intent(inout) :: work(lwork)

    real(8), allocatable   :: fao(:), fmo(:), tb_dmat(:)
    integer, allocatable   :: npos_array(:)
    integer :: size_tb_dmat, size_fao, size_fmo
!-----------------------------------------------------------------------

!   set interface for f77 work array:
    call di_set_f77_work(work, lwork)

!   decide how many contributions to calculate
!   nr_xyz_comp = 3 -> all three Bx, By, Bz
!   nr_xyz_comp = 1 -> only one of them (for visual)
    if (bxlao .or. bylao .or. bzlao) then
       nr_xyz_comp = 1
    else
       nr_xyz_comp = 3
    end if

    size_tb_dmat = n2bbasxq*nr_xyz_comp
    call alloc(tb_dmat, size_tb_dmat)
    tb_dmat = 0.0d0

!   get (C+)(T^B)C + c.c.
!   gosia: if (lao_lr_rearrange) then it is (C+)(S^B)C written on tb_dmat
    call get_tb_dmat(tb_dmat, cmo, ibeig)

    size_fao = n2bbasxq*nr_xyz_comp
    call alloc(fao, size_fao)
    fao = 0.0d0

!   calculate "G" part in AO basis
    call get_gtb(fao, tb_dmat, cmo, ibeig)

    call dealloc(tb_dmat)

!   do AO-to-MO transformation on "G" part
    size_fmo = n2orbxq*nr_xyz_comp
    call alloc(fmo, size_fmo)
    fmo = 0.0d0
    call do_ao_to_mo(fmo, fao, cmo, ibeig)
    call dealloc(fao)

!   calculate "F" part, add to "G" and write result to file
    call get_ftb(fmo, cmo, ibeig)

    call dealloc(fmo)

  end subroutine

  subroutine get_ftb(fmat, cmo, ibeig)
!-----------------------------------------------------------------------
    real(8), intent(in)    :: fmat(n2orbxq, *)
    !real(8), intent(inout)    :: fmat(n2orbxq, *)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork, kfree

    real(8), allocatable   :: fmo(:)
    real(8), allocatable   :: buf(:)
    real(8), allocatable   :: tb(:)

    integer :: isym(nr_xyz_comp), iopsy(nr_xyz_comp), irep(nr_xyz_comp)
    integer :: isymop(nr_xyz_comp), ihrmop(nr_xyz_comp), ifckop(nr_xyz_comp)
    integer :: indxpr
    integer :: icomp, i, j, t, isymh, iz
    logical :: ah, ee, ep
    real(8) :: a
!-----------------------------------------------------------------------

    call alloc(fmo, n2bbasxq)
    fmo = 0.0d0
!   get unperturbed Fock matrix in MO basis:
    lwork = len_f77_work
    kfree = 1
    call di_select_wrk(work, lwork)
    call getfck(fmo,         &
                iprprp,      &
                work,        &
                kfree,       &
                lwork)

    call di_deselect_wrk(work, lwork)

#ifdef MOD_LAO_REARRANGED
    if (lao_lr_rearrange) then
      call get_m_ij(fmo)
    end if
#endif

    call alloc(tb,  n2orbxq)
    tb = 0.0d0
    call alloc(buf, n2orbxq)
    buf = 0.0d0

    call opnfil(lu1int, 'AOPROPER', 'OLD', 'get_ftb')

    open(lufckl, file='FCKLON', form='unformatted',   &
         access='direct',recl=8*n2orbxq,                 &
         status='unknown')

    open(tbmo_unit, file='TBMO', form='unformatted',   &
         access='sequential', action='write',          &
         status='unknown')

    do icomp = 1, nr_xyz_comp
      indxpr       = ipcon(icomp)
      isym(icomp)  = iprpsym(indxpr)
      irep(icomp)  = isym(icomp) - 1
      iopsy(icomp) = jbtof(irep(icomp), 1)

      debug_me = .false.
      if (debug_me) then
        write(*, *) 'get_ftb, entering fmat, icomp = ', icomp
        call prqmat(fmat(1, icomp), &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, irep(icomp)), lupri)
#ifdef MOD_LAO_REARRANGED
        call check_symmetry(fmat(1, icomp), norbt, nz)
#endif
      end if

      call dzero(tb, n2orbxq)

      lwork = len_f77_work
      kfree = 1
      call di_select_wrk(work, lwork)

!     get T^B_icomp matrix to tb (from AOPROPER file)
      call prpmat(indxpr,        &
                  iopsy(icomp),  &
                  tb,            &
                  cmo,           &
                  ibeig,         &
                  icmoq,         &
                  norb,          &
                  work,          &
                  kfree,         &
                  lwork,         &
                  iprprp)

      call di_deselect_wrk(work, lwork)

!!!TEST
!      call test_new_routines(tb)

! write to TBMO file, to be reused in visual
! ------------------------------------------
! *) how many components:
      write(tbmo_unit) nr_xyz_comp
! *) which icomp:
      write(tbmo_unit) icomp
! *) irep of T^B_icomp:
      write(tbmo_unit) (isym(icomp) - 1)
! *) time-reversal symmetry:
!      write(tbmo_unit) 
! *) reorthonormalization over both (p+e) shells or only over e shells:
      write(tbmo_unit) doeprn 
! *) T^B_icomp matrix in MO basis:
      write(tbmo_unit) tb

      debug_me = .false.
      if (debug_me) then
         write(lupri, *) 'T^B in MO basis,icomp = ', icomp
         write(lupri, *) 'isym,iopsy = ', isym(icomp),iopsy(icomp)
         call prqmat(tb, &
                     norb, norb, norb, norb, &
                     nz, ipqtoq(1, irep(icomp)), lupri)
      end if

! in lao_lr_rearrange take only 'ai' blocks:
#ifdef MOD_LAO_REARRANGED
      if (lao_lr_rearrange) then
        ee = .true.
        ep = .true.
        if (xlr_skipep .or. .not. doeprn) ep = .false.
        if (xlr_skipee) ee = .false.
        call get_m_ai(tb, ee, ep)
      end if
#endif

! get final result to buf

      call dzero(buf, n2orbxq)

!     it is = 1, because unperturbed Fock matrix is symmetric (its real part)
!     it is set in getfck
      isymh = 1
!     ah = .true. because T is antihermitian
      ah = .true.

!     do 1-index transformation
      call ttra(doeprn,          &
                ah,              &
                isym(icomp),     &
                isymh,           &
                tb,              &
                fmo,             &
                buf,             &
                iprprp)

      debug_me = .false.
      if (debug_me) then
        write(*, *) 'buf in get_ftb after ttra, icomp = ', icomp
        call prqmat(buf, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, irep(icomp)), lupri)
      end if

!     add result to fmat
      call daxpy(n2orbxq,        &
                 1.0d0,          &
                 buf,            &
                 1,              & 
                 fmat(1, icomp), & 
                 1)

!     write to FCKLON file:
      call wrtdac(lufckl, norbt*norbt*nz,   &
                  fmat(1, icomp), icomp)

      debug_me = .false.
      if (debug_me) then
        write(*, *) 'fmat final in london_reorth, icomp = ', icomp
        call prqmat(fmat(1, icomp), &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, irep(icomp)), lupri)

#ifdef MOD_LAO_REARRANGED
        call check_symmetry(fmat, norbt, nz)
#endif
      end if

    end do


    close(tbmo_unit,   status='keep')
    close(lufckl, status='keep')
    close(lu1int,   status='keep')

    call dealloc(tb)
    call dealloc(buf)
    call dealloc(fmo)

  end subroutine


  subroutine get_gtb(fao, tb_dmat, cmo, ibeig)
    use xcint_main
!-----------------------------------------------------------------------
    real(8), intent(out)   :: fao(*)
    real(8), intent(in)    :: tb_dmat(*)

    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork

    integer, allocatable   :: npos_array(:)
    integer                :: npos
    
    real(8), allocatable   :: dmat(:)
    integer :: isym(nr_xyz_comp), iopsy(nr_xyz_comp), irep(nr_xyz_comp)
    integer :: isymop(nr_xyz_comp), ihrmop(nr_xyz_comp), ifckop(nr_xyz_comp)
    integer :: indxpr
    integer :: icomp, ioff
!-----------------------------------------------------------------------
#include "dcbgen.h"

!   set for twofck:
    do icomp = 1, nr_xyz_comp
      indxpr        = ipcon(icomp)
      isym(icomp)   = iprpsym(indxpr)
      iopsy(icomp)  = jbtof(isym(icomp) - 1, 1)
      irep(icomp)   = isym(icomp) - 1
      isymop(icomp) = isym(icomp)
      ihrmop(icomp) = -1
      ifckop(icomp) = jbtof(irep(icomp), 1)
    end do

!gosia fixme parallel run
    call my_get_npos(npos)
    call alloc(npos_array, npos)
    npos_array(1:npos) = 0

    lwork = len_f77_work
    call di_select_wrk(work, lwork)

    call twofck(isymop,      &
                ihrmop,      &
                ifckop,      &
                fao,         &
                tb_dmat,     &
                nr_xyz_comp, &
                npos_array,  &
                intnmr,      &
                iprprp,      &
                work,        &
                lwork)

    call di_deselect_wrk(work, lwork)

    debug_me = .false.
    if (debug_me) then
      do icomp = 1, nr_xyz_comp
        write(*, *) 'fao in get_gtb, icomp = ', icomp
        ioff = (icomp - 1)*nz*ntbas(0)*ntbas(0) + 1
        call prqmat(fao(ioff), &
                        ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
                        nz, ipqtoq(1, irep(icomp)), lupri)
      end do
    end if

    if (allocated(npos_array)) then
      call dealloc(npos_array)
    end if

!       construct kohn-sham contribution
!       ================================
!   in dft Bx/By/Bz contributions are added on top of each other!
    if (dirac_cfg_dft_calculation .and. .not. dft_cfg_no_sdft) then
    !if (dirac_cfg_dft_calculation) then

      debug_me = .false.
      if (debug_me) then
        do icomp = 1, nr_xyz_comp
          ioff = (icomp - 1)*nz*ntbas(0)**2 + 1
          write(lupri, *) 'london_reorth: fmat before integrate_xc, icomp = ', icomp
          call prqmat(fao(ioff), &
                      ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
                      nz, ipqtoq(1, isymax(icomp, 2)), lupri)
          write(lupri, *) 'london_reorth: tb_dmat before integrate_xc, icomp = ', icomp
          call prqmat(tb_dmat(ioff), &
                      ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
                      nz, ipqtoq(1, isymax(icomp, 2)), lupri)
        end do
      end if

      call alloc(dmat, n2bbasxq)
      dmat = 0.0d0
      call genden(dmat, cmo, 1, iprprp)

      call generate_num_grid(dmat)
#ifdef VAR_MPI
      if (parcal) call dirac_parctl(4)
#endif
      call integrate_xc(xc_mat_dim          = ntbas(0), &
                         xc_nz               = nz,       &
                         xc_dmat_0           = dmat,        &
                         xc_nr_dmat          = nr_xyz_comp, &
                         xc_nr_fmat          = nr_xyz_comp, &
                         xc_dmat             = tb_dmat,     &
                         xc_fmat             = fao,         &
                         xc_fmat_pg_sym      = isymop,      &
                         xc_dmat_pg_sym      = isymop,      &
                         xc_dmat_ih_sym      = ihrmop,      &
                         xc_do_london_rhs_ro = .true.)       

      call dealloc(dmat)

      if (debug_me) then
        do icomp = 1, nr_xyz_comp
          ioff = (icomp - 1)*nz*ntbas(0)**2 + 1
          write(lupri, *) 'london_reorth: fmat after integrate_xc, icomp = ', icomp
          call prqmat(fao(ioff), &
                      ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
                      nz, ipqtoq(1, isymax(icomp, 2)), lupri)
        end do
      end if

    end if

  end subroutine get_gtb


      SUBROUTINE my_GET_NPOS(NPOS)
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "cbihr2.h"
#include "blocks.h"
#include "dcbfir.h"
      integer npos
      call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. ,.TRUE. /))
      call SetIntTaskArrayDimension(NPOS,PARCAL)
      end subroutine



  subroutine do_ao_to_mo(fmo, fao, cmo, ibeig)
!-----------------------------------------------------------------------
    real(8), intent(inout) :: fmo(n2orbxq , *)
    real(8), intent(in)    :: fao(n2bbasxq, *)

    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    integer :: itim(nr_xyz_comp), isym(nr_xyz_comp), iopsy(nr_xyz_comp), irep(nr_xyz_comp)
    integer :: isymop(nr_xyz_comp), ihrmop(nr_xyz_comp), ifckop(nr_xyz_comp)
    integer :: indxpr
    integer :: icomp, i, j
!-----------------------------------------------------------------------

    do icomp = 1, nr_xyz_comp

      indxpr        = ipcon(icomp)
      isym(icomp)   = iprpsym(indxpr)
      irep(icomp)   = isym(icomp) - 1
      iopsy(icomp)  = jbtof(irep(icomp), 1)

      do i = 1, nfsym 
        j = mod(i + iopsy(icomp), 2) + 1 

        call qtrans90('AOMO',                        &
                      'S',                           &
                      0.0d0,                         &
                      nfbas(i, 0),                   &
                      nfbas(j, 0),                   &
                      norb(i),                       &
                      norb(j),                       &
!                     --------------------------------
                      fao(1 + i2basx(i, j), icomp),  &
                      ntbas(0),                      &
                      ntbas(0),                      &
                      nz,                            &
                      ipqtoq(1, irep(icomp)),        &
!                     --------------------------------
                      fmo(1 + i2orbx(i, j),icomp),   &
                      norbt,                         &
                      norbt,                         &
                      nz,                            &
                      ipqtoq(1, irep(icomp)),        &
!                     --------------------------------
                      cmo(1 + icmoq(i)),             &
                      nfbas(i, 0),                   &
                      norb(i),                       &
                      nz,                            &
                      ipqtoq(1, 0),                  &
!                     --------------------------------
                      cmo(1 + icmoq(j)),             &
                      nfbas(j, 0),                   & 
                      norb(j),                       &
                      nz,                            &
                      ipqtoq(1, 0),                  &
!                     --------------------------------
                      iprprp)

      debug_me = .false.
      if (debug_me) then
        write(*, *) 'fao in do_ao_to_mo after qtrans, icomp, i = ', icomp, i
        call prqmat(fao(1, icomp), &
                    ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
                    nz, ipqtoq(1, irep(icomp)), lupri)
        write(*, *) 'fmo in do_ao_to_mo after qtrans, icomp, i = ', icomp, i
        call prqmat(fmo(1, icomp), &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, irep(icomp)), lupri)
      end if

      end do
    end do

  end subroutine


  subroutine get_tb_dmat(tb_dmat, cmo, ibeig)
!-----------------------------------------------------------------------
    real(8), intent(out)   :: tb_dmat(n2bbasxq, *)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork, kfree

    real(8), allocatable   :: tb(:)
    integer                :: isym(nr_xyz_comp), iopsy(nr_xyz_comp)
    integer                :: irep(nr_xyz_comp)
    integer                :: indxpr, icomp
!-----------------------------------------------------------------------

    call opnfil(lu1int, 'AOPROPER', 'OLD', 'get_tb_dmat')

    call alloc(tb, n2orbxq)
    tb = 0.0d0

    do icomp = 1, nr_xyz_comp
      indxpr       = ipcon(icomp)
      isym(icomp)  = iprpsym(indxpr)
      irep(icomp)  = isym(icomp) - 1
      iopsy(icomp) = jbtof(irep(icomp), 1)

      call dzero(tb, n2orbxq)

      lwork = len_f77_work
      kfree = 1
      call di_select_wrk(work, lwork)

!     get T^B_icomp matrix into tb
!     gosia: if lao_lr_rearrange, then here we read in 'dS/dB' integrals
!     it was set in pammag.F in def_d1hblond subroutine
      call prpmat(indxpr,           &
                  iopsy(icomp),     &
                  tb,               &
                  cmo,              &
                  ibeig,            &
                  icmoq,            &
                  norb,             &
                  work,             &
                  kfree,            &
                  lwork,            &
                  iprprp)

      call di_deselect_wrk(work, lwork)

#ifdef MOD_LAO_REARRANGED
      if (lao_lr_rearrange) then
!       we only need occ-occ block of S^B matrix, so zero-out the rest:
        call get_m_ij(tb) 
      end if
#endif

      debug_me = .false.
      if (debug_me) then
        write(*, *) 'tb matrix in MO basis in get_tb_dmat'
        call prqmat(tb, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, irep(icomp)), lupri)
      end if


      lwork = len_f77_work
      kfree = 1
      call di_select_wrk(work, lwork)

      call getmdm(doeprn,            &
                  isym(icomp),       &
                  tb_dmat(1, icomp), &
                  cmo,               &
                  ibeig,             &
                  tb,                &
                  work,              &
                  kfree,             &
                  lwork,             &
                  iprprp)

      call di_deselect_wrk(work, lwork)

    end do

    close(lu1int)
    call dealloc(tb)

  end subroutine

      SUBROUTINE GET_NPOS(NPOS)
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "cbihr2.h"
#include "blocks.h"
#include "dcbfir.h"
      integer npos
      call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. ,.TRUE. /))
      call SetIntTaskArrayDimension(NPOS,PARCAL)
      end subroutine

end module
