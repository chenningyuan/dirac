!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

module london_direct_hb

  use memory_allocator 
  use dirac_cfg
  use dft_cfg
  use num_grid_gen
  use dirac_interface

  implicit none

  public london_twoel_hb

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

#include "infpar.h"



  integer, parameter :: max_nr_xyz_comp = 3
! there should be nr_mat=nzc1 (which is 4 except 2component spinfree... see dirac/dirrdn.F, grep for nzC1)
  integer, parameter :: nr_mat = 4

  integer :: isym_2el_lao(max_nr_xyz_comp)
  integer :: irep_2el_lao(max_nr_xyz_comp)
  integer :: iherm_2el_lao(max_nr_xyz_comp)
  integer :: ifckop_2el_lao(max_nr_xyz_comp)
  integer :: jbtof_2el_lao(max_nr_xyz_comp)

  integer :: ifock_type(nr_mat*max_nr_xyz_comp)
  integer :: irepdm(nr_mat*max_nr_xyz_comp)


  integer, parameter :: iqmult(4, 4, 4) = reshape((/1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1,   &
                                                    2, 1, 4, 3, 1, 2, 3, 4, 4, 3, 2, 1, 3, 4, 1, 2,   &
                                                    3, 4, 1, 2, 4, 3, 2, 1, 1, 2, 3, 4, 2, 1, 4, 3,   &
                                                    4, 3, 2, 1, 3, 4, 1, 2, 2, 1, 4, 3, 1, 2, 3, 4/), &
                                                  (/4, 4, 4/))

  integer, parameter :: IQPHASE(4, 4, 4) = reshape((/1, 1, 1, 1, 1,-1,-1, 1, 1, 1,-1,-1, 1,-1, 1,-1,   &
                                                     1,-1,-1, 1,-1,-1,-1,-1, 1,-1, 1,-1,-1,-1, 1, 1,   &
                                                     1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1, 1,   &
                                                     1,-1, 1,-1, 1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1/), &
                                                   (/4, 4, 4/))


contains


  subroutine london_twoel_hb(cmo, ibeig, work, lwork)
!-----------------------------------------------------------------------
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)
    integer, intent(in)    :: lwork
    real(8), intent(inout) :: work(lwork)

    real(8), allocatable :: dmat_ao(:, :)
    real(8), allocatable :: fmat(:, :, :)
!-----------------------------------------------------------------------

    call di_set_f77_work(work, lwork)

    call initialize_direct_lao()

!   density matrix
    call alloc(dmat_ao, n2bbasx, nr_mat)
    dmat_ao = 0.0d0
    call get_dmat_ao(cmo, ibeig, dmat_ao)

!   fock matrix
    call alloc(fmat, n2bbasx, nr_mat, max_nr_xyz_comp)
    fmat = 0.0d0

!   get 2-electron lao-integrals 
    call get_ao_fock_2el_lao(cmo, ibeig, fmat, dmat_ao)

    call dealloc(dmat_ao)

!   if calculations with symmetry we have to 'pack' fmat
!   before we add dft-lao contribution
    if (maxrep > 0) then
       call sym_antisym_fock_2el_lao(fmat)
       call ao_to_so_fock_2el_lao(fmat)
    end if

!   add xc contribution if dft calculations
    if (dirac_cfg_dft_calculation) then
       call xc_contrib_2el_lao(cmo, fmat)
    end if

!   transform from ao to mo basis
    call ao_to_mo_2el_lao(cmo, ibeig, fmat)

!   clean
    call dealloc(fmat)

  end subroutine


  subroutine sym_antisym_fock_2el_lao(fmat)
!   --------------------------------------------------------------------
    real(8), intent(inout) :: fmat(n2bbasx, nr_mat, max_nr_xyz_comp)
    integer :: n, h(4), icomp, i, j, ioff1, ioff2
    real(8) :: f
!   --------------------------------------------------------------------

    h(1)   = -1
    h(2:4) = 1

    do icomp = 1, max_nr_xyz_comp
      do n = 1, nr_mat
        do i = 1, ntbas(0)
          do j = 1, i
            ioff1 = (i - 1)*ntbas(0) + j
            ioff2 = (j - 1)*ntbas(0) + i
            f = 0.5d0*(fmat(ioff2, n, icomp) + h(n)*fmat(ioff1, n, icomp))
            fmat(ioff2, n, icomp) =      f
            fmat(ioff1, n, icomp) = h(n)*f
          end do
        end do
      end do
    end do


  end subroutine


  subroutine ao_to_so_fock_2el_lao(fmat)
!   --------------------------------------------------------------------
    real(8), intent(inout) :: fmat(n2bbasx, nr_mat, max_nr_xyz_comp)
    real(8), allocatable   :: fmat_help(:, :)
    integer :: i, icomp, iq, iz, ipq, irepd, irep, ipari, ioff, itim, ih
!   --------------------------------------------------------------------

    call alloc(fmat_help, n2bbasx, nz)
    fmat_help = 0.0d0

    do icomp = 1, max_nr_xyz_comp

      call dzero(fmat_help, n2bbasxq)

      irep  = isymax(icomp, 2)
      ipari = jbtof(irep, 1)

!     transform fmat from AO to SO basis
!     ----------------------------------
      do iz = 1, nr_mat
        irepd = irqmat(iz, irep)
        iq    = iqmult(1, jqbas(irepd, ipari), iz)
        ipq   = iqtopq(iq, irep)
        ioff  = (iz - 1)*ntbas(0)*ntbas(0) + 1
        call mtaoso(fmat(1, iz, icomp), fmat_help(1, ipq), ntbas(0), irepd, iprprp)
      end do

!     transform fmat from Hermit-sorted to DIRAC-sorted AO basis
!     ----------------------------------------------------------
      call butobs_no_work(fmat_help, nz)

!     insert half-phases
!     ------------------
      if (nz .lt. 4) then
        do iz = 1, nz
          iq   = ipqtoq(iz, irep)
          call q2bphase('F', iq, 1, fmat_help(1, iz))
        end do
      end if

      call dcopy(ntbas(0)*ntbas(0)*nz, fmat_help, 1, fmat(1, 1, icomp), 1)

    end do

    call dealloc(fmat_help)

  end subroutine



  subroutine initialize_direct_lao()
!   --------------------------------------------------------------------
    integer :: i, j

#include "dcbfir.h"
!   --------------------------------------------------------------------

!     symmetry info for 2-el LAO integrals
      do i = 1, max_nr_xyz_comp
      
!       old code: isymop
        isym_2el_lao(i)  = isymax(i, 2) + 1
        irep_2el_lao(i)  = isymax(i, 2)
      
!       old code: ihrmop
!                 0 -> no symmetry
!                 1 -> hermitian
!                -1 -> antihermitian
        iherm_2el_lao(i) = -1
      
!       old code: ifckop:
!                 1 -> only Coulomb
!                 2 -> only Exchange
!                 3 -> Coulomb+Exchange
        ifckop_2el_lao(i) = 1
      
        jbtof_2el_lao(i)  = jbtof(irep_2el_lao(i), 1)
      end do

  end subroutine


  subroutine get_ao_fock_2el_lao(cmo, ibeig, fmat, dmat)
!   --------------------------------------------------------------------
    real(8), intent(inout) :: fmat(n2bbasx, nr_mat, max_nr_xyz_comp)
    real(8), intent(inout) :: dmat(n2bbasx, nr_mat)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork

    integer, allocatable   :: npos(:)   
    integer :: i2typ
    integer :: itype, max_deriv, iatom
    logical :: nodv, nopv, nocont, retur, relcal, ddfock, fckddr, &
               london, first(3)

    logical :: debug_me
    integer :: iq, ipq, iz, irepd, i, ioff(3)
    integer :: mtottk(3)
    real(8) :: cpu1, cpu2, wall1, wall2, cputot, wlltot
!   --------------------------------------------------------------------

!   set variables used in this subroutine 
!   -------------------------------------
    itype     = -5
    max_deriv = 1
    iatom     = 0
    nodv      = .true.
    nopv      = .true.
    nocont    = .false.
    retur     = .false.
    relcal    = .true.
    ddfock    = .true.
    fckddr    = .true.
    london    = .true.
    first     = .true.

    call setfck(ifock_type, irepdm, max_nr_xyz_comp, nr_mat, & 
                isym_2el_lao, iherm_2el_lao, ifckop_2el_lao, &
                iprprp)

#ifdef VAR_MPI
      mtottk(1) = nlrgbl*(nlrgbl+1)/2
      mtottk(2) = nlrgbl*(nlrgbl+1)/2
      mtottk(3) = nsmlbl*(nsmlbl+1)/2
! from dcbfir.h,
! her_pardrv needs it and it is passed through common blocks
      first1 = .true.
      first2 = .true.
      first3 = .true.
#endif


!   get 2-electron London integrals
!   -------------------------------
    call gettim(cpu1,wall1)


    do i2typ = 1, 3
      if (btest(intnmr, (i2typ - 1))) then


#ifdef VAR_MPI
        if (parcal) then

           call alloc(npos, mtottk(i2typ))
           npos = 0
           call dirac_parctl(1)

           lwork = len_f77_work
           call di_select_wrk(work, lwork)

           call her_pardrv(work,lwork,fmat,dmat,nr_mat,irepdm,                   &
                           ifock_type,itype,max_deriv,iatom,nodv,nopv,nocont,   &
                           tktime,.false.,first, npos,mtottk(i2typ),          &
                           i2typ,(/0/),(/0.0d0/),(/0.0d0/),(/0.0d0/),(/0.0d0/),(/0.0d0/))

           call di_deselect_wrk(work, lwork)
           call dealloc(npos)

        else

          lwork = len_f77_work
          call di_select_wrk(work, lwork)
    
          call twoint(work, lwork, fmat, dmat, nr_mat, irepdm, ifock_type, &
                      (/0.0d0/), (/0/), (/0/), itype, max_deriv, iatom,   &
                      nodv, nopv, nocont, tktime, iprprp, iprnta, &
                      iprntb, iprntc, iprntd, retur, (/0/), i2typ,  &
                      (/0/), (/0.0d0/), (/0.0d0/), (/0.0d0/), (/0.0d0/),  &
                      (/0.0d0/), relcal, .false., (/0/), (/0.0d0/))
    
          call di_deselect_wrk(work, lwork)

        end if

#else

          lwork = len_f77_work
          call di_select_wrk(work, lwork)
    
          call twoint(work, lwork, fmat, dmat, nr_mat, irepdm, ifock_type, &
                      (/0.0d0/), (/0/), (/0/), itype, max_deriv, iatom,   &
                      nodv, nopv, nocont, tktime, iprprp, iprnta, &
                      iprntb, iprntc, iprntd, retur, (/0/), i2typ,  &
                      (/0/), (/0.0d0/), (/0.0d0/), (/0.0d0/), (/0.0d0/),  &
                      (/0.0d0/), relcal, .false., (/0/), (/0.0d0/))
    
          call di_deselect_wrk(work, lwork)

#endif

      end if
    end do

    call gettim(cpu2,wall2)
    cputot = cpu2 - cpu1
    wlltot = wall2 - wall1
    write(lupri,'(/a,f13.4,a,f13.4,a/)') 'cpu time:',cputot,' wall time:',wlltot, &
                                           ' used in call to twoint(in get_ao_fock_2el_lao).'

    if (allocated(npos)) then
      call dealloc(npos)
    end if

  end subroutine


  subroutine get_dmat_ao(cmo, ibeig, dmat_ao)
!   --------------------------------------------------------------------
    real(8), intent(inout) :: dmat_ao(n2bbasx, *)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork
    real(8), allocatable   :: dmat(:, :)
    integer                :: iz, irep, ipari, irepd, iq, ipq, ioff
!   --------------------------------------------------------------------

    call alloc(dmat, ntbas(0)*ntbas(0), nz)
    dmat = 0.0d0
    call genden(dmat, cmo, 1, iprprp)
    call dscal(ntbas(0)*ntbas(0)*nz, 2.0d0, dmat, 1)
    call bstobu_no_work(dmat, nz)

    isymop(1) = 1

    if (nz .lt. 4) then
      irep = isymop(1) - 1
      do iz = 1, nz
        iq = ipqtoq(iz, irep)
        call q2phase('D', iq, 1, dmat(1, iz))
      end do
    end if

    irep = isymop(1) - 1
    ipari = jbtof(irep, 1)
    do iz = 1, nr_mat
      irepd = irqmat(iz, irep)
      iq    = iqmult(1, jqbas(irepd, ipari), iz)
      ipq   = iqtopq(iq, irep)
      call dtsoao(dmat(1, ipq), dmat_ao(1, iz), ntbas(0), irepd, iprprp)
    end do

    call dealloc(dmat)


  end subroutine


  subroutine xc_contrib_2el_lao(cmo, fmat)
    use xcint_main
!   --------------------------------------------------------------------
    real(8), intent(in)    :: cmo(*)
    real(8), intent(inout) :: fmat(n2bbasx, nz, max_nr_xyz_comp)
    real(8), allocatable   :: dmat(:, :)

    integer :: icomp
    logical :: debug_me
!   --------------------------------------------------------------------

    call alloc(dmat, n2bbasx, nz)
    dmat = 0.0d0
    call genden(dmat, cmo, 1, iprprp)

!   gosia todo: Bx, By and Bz contributions are added on top, if we want just one
!   of these contributions (for instance for VISUAL), then routines in density_eval
!   need to be adapted...

    call generate_num_grid(dmat)
#ifdef VAR_MPI
            if (parcal) call dirac_parctl(4)
#endif
    call integrate_xc(xc_mat_dim              = ntbas(0),     &
                       xc_nz                   = nz,           &
                       xc_dmat_0               = dmat,            &
                       xc_nr_dmat              = 0,               &
                       xc_nr_fmat              = max_nr_xyz_comp, &
                       xc_fmat                 = fmat,            &
                       xc_fmat_pg_sym          = isym_2el_lao,    &
                       xc_do_london_rhs_direct = .true.)

    call dealloc(dmat)

  end subroutine


  subroutine ao_to_mo_2el_lao(cmo, ibeig, fmat)
!   --------------------------------------------------------------------
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)
    real(8), intent(inout) :: fmat(n2bbasx, nz, max_nr_xyz_comp)

    real(8), allocatable :: gmat(:)
    integer :: icomp, ii, jj, irep, iopsy
    logical :: debug_me
!   --------------------------------------------------------------------

    call alloc(gmat, n2orbxq)
    gmat = 0.0d0

    open(unit = lutwol, file = 'TWOLON', form = 'unformatted',  &
         access = 'direct', recl = 8*n2orbxq, action = 'write', &
         status = 'unknown')

    do icomp = 1, max_nr_xyz_comp
      irep = isymax(icomp, 2)

      call dzero(gmat, n2orbxq)
     
      do ii = 1, nfsym
        iopsy = jbtof(irep, 1)
        jj = mod(ii + iopsy, 2) + 1 

! todo add qbtrans for spin-free case

        call qtrans90('AOMO',                             &
                      'S',                                &
                      0.0d0,                              &
                      nfbas(ii, 0), nfbas(jj, 0),         &
                      norb(ii), norb(jj),                 &
!                     -------------------------------------
                      fmat(i2basx(ii, jj) + 1, 1, icomp), &
                      ntbas(0), ntbas(0), nz,             &
                      ipqtoq(1, irep_2el_lao(icomp)),     &
!                     -------------------------------------
                      gmat(i2orbx(ii, jj) + 1),           &
                      norbt,    norbt,    nz,             &
                      ipqtoq(1, irep_2el_lao(icomp)),     &
!                     -------------------------------------
                      cmo(icmoq(ii) + 1),                 &
                      nfbas(ii,0), norb(ii), nz,          &
                      ipqtoq(1, 0),                       &
!                     -------------------------------------
                      cmo(icmoq(jj) + 1),                 &
                      nfbas(jj, 0), norb(jj), nz,         &
                      ipqtoq(1, 0),                       &
!                     -------------------------------------
                      iprprp)

      end do

!     write the two-electron London contribution to file   
!     --------------------------------------------------
      call wrtdac(lutwol, n2orbxq, gmat(1), icomp)

      debug_me = .false.
      if (debug_me) then
        write(*, *) 'fmat in AO, icomp = ', icomp
        call prqmat(fmat(1, 1, icomp), &
                    ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
                    nz, ipqtoq(1, irep), lupri)
        write(*, *) 'gmat in MO, icomp = ', icomp
        call prqmat(gmat(1), &
                    norb, norb, norb, norb, &
                    nz, ipqtoq(1, irep), lupri)
      end if

    end do

    close(lutwol, status = 'keep')
    call dealloc(gmat)

  end subroutine

end module
