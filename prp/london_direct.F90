!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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

  use london_prp_cfg
  use memory_allocator 
  use dirac_cfg
  use dft_cfg
  use num_grid_gen

  implicit none

  public london_twoel_hb

  save

  private

  real(8), pointer :: f77_work(:)       ! Fortran 77 work array pointer
  integer          :: len_f77_work  = 0 ! Length of the work array
  integer          :: memory_usage  = 0 ! Amount of work array used
  integer          :: mem_wrk_count = 0 ! Counts

  integer, parameter :: max_nr_xyz_comp = 3
! there should be nr_mat=nzc1 (which is 4 except 2component spinfree... see dirac/dirrdn.F, grep for NZC1)
  integer, parameter :: nr_mat = 4

  integer :: isym_2el_lao(max_nr_xyz_comp)
  integer :: irep_2el_lao(max_nr_xyz_comp)
  integer :: iherm_2el_lao(max_nr_xyz_comp)
  integer :: ifckop_2el_lao(max_nr_xyz_comp)
  integer :: jbtof_2el_lao(max_nr_xyz_comp)

  integer :: ifock_type(nr_mat*max_nr_xyz_comp)
  integer :: irepdm(nr_mat*max_nr_xyz_comp)


  logical :: it_is_master, it_is_slave



  integer :: ierr
  logical :: test_sdft_lao


! parcal is in dcbgen.h
#include "dcbgen.h"
! first1, first2, first3 are in dcbfir.h, needed by her_pardrv
! (and passed by common blocks, not as arguments...)
#include "dcbfir.h"

contains


  subroutine london_twoel_hb(cmo, ibeig, work, lwork)
!-----------------------------------------------------------------------
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)
    integer, intent(in)    :: lwork
    real(8), intent(inout) :: work(lwork)

    integer :: icomp
    real(8), allocatable   :: dmat_work(:, :)
    real(8), allocatable   :: dmat_4_work(:, :)
    real(8), allocatable, target   :: fmat_work(:, :, :)
    real(8), allocatable   :: fmat_help_work(:, :)
    real(8), pointer               :: fmat_p(:, :, :)

! gosia:
! *) this is the main subroutine that calculates the 
!    2-electron lao contribution to property gradient (dH/dB|_0)
!-----------------------------------------------------------------------

    call set_f77_work(work, lwork)

    call initialize_direct_lao()

!   density matrix
!   --------------
    call alloc(dmat_4_work, nr_ao_total*nr_ao_total, nr_mat)
    call alloc(dmat_work, nr_ao_total*nr_ao_total, nr_qu)
    call call_gtcnvd(cmo, ibeig, dmat_work, dmat_4_work)
    call dealloc(dmat_work)


!   fock matrix
!   -----------
    call alloc(fmat_work, nr_ao_total*nr_ao_total, nr_mat, max_nr_xyz_comp)

!   get 2-electron lao-integrals 
!   ----------------------------
    call get_ao_fock_2el_lao(cmo, ibeig, fmat_work, dmat_4_work)

    call dealloc(dmat_4_work)

!   if calculations with symmetry we have to 'pack' fmat
!   before we add dft-lao contribution
!   --------------------------------------------------------------------
    if (p_maxrep .gt. 0) then

      call sym_antisym_fock_2el_lao(fmat_work)
      
      call alloc(fmat_help_work, nr_ao_total*nr_ao_total, nr_qu)
      
      call ao_to_so_fock_2el_lao(fmat_help_work, fmat_work)
      
      call dealloc(fmat_help_work)

      nullify(fmat_p)
      fmat_p => fmat_work(:, 1:nr_qu, :)

    else
      
      nullify(fmat_p)
      fmat_p => fmat_work
      
    end if


!   add xc contribution if dft calculations
!   ---------------------------------------
    if (dirac_cfg_dft_calculation) then
    
      call alloc(dmat_work, nr_ao_total*nr_ao_total, nr_qu)
      call xc_contrib_2el_lao(cmo, dmat_work, fmat_p)
    
      if (allocated(dmat_work))   call dealloc(dmat_work)

    end if


!   transform from ao to mo basis
!   -----------------------------
    call ao_to_mo_2el_lao(cmo, ibeig, fmat_p)


!   clean
!   -----
    if (associated(fmat_p))   nullify(fmat_p)
    if (allocated(fmat_work)) call dealloc(fmat_work)

  end subroutine


  subroutine sym_antisym_fock_2el_lao(fmat)
!   --------------------------------------------------------------------
    real(8), intent(inout) :: fmat(nr_ao_total*nr_ao_total, nr_mat, max_nr_xyz_comp)
    integer :: n, h(4), icomp, i, j, ioff1, ioff2
    real(8) :: f
!   --------------------------------------------------------------------

    h(1)   = -1
    h(2:4) = 1

    do icomp = 1, max_nr_xyz_comp
      do n = 1, nr_mat
        do i = 1, nr_ao_total
          do j = 1, i
            ioff1 = (i - 1)*nr_ao_total + j
            ioff2 = (j - 1)*nr_ao_total + i
            f = 0.5d0*(fmat(ioff2, n, icomp) + h(n)*fmat(ioff1, n, icomp))
            fmat(ioff2, n, icomp) =      f
            fmat(ioff1, n, icomp) = h(n)*f
          end do
        end do
      end do
    end do


  end subroutine


  subroutine ao_to_so_fock_2el_lao(fmat_help, fmat)
!   --------------------------------------------------------------------
    real(8), intent(inout) :: fmat_help(nr_ao_total*nr_ao_total, nr_qu)
    real(8), intent(inout) :: fmat(nr_ao_total*nr_ao_total, nr_mat, max_nr_xyz_comp)
    integer :: i, icomp, iq, iz, ipq, irepd, irep, ipari, ioff, itim, ih
!   --------------------------------------------------------------------

    do icomp = 1, max_nr_xyz_comp

      call dzero(fmat_help, nr_ao_total*nr_ao_total*nr_qu)

      irep  = p_isymax(icomp, 2)
      ipari = p_jbtof(irep, 1)

!     transform fmat from AO to SO basis
!     ----------------------------------
      do iz = 1, nr_mat
        irepd = p_irqmat(iz, irep)
        iq    = p_iqmult(1, p_jqbas(irepd, ipari), iz)
        ipq   = p_iqtopq(iq, irep)
        ioff  = (iz - 1)*nr_ao_total*nr_ao_total + 1
        call mtaoso(fmat(1, iz, icomp), fmat_help(1, ipq), nr_ao_total, irepd, prp_iprint)
      end do

!     transform fmat from Hermit-sorted to DIRAC-sorted AO basis
!     ----------------------------------------------------------
      call butobs_no_work(fmat_help, nr_qu)

!     insert half-phases
!     ------------------
      if (nr_qu .lt. 4) then
        do iz = 1, nr_qu
          iq   = p_ipqtoq(iz, irep)
          call q2bphase('F', iq, 1, fmat_help(1, iz))
        end do
      end if

      call dcopy(nr_ao_total*nr_ao_total*nr_qu, fmat_help, 1, fmat(1, 1, icomp), 1)

    end do

  end subroutine



  subroutine initialize_direct_lao()
!   --------------------------------------------------------------------
    integer :: i, j

#include "dcbfir.h"
!   --------------------------------------------------------------------

!   set variables necessary for london_direct_hb
!   --------------------------------------------
    call common_blocks_vars()

    it_is_master = .true.
    it_is_slave  = .false.

!     symmetry info for 2-el LAO integrals
      do i = 1, max_nr_xyz_comp
      
!       old code: isymop
        isym_2el_lao(i)  = p_isymax(i, 2) + 1
        irep_2el_lao(i)  = p_isymax(i, 2)
      
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
      
        jbtof_2el_lao(i)  = p_jbtof(irep_2el_lao(i), 1)
      end do

  end subroutine


  subroutine get_ao_fock_2el_lao(cmo, ibeig, fmat, dmat_4)
!   --------------------------------------------------------------------
    real(8), intent(inout)    :: fmat(nr_ao_total*nr_ao_total, nr_mat, max_nr_xyz_comp)
    real(8), intent(inout)    :: dmat_4(nr_ao_total*nr_ao_total, nr_mat)
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
                prp_iprint)

#ifdef VAR_MPI
      mtottk(1) = p_nlrgbl*(p_nlrgbl+1)/2
      mtottk(2) = p_nlrgbl*(p_nlrgbl+1)/2
      mtottk(3) = p_nsmlbl*(p_nsmlbl+1)/2
! from dcbfir.h,
! it is not nice, but her_pardrv needs it and it is passed through common blocks
      first1 = .true.
      first2 = .true.
      first3 = .true.
#endif


!   get 2-electron London integrals
!   -------------------------------
    call gettim(cpu1,wall1)


    do i2typ = 1, 3
      if (btest(nmr_int, (i2typ - 1))) then


#ifdef VAR_MPI
        if (parcal) then

           call alloc(npos, mtottk(i2typ))
           call dirac_parctl(1)

           lwork = len_f77_work
           !lwork = max_nr_xyz_comp*nr_mat*nr_ao_total*nr_ao_total*2
           call select_wrk(work, lwork)

           call her_pardrv(work,lwork,fmat,dmat_4,nr_mat,irepdm,                   &
                           ifock_type,itype,max_deriv,iatom,nodv,nopv,nocont,   &
                           p_tktime,.false.,first, npos,mtottk(i2typ),          &
                           i2typ,(/0/),(/0.0d0/),(/0.0d0/),(/0.0d0/),(/0.0d0/),(/0.0d0/))

           call deselect_wrk(work, lwork)
           call dealloc(npos)

        else

          lwork = len_f77_work
          call select_wrk(work, lwork)
    
          call twoint(work, lwork, fmat, dmat_4, nr_mat, irepdm, ifock_type, &
                      (/0.0d0/), (/0/), (/0/), itype, max_deriv, iatom,   &
                      nodv, nopv, nocont, p_tktime, prp_iprint, p_iprnta, &
                      p_iprntb, p_iprntc, p_iprntd, retur, (/0/), i2typ,  &
                      (/0/), (/0.0d0/), (/0.0d0/), (/0.0d0/), (/0.0d0/),  &
                      (/0.0d0/), relcal, .false., (/0/), (/0.0d0/))
    
          call deselect_wrk(work, lwork)

        end if

#else

          lwork = len_f77_work
          call select_wrk(work, lwork)
    
          call twoint(work, lwork, fmat, dmat_4, nr_mat, irepdm, ifock_type, &
                      (/0.0d0/), (/0/), (/0/), itype, max_deriv, iatom,   &
                      nodv, nopv, nocont, p_tktime, prp_iprint, p_iprnta, &
                      p_iprntb, p_iprntc, p_iprntd, retur, (/0/), i2typ,  &
                      (/0/), (/0.0d0/), (/0.0d0/), (/0.0d0/), (/0.0d0/),  &
                      (/0.0d0/), relcal, .false., (/0/), (/0.0d0/))
    
          call deselect_wrk(work, lwork)

#endif

      end if
    end do


    call gettim(cpu2,wall2)
    cputot = cpu2 - cpu1
    wlltot = wall2 - wall1
    write(p_lupri,'(/a,f13.4,a,f13.4,a/)') 'cpu time:',cputot,' wall time:',wlltot, &
                                           ' used in call to twoint(in get_ao_fock_2el_lao).'



    if (allocated(npos)) then
      call dealloc(npos)
    end if

  end subroutine


  subroutine call_gtcnvd(cmo, ibeig, dmat, dmat_4)
!   --------------------------------------------------------------------
    real(8), intent(inout)    :: dmat(nr_ao_total*nr_ao_total, nr_qu)
    real(8), intent(inout)    :: dmat_4(nr_ao_total*nr_ao_total, nr_mat)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)
    real(8), pointer       :: work(:)
    integer                :: lwork

    integer :: iz, irep, ipari, irepd, iq, ipq, ioff
! isymop is from dcbdhf.h:
#include "dcbdhf.h"
!   --------------------------------------------------------------------

!    lwork = nr_ao_total*nr_ao_total*nr_qu*2
!    call select_wrk(work, lwork)
!    call gtcnvd(dmat_4, cmo, ibeig, prp_iprint, work, lwork)
!    call deselect_wrk(work, lwork)

    call genden(dmat, cmo, 1, prp_iprint)
    call dscal(nr_ao_total*nr_ao_total*nr_qu, 2.0d0, dmat, 1)
    call bstobu_no_work(dmat, nr_qu)

    if (nr_qu .lt. 4) then
      irep  = irep_2el_lao(1)
      do iz = 1, nr_qu
        iq = p_ipqtoq(iz, irep)
        call q2phase('D', iq, 1, dmat(1, iz))
      end do
    end if

    irep  = irep_2el_lao(1)
    ipari = p_jbtof(irep, 1)
    do iz = 1, nr_mat
      irepd = p_irqmat(iz, irep)
      iq    = p_iqmult(1, p_jqbas(irepd, ipari), iz)
      ipq   = p_iqtopq(iq, irep)
      write(*, *) 'alaola irep,irepd,ipari,iq,ipq=', &
                          irep,irepd,ipari,iq,ipq
      CALL OUTPUT(dmat(1,ipq),1,nr_ao_total,1,nr_ao_total,&
                  nr_ao_total,nr_ao_total,-1,p_lupri)
      call dtsoao(dmat(1, ipq), dmat_4(1, iz), nr_ao_total, irepd, prp_iprint)
    end do

  end subroutine


  subroutine xc_contrib_2el_lao(cmo, dmat, fmat)
    use xc_driver
!   --------------------------------------------------------------------
    real(8), intent(in)    :: cmo(*)
    real(8), intent(inout)    :: dmat(nr_ao_total*nr_ao_total, nr_qu)
    real(8), intent(inout)    :: fmat(nr_ao_total*nr_ao_total, nr_qu, max_nr_xyz_comp)

    integer :: icomp
    logical :: debug_me
!   --------------------------------------------------------------------

    call genden(dmat, cmo, 1, prp_iprint)

    !debug_me = .false.
    !if (debug_me) then
    !  do icomp = 1, max_nr_xyz_comp
    !    write(p_lupri, *) 'london_direct: fmat before xc_integrator, icomp = ', icomp
    !    call prqmat(fmat_p(1, 1, icomp), &
    !                nr_ao_total, nr_ao_total, nr_ao_total, nr_ao_total,  &
    !                nr_qu, p_ipqtoq(1, p_isymax(icomp, 2)), p_lupri)
    !  end do
    !end if


!   gosia todo: Bx, By and Bz contributions are added on top, if we want just one
!   of these contributions (for instance for VISUAL), then routines in density_eval
!   need to be adapted...

    call generate_num_grid(dmat)
#ifdef VAR_MPI
            if (parcal) call dirac_parctl(4)
#endif
    if (dft_cfg_no_sdft) then
      call xc_integrator(xc_mat_dim              = nr_ao_total,   &
                         xc_nz                   = nr_qu,        &
                         xc_dmat_0               = dmat,         &
                         xc_nr_dmat              = 0,            &
                         xc_nr_fmat              = max_nr_xyz_comp,  &
                         xc_fmat                 = fmat,         &
                         xc_fmat_pg_sym          = isym_2el_lao, &
                         xc_do_london_rhs_direct = .true.)
    else
      call xc_integrator(xc_mat_dim              = nr_ao_total,      &
                         xc_nz                   = nr_qu,           &
                         xc_dmat_0               = dmat,            &
                         xc_nr_dmat              = 0,               &
                         xc_nr_fmat              = max_nr_xyz_comp,     &
                         xc_fmat                 = fmat,         &
                         xc_fmat_pg_sym          = isym_2el_lao,    &
                         xc_do_london_rhs_direct = .true.,       &
                         xc_response_order_mo    = 1)
    end if

    !if (debug_me) then
    !  do icomp = 1, max_nr_xyz_comp
    !    write(p_lupri, *) 'london_direct: fmat after xc_integrator, icomp = ', icomp
    !    call prqmat(fmat_p(1, 1, icomp), &
    !                nr_ao_total, nr_ao_total, nr_ao_total, nr_ao_total,  &
    !                nr_qu, p_ipqtoq(1, p_isymax(icomp, 2)), p_lupri)
    !  end do
    !end if


  end subroutine


  subroutine ao_to_mo_2el_lao(cmo, ibeig, fmat)
!   --------------------------------------------------------------------
    real(8), intent(in) :: cmo(*)
    integer, intent(in) :: ibeig(*)
    real(8), intent(inout) :: fmat(nr_ao_total*nr_ao_total, nr_mat, max_nr_xyz_comp)

    real(8), pointer       :: work(:)
    integer                :: lwork

    real(8), allocatable :: gmat(:)
    integer :: gmat_size, icomp, ii, jj, irep, iopsy
!   --------------------------------------------------------------------

    gmat_size = p_n2orbxq
    call alloc(gmat, gmat_size)

    open(twolon_unit, file = 'TWOLON', form = 'unformatted',        &
         access = 'direct', recl = 8*gmat_size, status = 'unknown')

    do icomp = 1, max_nr_xyz_comp
      irep = p_isymax(icomp, 2)

      call dzero(gmat, gmat_size)
     
      do ii = 1, p_nfsym
        iopsy = p_jbtof(irep, 1)
        jj = mod(ii + iopsy, 2) + 1

        lwork = len_f77_work
        call select_wrk(work, lwork)

! todo add qbtrans for spin-free case
! qtrans is not parallelized, so collect data from nodes before calling this

        call qtrans('AOMO',                            &
                    'S',                               &
                    0.0d0,                             &
                    p_nfbas(ii, 0),                    &
                    p_nfbas(jj, 0),                    &
                    p_norb(ii),                        &
                    p_norb(jj),                        &
                    fmat(p_i2basx(ii, jj) + 1, 1, icomp), &
                    nr_ao_total,                        &
                    nr_ao_total,                        &
                    nr_qu,                             &
                    p_ipqtoq(1, irep_2el_lao(icomp)),  &
                    gmat(p_i2orbx(ii, jj) + 1),        &
                    p_norbt,                           &
                    p_norbt,                           &
                    nr_qu,                             &
                    p_ipqtoq(1, irep_2el_lao(icomp)),  &
                    cmo(p_icmoq(ii) + 1),              &
                    p_nfbas(ii,0),                     &
                    p_norb(ii),                        &
                    nr_qu,                             &
                    p_ipqtoq(1, 0),                    &
                    cmo(p_icmoq(jj) + 1),              &
                    p_nfbas(jj, 0),                    &
                    p_norb(jj),                        &
                    nr_qu,                             &
                    p_ipqtoq(1, 0),                    &
                    work,                              &
                    lwork,                             &
                    prp_iprint)

        call deselect_wrk(work, lwork)

!       write the two-electron London contribution to file   
!       --------------------------------------------------
        call wrtdac(twolon_unit, gmat_size, gmat, icomp)

      end do
    end do

    close(twolon_unit, status = 'keep')
    call dealloc(gmat)

  end subroutine


! ******************************************************************************

  SUBROUTINE set_f77_work(wrk,lwork)
     real(8),target     :: wrk(:)
     integer            :: lwork, ltempwork
     logical            :: ldummy1,ldummy2
     len_f77_work = lwork
     f77_work     => wrk
  END SUBROUTINE

  SUBROUTINE select_wrk(wrk, mem_req)
     real(8),pointer     :: wrk(:)
     integer, intent(in) :: mem_req

     wrk => f77_work

     memory_usage  = mem_req
     mem_wrk_count = mem_wrk_count + 1
     if (mem_wrk_count > 1) call quit('counter > 1 not implemented')
  end subroutine

  SUBROUTINE deselect_wrk(wrk, mem_req)
     real(8),pointer     :: wrk(:)
     integer, intent(in) :: mem_req
     nullify(wrk)
     memory_usage  = 0
     mem_wrk_count = mem_wrk_count - 1
  end subroutine

! ******************************************************************************

end module
