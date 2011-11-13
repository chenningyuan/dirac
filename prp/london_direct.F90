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

module london_direct_hb

  use london_prp_cfg
  use memory_allocator 

  use dirac_cfg
  use dft_cfg
  use xc_driver


  implicit none

  public london_twoel_hb

  save

  private

  real(8), pointer :: f77_work(:)       ! Fortran 77 work array pointer
  integer          :: len_f77_work  = 0 ! Length of the work array
  integer          :: memory_usage  = 0 ! Amount of work array used
  integer          :: mem_wrk_count = 0 ! Counts

  integer :: mat_dim_ao
  integer :: mat_dim_mo
  logical :: get_nmr_integrals_ll
  logical :: get_nmr_integrals_ls
  logical :: get_nmr_integrals_ss

  integer, parameter :: max_nr_xyz_comp = 3
  integer :: nr_shells_temp
  logical :: only_large_shells
  integer :: isym_2el_lao(max_nr_xyz_comp)
  integer :: irep_2el_lao(max_nr_xyz_comp)
  integer :: iherm_2el_lao(max_nr_xyz_comp)
  integer :: ifckop_2el_lao(max_nr_xyz_comp)
  integer :: jbtof_2el_lao(max_nr_xyz_comp)

  integer :: ifock_type(4, max_nr_xyz_comp)
  integer :: irepdm(4, max_nr_xyz_comp)

! there should be nr_mat=nzc1 (which is 4 except 2component spinfree... see dirac/dirrdn.F, grep for NZC1)
  integer, parameter :: nr_mat = 4

contains


  subroutine london_twoel_hb(cmo, ibeig, work, lwork)
!-----------------------------------------------------------------------
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)
    integer, intent(in)    :: lwork
    real(8), intent(inout) :: work(lwork)

    real(8), allocatable   :: fmat(:)
    integer                :: fmat_size

!-----------------------------------------------------------------------

!   set interface for f77 work array
!   --------------------------------
    call set_f77_work(work, lwork)

!   set variables used in this module
!   ---------------------------------
    call initialize_direct_lao()

    !fmat_size = mat_dim_ao*mat_dim_ao*nr_qu*nr_xyz_comp
    fmat_size = mat_dim_ao*mat_dim_ao*nr_mat*nr_xyz_comp
    call alloc(fmat, fmat_size)

!   twolo1a:
!    call bla1(fmat, cmo, ibeig)
    call get_ao_fock_2el_lao(fmat, cmo, ibeig)

!   rsklond
    call bla2(fmat, cmo, ibeig)
!    call symmetrize_antisymmetrize_fock_2el_lao(fmat)

!   if symmetry:
!   transform to SA-AO basis
!    call build_fock_ao(fmat, cmo, ibeig)

!   twolo1b
    call bla3(fmat, cmo, ibeig)

!   add xc contribution if dft calculations
!   ---------------------------------------
    if (dirac_cfg_dft_calculation) then
      call xc_contrib_2el_lao(fmat, cmo)
    end if


!   twolo1c
!    call bla5(fmat, cmo, ibeig)
    call ao_to_mo_2el_lao(fmat_size, fmat, cmo, ibeig)


    call dealloc(fmat)

  end subroutine london_twoel_hb








  subroutine bla1(fmat, cmo, ibeig)
!   --------------------------------------------------------------------
    real(8), intent(inout) :: fmat(*)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork, kfree

    real(8), allocatable :: dmat(:)
    integer :: irep(nr_qu*nr_xyz_comp)
    integer :: ifct(nr_qu*nr_xyz_comp)
    integer :: dmat_size
!   --------------------------------------------------------------------

    dmat_size = p_n2bbasxq 
    call alloc(dmat, dmat_size)

    lwork = len_f77_work
    kfree = 1
    call select_wrk(work, lwork)

    call twolo1a(fmat, dmat, cmo, ibeig, nr_qu, irep, ifct, &
                 work, lwork, prp_iprint)

    call deselect_wrk(work, lwork)

    call dealloc(dmat)

  end subroutine bla1


  subroutine bla2(fmat, cmo, ibeig)
!   --------------------------------------------------------------------
    !real(8), intent(inout) :: fmat(*)
    real(8), intent(inout) :: fmat(p_n2bbasx, nr_qu, *)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork, kfree

    real(8), allocatable :: skeleton(:)
    real(8), allocatable :: fac(:)
    real(8), allocatable :: kindex(:)
    real(8), allocatable :: point(:)

    integer :: i2typ
!   --------------------------------------------------------------------

    call alloc(skeleton, p_n2bbasxq)
    call alloc(fac     , mat_dim_ao*(p_maxrep + 1))
    call alloc(kindex  , (p_kmax*(p_maxrep + 1)*p_khk(p_nhtyp) + 1)/p_irat)
    call alloc(point   , (mat_dim_ao*(p_maxrep + 1) + 1)/p_irat)

    !call rsklond(fmat, skeleton, fac, kindex, point, nr_qu, irep_2el_lao, prp_iprint)
    call rsklond(fmat, skeleton, fac, kindex, point, nr_qu, irepdm, prp_iprint)

    call dealloc(skeleton)
    call dealloc(fac     )
    call dealloc(kindex  )
    call dealloc(point   )


        open(77, file='after_rsklond', status='unknown',  & 
             form='formatted', access = 'append')

        do i2typ = 1, 3
          write(77, *) 'gosia: i2typ, p_isymax(i2typ,2) = ', & 
                               i2typ, p_isymax(i2typ,2)
          call prqmat(fmat(1, 1, i2typ), & 
                      mat_dim_ao, mat_dim_ao, mat_dim_ao, mat_dim_ao, &
                      nr_qu, p_ipqtoq(1, irep_2el_lao(i2typ)), 77)
        end do
        close(77)



  end subroutine bla2



  subroutine bla3(fmat, cmo, ibeig)
!   --------------------------------------------------------------------
    !real(8), intent(inout) :: fmat(p_n2bbasx, nr_qu, *)
    real(8), intent(inout) :: fmat(p_n2bbasx, 4, *)
    real(8), allocatable   :: fock(:, :)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork, kfree

    integer :: icomp, iq, iz, ipq, irepd, ih, irep, ipari, itim
!   --------------------------------------------------------------------

        open(66, file='twolo1b', status='unknown', form='formatted', &
             access = 'append')



    call alloc(fock, p_n2bbasx, nr_qu)

    do icomp = 1, nr_xyz_comp

      call dzero(fock, p_n2bbasxq)
      irep  = p_isymax(icomp, 2)
      ipari = p_jbtof(irep, 1)

      write(p_lupri, *) 'twolo1b, icomp = ', icomp
      write(p_lupri, *) 'twolo1b, irep  = ', irep
      write(p_lupri, *) 'twolo1b, ipari = ', ipari

      !do iz = 1, nr_qu
      do iz = 1, 4
        irepd = p_irqmat(iz, irep)
        iq    = p_iqmult(1, p_jqbas(irepd, ipari), iz)
        ipq   = p_iqtopq(iq, irep)

        write(p_lupri, *) 'twolo1b, iz    = ', iz
        write(p_lupri, *) 'twolo1b, irepd = ', irepd
        write(p_lupri, *) 'twolo1b, iq    = ', iq
        write(p_lupri, *) 'twolo1b, ipq   = ', ipq

        call mtaoso(fmat(1, iz, icomp), fock(1, ipq), nr_ao_total, irepd, prp_iprint)

        !irepd = p_irqmat(iz, irep_2el_lao(icomp))
        !iq    = p_iqmult(1, p_jqbas(irepd, jbtof_2el_lao(icomp)), iz)
        !ipq   = p_iqtopq(iq, irep_2el_lao(icomp))
        !call mtaoso(fmat(1, iz, icomp), fock(1, ipq), nr_ao_total, irepd, prp_iprint)
      end do


!     insert half-phases:
      if (nr_qu .lt. 4) then
        do iz = 1, nr_qu
          !iq   = p_ipqtoq(iz, irep_2el_lao(icomp))
          iq   = p_ipqtoq(iz, irep)
          call q2phase('D', iq, 1, fock(1, iz))
        end do
      end if

      do iz = 1, nr_qu
        call dsltsp(mat_dim_ao, fock(1,iz), fmat(1,iz,icomp))
      end do

      call dzero(fock, p_n2bbasxq)
      itim = -1
      do iz = 1, nr_qu
        !iq = p_ipqtoq(iz, irep_2el_lao(icomp))
        !ih = p_ihqmat(iq, iherm_2el_lao(icomp))
        iq = p_ipqtoq(iz, irep)
        ih = p_ihqmat(iq, itim)
        if (ih .eq. 1) then
          call dsptsi(mat_dim_ao, fmat(1, iz, icomp), fock(1,iz))
        else if(ih .eq. 2) then
          call daptge(mat_dim_ao, fmat(1, iz, icomp), fock(1,iz))
        else
          call quit('twolo1: not defined ih!')
        endif
      enddo

      call butobs_no_work(fock, nr_qu)

      call dcopy(p_n2bbasxq, fock, 1, fmat(1, 1, icomp), 1)

        write(66, *) 'gosia: icomp = ', icomp
        call prqmat(fmat(1, 1, icomp),  &
                    mat_dim_ao, mat_dim_ao, mat_dim_ao, mat_dim_ao, &
                    !nr_qu, p_ipqtoq(1, irep_2el_lao(icomp)), 66)
                    nr_qu, p_ipqtoq(1, irep), 66)
!
    end do

    call dealloc(fock)
     close(66)

  end subroutine bla3


  subroutine bla4(fmat, cmo, ibeig)
!   --------------------------------------------------------------------
    real(8), intent(inout) :: fmat(*)
    real(8), intent(in)  :: cmo(*)
    integer, intent(in)  :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork, kfree

    real(8), allocatable :: dmat0(:)
    integer :: dmat0_size
!   --------------------------------------------------------------------


!   add xc contribution if dft calculations
!   ---------------------------------------
    if (dirac_cfg_dft_calculation) then
!      call xc_contrib_2el_lao(fmat, cmo)
      dmat0_size = p_n2bbasxq
      call alloc(dmat0, dmat0_size)

      lwork = len_f77_work
      kfree = 1
      call select_wrk(work, lwork)
      call add_dft_to_twolon(fmat, dmat0, cmo, &
                             work, lwork, prp_iprint)
      call deselect_wrk(work, lwork)
      call dealloc(dmat0)

    end if

  end subroutine bla4



  subroutine bla5(fmat, cmo, ibeig)
!   --------------------------------------------------------------------
    real(8), intent(inout) :: fmat(*)
    real(8), intent(in)  :: cmo(*)
    integer, intent(in)  :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork, kfree

    real(8), allocatable :: gbmo(:)
    integer :: gbmo_size
!   --------------------------------------------------------------------

    gbmo_size = p_n2orbxq 
    call alloc(gbmo, gbmo_size)

    lwork = len_f77_work
    call select_wrk(work, lwork)

    call twolo1c(gbmo, fmat, cmo, nr_xyz_comp, work, lwork, prp_iprint)

    call deselect_wrk(work, lwork)
    call dealloc(gbmo)

  end subroutine bla5







  subroutine initialize_direct_lao()
!   --------------------------------------------------------------------
    integer :: i, j
    logical :: print_info
!   --------------------------------------------------------------------

!   read basic variables from basis_info:
!   for a moment - no parallel calculations
!    call basis_info_read(.true., .false.)

!   set other variables necessary for london_direct_hb
    call delete_direct_info()

!   variables from common blocks
    call common_blocks_vars()


!   matrix dimensions
    mat_dim_ao = nr_ao_total
    mat_dim_mo = p_norbt


!   which type of integrals we include (LL, LS, SS)
    if (btest(nmr_int, 0)) get_nmr_integrals_ll = .true.
    if (btest(nmr_int, 1)) get_nmr_integrals_ls = .true.
    if (btest(nmr_int, 2)) get_nmr_integrals_ss = .true.


!   symmetry info for 2-el LAO integrals
    do i = 1, nr_xyz_comp

!     old code: isymop
      isym_2el_lao(i)  = p_isymax(i, 2) + 1
      irep_2el_lao(i)  = p_isymax(i, 2)

!     old code: ihrmop
!               0 -> no symmetry
!               1 -> hermitian
!              -1 -> antihermitian
      iherm_2el_lao(i) = -1

!     old code: ifckop:
!               1 -> only Coulomb
!               2 -> only Exchange
!               3 -> Coulomb+Exchange
      ifckop_2el_lao(i) = 1

      jbtof_2el_lao(i)  = p_jbtof(irep_2el_lao(i), 1)
    end do


    !if (levyle .or. x2c) then
    if (p_levyle) then
      only_large_shells = .true.
      nr_shells_temp    = p_nlrgsh
    else
      nr_shells_temp    = p_kmax
    end if



    print_info = .true.

    if (print_info) then
      write(p_lupri, *) 'info from initialize_direct_lao: '
      write(p_lupri, *) 'x ------------------------------------- x'
      write(p_lupri, *) 'mat_dim_ao  = ', mat_dim_ao  
      write(p_lupri, *) 'mat_dim_mo  = ', mat_dim_mo 
      write(p_lupri, *) 'nr_xyz_comp = ', nr_xyz_comp
      do i = 1, nr_xyz_comp
        write(p_lupri, *) 'i, isym_2el_lao(i)   = ', i, isym_2el_lao(i)
        write(p_lupri, *) 'i, irep_2el_lao(i)   = ', i, irep_2el_lao(i)
        write(p_lupri, *) 'i, iherm_2el_lao(i)  = ', i, iherm_2el_lao(i)
        write(p_lupri, *) 'i, ifckop_2el_lao(i) = ', i, ifckop_2el_lao(i)
        write(p_lupri, *) 'i, jbtof_2el_lao(i)  = ', i, jbtof_2el_lao(i)
      end do
      write(p_lupri, *) 'nr_shells_temp       = ', nr_shells_temp
      write(p_lupri, *) 'get_nmr_integrals_ll = ', get_nmr_integrals_ll
      write(p_lupri, *) 'get_nmr_integrals_ls = ', get_nmr_integrals_ls
      write(p_lupri, *) 'get_nmr_integrals_ss = ', get_nmr_integrals_ss
      write(p_lupri, *) 'x ------------------------------------- x'
    end if

  end subroutine initialize_direct_lao



  subroutine delete_direct_info()
!   --------------------------------------------------------------------
!   --------------------------------------------------------------------

    get_nmr_integrals_ll = .false.
    get_nmr_integrals_ls = .false.
    get_nmr_integrals_ss = .false.

    isym_2el_lao   = 0
    irep_2el_lao   = 0
    iherm_2el_lao  = 0
    ifckop_2el_lao = 0
    jbtof_2el_lao  = 0

    only_large_shells = .false.
    nr_shells_temp    = 0

!
    irepdm = 0
    ifock_type = 0

  end subroutine delete_direct_info



  subroutine get_ao_fock_2el_lao(fmat, cmo, ibeig)
!   --------------------------------------------------------------------
    !real(8), intent(out)   :: fmat(*)
    real(8), intent(out)   :: fmat(p_n2bbasx, nr_mat, *)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork, kfree

    real(8), allocatable   :: dmat(:, :), temp(:, :)
    integer :: i2typ
    integer :: itype, max_deriv, iatom
    logical :: nodv, nopv, nocont, retur, relcal, ddfock, fckddr, &
               london, first
    integer :: dmat_size

    logical :: debug_me
    integer :: i1(nr_xyz_comp), i2(nr_xyz_comp), i3(nr_xyz_comp), i
    integer :: iq, ipq, iz, irepd

!   gosia: on output we get fmat
!   - nonsymmetric, square (ntbas(0) x ntbas(0))
!   - Hermit sorted ('boson unsorted') 
!   - in SA-AO basis
!   - unpacked (nr_qu matrices)
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

!   set ifctyp and irepdm
!   ---------------------
    do i = 1, nr_xyz_comp
      i1(i) = 1 
      i2(i) = 1 
      i3(i) = 1 
    end do

    call setfck(ifock_type,     &
                irepdm,         &
                nr_xyz_comp,    &
                nr_qu,          &
                i1,             &
                i2,             &
                i3,             &
                prp_iprint)



    call alloc(dmat, p_n2bbasx, nr_qu)
    call alloc(temp, p_n2bbasx, nr_qu)

!   generate (unperturbed) density matrix
!    lwork = len_f77_work
!    call select_wrk(work, lwork)
!    call gtcnvd(dmat, cmo, ibeig, prp_iprint, work, lwork)
!    call deselect_wrk(work, lwork)

!   ----------------------------------------------------
    call genden(temp, cmo, 1, prp_iprint)
    call dscal(p_n2bbasxq, 2.0d0, temp, 1)



    call bstobu_no_work(temp, nr_qu)


    if(nr_qu .lt. 4) then
      do iz = 1, nr_qu
        iq = p_ipqtoq(iz, 0)
        call q2phase('D', iq, 1, temp(1, iz))
      enddo
    endif


    do iz = 1, nr_qu
      irepd = p_irqmat(iz, 0)
      iq    = p_iqmult(1, p_jqbas(irepd, p_jbtof(0, 1)), iz)
      ipq   = p_iqtopq(iq, 0)
      if(ipq.gt.nr_qu) call quit('gtcnvd.density: ipq.gt.nz !')
      if(ipq.le.0 ) call quit('gtcnvd.density: ipq.le.0 !')

      call dtsoao(temp(1, ipq), dmat(1, iz), nr_ao_total, irepd, prp_iprint)
    enddo

    write(p_lupri, *) 'density matrix before twoint' 
    call prqmat(dmat(1, 1),                               & 
                mat_dim_ao, mat_dim_ao, mat_dim_ao, mat_dim_ao,  &
                nr_qu, p_ipqtoq(1, 0), p_lupri)



    write(p_lupri, *) 'new twolo1a, get_nmr_integrals_ll =', get_nmr_integrals_ll
    write(p_lupri, *) 'new twolo1a, get_nmr_integrals_ls =', get_nmr_integrals_ls
    write(p_lupri, *) 'new twolo1a, get_nmr_integrals_ss =', get_nmr_integrals_ss



    open(88, file='twolo1a', status='unknown', form='formatted', &
         access = 'append')



!   get 2-electron London integrals
!   -------------------------------
    if (get_nmr_integrals_ll) then
      i2typ = 1

      lwork = len_f77_work
      call select_wrk(work, lwork)

      call twoint(work, lwork, fmat, dmat, nr_qu, irepdm, ifock_type, &
                  (/0.0d0/), (/0/), (/0/), itype, max_deriv, iatom,   &
                  nodv, nopv, nocont, p_tktime, prp_iprint, p_iprnta, &
                  p_iprntb, p_iprntc, p_iprntd, retur, (/0/), i2typ,  &
                  (/0/), (/0.0d0/), (/0.0d0/), (/0.0d0/), (/0.0d0/),  &
                  (/0.0d0/), relcal, .false., (/0/), (/0.0d0/))

      call deselect_wrk(work, lwork)

        write(88, *) 'gosia: i2typ = ', i2typ
        call prqmat(fmat(1, 1, i2typ),                               & 
                    mat_dim_ao, mat_dim_ao, mat_dim_ao, mat_dim_ao,  &
                    nr_qu, p_ipqtoq(1, 0), 88)


    end if

    if (get_nmr_integrals_ls) then
      i2typ = 2

      lwork = len_f77_work
      call select_wrk(work, lwork)

      call twoint(work, lwork, fmat, dmat, nr_qu, irepdm, ifock_type, &
                  (/0.0d0/), (/0/), (/0/), itype, max_deriv, iatom,   &
                  nodv, nopv, nocont, p_tktime, prp_iprint, p_iprnta, &
                  p_iprntb, p_iprntc, p_iprntd, retur, (/0/), i2typ,  &
                  (/0/), (/0.0d0/), (/0.0d0/), (/0.0d0/), (/0.0d0/),  &
                  (/0.0d0/), relcal, .false., (/0/), (/0.0d0/))

      call deselect_wrk(work, lwork)

        write(88, *) 'gosia: i2typ = ', i2typ
        call prqmat(fmat(1, 1, i2typ),                               & 
                    mat_dim_ao, mat_dim_ao, mat_dim_ao, mat_dim_ao,  &
                    nr_qu, p_ipqtoq(1, 0), 88)


    end if

    if (get_nmr_integrals_ss) then
      i2typ = 3

      lwork = len_f77_work
      call select_wrk(work, lwork)

      call twoint(work, lwork, fmat, dmat, nr_qu, irepdm, ifock_type, &
                  (/0.0d0/), (/0/), (/0/), itype, max_deriv, iatom,   &
                  nodv, nopv, nocont, p_tktime, prp_iprint, p_iprnta, &
                  p_iprntb, p_iprntc, p_iprntd, retur, (/0/), i2typ,  &
                  (/0/), (/0.0d0/), (/0.0d0/), (/0.0d0/), (/0.0d0/),  &
                  (/0.0d0/), relcal, .false., (/0/), (/0.0d0/))

      call deselect_wrk(work, lwork)

        write(88, *) 'gosia: i2typ = ', i2typ
        call prqmat(fmat(1, 1, i2typ),                               & 
                    mat_dim_ao, mat_dim_ao, mat_dim_ao, mat_dim_ao,  &
                    nr_qu, p_ipqtoq(1, 0), 88)


    end if

    call dealloc(dmat)
    call dealloc(temp)

    close(88)


  end subroutine get_ao_fock_2el_lao


  subroutine symmetrize_antisymmetrize_fock_2el_lao(fmat)
!   --------------------------------------------------------------------
    real(8), intent(inout) :: fmat(mat_dim_ao, mat_dim_ao, nr_qu, *)

    real(8) :: f
    integer :: h(nr_qu), i, j, k, l
!   --------------------------------------------------------------------

!   antisymmetrize/symmetrize fmat 
!   (Areal: antisymmetrize, Aim, Breal, Bim: symmetrize, if present) 

    h(1) = -1
    if (nr_qu .gt. 1) then
      do i = 2, nr_qu
        h(i) = 1
      end do
    end if

    do i = 1, nr_xyz_comp
      do j = 1, nr_qu
        do k = 1, mat_dim_ao
          do l = 1, k
            f = 0.5d0*(fmat(k, l, j, i) + h(j)*fmat(l, k, j, i))
            fmat(k, l, j, i) =      f
            fmat(l, k, j, i) = h(j)*f
          end do
        end do
      end do
    end do


  end subroutine symmetrize_antisymmetrize_fock_2el_lao




  subroutine build_fock_ao(fmat, cmo, ibeig)
!   --------------------------------------------------------------------
    real(8), intent(inout) :: fmat(mat_dim_ao, mat_dim_ao, nr_qu, *)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)


    real(8), allocatable :: skeleton(:, :)
    real(8), allocatable :: factor(:, :)
    real(8), allocatable :: kindex(:, :, :)
    real(8), allocatable :: ipoint(:, :)
!   --------------------------------------------------------------------

!   1. part directly from old rsklond:
    if (p_maxrep .gt. 0) then

      call alloc(skeleton, mat_dim_ao, mat_dim_ao)
      call alloc(factor  , mat_dim_ao, p_maxrep + 1)
      call alloc(kindex  , p_kmax, p_maxrep + 1, (p_khk(p_nhtyp) + 1)/p_irat)
      call alloc(ipoint  , mat_dim_ao, ((p_maxrep + 1) + 1)/p_irat)







      call dealloc(skeleton)
      call dealloc(factor  )
      call dealloc(kindex  )
      call dealloc(ipoint  )

    end if



!g!!!!!!!!!!!!!!!!!!!
!g!   1. part from old rsklond:
!g!   --------------------------
!g    if (maxrep .gt. 0) then
!g
!g!     call alloc(skeleton, p_n2bbasxq)
!g      call alloc(skeleton, mat_dim_ao, mat_dim_ao)
!g      !call alloc(fac     , mat_dim_ao*(p_maxrep + 1))
!g      call alloc(fac     , mat_dim_ao, p_maxrep + 1)
!g      !call alloc(kindex  , (p_kmax*(p_maxrep + 1)*p_khk(p_nhtyp) + 1)/p_irat)
!g      call alloc(kindex  , p_kmax, p_maxrep + 1, (p_khk(p_nhtyp) + 1)/p_irat)
!g      !call alloc(point   , (mat_dim_ao*(p_maxrep + 1) + 1)/p_irat)
!g      call alloc(point   , mat_dim_ao, ((p_maxrep + 1) + 1)/p_irat)
!g
!g!     set kindex:
!g      iorb = 0
!g      do ishell = 1, nr_shells_temp
!g        do irep = 0, maxrep
!g          if (iand(irep, istbao(ishell)) == 0) then
!g!           symmetry independent shells
!g            do ideg = 1, p_khkt(ishell)
!g              iorb = iorb + 1
!g              kindex(ishell, irep, ideg) = iorb
!g            end do
!g          else
!g!           shells 'copied' by symmetry operations
!g            isymf = ibtxor(irep, ibtand(irep, istbao(ishell)))
!g            do ideg = 1, p_khkt(ishell)
!g              kindex(ishell, irep, ideg) = kindex(ishell, isymf, ideg)
!g            end do
!g          end if
!g        end do 
!g      end do
!g
!g!   set point and factor
!g      iorb = 0
!g      do ishell = 1, nr_shells_temp
!g        nhkta = p_nhkt(ishell)
!g        do irep = 0, maxrep
!g          if (iand(irep, istbao(ishell)) == 0) then
!g            do ideg = 1, p_khkt(ishell)
!g              iorb = iorb + 1
!g              do irep2 = 0, maxrep
!g                isympr = ibtxor(irep, irep2)
!g                point(iorb, irep2) = index(ishell, isympr, ideg)
!g                if (ishell .gt. p_nlrgsh) then
!g                  factor(iorb, irep2) = p_ipar(irep2)*pt(iand(irep2, isymao(nhkta, ideg)))
!g                else
!g                  factor(iorb, irep2) = pt(iand(irep2, isymao(nhkta, ideg)))
!g                end if
!g              end do
!g            end do
!g          end if
!g        end do
!g      end do
!g
!g
!g      do n = 1, ndmat
!g        do  x = 1, 3
!g          call dcopy(mat_dim_ao*mat_dim_ao, fmat(1, 1, n, x), 1, sklton, 1)
!g          do isymop = 1, maxrep
!g            fp = pt(iand(isymop, ibtxor(isymdm(n), p_isymax(x, 2))))
!g            do a = 1, mat_dim_ao
!g              ap = ipoint(a,isymop)
!g              fa = factor(a,isymop)
!g              do b = 1, mat_dim_ao
!g                bp = ipoint(b,isymop)
!g                fb = factor(b,isymop)
!g                fmat(ap,bp,n,x)=fmat(ap,bp,n,x) +fa*fb*fp*sklton(a,b)
!g              end do    
!g            end do   
!g          end do
!g        end do
!g      end do
!g
!g         fac = d1/sdpre(maxrep + 1)
!g         nn = ndmat*3*nbasis*nbasis
!g         call dscal(nn,fac,fmat,1)
!g
!g
!g
!g    end if ! if (maxrep .gt. 0)

! the part that was in old twolo1b:

!
!    do i = 1, nr_xyz_comp
!
!!     transform AO to SO
!!     call mtaoso
!
!!     insert half-phases
!      if (nr_qu .lt. 4) then
!        do j = 1, nr_qu
!          iq = p_ipqtoq(j, irep_2el_lao(i))
!          call q2phase()
!        end do
!      end if
!
!    end do


  end subroutine build_fock_ao




  subroutine xc_contrib_2el_lao(fmat, cmo)
!   --------------------------------------------------------------------
    real(8), intent(inout) :: fmat(*)
    real(8), intent(in)    :: cmo(*)

    real(8), allocatable   :: dmat(:)
    integer :: dmat_size, i
!   --------------------------------------------------------------------

!   get unperturbed density matrix to dmat
!   --------------------------------------
    dmat_size = mat_dim_ao*mat_dim_ao*nr_qu
    call alloc(dmat, dmat_size)
    call genden(dmat, cmo, 1, prp_iprint)

!   add xc contribution
!   gosia fixme: it won't work now with visual+london now

    call xc_mpi_wake_up_coworkers()
    if (dft_cfg_no_sdft) then
      call xc_integrator(xc_mat_dim              = mat_dim_ao,   &
                         xc_nz                   = nr_qu,        &
                         xc_dmat_0               = dmat,         &
                         xc_nr_dmat              = 0,            &
                         xc_nr_fmat              = nr_xyz_comp,  &
                         xc_fmat                 = fmat,         &
                         xc_fmat_pg_sym          = isym_2el_lao, &
                         xc_do_london_rhs_direct = .true.)
    else
      call xc_integrator(xc_mat_dim           = mat_dim_ao,      &
                         xc_nz                = nr_qu,           &
                         xc_dmat_0            = dmat,            &
                         xc_nr_dmat           = 0,               &
                         xc_nr_fmat           = nr_xyz_comp,     &
                         xc_fmat              = fmat,            &
                         xc_fmat_pg_sym       = isym_2el_lao,    &
                         xc_response_order_mo = 1,               &
                         xc_do_london_rhs_direct = .true.)
    end if
    call xc_mpi_release_coworkers()

  end subroutine xc_contrib_2el_lao








  subroutine ao_to_mo_2el_lao(fmat_size, fmat, cmo, ibeig)
!   --------------------------------------------------------------------
    integer, intent(in) :: fmat_size
    !real(8), intent(in) :: fmat(fmat_size, *)
    real(8), intent(in) :: fmat(nr_ao_total*nr_ao_total, nr_qu, *)
    real(8), intent(in) :: cmo(*)
    integer, intent(in) :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork

    real(8), allocatable :: gmat(:)
    integer :: gmat_size, icomp, ii, jj
!   --------------------------------------------------------------------

    !gmat_size = mat_dim_mo*mat_dim_mo*nr_qu
    gmat_size = p_n2orbxq
    call alloc(gmat, gmat_size)


!    lwork = len_f77_work
!    call select_wrk(work, lwork)
!    call twolo1c(gmat, fmat, cmo, nr_xyz_comp, work, lwork, prp_iprint)
!    call deselect_wrk(work, lwork)

    open(twolon_unit, file = 'TWOLON', form = 'unformatted',        &
         access = 'direct', recl = 8*gmat_size, status = 'unknown')

    do icomp = 1, nr_xyz_comp

      call dzero(gmat, gmat_size)
     
      do ii = 1, p_nfsym
        jj = mod(ii + jbtof_2el_lao(icomp), 2) + 1

        lwork = len_f77_work
        call select_wrk(work, lwork)


        call qtrans('AOMO',                            &
                    'S',                               &
                    0.0d0,                             &
                    p_nfbas(ii, 0),                    &
                    p_nfbas(jj, 0),                    &
                    p_norb(ii),                        &
                    p_norb(jj),                        &
                    !fmat(p_i2basx(ii, jj) + 1, icomp), &
                    fmat(p_i2basx(ii, jj) + 1, 1, icomp), &
                    nr_ao_total,                       &
                    nr_ao_total,                       &
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


  end subroutine ao_to_mo_2el_lao












  SUBROUTINE set_f77_work(wrk,lwork)
     real(8),target :: wrk(:)
     integer            :: lwork, ltempwork
     logical            :: ldummy1,ldummy2
     len_f77_work = lwork
     f77_work     => wrk
  END SUBROUTINE set_f77_work

  SUBROUTINE select_wrk(wrk, mem_req)
     real(8),pointer :: wrk(:)
     integer, intent(in) :: mem_req

        wrk => f77_work

        memory_usage  = mem_req
        mem_wrk_count = mem_wrk_count + 1
        if (mem_wrk_count > 1) call quit('counter > 1 not implemented')
  end subroutine

  SUBROUTINE deselect_wrk(wrk, mem_req)
     real(8),pointer :: wrk(:)
     integer, intent(in) :: mem_req
        nullify(wrk)
        memory_usage  = 0
        mem_wrk_count = mem_wrk_count - 1
  end subroutine




end module
