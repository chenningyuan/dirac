module london_response

  use interface_for_pammag
  use pammag_info
  use memory_allocator 

  use dirac_cfg
#ifdef MOD_DFT
  use dft_cfg
  use dft_main
#endif

  implicit none

  public rntlon_new

  save
 
! copied from aoresponse
  private
  real(8), pointer :: f77_work(:)       ! Fortran 77 work array pointer
  integer          :: len_f77_work  = 0 ! Length of the work array
  integer          :: memory_usage  = 0 ! Amount of work array used
  integer          :: mem_wrk_count = 0 ! Counts
 
contains

  subroutine rntlon_new(cmo, ibeig, work, lwork)
!-----------------------------------------------------------------------
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    integer, intent(in)    :: lwork
    real(8), intent(inout) :: work(lwork)

    real(8), allocatable   :: mat1(:), mat2(:)
    real(8), allocatable   :: tb(:), buf(:)
    integer, allocatable   :: npos_array(:)


! what we do now with memory:
!============================
! 1. allocate 3*t_n2bbasxq for mat1 and 3*t_n2bbasxq for mat2
! 2. inside get_tb_dbmat:
!    *) loop: do icomp =1, 3
!    *) get T^B_icomp -> use 1*t_n2orbxq and store it in mat1(1, icomp)
!    *) call getmdm and get tb_dmat -> use 1*t_n2bbasxq (for 1 icomp) and write tb_dmat to mat2(1, icomp)
! 3. zero out mat1
! 4. call twofck - result write to mat1 (used all 3*t_n2bbasxq)
!    *) add dft
! 5. zero out mat2
! 6. do AO-to-MO transformation - result written to mat2 (it is "G" - uses 3*t_n2orbxq, we have to keep it)
! 7. zero out mat1 
! 8. get Fock in MO basis - store it in mat1
! 9(!) allocate some additional memory: t_n2orbxq for T^B_icomp (tb) and t_n2orbxq for buf
! 10. and in a do icomp = 1, 3 loop do:
!    *) get T^B_icomp to tb
!    *) call ttra (needs tb and mat1)- store result in buf
!    *) add buf+mat2(1, icomp) and write to mat2(1, icomp)
!    *) zero out buf
!
!-----------------------------------------------------------------------

! set interface for f77 work array:
    call set_f77_work(work, lwork)

! set common blocks variables used (interface_for_pammag.F90):
    call common_blocks_vars()

! allocate memory:
    call alloc(mat1, 3*t_n2bbasxq)
    call alloc(mat2, 3*t_n2bbasxq)


! get (C+)(T^B)C into mat2 (use mat1 for T^B):
    call get_tb_dmat(mat1, mat2, cmo, ibeig)

! calculate "G" 
    call get_gtb(mat1, mat2, cmo, ibeig)

! calculate "F"
    call dzero(mat1, 3*t_n2bbasxq)
    call get_ftb(mat1, mat2, cmo, ibeig)

! dealloc:
    call dealloc(mat1)
    call dealloc(mat2)

  end subroutine


  subroutine get_ftb(mat1, mat2, cmo, ibeig)
!-----------------------------------------------------------------------
    real(8) :: mat1(*)
    !real(8), intent(inout) :: mat1(*)
    real(8), intent(inout) :: mat2(*)

    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork

    real(8), allocatable   :: buf(:)
    real(8), allocatable   :: tb(:)

    integer :: isym(nr_fock), iopsy(nr_fock)
    integer :: isymop(nr_fock), ihrmop(nr_fock), ifckop(nr_fock)
    integer :: indxpr
    integer :: icomp, i, j, isymh, kfree
    logical :: ah
!-----------------------------------------------------------------------

!get unperturbed Fock matrix in MO basis:
!todo: reduce memory and get it from DFFCK1 and DFFCK2 files:
! (but there we have Fock in AO basis...)

    !call dzero(mat1, 3*t_n2bbasxq)

    lwork = len_f77_work
    kfree = 1
    call select_wrk(work, lwork)

    call getfck(mat1,             & !fmo (t_n2bbasxq)
                print_prp_level,  &
                work,             &
                kfree,            &
                lwork)

    call deselect_wrk(work, lwork)

! allocate some more memory...
    call alloc(tb,  t_n2orbxq)
    call alloc(buf, t_n2orbxq)

! we need AOPROPER and FCKLON files:

    call opnfil(int1_unit, 'AOPROPER', 'OLD', 'get_ftb')

    open(fcklon_unit, file='FCKLON', form='unformatted',   &
         access='direct',recl=8*t_n2orbxq,                 &
         status='unknown')

    do icomp = 1, 3
! get T^B to tb
      indxpr       = ip_conmat(icomp)
      isym(icomp)  = prp_boson_sym(indxpr)
      iopsy(icomp) = t_jbtof(isym(icomp)-1, 1)

      call dzero(tb, t_n2orbxq)

      lwork = len_f77_work
      kfree = 1
      call select_wrk(work, lwork)

      call prpmat(indxpr,          &
                  iopsy(icomp),    &
                  tb,              & ! T^B matrix
                  cmo,             &
                  ibeig,           &
                  cmoq_off,        &
                  t_norb,          &
                  work,            &
                  kfree,           &
                  lwork,           &
                  print_prp_level)

      call deselect_wrk(work, lwork)

! get final result to buf

      ah = .true.
      isymh = 1

      call dzero(buf, t_n2orbxq)

      call ttra(sum_ep_rn,          & !epmix
                ah,                 & !ah
                isym(icomp),        & !isymt
                isymh,              & !isymh
                tb,                 & !tmat,  n2orbxq
                mat1(1:t_n2bbasxq),               & !hmat,  n2bbasxq
                buf,                & !thmat, n2orbxq (for 1 icomp)
                print_prp_level)

! add buf to mat2 and store in mat2
      call daxpy(t_n2orbxq,         &
                 1.0d0,             &
                 buf,               &
                 1,                 & 
                 mat2(1 + (icomp - 1)*t_n2orbxq),    & 
                 1)

! finally write to file:

      call wrtdac(fcklon_unit, t_n2orbxq,   &
                  mat2(1 + (icomp - 1)*t_n2orbxq), icomp)

    end do ! do i = 1, 3

    close(fcklon_unit, status='keep')
    close(int1_unit,   status='keep')

    call dealloc(tb)
    call dealloc(buf)

  end subroutine


  subroutine get_gtb(mat1, mat2, cmo, ibeig)
!-----------------------------------------------------------------------
    real(8), intent(inout) :: mat1(*)
    real(8), intent(inout) :: mat2(*)

    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork

    integer, allocatable   :: npos_array(:)
    integer                :: npos
    
    real(8), allocatable :: dmat(:)
    integer :: isym(nr_fock), iopsy(nr_fock), irep(nr_fock)
    integer :: isymop(nr_fock), ihrmop(nr_fock), ifckop(nr_fock)
    integer :: indxpr
    integer :: i
!-----------------------------------------------------------------------

! set for twofck:
    do i = 1, 3
      indxpr    = ip_conmat(i)
      isym(i)   = prp_boson_sym(indxpr)
      iopsy(i)  = t_jbtof(isym(i) - 1, 1)
      irep(i)   = isym(i) - 1
      isymop(i) = isym(i)
      ihrmop(i) = -1
      ifckop(i) = t_jbtof(isym(i) - 1, 1)
    end do

! get Fock matrix in AO basis to mat1
    call stex_get_npos(npos)
    call alloc(npos_array, npos)

    call dzero(mat1, 3*t_n2bbasxq)

    lwork = len_f77_work
    call select_wrk(work, lwork)

    call twofck(isymop,              &
                ihrmop,              &
                ifckop,              &
                mat1,                & !fao
                mat2,                & !tb_dmat
                nr_fock,             & ! = 3
                npos_array,          &
                nmr_int,             &
                print_prp_level,     &
                work,                &
                lwork)

    call deselect_wrk(work, lwork)

    if (allocated(npos_array)) then
      call dealloc(npos_array)
    end if

#ifdef MOD_DFT
!       construct kohn-sham contribution
!       ================================
! (stolen from pamrvc.F)

    if (dirac_cfg_dft_calculation) then

!generate standard density matrix
!need some more memory...
      call alloc(dmat, t_n2bbasxq)
      call genden(dmat, cmo, 1, print_prp_level)
!fixme should be response_order=1
! tsaue+gosia: This contribution is wrong, in fact we get half the direct contribution !!!
      call dftexc(nr_ao_total,   &
                  nr_qu,         &
                  nr_fock,       &
                  mat1,          & !fmat
                  dmat,          & !dmat_0
                  nr_fock,       & !nr_dmat_p
                  mat2,          & !dmat_p
                  isymop,        &
                  ihrmop,        &
                  0,             & !response_order
                  0,             & !response_order_ao
                  0,             & !g_order
                  1)               !m_order

      call dealloc(dmat)

    end if
#endif /* MOD_DFT */


! do AO-to-MO transformation (and store in mat2):

    call dzero(mat2, 3*t_n2bbasxq)

    call do_ao_to_mo(mat2,        & ! result in MO basis (3*t_n2orbxq)
                     mat1,        & ! result in AO basis 
                     cmo,         &
                     ibeig)

  end subroutine


  subroutine do_ao_to_mo(fmo, fao, cmo, ibeig)
!-----------------------------------------------------------------------
    real(8), intent(inout) :: fmo(*)
    real(8), intent(in)    :: fao(*)

    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork

    integer :: itim(nr_fock), isym(nr_fock), iopsy(nr_fock), irep(nr_fock)
    integer :: isymop(nr_fock), ihrmop(nr_fock), ifckop(nr_fock)
    integer :: indxpr
    integer :: icomp, i, j
!-----------------------------------------------------------------------

!#define test_gmolitx2
!remove this - it gives wrong result...
#ifdef test_gmolitx2

    do i = 1, 3

      lwork = len_f77_work
      call select_wrk(work, lwork)

      call gmolitx2(fmo(1+(i-1)*t_n2orbxq),  &
                    fao(1+(i-1)*t_n2bbasxq),  &
                    cmo,        &
                    ibeig,      &
                    work,       &
                    lwork)

      call deselect_wrk(work, lwork)

    end do

#else 

    do icomp = 1, 3

      indxpr        = ip_conmat(icomp)
      isym(icomp)   = prp_boson_sym(indxpr)
      iopsy(icomp)  = t_jbtof(isym(icomp) - 1, 1)
      irep(icomp)   = isym(icomp) - 1

      do i = 1, nr_fermion_irr
        j = mod(i + iopsy(icomp), 2) + 1

!fixme: len_f77_work is too much, maybe try to evaluate: lwork=n2bbasxq*some_number?

        lwork = len_f77_work
        call select_wrk(work, lwork)

        call qtrans('AOMO',                               &
                    'S',                                  &
                    0.0d0,                                &
                    nfbas_fi(i, 0),                       &
                    nfbas_fi(j, 0),                       &
                    t_norb(i),                            &
                    t_norb(j),                            &
                    fao(1 + basx_off(i, j) + (icomp-1)*t_n2bbasxq), &
                    nr_ao_total,                           &
                    nr_ao_total,                           &
                    nr_qu,                                 &
                    prp_cb_pq_to_uq(1, irep(icomp)),       &
                    fmo(1 + orbx_off(i, j) + (icomp-1)*t_n2orbxq),  &
                    t_norbt,                               &
                    t_norbt,                               &
                    nr_qu,                                 &
                    prp_cb_pq_to_uq(1, irep(icomp)),       &
                    cmo(1 + cmoq_off(i)),                  &
                    nfbas_fi(i, 0),                        &
                    t_norb(i),                             &
                    nr_qu,                                 &
                    prp_cb_pq_to_uq(1, 0),                 &
                    cmo(1 + cmoq_off(j)),                  &
                    nfbas_fi(j, 0),                        & 
                    t_norb(j),                             &
                    nr_qu,                                 &
                    prp_cb_pq_to_uq(1, 0),                 &
                    work,                                  &
                    lwork,                                 &
                    print_prp_level)

        call deselect_wrk(work, lwork)

      end do
    end do

#endif /* test_gmolitx2 */

  end subroutine


  subroutine get_tb_dmat(tb, tb_dmat, cmo, ibeig)
!-----------------------------------------------------------------------
    real(8), intent(out)   :: tb_dmat(*)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)
    real(8), intent(in)    :: tb(t_n2orbxq)

    real(8), pointer       :: work(:)
    integer                :: lwork

    integer :: isym(nr_fock), iopsy(nr_fock)
    integer :: indxpr, icomp 
    integer :: kfree
!-----------------------------------------------------------------------

    call opnfil(int1_unit, 'AOPROPER', 'OLD', 'get_mod_den')

    do icomp = 1, 3
      indxpr       = ip_conmat(icomp)
      isym(icomp)  = prp_boson_sym(indxpr)
      iopsy(icomp) = t_jbtof(isym(icomp)-1, 1)

      call dzero(tb, t_n2orbxq)

      lwork = len_f77_work
      kfree = 1
      call select_wrk(work, lwork)

      call prpmat(indxpr,           &
                  iopsy(icomp),     &
                  tb,               &
                  cmo,              &
                  ibeig,            &
                  cmoq_off,         &
                  t_norb,           &
                  work,             &
                  kfree,            &
                  lwork,            &
                  print_prp_level)

      call deselect_wrk(work, lwork)

      lwork = len_f77_work
      kfree = 1
      call select_wrk(work, lwork)

      call getmdm(sum_ep_rn,        &
                  isym(icomp),      &
                  tb_dmat(1+(icomp-1)*t_n2bbasxq),  &
                  cmo,              &
                  ibeig,            &
                  tb,               &
                  work,             &
                  kfree,            &
                  lwork,            &
                  print_prp_level)

      call deselect_wrk(work, lwork)

    end do

    close(int1_unit)

  end subroutine



!-----------------------------------------------------------------------
! three subroutines from aoresponse/dirac_interface.F90
! i put them here as for testing 
! it's simpler then 'use dirac_interface'...

  SUBROUTINE set_f77_work(wrk,lwork)
     real(8),target :: wrk(:)
     integer            :: lwork, ltempwork
     logical            :: ldummy1,ldummy2
     len_f77_work = lwork
     f77_work     => wrk
!    call get_linsca_work(ldummy1,ldummy2,ltempwork,lwork)
!    call di_set_coulomb_mem_req(ltempwork)
!    ! Simen: Refine!
!    call di_set_dft_mem_req(ltempwork)
  END SUBROUTINE set_f77_work

  SUBROUTINE select_wrk(wrk, mem_req)
     real(8),pointer :: wrk(:)
     integer, intent(in) :: mem_req
!     if (mem_req .GT. len_f77_work - memory_usage) then
!    if (mem_req .GT. len_f77_work) then
!       call quit('code di_create_wrk first')
!       call di_create_wrk(wrk,mem_req)
!    else
!        wrk => f77_work((memory_usage+1):len_f77_work)
!        memory_usage = memory_usage + mem_req

        wrk => f77_work

        memory_usage  = mem_req
        mem_wrk_count = mem_wrk_count + 1
        if (mem_wrk_count > 1) call quit('counter > 1 not implemented')
!    endif
  end subroutine

  SUBROUTINE deselect_wrk(wrk, mem_req)
     real(8),pointer :: wrk(:)
     integer, intent(in) :: mem_req
!    if (mem_req .GT. len_f77_work) then
!       call di_delete_wrk(wrk)
!       call quit('code di_delete_wrk first')
!    else
        nullify(wrk)
!        memory_usage = memory_usage - mem_req
        memory_usage  = 0
        mem_wrk_count = mem_wrk_count - 1
!    endif
  end subroutine

end module
