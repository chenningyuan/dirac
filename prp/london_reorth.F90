module london_reorth_hb

  use interface_for_pammag
  use memory_allocator 

  use dirac_cfg
  use dft_cfg
  use xc_driver


  implicit none

  public london_reort_hb

  save

  private
  real(8), pointer :: f77_work(:)       ! Fortran 77 work array pointer
  integer          :: len_f77_work  = 0 ! Length of the work array
  integer          :: memory_usage  = 0 ! Amount of work array used
  integer          :: mem_wrk_count = 0 ! Counts

  integer          :: nr_fock       = 3
  integer          :: ncomp         = 3 ! x,y,z

contains


  subroutine london_reort_hb(cmo, ibeig, work, lwork)

! calculates the reorthonormalization contribution to dH(0)/dB|B=0
! when London atomic orbitals are used:
! {T^B, H(0)} = {T^B, h(0)} + {T^B, L}
! please see two last terms in eq. (31) in JCP 131, 124119 (2009)

!-----------------------------------------------------------------------
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    integer, intent(in)    :: lwork
    real(8), intent(inout) :: work(lwork)

    real(8), allocatable   :: fao(:), fmo(:), tb_dmat(:)
    integer, allocatable   :: npos_array(:)
!-----------------------------------------------------------------------

!   set interface for f77 work array:
    call set_f77_work(work, lwork)

!   set common blocks variables used (interface_for_pammag.F90):
    call common_blocks_vars()

    call alloc(tb_dmat, p_n2bbasxq*nr_fock)

!   get (C+)(T^B)C
    call get_tb_dmat(tb_dmat, cmo, ibeig)

    call alloc(fao, p_n2bbasxq*nr_fock)

!   calculate "G" part in AO basis
    call get_gtb(fao, tb_dmat, cmo, ibeig)

    call dealloc(tb_dmat)

!   do AO-to-MO transformation on "G" part
    call alloc(fmo, p_n2orbxq*nr_fock)
    call do_ao_to_mo(fmo, fao, cmo, ibeig)

!   calculate "F" part, add to "G" and write result to file
    call get_ftb(fmo, cmo, ibeig)

    call dealloc(fmo)

  end subroutine


  subroutine get_ftb(fmat, cmo, ibeig)
!-----------------------------------------------------------------------
    real(8), intent(in)    :: fmat(p_n2orbxq, *)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork, kfree

    real(8), allocatable   :: fmo(:)
    real(8), allocatable   :: buf(:)
    real(8), allocatable   :: tb(:)

    integer :: isym(nr_fock), iopsy(nr_fock)
    integer :: isymop(nr_fock), ihrmop(nr_fock), ifckop(nr_fock)
    integer :: indxpr
    integer :: icomp, i, j, isymh
    logical :: ah
!-----------------------------------------------------------------------

  
    call alloc(fmo, p_n2bbasxq)

    lwork = len_f77_work
    kfree = 1
    call select_wrk(work, lwork)

!   get unperturbed Fock matrix in MO basis:
    call getfck(fmo,         &
                prp_iprint,  &
                work,        &
                kfree,       &
                lwork)

    call deselect_wrk(work, lwork)

    call alloc(tb,  p_n2orbxq)
    call alloc(buf, p_n2orbxq)

    call opnfil(int1_unit, 'AOPROPER', 'OLD', 'get_ftb')

!   gosia:why 8* and not 3* (now - copied from RNTLON)?
    open(fcklon_unit, file='FCKLON', form='unformatted',   &
         access='direct',recl=8*p_n2orbxq,                 &
         status='unknown')

    do icomp = 1, ncomp
      indxpr       = ip_conmat(icomp)
      isym(icomp)  = p_iprpsym(indxpr)
      iopsy(icomp) = p_jbtof(isym(icomp)-1, 1)

      call dzero(tb, p_n2orbxq)

      lwork = len_f77_work
      kfree = 1
      call select_wrk(work, lwork)

!     get T^B_icomp matrix to tb (from AOPROPER file)
      call prpmat(indxpr,        &
                  iopsy(icomp),  &
                  tb,            &
                  cmo,           &
                  ibeig,         &
                  p_icmoq,       &
                  p_norb,        &
                  work,          &
                  kfree,         &
                  lwork,         &
                  prp_iprint)

      call deselect_wrk(work, lwork)

! get final result to buf

      ah = .true.
      isymh = 1

      call dzero(buf, p_n2orbxq)

!     do 1-index transformation
      call ttra(sum_ep_rn,       &
                ah,          &
                isym(icomp),     &
                isymh,           &
                tb,              &
                fmo,             &
                buf,             &
                prp_iprint)

!     add result to fmat
      call daxpy(p_n2orbxq,      &
                 1.0d0,          &
                 buf,            &
                 1,              & 
                 fmat(1, icomp), & 
                 1)

!     write to FCKLON file:
      call wrtdac(fcklon_unit, p_n2orbxq,   &
                  fmat(1, icomp), icomp)

    end do

    close(fcklon_unit, status='keep')
    close(int1_unit,   status='keep')

    call dealloc(tb)
    call dealloc(buf)
    call dealloc(fmo)

  end subroutine


  subroutine get_gtb(fao, tb_dmat, cmo, ibeig)
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
    integer :: isym(nr_fock), iopsy(nr_fock), irep(nr_fock)
    integer :: isymop(nr_fock), ihrmop(nr_fock), ifckop(nr_fock)
    integer :: indxpr
    integer :: icomp
!-----------------------------------------------------------------------

!   set for twofck:
    do icomp = 1, ncomp
      indxpr        = ip_conmat(icomp)
      isym(icomp)   = p_iprpsym(indxpr)
      iopsy(icomp)  = p_jbtof(isym(icomp) - 1, 1)
      irep(icomp)   = isym(icomp) - 1
      isymop(icomp) = isym(icomp)
      ihrmop(icomp) = -1
      ifckop(icomp) = p_jbtof(isym(icomp) - 1, 1)
    end do

!gosia fixme parallel run
    call stex_get_npos(npos)
    call alloc(npos_array, npos)

    lwork = len_f77_work
    call select_wrk(work, lwork)

    call twofck(isymop,      &
                ihrmop,      &
                ifckop,      &
                fao,         &
                tb_dmat,     &
                nr_fock,     &
                npos_array,  &
                nmr_int,     &
                prp_iprint,  &
                work,        &
                lwork)

    call deselect_wrk(work, lwork)

    if (allocated(npos_array)) then
      call dealloc(npos_array)
    end if

!       construct kohn-sham contribution
!       ================================

! check if sdft, as only in this case this contribution is != 0
    if (dirac_cfg_dft_calculation .and. .not. dft_cfg_no_sdft) then

      call alloc(dmat, p_n2bbasxq)
      call genden(dmat, cmo, 1, prp_iprint)

         call xc_mpi_wake_up_coworkers()
         call xc_integrator(xc_mat_dim           = nr_ao_total,   &
                            xc_nz                = nr_qu,         &
                            xc_dmat_0            = dmat,          &
                            xc_nr_dmat           = nr_fock,       &
                            xc_nr_fmat           = nr_fock,       &
                            xc_dmat              = tb_dmat,       &
                            xc_fmat              = fao,           &
                            xc_fmat_pg_sym       = isymop,        &
                            xc_dmat_pg_sym       = isymop,        &
                            xc_dmat_ih_sym       = ihrmop,        &
                            xc_response_order_mo = 1,             &
                            xc_do_london_rhs_ro  = .true.)       
         call xc_mpi_release_coworkers()

      call dealloc(dmat)

    end if

  end subroutine get_gtb



  subroutine do_ao_to_mo(fmo, fao, cmo, ibeig)
!-----------------------------------------------------------------------
    real(8), intent(inout) :: fmo(p_n2orbxq , *)
    real(8), intent(in)    :: fao(p_n2bbasxq, *)

    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork

    integer :: itim(nr_fock), isym(nr_fock), iopsy(nr_fock), irep(nr_fock)
    integer :: isymop(nr_fock), ihrmop(nr_fock), ifckop(nr_fock)
    integer :: indxpr
    integer :: icomp, i, j
!-----------------------------------------------------------------------

    do icomp = 1, ncomp

      indxpr        = ip_conmat(icomp)
      isym(icomp)   = p_iprpsym(indxpr)
      iopsy(icomp)  = p_jbtof(isym(icomp) - 1, 1)
      irep(icomp)   = isym(icomp) - 1

      do i = 1, p_nfsym 
        j = mod(i + iopsy(icomp), 2) + 1

        lwork = len_f77_work
        call select_wrk(work, lwork)

        call qtrans('AOMO',                          &
                    'S',                             &
                    0.0d0,                           &
                    p_nfbas(i, 0),                   &
                    p_nfbas(j, 0),                   &
                    p_norb(i),                       &
                    p_norb(j),                       &
                    fao(1 + p_i2basx(i, j), icomp),  &
                    nr_ao_total,                     &
                    nr_ao_total,                     &
                    nr_qu,                           &
                    p_ipqtoq(1, irep(icomp)),        &
                    fmo(1 + p_i2orbx(i, j),icomp),   &
                    p_norbt,                         &
                    p_norbt,                         &
                    nr_qu,                           &
                    p_ipqtoq(1, irep(icomp)),        &
                    cmo(1 + p_icmoq(i)),             &
                    p_nfbas(i, 0),                   &
                    p_norb(i),                       &
                    nr_qu,                           &
                    p_ipqtoq(1, 0),                  &
                    cmo(1 + p_icmoq(j)),             &
                    p_nfbas(j, 0),                   & 
                    p_norb(j),                       &
                    nr_qu,                           &
                    p_ipqtoq(1, 0),                  &
                    work,                            &
                    lwork,                           &
                    prp_iprint)

        call deselect_wrk(work, lwork)

      end do
    end do

  end subroutine


  subroutine get_tb_dmat(tb_dmat, cmo, ibeig)
!-----------------------------------------------------------------------
    real(8), intent(out)   :: tb_dmat(p_n2bbasxq, *)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    real(8), pointer       :: work(:)
    integer                :: lwork, kfree

    real(8), allocatable   :: tb(:)
    integer                :: isym(nr_fock), iopsy(nr_fock)
    integer                :: indxpr, icomp 
!-----------------------------------------------------------------------

    call opnfil(int1_unit, 'AOPROPER', 'OLD', 'get_tb_dmat')
    call alloc(tb, p_n2orbxq)

    do icomp = 1, ncomp
      indxpr       = ip_conmat(icomp)
      isym(icomp)  = p_iprpsym(indxpr)
      iopsy(icomp) = p_jbtof(isym(icomp)-1, 1)
      !write(*, *) 'alamakota icomp, indxpr, isym, jbtof(irep,1) = ', &
      !             icomp, indxpr, isym(icomp), p_jbtof(isym(icomp)-1, 1)

      call dzero(tb, p_n2orbxq)

      lwork = len_f77_work
      kfree = 1
      call select_wrk(work, lwork)

!     get T^B_icomp matrix into tb
      call prpmat(indxpr,           &
                  iopsy(icomp),     &
                  tb,               &
                  cmo,              &
                  ibeig,            &
                  p_icmoq,          &
                  p_norb,           &
                  work,             &
                  kfree,            &
                  lwork,            &
                  prp_iprint)

      call deselect_wrk(work, lwork)

      lwork = len_f77_work
      kfree = 1
      call select_wrk(work, lwork)

      call getmdm(sum_ep_rn,         &
                  isym(icomp),       &
                  tb_dmat(1, icomp), &
                  cmo,               &
                  ibeig,             &
                  tb,                &
                  work,              &
                  kfree,             &
                  lwork,             &
                  prp_iprint)

      call deselect_wrk(work, lwork)

    end do

    close(int1_unit)

  end subroutine get_tb_dmat



!-----------------------------------------------------------------------
! three subroutines from aoresponse/dirac_interface.F90

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
