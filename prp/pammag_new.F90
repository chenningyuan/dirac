module london_response

  use interface_for_pammag
  use memory_allocator 

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

! ------------------ test part... --------------------------------------

  subroutine rntlon_new(cmo, ibeig, work, lwork)
!-----------------------------------------------------------------------
! test routine
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)

    integer, intent(in)    :: lwork
    real(8), intent(inout) :: work(lwork)

    real(8), allocatable   :: h0tb1(:)
    integer                :: i
!-----------------------------------------------------------------------
!initialize: len_f77_work=lwork, f77_work points to work array(work(kfree))
    call set_f77_work(work, lwork)

    call common_blocks_vars()
    call alloc(h0tb1, nr_ao_total*nr_ao_total*nr_qu)
    call test_mem_problems_g(h0tb1, cmo, ibeig)

!needed for getdhb:
    open(fcklon_unit,file='FCKLON',form='unformatted',         &
         access='direct',recl=8*nr_ao_large*nr_ao_large*nr_qu, &
         status='unknown')
    do i = 1, 3
      call wrtdac(fcklon_unit,nr_ao_large*nr_ao_large*nr_qu,   &
                  (/0.0d0/),i)
    enddo
    close(fcklon_unit, status='keep')

    call dealloc(h0tb1)

  end subroutine


  subroutine test_mem_problems_g(h0tb1, cmo, ibeig)
    real(8), intent(inout) :: h0tb1(*)
    real(8), intent(in)    :: cmo(*)
    integer, intent(in)    :: ibeig(*)
    real(8), pointer       :: work(:)
    integer                :: lwork, i
    real(8), allocatable   :: tab(:)

    call alloc(tab, nr_ao_total*nr_ao_total*nr_qu)

    lwork = len_f77_work
    call select_wrk(work, lwork)

!   call getfck(tab,print_prp_level,work,1,lwork)
!for line above i get valgrind msg:
! "Bad permissions for mapped region ..."
! which then points to prp/pamset.F line 469 (call memget)
! and then to gp/mempkg.F line 146
! i'll write details in mail

!radovan: this one works for me:
    call getfck_wrap(tab,print_prp_level,work,lwork)

    call deselect_wrk(work, lwork)

    call dcopy(nr_ao_total*nr_ao_total*nr_qu,    &
               tab,1,h0tb1,1)

    call dealloc(tab)

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
