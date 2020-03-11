module exacorr_utils

!  Collection of simple utility functions used in exacorr

  implicit none
  private

  public relative_error
  public print_date
  public get_free_fileunit
  public print_exacorr_logo
  public print_iteration
  public write_tensor

  contains

    real(8) function relative_error(quantity_1,quantity_2)
    ! Gives percentage error for values > 1, otherwise absolute error

      real(8),intent(in):: quantity_1,quantity_2
      real(8)           :: error

      error = quantity_1 - quantity_2

      if (abs(quantity_1).gt.1.0) then
        relative_error = abs (error / quantity_1)
      else
        relative_error = abs (error)
      end if

    end function relative_error
       
    subroutine print_date(text)
    ! Prints time and date with a string
      character(len=10) :: datex
      character(len=8)  :: timex
      character(len=*)  :: text
      integer           :: lentext
      lentext = len(text)
      call daytime(datex,timex)
      write (*,'(1x,a10,1x,a8,1x,a)') datex,timex,text(1:lentext)
    end subroutine print_date

    subroutine get_free_fileunit(unit)
    ! idea from http://fortranwiki.org, will become obsolete once fortran2008 is supported everywhere
    ! gets free file unit
      integer, intent(out) :: unit
      integer, parameter   :: LUN_MIN=100, LUN_MAX=1000
      logical :: opened
      integer :: lun
      unit=-1
      do lun=LUN_MIN,LUN_MAX
        inquire(unit=lun,opened=opened)
        if (.not. opened) then
          unit=lun
          exit
        end if
      end do
    end subroutine get_free_fileunit

    subroutine print_exacorr_logo()
    ! Prints logo and authors of exacorr

      print*, ""
      print*, ""
      print*, "                   *****************************************************"
      print*, "                   ***     Entering the ExaCorr module in DIRAC      ***"
      print*, "                   ***                                               ***"
      print*, "                   ***  Authors:         - Lucas Visscher            ***"
      print*, "                   ***                   - Anastasios Papadopoulos   ***"
      print*, "                   ***  Contributors:    - Johann Potoschnig         ***"
      print*, "                   ***                   - Andre S. P. Gomes         ***"
      print*, "                   ***  Features:        - CCD/CCSD                  ***"
      print*, "                   ***                   - 1DM                       ***"
      print*, "                   ***                   - ReSpect interface         ***"
      print*, "                   ***                                               ***"
      print*, "                   ***  Relativistic Coupled Cluster code using the  ***"
      print*, "                   ***  TALSH and ExaTensor libraries developed by   ***"
      print*, "                   ***  Dmitry Lyakh                                 ***"
      print*, "                   *****************************************************"
      print*, ""
      print*, ""

    end subroutine print_exacorr_logo

    subroutine print_iteration(iteration,convergence,energy)
    ! Formatted print for iterations

      ! input variables
      integer, intent(in)           :: iteration
      real(8), intent(in)           :: convergence
      real(8), intent(in), optional :: energy

      if (iteration.eq.1 .and. present(energy)) then
        write(*,*) ""
        write(*,*) "-----------------------------------------------"
        write(*,*) "  It.     Energy                  Convergence  "
        write(*,*) "-----------------------------------------------"
        write(*,'(A3,I2,A4,F20.15,A3,E10.1)') "   ", iteration, "    ", energy, "   ", convergence
      else if (iteration.eq.1) then
        write(*,*) ""
        write(*,*) "------------------------"
        write(*,*) "  It.      Convergence  "
        write(*,*) "------------------------"
        write(*,'(A3,I2,A4,E10.1)') "   ", iteration, "    ", convergence
      else if (iteration.ne.1 .and. present(energy)) then
        write(*,'(A3,I2,A4,F20.15,A3,E10.1)') "   ", iteration, "    ", energy, "   ", convergence
      else if (iteration.ne.1) then
        write(*,'(A3,I2,A4,E10.1)') "   ", iteration, "    ", convergence
      end if

    end subroutine print_iteration

    subroutine write_tensor(exa_tensor, i_dims)

      use exatensor
      use talsh
      use tensor_algebra

      type(tens_rcrsv_t)              :: exa_tensor
      integer(INTL), dimension(2)     :: i_dims
      integer(INTD), dimension(2)     :: h_dims
      type(talsh_tens_t)              :: talsh_tensor
      complex(8), pointer, contiguous :: talsh_tens(:,:)
      type(C_PTR)                     :: body_p
      integer                         :: i,j
      integer(INTD)                   :: ierr
      integer(C_SIZE_T)               :: buf_size=1_8*1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
      integer(C_INT)                  :: host_arg_max

      ierr=talsh_init(buf_size,host_arg_max)
      call print_date('Initialized talsh library')
      write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max

      h_dims=i_dims
      print *, '+++ write to file'
      ierr=talsh_tensor_construct(talsh_tensor,C8,h_dims,init_val=(0.D0,0.D0))
      print *, 'construct err:', ierr
      ierr=exatns_tensor_get_slice(exa_tensor,talsh_tensor)
      print *, 'get_slice err:', ierr
      ierr=talsh_tensor_get_body_access(talsh_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(body_p,talsh_tens,h_dims) 
      do i=1,h_dims(2)
          print*, (talsh_tens(i,j), j=1,h_dims(1))
      enddo
      ierr=talsh_tensor_destruct(talsh_tensor)

      ierr = talsh_shutdown()
      
    end subroutine write_tensor

end module exacorr_utils