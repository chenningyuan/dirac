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
  public quicksort
  public print_exainput
  public print_exaoutput

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

    subroutine print_exainput(ccd, lambda, moint_scheme, ncycles, t_econv, t_cholesky, level_shift, &
      talsh, exa_blocksize, talsh_buff, nao, nocc, nvir, print_level)

      logical, intent(in) :: ccd, lambda, talsh
      integer, intent(in) :: moint_scheme, ncycles, exa_blocksize, talsh_buff, nao, nocc, nvir, print_level
      real(8), intent(in) :: t_econv, t_cholesky, level_shift

      write(*,*) ""
      write(*,'(A)') "      Characteristics of this calculation      "
      write(*,'(A)') "-----------------------------------------------"
      write(*,'(A,L1)') " - Solve CCD equations :                ", ccd
      write(*,'(A,L1)') " - Solve Lambda equations :             ", lambda
      write(*,'(A,I5)') " - Integral transformation scheme : ", moint_scheme
      write(*,'(A,I5)') " - Maximum number of cycles :       ", ncycles
      write(*,'(A,ES10.1)') " - Target precision of amplitudes :  ", t_econv
      if (moint_scheme.eq.42) then
        write(*,'(A,ES10.1)') " - Cholesky decomposition threshold :",t_cholesky
      end if
      if (level_shift.gt.1.0D-14) then
        write(*,'(A,ES10.1)') " - Level Shift :                      ",level_shift
      end if
      write(*,*) "" 
      write(*,'(A,I5)') " - Number of atomic orbitals :      ", nao
      write(*,'(A,I5)') " - Number of occupied spinors :     ", nocc
      write(*,'(A,I5)') " - Number of virtual spinors :      ", nvir
      write(*,*) ""
      if (talsh) then 
        write(*,'(A)') " Using the single-node TALSH library "
        write(*,'(A,I5,A)') " - TALSH buffer size:               ",talsh_buff," GB"
      else 
        write(*,'(A)') " Using the multi-node ExaTensor library (MPI) "
        write(*,'(A,I5)') " - ExaTensor block size :           ", exa_blocksize
        write(*,'(A)') "!! ATTENTION !! Properties require TALSH, density matrix must be smaller than: "
        write(*,'(A,I5,A)') " - TALSH buffer size:              ",talsh_buff," GB"
      end if
      write(*,*) ""
      write(*,'(A,I5)') " - Print Level :                    ", print_level
      write(*,'(A)') "-----------------------------------------------"
      write(*,*) ""

    end subroutine print_exainput

    subroutine print_exaoutput(talsh, ccd, scf_energy, mp2_energy, ccsd_energy, t_energy)

      logical, intent(in)           :: talsh, ccd
      real(8), intent(in)           :: scf_energy, mp2_energy, ccsd_energy
      real(8), intent(in), optional :: t_energy(3)

      write(*,'(A)') "     -     Final results from EXACORR     -    "
      write(*,'(A)') "-----------------------------------------------"
      if (talsh) then 
        write(*,'(A)') " Using the single-node TALSH library "
      else 
        write(*,'(A)') " Using the multi-node ExaTensor library (MPI) "
      end if
      write(*,*) ""
      write (*,'(A,F20.15)') "  MP2 energy :                      ",mp2_energy
      if (ccd) then 
        write (*,'(A,F20.15)') "  Final  CCD energy :               ",ccsd_energy
      else
        write (*,'(A,F20.15)') "  Final CCSD energy :               ",ccsd_energy
      end if
      if (present(t_energy)) then 
        write (*,'(A,F20.15)') "  4th order triples correction :    ",t_energy(1)
        write (*,'(A,F20.15)') "  5th order triples (T) correction :",t_energy(2)
        write (*,'(A,F20.15)') "  5th order triples -T  correction :",t_energy(3)
        write (*,'(A,F20.15)') "  CCSD+T  correlation energy :      ",ccsd_energy+t_energy(1)
        write (*,'(A,F20.15)') "  CCSD(T) correlation energy :      ",ccsd_energy+t_energy(1)+t_energy(2)
        write (*,'(A,F20.15)') "  CCSD-T  correlation energy :      ",ccsd_energy+t_energy(1)+t_energy(3)
      end if
      write(*,*) ""
      write (*,'(A,F25.15)') "  SCF energy :                 ",scf_energy
      write (*,'(A,F25.15)') "  Total MP2 energy :           ",scf_energy+mp2_energy
      if (ccd) then 
        write (*,'(A,F25.15)') "  Total  CCD energy :          ",scf_energy+ccsd_energy
      else
        write (*,'(A,F25.15)') "  Total CCSD energy :          ",scf_energy+ccsd_energy
      end if
      if (present(t_energy)) then 
        write (*,'(A,F25.15)') "  Total CCSD+T  energy :       ",scf_energy+ccsd_energy+t_energy(1)
        write (*,'(A,F25.15)') "  Total CCSD(T) energy :       ",scf_energy+ccsd_energy+t_energy(1)+t_energy(2)
        write (*,'(A,F25.15)') "  Total CCSD-T  energy :       ",scf_energy+ccsd_energy+t_energy(1)+t_energy(3)
      end if
      write(*,*) ""
      write(*,'(A)') "-----------------------------------------------"
      write(*,*) ""

    end subroutine print_exaoutput

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
      write(*,*) '+++ write to file'
      ierr=talsh_tensor_construct(talsh_tensor,C8,h_dims,init_val=(0.D0,0.D0))
      write(*,*) 'construct err:', ierr
      ierr=exatns_tensor_get_slice(exa_tensor,talsh_tensor)
      write(*,*) 'get_slice err:', ierr
      ierr=talsh_tensor_get_body_access(talsh_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(body_p,talsh_tens,h_dims) 
      do i=1,h_dims(2)
          write(*,*) (talsh_tens(i,j), j=1,h_dims(1))
      enddo
      ierr=talsh_tensor_destruct(talsh_tensor)

      ierr = talsh_shutdown()
      
    end subroutine write_tensor

    recursive subroutine quicksort(mo_list, i_beg, i_end)

        implicit none

        integer,intent(inout) :: mo_list(:)
        integer, intent(in) :: i_beg, i_end
        integer :: i,j
        integer :: ref, tem

        i=floor((float(i_end)+float(i_beg))/2)
        ref=mo_list(i)
        i=i_beg
        j=i_end

        do
          do while (mo_list(i)<ref)
            i=i+1
          end do
          do while (mo_list(j)>ref)
            j=j-1
          end do
          if (i<j) then
            tem=mo_list(i)
            mo_list(i)=mo_list(j)
            mo_list(j)=tem
            i=i+1
            j=j-1
          else if(i.eq.j) then
            i=i+1
            EXIT
          else
            EXIT
          end if
        end do

        if (i_beg<j) call quicksort(mo_list, i_beg, j)
        if (i<i_end) call quicksort(mo_list, i, i_end)

    end subroutine quicksort

end module exacorr_utils