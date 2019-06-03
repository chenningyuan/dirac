module exacorr_utils

!  Collection of simple utility functions used in exacorr

  implicit none
  private

  public relative_error
  public print_date
  public get_free_fileunit
  public print_exacorr_logo
  public print_iteration
  public make_spinor_list

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

      subroutine make_spinor_list ( nesh, exa_input, string_occupied, & 
        string_occ_beta, string_virtual, string_vir_beta)

        use exacorr_datatypes

        implicit none

        integer, intent(in)              :: nesh                     ! number of positive energy solutions
        type(exacc_input), intent(inout) :: exa_input
        character(len = *), intent(in)   :: string_occupied
        character(len = *), intent(in)   :: string_virtual
        character(len = *), intent(in)   :: string_occ_beta
        character(len = *), intent(in)   :: string_vir_beta

        integer              :: kvec                     ! Dummy variable needed in first call to numlst
        integer, allocatable :: mobe_occ(:), mobe_vir(:)
        integer              :: i
        integer              :: nmo(4)
        integer              :: imo, ialpha, ibeta

        if (exa_input%beta_occ) then
          if (exa_input%beta_vir) then
            nmo(1:4)=0
            call numlst (string_occupied,kvec,1,1,nesh,1,nmo(1))
            call numlst (string_occ_beta,kvec,1,1,nesh,1,nmo(2))
            call numlst (string_virtual,kvec,1,1,nesh,1,nmo(3))
            call numlst (string_vir_beta,kvec,1,1,nesh,1,nmo(4))

            exa_input%nkr_occ=nmo(1)
            allocate (exa_input%mokr_occ(exa_input%nkr_occ))
            allocate (mobe_occ(nmo(2)))
            exa_input%nkr_vir=nmo(3)
            allocate (exa_input%mokr_vir(exa_input%nkr_vir))
            allocate (mobe_vir(nmo(4)))

            call numlst(string_occupied,exa_input%mokr_occ, &
                  nmo(1),1,nesh,1,nmo(1))
            call numlst(string_occ_beta,mobe_occ, &
                  nmo(2),1,nesh,1,nmo(2))
            call numlst(string_virtual,exa_input%mokr_vir, &
                  nmo(3),1,nesh,1,nmo(3))
            call numlst(string_vir_beta,mobe_vir, &
                  nmo(4),1,nesh,1,nmo(4))

            if(exa_input%mokr_occ(1).ne.mobe_occ(1)) then
              print *, "core orbital need to be closed shell"
              call quit ("change .OCC_BETA in EXACORR")
            end if
            
            exa_input%nocc=nmo(1)+nmo(2)
            exa_input%nvir=nmo(3)+nmo(4)
            allocate (exa_input%mo_occ(exa_input%nocc))
            allocate (exa_input%mo_vir(exa_input%nvir))

            imo=0
            do i=1, max(nmo(1), nmo(2))
              if(i.le.nmo(1)) then
                imo=imo+1
                exa_input%mo_occ(imo)=exa_input%mokr_occ(i)*2-1
              end if
              if(i.le.nmo(2)) then
                imo=imo+1
                exa_input%mo_occ(imo)=mobe_occ(i)*2
              end if
            end do
         
            ialpha=1
            ibeta=1
            imo=0
            do i=1, nmo(3)+nmo(4)
              if(ibeta.gt.nmo(4)) then
                if (ialpha.le.nmo(3)) then
                  imo=imo+1
                  exa_input%mo_vir(imo)=exa_input%mokr_vir(ialpha)*2-1
                  ialpha=ialpha+1
                end if
              else if (ialpha.gt.nmo(3)) then
                if (ibeta.le.nmo(4)) then
                  imo=imo+1
                  exa_input%mo_vir(imo)=mobe_vir(ibeta)*2
                  ibeta=ibeta+1
                end if
           else if(exa_input%mokr_vir(ialpha).eq.mobe_vir(ibeta)) then
                imo=imo+1
                exa_input%mo_vir(imo)=exa_input%mokr_vir(ialpha)*2-1
                ialpha=ialpha+1
                imo=imo+1
                exa_input%mo_vir(imo)=mobe_vir(ibeta)*2
                ibeta=ibeta+1
           else if (exa_input%mokr_vir(ialpha).lt.mobe_vir(ibeta)) then
                imo=imo+1
                exa_input%mo_vir(imo)=exa_input%mokr_vir(ialpha)*2-1
                ialpha=ialpha+1
              else                   
                imo=imo+1
                exa_input%mo_vir(imo)=mobe_vir(ibeta)*2
                ibeta=ibeta+1
              end if
            end do

            deallocate (mobe_occ)
            deallocate (mobe_vir)
          else
            call quit ("also define .OCC_BETA in EXACORR")
          end if
        else if (exa_input%beta_vir) then
          call quit ("also define .VIR_BETA in EXACORR")
        else
          exa_input%nkr_occ = 0
          call numlst (string_occupied,kvec,1,1,nesh,1, &
              exa_input%nkr_occ)
          allocate (exa_input%mokr_occ(exa_input%nkr_occ))
          call numlst (string_occupied,exa_input%mokr_occ, &
              exa_input%nkr_occ, 1,nesh,1,exa_input%nkr_occ)
          exa_input%nocc = 2 * exa_input%nkr_occ
          allocate (exa_input%mo_occ(exa_input%nocc))
          do i = 1, exa_input%nocc
            exa_input%mo_occ(i) = (exa_input%mokr_occ(1)-1)*2+i
          end do

          exa_input%nkr_vir = 0
          call numlst (string_virtual,kvec,1,1,nesh,1, &
              exa_input%nkr_vir)
          allocate (exa_input%mokr_vir(exa_input%nkr_vir))
          call numlst (string_virtual,exa_input%mokr_vir, &
              exa_input%nkr_vir,1,nesh,1,exa_input%nkr_vir)
          exa_input%nvir = 2 * exa_input%nkr_vir
          allocate (exa_input%mo_vir(exa_input%nvir))
          do i = 1, exa_input%nvir
            exa_input%mo_vir(i) = &
                  exa_input%mo_occ(exa_input%nocc)+i
          end do
        end if

      end subroutine make_spinor_list

end module exacorr_utils